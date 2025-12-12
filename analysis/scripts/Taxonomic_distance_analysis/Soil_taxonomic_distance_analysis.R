###taxonomic distance comparison
setwd("~/Dropbox_Harvard/hutlab/Jacob/Projects/Enviro_revs/environmental_methods_comparison/")
source("analysis/scripts/helpers.R")

library(dplyr)
#make tree for gut data
dataset_soil <- "soil"
#set level to species
level <- 7
#set abundance threshold to 0.1%
threshold <- .1

#generate the profiles and only keep NCBI taxa
profiles <- remove_non_ncbi(preprocess(dataset_soil, level, FALSE))
#grab the truth samples and remove the unclassified portion
truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset_soil, level)))
#normalize them to TSS
profiles = lapply(profiles, renormalize)
## apply a threshold if there is one, threshold is a global variable
profiles = lapply(profiles, threshold_sample, threshold)
## remove rows that sum to 0
profiles = lapply(profiles, function(x){
  if(length(which(rowSums(x[,-1])==0))!=0)
    return(x[-which(rowSums(x[,-1])==0),])
  else
    return(x)
})
## re-nomarlize after applying the threshold
profiles = lapply(profiles, renormalize)

label_converts <-  c('2599805'=931866, '1404649'=2823807, '2735433'=3014751, '2829818'=2952571, '2840472'=2840469, 
                     '2884447'=2792859, '80870'=80867, '1532'=33035, '2778071'=29523, '285567'=1942, '290385'=104623, '37480'=67362,
                     '423539'=3074428, '50340'=53407, '84292'=162393, '88075'=88074, '96101'=46163, '1497613'=1505087, 
                     '1915400'=67332, '335659'=1404864, '90270'=2754056, '80870'=80867, '1532'=33035, 
                     '285567'=1942, '290385'=104623, '37480'=67362, '423539'=3074428, '50340'=53407,
                     '84292'=162393, '88075'=88074, '96101'=46163, '335659'=1404864, '90270'=2754056)
###update labels needs to search if they convert to the same and sum them
truth <- update_labels(truth, label_converts)
profiles <- lapply(profiles,function(x) update_labels(x,label_converts))

taxIDs <- c()
for(i in profiles){
  taxIDs <- c(taxIDs, i$TaxIDs)
}
all_detected_IDs <- unique(taxIDs)


##load in the tree
library(ape)
tree_soil <- read.tree("analysis/databases/Taxonomic_trees/soil_simulation_tree.nwk")

tree_soil$edge.length <- rep(1, nrow(tree_soil$edge))

tree_soil$tip.label <- gsub("[[:punct:] ]", "", tree_soil$tip.label)

missing_IDs <- c()
for(i in profiles){
  missing <- which(is.na(match(i$TaxIDs, tree_soil$tip.label)))
  missing_names <- i$TaxIDs[missing]
  missing_IDs <- c(missing_IDs, missing_names)
}

##double check that these all should be removed.
library(taxonomizr)
##they all viruses so they need to be remove and then renormalize the profile again.
getTaxonomy(missing_IDs, desiredTaxa = "acellular root")
profiles <- lapply(profiles, function(x) remove_specific(x, missing_IDs))
profiles <- lapply(profiles, renormalize)

colnames(truth)[1] <- "TaxIDs"

##okay profiles and tree is now preppred we can now calculate the distances.
##remove truth samples that sum to 0
remove_samples <- names(which(colSums(truth[,-1])==0))
##all of these samples have an expected 100% unclassified porition 
truth <- truth %>% dplyr::select(-remove_samples)
profiles <- lapply(profiles, function(x) dplyr::select(x,-remove_samples))

#calculate the distances for each tool
distances <- lapply(profiles, function(x) calc_unifrac(truth, x, tree=tree_soil))
names(distances) <- c("centrifuge", "kraken", "MPA2", "MPA3", "MPA4", "metaxa", "motus", "gtdbtk_megahit",
                      "gtdbtk_metaspades", "phylophlan_megahit", "phylophlan_metaspades")

distances_df <- do.call(cbind, distances)
distances_df[which(is.na(distances_df))] <- 1
distances_df <- data.frame(distances_df, check.name=F)
distances_df[,"Sample"] <- rownames(distances_df)
distances_df_melt <- reshape2::melt(distances_df)

library(ggplot2)
distances_df_melt %>% ggplot(aes(x=variable, y=value)) + geom_boxplot() +
  ylab("Taxonomic - Weighted Unifrac Distance") + xlab("Method") +
  theme_bw(base_size = 12)

### only for "core" samples





##divive up by sample type

distances_df_melt$Species_count <- str_extract(distances_df_melt$Sample, "\\.n\\.[0-9]*")
distances_df_melt$Sample_size <- str_extract(distances_df_melt$Sample, "\\.size\\..*k")
distances_df_melt$Unknown_amount <- str_extract(distances_df_melt$Sample, "\\.up\\..*mut")
distances_df_melt$mut_rate <- str_extract(distances_df_melt$Sample, "\\_rate.*_")
distances_df_melt$genome_skew <- str_extract(distances_df_melt$Sample, "k\\..*sig")
distances_df_melt$sigma <- str_extract(distances_df_melt$Sample, "ma\\..*up")


distances_df_melt_longer <- tidyr::pivot_longer(distances_df_melt,
                                                cols=c("Species_count", "Sample_size", "Unknown_amount", "mut_rate"),
                                                names_to="facet",
                                                values_to="facet_val")


### filter to core samples and plot

core_samples <- distances_df_melt %>% filter(Species_count==".n.300") %>% filter(Sample_size==".size.7.5.k") %>%
  filter(mut_rate=="_rate.0.0_sample_") %>% filter(genome_skew=="k.100.sig") %>% filter(sigma=="ma.1.0.up")
core_samples <- core_samples %>% filter(Unknown_amount==".up.0.75.mut")



core_samples %>% filter(variable %in% c("kraken", "MPA4", "motus", "gtdbtk_megahit", "phylophlan_megahit")) %>% 
  ggplot(aes(x=variable, y=1-value, fill=variable)) + geom_point(size=3, pch=21) + ylab("1 - Weighted Taxonomic distance") +
  xlab("Method") +
  theme_linedraw() +
  theme(text = element_text(size=15)) +
  scale_fill_manual(values=c("kraken"="#A15BE4", "MPA4"="#142755", "motus"="#81C784", 
                              "gtdbtk_megahit"="#FF7300", "phylophlan_megahit"="#AC0911")) +
  scale_x_discrete(labels = c("kraken" = "Kraken 2 / Braken 2",
                              "MP4" = "MetaPhlAn 4",
                              "motus" = "mOTUs 3",
                              "gtdbtk_megahit" = "GTDB-Tk MEGAHIT",
                              "phylophlan_megahit"= "PhyloPhlAn MEGAHIT")) +
  ggtitle("Taxonomic Distances Soil") +
  theme(axis.text.x=element_text(angle=90)) +
  labs(fill="Method") +
  ylim(0,1)
  

ggsave("analysis/figures/taxonomic_distances/soil_2B.pdf", width=5, height=5)



distances_df_melt_longer$facet_val

distances_df_melt_longer %>% ggplot(aes(x=facet_val, y=value, group=variable, color=variable)) + geom_point() +
  facet_grid(cols=vars(facet), scales="free_x") + geom_path()

keep_tools <- c("kraken", "MPA4", "motus", "gtdbtk_megahit", "phylophlan_megahit")

distances_df_melt_longer %>% filter(variable%in%keep_tools) %>% ggplot(aes(x=facet_val, y=value, group=variable, color=variable)) + geom_point() +
  facet_grid(cols=vars(facet), scales="free_x") + geom_line() + theme_bw(base_size=12) + geom_line()


#plot for mut_rate
mut_count_df <- distances_df_melt %>% filter(Species_count==".n.300") %>% 
  filter(Sample_size==".size.7.5.k") %>% filter(Unknown_amount==".up.0.5.mut") %>%
  filter(genome_skew=="k.100.sig") %>% filter(sigma=="ma.1.0.up")

mut_count_df$mut_rate_fix <- gsub("[^0-9.]", "", mut_count_df$mut_rate)
mut_count_df$mut_rate_fix <- gsub("^\\.", "", mut_count_df$mut_rate_fix)
mut_count_df$mut_rate_fix <- as.numeric(mut_count_df$mut_rate_fix)

p1 <- mut_count_df %>% filter(variable%in%keep_tools) %>% ggplot(aes(x=mut_rate_fix, y=1-value, group=variable, color=variable)) + 
 theme_bw() + geom_smooth() + ylab("1 - Weighted Taxonomic distance") +
  xlab("Mutation rate") +
  geom_point()




### Now do it for species count

spec_count_df <- distances_df_melt %>% filter(mut_rate=="_rate.0.0_sample_") %>% 
  filter(Sample_size==".size.7.5.k") %>% filter(Unknown_amount==".up.0.5.mut") %>%
  filter(genome_skew=="k.100.sig") %>% filter(sigma=="ma.1.0.up")


spec_count_df$Species_count <- as.numeric(gsub("\\.n\\.", "", spec_count_df$Species_count))

p2 <- spec_count_df %>% filter(variable%in%keep_tools) %>% ggplot(aes(x=Species_count, y=1-value, group=variable, color=variable)) + 
  theme_bw() + geom_smooth() + ylab("1 -  Weighted Taxonomic distance") +
  xlab("Species Count") +
  geom_point()



## sample size
sample_count_df <- distances_df_melt %>% filter(mut_rate=="_rate.0.0_sample_") %>% 
  filter(Species_count==".n.300") %>% filter(Unknown_amount==".up.0.5.mut") %>%
  filter(genome_skew=="k.100.sig") %>% filter(sigma=="ma.1.0.up")


sample_count_df$Sample_size <- gsub("\\.size\\.", "", sample_count_df$Sample_size)
sample_count_df$Sample_size <- factor(sample_count_df$Sample_size, levels=c("0.05.k", "0.5.k",  "1.5.k",  "7.5.k", "30.0.k"))

p3 <- sample_count_df %>% filter(variable%in%keep_tools) %>% ggplot(aes(x=Sample_size, y=1-value, group=variable, color=variable)) + 
  theme_bw() + geom_smooth() + ylab("1 - Weighted Taxonomic Distance") +
  xlab("Sample Size") +
  geom_point()

## unknown %

unknown_count_df <- distances_df_melt %>% filter(mut_rate=="_rate.0.0_sample_") %>% 
  filter(Species_count==".n.300") %>% filter(Sample_size==".size.7.5.k") %>%
  filter(genome_skew=="k.100.sig") %>% filter(sigma=="ma.1.0.up")


unknown_count_df$Unknown_amount <- gsub("\\.up\\.", "", unknown_count_df$Unknown_amount)
unknown_count_df$Unknown_amount <- gsub("\\.mut", "", unknown_count_df$Unknown_amount)
unknown_count_df$Unknown_amount <- as.numeric(unknown_count_df$Unknown_amount)



p4 <- unknown_count_df %>% filter(variable%in%keep_tools) %>% ggplot(aes(x=Unknown_amount, y=1-value, group=variable, color=variable)) + 
  theme_bw() + geom_smooth() + ylab("1 - Weighted Taxonomic Distance") +
  xlab("Unknown %") +
  geom_point()

library(patchwork)
p1 + p2 + p3 + p4 + plot_layout(guides="collect")


distances <- lapply(profiles, function(x) calc_bc(truth, x))
names(distances) <- c("centrifuge", "kraken", "MPA2", "MPA3", "MPA4", "metaxa", "motus", "gtdbtk_megahit",
                      "gtdbtk_metaspades", "phylophlan_megahit", "phylophlan_metaspades")

distances_df <- do.call(cbind, distances)
distances_df[which(is.na(distances_df))] <- 1

distances_df_melt <- reshape2::melt(distances_df)
library(ggplot2)
distances_df_melt %>% ggplot(aes(x=Var2, y=value)) + geom_boxplot() +
  ylab("Taxonomic - Weighted Unifrac Distance") + xlab("Method") +
  theme_bw(base_size = 12)
