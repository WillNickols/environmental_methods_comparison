rm(list = ls())
library(reshape2)
library(ggplot2)
library(Maaslin2)
library(data.table)
library(glue)
library(Rtsne)
library(ggbeeswarm)
library(ggtext)

source("analysis/scripts/helpers.R")

preprocess_all_real()

threshold = 0.05
tool_names = c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
tool_core = c("Kraken 2 / Bracken 2", "MetaPhlAn 4", "mOTUs 3", "GTDB-Tk MEGAHIT", "PhyloPhlAn MEGAHIT")
colAdd <- c("#CF9FFF", "#A15BE4", "#00FFFF", "#0061FE", "#1434A4", "#81C784", "#2E7D32", "#FF7300", "#FBB15B", "#AC0911", "#E95420")
col = scale_color_manual(values = colAdd)
fil = scale_fill_manual(values = colAdd)

metadata_overview <- function() {
  read_info <- read.csv("metadata/sample_info/read_info.tsv", sep="\t")
  df_to_plot <- data.frame(table(read_info$Study))
  colnames(df_to_plot) <- c("Study", "Sample number")
  df_to_plot$Study <- sapply(tolower(df_to_plot$Study), simpleCap)
  df_to_plot$`Metadata variables` <- 0
  for (file in list.files("metadata")) {
    if (file.exists(paste0("metadata/", file, "/merged_metadata.tsv"))) {
      study_tmp <- case_when(file == "forest_soil" ~ "Forest soil",
                             file == "human" ~ "Human gut",
                             file == "acid_mine" ~ "Acid mine drainage",
                             file == "animal_gut" ~ "Wild animal gut",
                             file == "gator_soil" ~ "Gator nest",
                             file == "saltmarsh" ~ "Salt marsh",
                             file == "tara_polar" ~ "Tara polar")
      df_to_plot$`Metadata variables`[df_to_plot$Study == study_tmp] <- length(colnames(read.csv(paste0("metadata/", file, "/merged_metadata.tsv"), sep="\t"))) - 1
    }
  }
  
  df_to_plot$Study[df_to_plot$Study == "Tara polar"] <- "Polar ocean"
  
  sf <- max(df_to_plot$`Sample number`) / max(df_to_plot$`Metadata variables`)
  df_to_plot$`Metadata variables` <- df_to_plot$`Metadata variables` * sf
  df_to_plot <- reshape2::melt(df_to_plot)
  df_to_plot$value <- df_to_plot$value / sf

  ggplot(df_to_plot, aes(x=factor(Study, rev(c("Human gut", "Cat gut", "Dog gut", "Wild animal gut", "Forest soil", "Gator nest", "Acid mine drainage", "Salt marsh", "Coastal sediment", "Polar ocean"))))) +
    geom_bar( aes(y = value, fill = factor(variable, levels = c("Metadata variables", "Sample number"))),
              stat="identity", position=position_dodge(),
              color="black", alpha=.6)  +
    scale_fill_manual(values = c("blue", "red")) +
    scale_y_continuous(name = "Metadata variables",labels = scales::comma,sec.axis = sec_axis(~.*sf, name="Sample number",
                                                                              labels = scales::comma))+
    labs(fill='')+
    theme_bw()+
    theme(legend.position = 'bottom',
          legend.direction = 'vertical',
          text = element_text(size = 20)
          ) + 
    xlab("") +
    coord_flip()
  
  dir.create(file.path("figures/overview_figures/"), showWarnings = FALSE)
  ggsave(paste0("figures/overview_figures/sample_number_and_metadata.png"), width=13, height=20, units = 'cm', dpi=1000)
}
metadata_overview()

PCoABray_all <- function() {
  set.seed(1)
  
  dist_mat <- read.csv("real_data_outputs/inputs/overall_dists/mat_abundance_jaccard.csv", sep = ";")
  id_names <- scan("real_data_outputs/inputs/overall_dists/in.txt", character(), quote = "", sep = "\n")
  conversion_df <- data.frame("IDs" = gsub("\\:.*", "", id_names), 
                             "Sample" = gsub(".* ", "", id_names) %>% 
                               gsub(pattern=".*\\/", replacement="") %>% 
                               gsub(pattern="_.*", replacement="") %>% 
                               gsub(pattern="\\..*", replacement="") %>% 
                               gsub(pattern="^re", replacement=""))
  dist_mat$X <- NULL
  colnames(dist_mat) <- conversion_df$Sample[match(colnames(dist_mat), conversion_df$IDs)]
  
  pc <- capscale(dist_mat~1, comm = dist_mat, na.action = na.omit)
  cap = data.frame(pc$CA$u)
  
  dataset_mapping <- read.csv("metadata/sample_info/read_info.tsv", sep="\t")
  cap$dataset <- dataset_mapping$Study[match(rownames(cap), dataset_mapping$ID)]

  s = summary(pc)
  
  colAdd <- c("#FC6C85", "#EE0911", "#FBB15B", "#964B00", "#2E7D32", "#81C784", "#FF7300", "#00FFFF", "#0061FE", "#1434A4")
  col = scale_color_manual(values = colAdd)
  
  tsne <- Rtsne(dist_mat, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
  cap[,1:2] <- tsne$Y

  # Plot ordination
  ggplot(cap, aes(MDS1, MDS2)) + 
    geom_point(aes(color=factor(dataset, levels = c("Human Gut", "Cat Gut", "Dog Gut", "Wild Animal Gut", "Forest Soil", "Gator Nest", "Acid Mine Drainage", "Salt Marsh", "Coastal Sediment", "Tara Polar"))), size = 3, alpha = 0.7) + 
    theme_bw() +
    labs(color = "Dataset") + 
    labs(x = "tSNE Axis 1",
         y = "tSNE Axis 2") +
    col + 
    theme(legend.position = "none",
          text = element_text(size = 16)
          )

  dir.create(file.path("figures/PCoABray/all_datasets/"), showWarnings = FALSE)
  ggsave(paste0("figures/PCoABray/all_datasets/all_datasets_tSNE.png"), width=14, height=14, units = 'cm', dpi=1000)
}
PCoABray_all()

read_num_all <- function() {
  sample_sizes <- read.csv("metadata/sample_info/read_info.tsv", sep="\t", check.names = F)
  sample_sizes$Study <- sapply(tolower(sample_sizes$Study), simpleCap)
  sample_sizes$Study[sample_sizes$Study == "Tara polar"] <- "Polar ocean"
  
  colAdd <- c("#FC6C85", "#EE0911", "#FBB15B", "#964B00", "#2E7D32", "#81C784", "#FF7300", "#00FFFF", "#0061FE", "#1434A4")
  col = scale_color_manual(values = colAdd)
  fil = scale_fill_manual(values = colAdd)
  
  p1 <- ggplot(sample_sizes, aes(`Raw read number`/10^6, 
                                 fill = factor(Study, levels = c("Human gut", "Cat gut", "Dog gut", "Wild animal gut", "Forest soil", "Gator nest", "Acid mine drainage", "Salt marsh", "Coastal sediment", "Polar ocean")), 
                                 color = factor(Study, levels = c("Human gut", "Cat gut", "Dog gut", "Wild animal gut", "Forest soil", "Gator nest", "Acid mine drainage", "Salt marsh", "Coastal sediment", "Polar ocean")))) +
    geom_histogram(alpha = 0.7, bins = 30) +
    scale_x_continuous(trans="log10") +
    theme_bw() + 
    ylab("Sample number") + 
    xlab("Million raw reads") + 
    labs(color = "Dataset", fill = "Dataset") + 
    col + 
    fil + 
    theme(text = element_text(size = 16),
          legend.position = "bottom"
    ) + 
    guides(colour = guide_legend(nrow = 4))
  
  dir.create(file.path("figures/overview_figures/"), showWarnings = FALSE)
  ggsave(paste0("figures/overview_figures/read_counts.png"), width=18, height=18, units = 'cm', dpi=1000)
}
read_num_all()

PCoABray <- function(dataset, level, ncbi_only) {
  set.seed(1)
  
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
    tool_names <- tool_names[-c(9, 11)]
    colAdd <- colAdd[-c(9, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  
  data <- preprocess(dataset, level)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  for (i in 1:length(profiles)) {
    tmp <- profiles[[i]]
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    profiles[[i]] <- tmp
  }
  
  meta = read.csv(paste0("metadata/sample_info/read_info.tsv"), sep="\t")
  sample_ids = read.csv(paste0("real_data_outputs/inputs/", dataset, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)[1,] %>% 
    gsub(pattern="_taxonomic_profile", replacement = "") %>% gsub(pattern="^re", replacement="")
  sample_ids <- sort(sample_ids[-1])
  
  place_holders = rownames(profiles[[1]])
  for (i in 1:length(profiles)) {
    rownames(profiles[[i]]) <- mapvalues(rownames(profiles[[i]]), place_holders, sample_ids, warn_missing = F)
    if (ncol(profiles[[i]]) == 1) {
      colname_tmp <- colnames(profiles[[i]])
      rownames_tmp <- rownames(profiles[[i]])
      profiles[[i]] <- data.frame(profiles[[i]][rownames(profiles[[i]]) %in% meta$ID,])
      colnames(profiles[[i]]) <- colname_tmp
      rownames(profiles[[i]]) <- rownames_tmp
      subset_vec <- rowSums(profiles[[i]]) != 0
      profiles[[i]] <- data.frame(profiles[[i]][subset_vec,])
      colnames(profiles[[i]]) <- colname_tmp
      rownames(profiles[[i]]) <- rownames_tmp[subset_vec]
    } else {
      profiles[[i]] <- profiles[[i]][rownames(profiles[[i]]) %in% meta$ID,]
      profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
    }
    if (nrow(profiles[[i]]) > 0) {
      rownames(profiles[[i]]) <- paste0(rownames(profiles[[i]]), "_", i)
    }
  }
  
  dat_taxa_species <- bind_rows(profiles[[1]], profiles[[2]])
  for (i in 3:length(profiles)) {
    dat_taxa_species <- bind_rows(dat_taxa_species, profiles[[i]])
  }
  dat_taxa_species[is.na(dat_taxa_species)] <- 0
  
  meta <- meta[meta$ID %in% sample_ids,]
  meta <- meta[order(meta$ID),]
  rownames(meta) <- meta$ID
  
  stat_meta <- data.frame(matrix(nrow = 0, ncol = nrow(profiles[[1]])))
  for (i in 1:length(profiles)) {
    tmp <- meta
    tmp$ID <- paste0(tmp$ID, "_", i)
    tmp$method <- tool_names[i]
    tmp <- tmp[tmp$ID %in% rownames(profiles[[i]]),]
    stat_meta <- rbind(stat_meta, tmp)
  }
  
  meta_wo_na = stat_meta[, which(names(stat_meta) != "ID")]
  rownames(meta_wo_na) <- stat_meta$ID
  
  meta_wo_na <- meta_wo_na[rownames(meta_wo_na) %in% rownames(dat_taxa_species), ]
  meta_wo_na$method_type <- case_when(meta_wo_na$method=="Centrifuge" | meta_wo_na$method=="Kraken 2 / Bracken 2" ~ "kmer",
                                      grepl("MetaPhlAn", meta_wo_na$method) ~ "unique marker",
                                      grepl("mOTUs|Metaxa", meta_wo_na$method) ~ "conserved marker",
                                      TRUE ~ "assembly")
  
  distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
  
  adonis_out <- adonis2(as.formula(paste("distances ~ method_type")), data = meta_wo_na[!is.na(meta_wo_na[,"method_type"]),], permutations = 9999)

  pc <- capscale(distances~1, comm = distances, na.action = na.omit)
  cap = data.frame(pc$CA$u)
  cap$Method <- tool_names[gsub(".*_", "", rownames(cap)) %>% as.numeric()]

  cap$Sample <- gsub("_[0-9]+$", "", rownames(cap))
  s = summary(pc)
  
  # Swap study name for plotting
  study = case_when(dataset == "acid_mine" ~ "Acid Mine Drainage",
                     dataset == "animal_gut" ~ "Wild Animal Gut",
                     dataset == "cat_gut" ~ "Cat Gut",
                     dataset == "coastal_sediment" ~ "Coastal Sediment",
                     dataset == "dog_gut" ~ "Dog Gut",
                     dataset == "forest_soil" ~ "Forest Soil",
                     dataset == "gator_soil" ~ "Gator Nest",
                     dataset == "human" ~ "Human Gut",
                     dataset == "saltmarsh" ~ "Salt Marsh",
                     dataset == "tara_polar" ~ "Polar ocean",)
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  study <- paste0(study, ", ", level_name, " level")
  
  sample_indexes = sample(unique(cap$Sample), 3)
  cap$Sample <- ifelse(cap$Sample %in% sample_indexes, cap$Sample, NA)
  
  colAdd <- colAdd[tool_names %in% cap$Method]
  col = scale_color_manual(values = colAdd)
  
  # Plot ordination
  ggplot(cap, aes(MDS1, MDS2)) + 
    geom_point(aes(color=factor(Method, levels=tool_names)), size = 2, alpha = 0.7) + 
    theme_classic() +
    ggtitle(study) +
    labs( color = "Tool") + 
    labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
         y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = '')) +
    col + guides(color=guide_legend(title="Method")) +
    geom_line(subset(cap,!is.na(Sample)), mapping=aes(MDS1, MDS2, linetype=Sample), alpha=0.5, linewidth=0.5) + 
    guides(linetype = "none") +
    annotate("text", x=mean(c(rep(min(cap$MDS1), 2), max(cap$MDS1))), y=max(cap$MDS2), label= paste0("PERMANOVA: R-sq ", round(adonis_out$R2[1], 3), ", P-value ", sprintf("%.4f",adonis_out$`Pr(>F)`[1])))

  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/PCoABray/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_", ncbi_only,"_PCoA.png"), width=18, height=12, units = 'cm', dpi=1000)
}
run_all_PCoABray <- function() {
  for (dataset in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")]) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_NCBI_only_PCoA.png")) | !file.exists(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_all_taxa_PCoA.png"))) {
        PCoABray(dataset, level, TRUE)
        PCoABray(dataset, level, FALSE)
      }
    }
  }
}
run_all_PCoABray()

PCoAJaccard <- function(dataset, level, ncbi_only) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
    tool_names <- tool_names[-c(9, 11)]
    colAdd <- colAdd[-c(9, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  data <- preprocess(dataset, level)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  for (i in 1:length(profiles)) {
    tmp <- profiles[[i]]
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    profiles[[i]] <- tmp
  }
  
  dat_taxa_species <- bind_rows(profiles[[1]], profiles[[2]])
  for (i in 3:length(profiles)) {
    dat_taxa_species <- bind_rows(dat_taxa_species, profiles[[i]])
  }
  dat_taxa_species[is.na(dat_taxa_species)] <- 0
  dat_taxa_species <- dat_taxa_species[rowSums(dat_taxa_species) != 0,]
  distances <- vegdist(dat_taxa_species, method="jaccard") %>% as.matrix(labels=TRUE)
  
  pc <- capscale(distances~1, comm = distances, na.action = na.omit)
  cap = data.frame(pc$CA$u)
  cap$Method <- tool_names[gsub(".*_", "", rownames(cap)) %>% as.numeric()]
  
  cap$Sample <- gsub("_[0-9]+$", "", rownames(cap))
  s = summary(pc)
  
  # Swap study name for plotting
  study = case_when(dataset == "acid_mine" ~ "Acid Mine Drainage",
                    dataset == "animal_gut" ~ "Wild Animal Gut",
                    dataset == "cat_gut" ~ "Cat Gut",
                    dataset == "coastal_sediment" ~ "Coastal Sediment",
                    dataset == "dog_gut" ~ "Dog Gut",
                    dataset == "forest_soil" ~ "Forest Soil",
                    dataset == "gator_soil" ~ "Gator Nest",
                    dataset == "human" ~ "Human Gut",
                    dataset == "saltmarsh" ~ "Salt Marsh",
                    dataset == "tara_polar" ~ "Polar ocean",)
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  study <- paste0(study, ", ", level_name, " level")
  
  sample_indexes = sample(unique(cap$Sample), 1)
  cap$Sample <- ifelse(cap$Sample %in% sample_indexes, cap$Sample, NA)
  
  colAdd <- colAdd[tool_names %in% cap$Method]
  col = scale_color_manual(values = colAdd)
  
  # Plot ordination
  ggplot(cap, aes(MDS1, MDS2)) + 
    geom_point(aes(color=factor(Method, levels=tool_names)), size = 1, alpha = 0.7) + 
    theme_classic() +
    ggtitle(study) +
    labs( color = "Tool") + 
    labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
         y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = '')) +
    col + guides(color=guide_legend(title="Method")) +
    geom_line(subset(cap,!is.na(Sample)), mapping=aes(MDS1, MDS2))
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/PCoAJaccard/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/PCoAJaccard/",dataset,"/",dataset,"_",level_name, "_", ncbi_only,"_PCoA.png"), width=14, height=12, units = 'cm', dpi=1000)
}
run_all_PCoAJaccard <- function() {
  for (dataset in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")]) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/PCoAJaccard/",dataset,"/",dataset,"_",level_name, "_NCBI_only_PCoA.png")) | !file.exists(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_all_taxa_PCoA.png"))) {
        PCoAJaccard(dataset, level, TRUE)
        PCoAJaccard(dataset, level, FALSE)
      }
    }
  }
}
run_all_PCoAJaccard()

grid_iou <- function(dataset, level, ncbi_only, core_only = TRUE) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(9, 11)]
    colAdd <- colAdd[-c(9, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    melted_cormat = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(melted_cormat) <- c("Var1", "Var2")
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    for (dataset_tmp in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("human", "overall_dists")]) {
      if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names_tmp <- tool_names[-c(3, 11)]
      } else {
        tool_names_tmp <- tool_names
      }
      
      data <- preprocess(dataset_tmp, level)
      
      if (ncbi_only) {
        profiles <- remove_non_ncbi(data)
      } else {
        profiles <- remove_unknown(data)
      }
      
      profiles = lapply(profiles, renormalize)
      profiles = lapply(profiles, threshold_sample, threshold)
      profiles = lapply(profiles, renormalize)
      
      cormat = matrix(nrow = length(data), ncol = length(data))
      for (i in 1:nrow(cormat)) {
        for (j in 1:ncol(cormat)) {
          cormat[i,j] <- mean(calc_iou(list_taxa_by_sample(profiles[[i]], threshold), list_taxa_by_sample(profiles[[j]], threshold)), na.rm=T)
        }
      }
      
      colnames(cormat) <- tool_names_tmp
      rownames(cormat) <- tool_names_tmp
      melted_cormat_tmp <- melt(get_upper_tri(cormat), na.rm = T)
      melted_cormat_tmp$Var1 <- as.character(melted_cormat_tmp$Var1)
      melted_cormat_tmp$Var2 <- as.character(melted_cormat_tmp$Var2)
      melted_cormat <- full_join(melted_cormat, melted_cormat_tmp, by=c("Var1", "Var2"))
    }
    melted_cormat <- data.frame("Var1" = melted_cormat$Var1, 
                                "Var2" = melted_cormat$Var2, 
                                "value" = rowMeans(melted_cormat[,-c(1,2)], na.rm = T))

  } else {
    data <- preprocess(dataset, level)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    cormat = matrix(nrow = length(data), ncol = length(data))
    for (i in 1:nrow(cormat)) {
      for (j in 1:ncol(cormat)) {
        cormat[i,j] <- mean(calc_iou(list_taxa_by_sample(profiles[[i]], threshold), list_taxa_by_sample(profiles[[j]], threshold)), na.rm=T)
      }
    }
    
    colnames(cormat) <- tool_names
    rownames(cormat) <- tool_names
    melted_cormat <- melt(get_upper_tri(cormat), na.rm = T)
  }
  
  melted_cormat <- melted_cormat[order(melted_cormat$Var2),]
  melted_cormat <- melted_cormat[order(melted_cormat$Var1),]
  
  if (core_only) {
    melted_cormat <- melted_cormat[melted_cormat$Var1 %in% tool_core & melted_cormat$Var2 %in% tool_core,]
  }
  
  melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2,]
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  melted_cormat$value <- round(melted_cormat$value, 2)
  ggheatmap <- ggplot(melted_cormat, aes(factor(Var2, tool_names), factor(Var1, tool_names), fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max(melted_cormat$value)/2,
                         limit = c(0,max(melted_cormat$value)), space = "Lab", 
                         name="Intersection\n over union") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 24, hjust = 1),
          axis.text.y = element_text(size = 24),
          text = element_text(size=24))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(Var2, tool_names), factor(Var1, tool_names), label = value), color = "black", size = 9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.title=element_text(size=24),
      legend.text = element_text(size=24))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/iou/",dataset,"/"), showWarnings = FALSE)
  if (core_only) {
    ggsave(paste0("figures/iou/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=24, height=18, units = 'cm', dpi=1000, bg='#ffffff')
  } else {
    ggsave(paste0("figures/iou/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_all_tools.png"), width=48, height=36, units = 'cm', dpi=1000, bg='#ffffff')
  }
}
make_all_grid_iou <- function() {
  for (dataset in c(list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")], "all")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/iou/",dataset,"/",dataset,"_",level_name,"_", "NCBI_only", ".png"))) {
        grid_iou(dataset, level, TRUE)
      }
      if (!file.exists(paste0("figures/iou/",dataset,"/",dataset,"_",level_name,"_", "NCBI_only", "_all_tools.png"))) {
        grid_iou(dataset, level, TRUE, FALSE)
      }
      if (!file.exists(paste0("figures/iou/",dataset,"/",dataset,"_",level_name,"_", "all_taxa", ".png"))) {
        grid_iou(dataset, level, FALSE)
      }
    }
  }
}
make_all_grid_iou()

grid_iou_by_dataset <- function(level, ncbi_only) {
  out_df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(out_df) <- c("Var1", "Var2", "value", "dataset")
  
  for (dataset_tmp in list.files("real_data_outputs/inputs")[list.files("real_data_outputs/inputs")!="overall_dists"]) {
    if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
      tool_names_tmp <- tool_names[-c(3, 11)]
    } else {
      tool_names_tmp <- tool_names
    }
    set.seed(1)
    
    data <- preprocess(dataset_tmp, level)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    cormat = matrix(nrow = length(data), ncol = length(data))
    for (i in 1:nrow(cormat)) {
      for (j in 1:ncol(cormat)) {
        cormat[i,j] <- mean(calc_iou(list_taxa_by_sample(profiles[[i]], threshold), list_taxa_by_sample(profiles[[j]], threshold)), na.rm=T)
        if (is.na(cormat[i,j])) {
          cormat[i,j] <- 0
        }
      }
    }
    
    colnames(cormat) <- tool_names_tmp
    rownames(cormat) <- tool_names_tmp
    melted_cormat <- melt(get_upper_tri(cormat), na.rm = T)
    melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2,]
    melted_cormat$Var1 <- case_when(grepl("Centrifuge|Kraken", melted_cormat$Var1) ~ "K-mer",
                                    grepl("MetaPhlAn", melted_cormat$Var1) ~ "Unique marker",
                                    grepl("Metaxa|mOTUs", melted_cormat$Var1) ~ "Universal marker",
                                    TRUE ~ "Assembly")
    melted_cormat$Var2 <- case_when(grepl("Centrifuge|Kraken", melted_cormat$Var2) ~ "K-mer",
                                    grepl("MetaPhlAn", melted_cormat$Var2) ~ "Unique marker",
                                    grepl("Metaxa|mOTUs", melted_cormat$Var2) ~ "Universal marker",
                                    TRUE ~ "Assembly")
    melted_cormat$dataset <- dataset_tmp
    out_df <- rbind(out_df, melted_cormat)
  }
  
  out_df$dataset <- case_when(out_df$dataset == "acid_mine" ~ "Acid mine drainage",
                              out_df$dataset == "animal_gut" ~ "Wild animal gut",
                              out_df$dataset == "cat_gut" ~ "Cat gut",
                              out_df$dataset == "coastal_sediment" ~ "Coastal sediment",
                              out_df$dataset == "dog_gut" ~ "Dog gut",
                              out_df$dataset == "forest_soil" ~ "Forest soil",
                              out_df$dataset == "gator_soil" ~ "Gator nest",
                              out_df$dataset == "human" ~ "Human gut",
                              out_df$dataset == "saltmarsh" ~ "Salt marsh",
                              out_df$dataset == "tara_polar" ~ "Polar ocean",
                              out_df$dataset == "all" ~ "all")
  
  out_df <- aggregate(value~Var1 + Var2 + dataset, data=out_df, FUN=mean)
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  out_df$value <- round(out_df$value, 2)
  tool_type_names <- c("K-mer", "Unique marker", "Universal marker", "Assembly")
  ggheatmap <- ggplot(out_df, aes(factor(Var2, tool_type_names), factor(Var1, tool_type_names), fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), space = "Lab", 
                         name="Intersection\n over union") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 14, hjust = 1),
          axis.text.y = element_text(size = 14))+
    coord_fixed() + 
    facet_wrap(~factor(dataset, c("Human gut", "Cat gut", "Dog gut", "Wild animal gut", "Forest soil", "Gator nest", "Acid mine drainage", "Salt marsh", "Coastal sediment", "Polar ocean")), nrow=1)
  
  ggheatmap + 
    geom_text(aes(factor(Var2, tool_type_names), factor(Var1, tool_type_names), label = value), color = "black", size = 4.5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      text = element_text(size = 20),
      legend.title = element_text(size=14))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/iou_all_datasets/"), showWarnings = FALSE)
  ggsave(paste0("figures/iou_all_datasets/",level_name,"_", ncbi_only, ".png"), width=60, height=10, units = 'cm', dpi=1000, bg='#ffffff')
  
}
make_all_grid_iou_by_dataset <- function() {
  for (level in 1:7) {
    level_name <- case_when(level == 1 ~ "kingdom",
                            level == 2 ~ "phylum", 
                            level == 3 ~ "class",
                            level == 4 ~ "order",
                            level == 5 ~ "family",
                            level == 6 ~ "genus",
                            level == 7 ~ "species")
    if (!file.exists(paste0("figures/iou_all_datasets/",level_name,"_", "NCBI_only", ".png"))) {
      grid_iou_by_dataset(level, TRUE)
    }
    if (!file.exists(paste0("figures/iou_all_datasets/",level_name,"_", "all_taxa", ".png"))) {
      grid_iou_by_dataset(level, FALSE)
    }
  }
}
make_all_grid_iou_by_dataset()

grid_bc <- function(dataset, level, ncbi_only, core_only = TRUE) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(9, 11)]
    colAdd <- colAdd[-c(9, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    melted_cormat = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(melted_cormat) <- c("Var1", "Var2")
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    for (dataset_tmp in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("human", "overall_dists")]) {
      if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names_tmp <- tool_names[-c(3, 11)]
      } else {
        tool_names_tmp <- tool_names
      }
      
      data <- preprocess(dataset_tmp, level)
      
      if (ncbi_only) {
        profiles <- remove_non_ncbi(data)
      } else {
        profiles <- remove_unknown(data)
      }
      
      profiles = lapply(profiles, renormalize)
      profiles = lapply(profiles, threshold_sample, threshold)
      profiles = lapply(profiles, renormalize)
      
      cormat = matrix(nrow = length(data), ncol = length(data))
      for (i in 1:nrow(cormat)) {
        for (j in 1:ncol(cormat)) {
          cormat[i,j] <- mean(calc_bc(profiles[[i]], profiles[[j]]), na.rm=T)
        }
      }
      
      colnames(cormat) <- tool_names_tmp
      rownames(cormat) <- tool_names_tmp
      melted_cormat_tmp <- melt(get_upper_tri(cormat), na.rm = T)
      melted_cormat_tmp$Var1 <- as.character(melted_cormat_tmp$Var1)
      melted_cormat_tmp$Var2 <- as.character(melted_cormat_tmp$Var2)
      melted_cormat <- full_join(melted_cormat, melted_cormat_tmp, by=c("Var1", "Var2"))
    }
    melted_cormat <- data.frame("Var1" = melted_cormat$Var1, 
                                "Var2" = melted_cormat$Var2, 
                                "value" = rowMeans(melted_cormat[,-c(1,2)], na.rm = T))
    
  } else {
    data <- preprocess(dataset, level)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    cormat = matrix(nrow = length(data), ncol = length(data))
    for (i in 1:nrow(cormat)) {
      for (j in 1:ncol(cormat)) {
        cormat[i,j] <- mean(calc_bc(profiles[[i]], profiles[[j]]), na.rm=T)
      }
    }
    
    colnames(cormat) <- tool_names
    rownames(cormat) <- tool_names
    melted_cormat <- melt(get_upper_tri(cormat), na.rm = T)
  }
  
  melted_cormat <- melted_cormat[order(melted_cormat$Var2),]
  melted_cormat <- melted_cormat[order(melted_cormat$Var1),]
  
  if (core_only) {
    melted_cormat <- melted_cormat[melted_cormat$Var1 %in% tool_core & melted_cormat$Var2 %in% tool_core,]
  }
  melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2,]
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  melted_cormat$value <- 1 - round(melted_cormat$value, 2)
  ggheatmap <- ggplot(melted_cormat, aes(factor(Var2, tool_names), factor(Var1, tool_names), fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max(melted_cormat$value)/2,
                         limit = c(0,max(melted_cormat$value)), space = "Lab", 
                         name="1 - Bray Curtis\ndissimilarity") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 24, hjust = 1),
          axis.text.y = element_text(size = 24),
          text = element_text(size=20),
          legend.text = element_text(size=24))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(Var2, tool_names), factor(Var1, tool_names), label = value), color = "black", size = 9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.title=element_text(size=24),
      legend.text = element_text(size=24))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/bray_grid/",dataset,"/"), showWarnings = FALSE)
  
  if (core_only) {
    ggsave(paste0("figures/bray_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=24, height=18, units = 'cm', dpi=1000, bg='#ffffff')
  } else {
    ggsave(paste0("figures/bray_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_all_tools.png"), width=48, height=36, units = 'cm', dpi=1000, bg='#ffffff')
  }
  
}
make_all_grid_bc <- function() {
  for (dataset in c(list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")], "all")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/bray_grid/",dataset,"/",dataset,"_",level_name,"_NCBI_only.png"))) {
        grid_bc(dataset, level, TRUE)
      }
      if (!file.exists(paste0("figures/bray_grid/",dataset,"/",dataset,"_",level_name,"_NCBI_only_all_taxa.png"))) {
        grid_bc(dataset, level, TRUE, FALSE)
      }
      if (!file.exists(paste0("figures/bray_grid/",dataset,"/",dataset,"_",level_name,"_all_taxa.png"))) {
        grid_bc(dataset, level, FALSE)
      }
    }
  }
}
make_all_grid_bc()

grid_bc_by_dataset <- function(level, ncbi_only) {
  out_df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(out_df) <- c("Var1", "Var2", "value", "dataset")
  
  for (dataset_tmp in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")]) {
    if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
      tool_names_tmp <- tool_names[-c(3, 11)]
    } else {
      tool_names_tmp <- tool_names
    }
    set.seed(1)
    
    data <- preprocess(dataset_tmp, level)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    cormat = matrix(nrow = length(data), ncol = length(data))
    for (i in 1:nrow(cormat)) {
      for (j in 1:ncol(cormat)) {
        cormat[i,j] <- mean(calc_bc(profiles[[i]], profiles[[j]]), na.rm=T)
        if (is.na(cormat[i,j])) {
          cormat[i,j] <- 0
        }
      }
    }
    
    colnames(cormat) <- tool_names_tmp
    rownames(cormat) <- tool_names_tmp
    melted_cormat <- melt(get_upper_tri(cormat), na.rm = T)
    melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2,]
    melted_cormat$Var1 <- case_when(grepl("Centrifuge|Kraken", melted_cormat$Var1) ~ "K-mer",
                                    grepl("MetaPhlAn", melted_cormat$Var1) ~ "Unique marker",
                                    grepl("Metaxa|mOTUs", melted_cormat$Var1) ~ "Universal marker",
                                    TRUE ~ "Assembly")
    melted_cormat$Var2 <- case_when(grepl("Centrifuge|Kraken", melted_cormat$Var2) ~ "K-mer",
                                    grepl("MetaPhlAn", melted_cormat$Var2) ~ "Unique marker",
                                    grepl("Metaxa|mOTUs", melted_cormat$Var2) ~ "Universal marker",
                                    TRUE ~ "Assembly")
    melted_cormat$dataset <- dataset_tmp
    out_df <- rbind(out_df, melted_cormat)
  }
  
  out_df$dataset <- case_when(out_df$dataset == "acid_mine" ~ "Acid mine drainage",
                              out_df$dataset == "animal_gut" ~ "Wild animal gut",
                              out_df$dataset == "cat_gut" ~ "Cat gut",
                              out_df$dataset == "coastal_sediment" ~ "Coastal sediment",
                              out_df$dataset == "dog_gut" ~ "Dog gut",
                              out_df$dataset == "forest_soil" ~ "Forest soil",
                              out_df$dataset == "gator_soil" ~ "Gator nest",
                              out_df$dataset == "human" ~ "Human gut",
                              out_df$dataset == "saltmarsh" ~ "Salt marsh",
                              out_df$dataset == "tara_polar" ~ "Polar ocean",
                              out_df$dataset == "all" ~ "all")
  
  out_df <- aggregate(value~Var1 + Var2 + dataset, data=out_df, FUN=mean)
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  out_df$value <- 1 - round(out_df$value, 2)
  tool_type_names <- c("K-mer", "Unique marker", "Universal marker", "Assembly")
  ggheatmap <- ggplot(out_df, aes(factor(Var2, tool_type_names), factor(Var1, tool_type_names), fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), space = "Lab", 
                         name="1 - Bray Curtis\ndissimilarity") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 14, hjust = 1),
          axis.text.y = element_text(size = 14))+
    coord_fixed() + 
    facet_wrap(~factor(dataset, c("Human gut", "Cat gut", "Dog gut", "Wild animal gut", "Forest soil", "Gator nest", "Acid mine drainage", "Salt marsh", "Coastal sediment", "Polar ocean")), nrow=1)
  
  ggheatmap + 
    geom_text(aes(factor(Var2, tool_type_names), factor(Var1, tool_type_names), label = value), color = "black", size = 4.5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      text = element_text(size = 20),
      legend.title = element_text(size=14))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/bray_grid_by_dataset/"), showWarnings = FALSE)
  ggsave(paste0("figures/bray_grid_by_dataset/",level_name,"_", ncbi_only, ".png"), width=60, height=10, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_grid_bc_by_dataset <- function() {
  for (level in 1:7) {
    level_name <- case_when(level == 1 ~ "kingdom",
                            level == 2 ~ "phylum", 
                            level == 3 ~ "class",
                            level == 4 ~ "order",
                            level == 5 ~ "family",
                            level == 6 ~ "genus",
                            level == 7 ~ "species")
    if (!file.exists(paste0("figures/bray_grid_by_dataset/",level_name,"_", "NCBI_only", ".png"))) {
      grid_bc_by_dataset(level, TRUE)
    }
    if (!file.exists(paste0("figures/bray_grid_by_dataset/",level_name,"_", "all_taxa", ".png"))) {
      grid_bc_by_dataset(level, FALSE)
    }
  }
}
make_all_grid_bc_by_dataset()

# Plot an ECDF of unknown percentage by sample
unknownECDF <- function(dataset, level) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(9, 11)]
    colAdd <- colAdd[-c(9, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    df <- data.frame(matrix(ncol = 3, nrow = 0))
    for (dataset_tmp in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("human", "overall_dists")]) {
      if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names_tmp <- tool_names[-c(3, 11)]
      } else {
        tool_names_tmp <- tool_names
      }
      
      data <- preprocess(dataset_tmp, level, real = TRUE, by_dataset = TRUE)
      data = lapply(data, renormalize)
      data = lapply(data, threshold_sample, threshold)
      data = lapply(data, renormalize)
      
      df_tmp <- data.frame(matrix(ncol = ncol(data[[1]]) - 1, nrow = 0))
      for (i in 1:length(data)) {
        tmp <- data[[i]]
        df_tmp[nrow(df_tmp) + 1,] <- c(as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1]), rep(NA, ncol(df_tmp) - ncol(tmp) + 1))
      }
      df_tmp$Methods <- tool_names_tmp
      
      df_tmp <- melt(df_tmp, id.vars = c("Methods"))
      
      df <- rbind(df, df_tmp)
    }
  } else {
    data <- preprocess(dataset, level, real = TRUE, by_dataset = TRUE)
    data = lapply(data, renormalize)
    data = lapply(data, threshold_sample, threshold)
    data = lapply(data, renormalize)
    
    df <- data.frame(matrix(ncol = ncol(data[[1]]) - 1, nrow = 0))
    for (i in 1:length(data)) {
      tmp <- data[[i]]
      df[nrow(df) + 1,] <- c(as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1]), rep(NA, ncol(df) - ncol(tmp) + 1))
    }
    df$Methods <- tool_names
    
    df <- melt(df, id.vars = c("Methods"))
  }
  
  study = case_when(dataset == "acid_mine" ~ "acid mine drainage",
                    dataset == "animal_gut" ~ "wild animal gut",
                    dataset == "cat_gut" ~ "cat gut",
                    dataset == "coastal_sediment" ~ "coastal sediment",
                    dataset == "dog_gut" ~ "dog gut",
                    dataset == "forest_soil" ~ "forest soil",
                    dataset == "gator_soil" ~ "gator nest",
                    dataset == "human" ~ "human gut",
                    dataset == "saltmarsh" ~ "salt marsh",
                    dataset == "tara_polar" ~ "Polar ocean",
                    dataset == "all" ~ "all")
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  colAdd <- colAdd[tool_names %in% df$Methods]
  col = scale_color_manual(values = colAdd)
  
  df$Methods <- factor(df$Methods, c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes"))
  
  ggplot(df, aes(x=value, color=Methods)) + 
    stat_ecdf(geom = "step", size = 1) +
    theme_bw() + ggtitle(paste0("Empirical cumulative density over ", study, " samples")) +
    ylab("Cumulative density") + xlab(paste0("Percent unknown at the ", level_name, " level")) +
    col +
    guides(color=guide_legend(title="Method"))
  dir.create(file.path("figures/unknownECDF/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/unknownECDF/",dataset,"/",dataset,"_",level_name,"_ECDFUnknowns.png"), width=18, height=10, units = 'cm', dpi=1000)
}
make_all_unknownECDF_plots <- function() {
  for (dataset in c(list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")], "all")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/unknownECDF/",dataset,"/",dataset,"_",level_name,"_ECDFUnknowns.png"))) {
        unknownECDF(dataset, level)
      }
    }
  }
}
make_all_unknownECDF_plots()

grid_bc_mantel <- function(dataset, level, ncbi_only, core_only = TRUE) {
  if (level == 7) {
    tool_names <- c(tool_names, "Simka")
  }
  
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(9, 11)]
    colAdd <- colAdd[-c(9, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    melted_cormat = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(melted_cormat) <- c("Var1", "Var2")
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    for (dataset_tmp in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("human", "overall_dists")]) {
      if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names_tmp <- tool_names[-c(3, 11)]
      } else {
        tool_names_tmp <- tool_names
      }
      
      data <- preprocess(dataset_tmp, level, real = TRUE, by_dataset = TRUE)
      
      if (ncbi_only) {
        profiles <- remove_non_ncbi(data)
      } else {
        profiles <- remove_unknown(data)
      }
      
      profiles = lapply(profiles, renormalize)
      profiles = lapply(profiles, threshold_sample, threshold)
      profiles = lapply(profiles, renormalize)
      
      if (level == 7) {
        cormat = matrix(nrow = length(data) + 1, ncol = length(data) + 1)
      } else {
        cormat = matrix(nrow = length(data), ncol = length(data))
      }
      for (i in 1:length(data)) {
        for (j in 1:length(data)) {
          if (i <= j) {
            cormat[i,j] <- run_bc_mantel(profiles[[i]], profiles[[j]])[1]
          }
        }
      }
      
      if (level == 7) {
        simka = read.csv(paste0("real_data_outputs/inputs/", dataset_tmp, "/Simka/mat_abundance_braycurtis.csv"), sep = ";")
        # Simka columns are ordered by accession number and the profiles are too.  Centrifuge will always have all columns
        colnames(simka) <- colnames(profiles[[1]])
        for (i in 1:(ncol(cormat) - 1)) {
          cormat[i,ncol(cormat)] <- run_bc_mantel(profiles[[i]], simka)[1]
        }
        cormat[nrow(cormat),ncol(cormat)] <- run_bc_mantel(simka, simka)[1]
      }
      
      colnames(cormat) <- tool_names_tmp
      rownames(cormat) <- tool_names_tmp
      melted_cormat_tmp <- melt(get_upper_tri(cormat), na.rm = T)
      melted_cormat_tmp$Var1 <- as.character(melted_cormat_tmp$Var1)
      melted_cormat_tmp$Var2 <- as.character(melted_cormat_tmp$Var2)
      melted_cormat <- full_join(melted_cormat, melted_cormat_tmp, by=c("Var1", "Var2"))
    }
    sds = apply(melted_cormat[,-c(1,2)], 1, sd, na.rm = TRUE)
    ns = rowSums(!is.na(melted_cormat[,-c(1,2)]))
    ses = sds/ns
    means = rowMeans(melted_cormat[,-c(1,2)], na.rm = T)
    tstats = abs(means)/ses
    pvals <- pt(tstats, ns - 1, lower.tail = F)
    melted_cormat <- data.frame("Var1" = melted_cormat$Var1, 
                                "Var2" = melted_cormat$Var2, 
                                "value" = means,
                                "stars" = case_when(pvals < 0.001 ~ "***",
                                                    pvals < 0.01 ~ "**",
                                                    pvals < 0.05 ~ "*",
                                                    TRUE ~ ""))
    
  } else {
    data <- preprocess(dataset, level, real = TRUE, by_dataset = TRUE)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    if (level == 7) {
      cormat = matrix(nrow = length(data) + 1, ncol = length(data) + 1)
      signif_mat = matrix(nrow = length(data) + 1, ncol = length(data) + 1)
    } else {
      cormat = matrix(nrow = length(data), ncol = length(data))
      signif_mat = matrix(nrow = length(data), ncol = length(data))
    }
    for (i in 1:length(data)) {
      for (j in 1:length(data)) {
        if (i <= j) {
          mantel_out <- run_bc_mantel(profiles[[i]], profiles[[j]])
          cormat[i,j] <- mantel_out[1]
          signif_mat[i,j] <- mantel_out[2]
        }
      }
    }
    
    if (level == 7) {
      simka = read.csv(paste0("real_data_outputs/inputs/", dataset, "/Simka/mat_abundance_braycurtis.csv"), sep = ";")
      # Simka columns are ordered by accession number and the profiles are too.  Centrifuge will always have all columns
      colnames(simka) <- colnames(profiles[[1]])
      for (i in 1:(ncol(cormat) - 1)) {
        mantel_out <- run_bc_mantel(simka, profiles[[i]])
        cormat[i,ncol(cormat)] <- mantel_out[1]
        signif_mat[i,ncol(cormat)] <- mantel_out[2]
      }
      
      mantel_out <- run_bc_mantel(simka, simka)
      cormat[nrow(cormat),ncol(cormat)] <- mantel_out[1]
      signif_mat[nrow(cormat),ncol(cormat)] <- mantel_out[2]
    }
    
    colnames(cormat) <- tool_names
    rownames(cormat) <- tool_names
    colnames(signif_mat) <- tool_names
    rownames(signif_mat) <- tool_names
    melted_cormat <- melt(get_upper_tri(cormat), na.rm = T)
    melted_signif <- melt(get_upper_tri(signif_mat), na.rm = T)
    melted_signif$value <- case_when(melted_signif$value < 0.001 ~ "***",
                                     melted_signif$value < 0.01 ~ "**",
                                     melted_signif$value < 0.05 ~ "*",
                                     TRUE ~ "")
    
    melted_cormat$stars <- melted_signif$value
  }
  
  melted_cormat <- melted_cormat[order(melted_cormat$Var2),]
  melted_cormat <- melted_cormat[order(melted_cormat$Var1),]
  
  if (core_only) {
    melted_cormat <- melted_cormat[melted_cormat$Var1 %in% c(tool_core, "Simka") & melted_cormat$Var2 %in% c(tool_core, "Simka"),]
  }
  melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2,]
  

  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (sum(melted_cormat$value < 0, na.rm = T) > 0) {
    scale_fill_special <- scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                                               limit = c(min(melted_cormat$value) * 1.1 ,max(melted_cormat$value) * 1.1), space = "Lab", 
                                               name="Mantel r")
  } else {
    scale_fill_special <- scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = mean(c(min(melted_cormat$value),max(melted_cormat$value))),
                                               limit = c(min(melted_cormat$value) * 0.9,max(melted_cormat$value) * 1.1), space = "Lab", 
                                               name="Mantel r")
  }
  
  melted_cormat$value <- round(melted_cormat$value, 2)
  
  melted_cormat$Var2 <- as.character(melted_cormat$Var2)
  # melted_cormat$Var2[melted_cormat$Var2 == "Simka"] <- glue("*Simka*")
  
  ggheatmap <- ggplot(melted_cormat, aes(factor(Var2, c(tool_names, glue("*Simka*"))), factor(Var1, c(tool_names, glue("*Simka*"))), fill = value))+
    geom_tile(color = "white")+
    scale_fill_special +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_markdown(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(Var2, c(tool_names, glue("*Simka*"))), factor(Var1, c(tool_names, glue("*Simka*"))), label = paste0(value, "\n", stars)), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.57, 0.75),
      legend.direction = "horizontal",
      legend.title=element_text(size=14)) +
    guides(fill = guide_colorbar(barwidth = 5.5, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/beta_mantel/",dataset,"/"), showWarnings = FALSE)
  if (core_only) {
    ggsave(paste0("figures/beta_mantel/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=13, height=11, units = 'cm', dpi=1000, bg='#ffffff')
  } else {
    ggsave(paste0("figures/beta_mantel/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_all_tools.png"), width=19.5, height=17.5, units = 'cm', dpi=1000, bg='#ffffff')
  }
}
make_all_grid_bc_mantel <- function() {
  for (dataset in c(list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")], "all")) {
    for (level in c(2, 6, 7)) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      dir.create(file.path(paste0("figures/beta_mantel/",dataset,"/")), showWarnings = FALSE)
      if (!file.exists(paste0("figures/beta_mantel/",dataset,"/",dataset,"_",level_name,"_all_taxa_all_tools.png"))) {
        grid_bc_mantel(dataset, level, FALSE, FALSE)
      }
      if (!file.exists(paste0("figures/beta_mantel/",dataset,"/",dataset,"_",level_name,"_all_taxa.png"))) {
        grid_bc_mantel(dataset, level, FALSE)
      }
    }
  }
}
make_all_grid_bc_mantel()

grid_alpha_mantel <- function(dataset, level, ncbi_only, core_only = TRUE) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(9, 11)]
    colAdd <- colAdd[-c(9, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    melted_cormat = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(melted_cormat) <- c("Var1", "Var2")
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    for (dataset_tmp in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("human", "overall_dists")]) {
      if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names_tmp <- tool_names[-c(3, 11)]
      } else {
        tool_names_tmp <- tool_names
      }
      
      data <- preprocess(dataset_tmp, level, real = TRUE, by_dataset = TRUE)
      
      if (ncbi_only) {
        profiles <- remove_non_ncbi(data)
      } else {
        profiles <- remove_unknown(data)
      }
      
      profiles = lapply(profiles, renormalize)
      profiles = lapply(profiles, threshold_sample, threshold)
      profiles = lapply(profiles, renormalize)
      
      cormat = matrix(nrow = length(data), ncol = length(data))
      
      for (i in 1:length(data)) {
        for (j in 1:length(data)) {
          cormat[i,j] <- run_alpha_cor(profiles[[i]], profiles[[j]])[1]
        }
      }
      
      colnames(cormat) <- tool_names_tmp
      rownames(cormat) <- tool_names_tmp
      melted_cormat_tmp <- melt(get_upper_tri(cormat), na.rm = T)
      melted_cormat_tmp$Var1 <- as.character(melted_cormat_tmp$Var1)
      melted_cormat_tmp$Var2 <- as.character(melted_cormat_tmp$Var2)
      melted_cormat <- full_join(melted_cormat, melted_cormat_tmp, by=c("Var1", "Var2"))
    }
    
    sds = apply(melted_cormat[,-c(1,2)], 1, sd, na.rm = TRUE)
    ns = rowSums(!is.na(melted_cormat[,-c(1,2)]))
    ses = sds/ns
    means = rowMeans(melted_cormat[,-c(1,2)], na.rm = T)
    tstats = abs(means)/ses
    pvals <- pt(tstats, ns - 1, lower.tail = F)
    melted_cormat <- data.frame("Var1" = melted_cormat$Var1, 
                                "Var2" = melted_cormat$Var2, 
                                "value" = means,
                                "stars" = case_when(pvals < 0.001 ~ "***",
                                                    pvals < 0.01 ~ "**",
                                                    pvals < 0.05 ~ "*",
                                                    TRUE ~ ""))
    
  } else {
    data <- preprocess(dataset, level, real = TRUE, by_dataset = TRUE)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    cormat = matrix(nrow = length(data), ncol = length(data))
    signif_mat = matrix(nrow = length(data), ncol = length(data))
    
    for (i in 1:length(data)) {
      for (j in 1:length(data)) {
        alpha_out <- run_alpha_cor(profiles[[i]], profiles[[j]])
        cormat[i,j] <- alpha_out[1]
        signif_mat[i,j] = alpha_out[2]
      }
    }
    
    colnames(cormat) <- tool_names
    rownames(cormat) <- tool_names
    melted_cormat <- melt(get_upper_tri(cormat), na.rm = T)

    colnames(signif_mat) <- tool_names
    rownames(signif_mat) <- tool_names
    melted_signif <- melt(get_upper_tri(signif_mat), na.rm = T)
    melted_signif$value <- case_when(melted_signif$value < 0.001 ~ "***",
                                     melted_signif$value < 0.01 ~ "**",
                                     melted_signif$value < 0.05 ~ "*",
                                     TRUE ~ "")
    
    melted_cormat$stars <- melted_signif$value
  }
  
  melted_cormat <- melted_cormat[order(melted_cormat$Var2),]
  melted_cormat <- melted_cormat[order(melted_cormat$Var1),]
  
  if (core_only) {
    melted_cormat <- melted_cormat[melted_cormat$Var1 %in% c(tool_core, "Simka") & melted_cormat$Var2 %in% c(tool_core, "Simka"),]
  }
  melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2,]
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (sum(melted_cormat$value < 0, na.rm = T) > 0) {
    scale_fill_special <- scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                                               limit = c(min(melted_cormat$value) * 1.1 ,max(melted_cormat$value) * 1.1), space = "Lab", 
                                               name="Spearman r")
  } else {
    scale_fill_special <- scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = mean(c(min(melted_cormat$value),max(melted_cormat$value))),
                                               limit = c(min(melted_cormat$value) * 0.9,max(melted_cormat$value) * 1.1), space = "Lab", 
                                               name="Spearman r")
  }
  
  
  melted_cormat$value <- round(melted_cormat$value, 2)
  melted_cormat$Var2 <- as.character(melted_cormat$Var2)
  melted_cormat$Var2[melted_cormat$Var2 == "Simka"] <- glue("*Simka*")
  
  ggheatmap <- ggplot(melted_cormat, aes(factor(Var2, c(tool_names, glue("*Simka*"))), factor(Var1, c(tool_names, glue("*Simka*"))), fill = value))+
    geom_tile(color = "white")+
    scale_fill_special +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_markdown(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1),
          axis.text.y = element_text(size = 12))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(Var2, c(tool_names, glue("*Simka*"))), factor(Var1, c(tool_names, glue("*Simka*"))), label = paste0(value, "\n", stars)), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.75),
      legend.direction = "horizontal",
      legend.title=element_text(size=14)) +
    guides(fill = guide_colorbar(barwidth = 5.5, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/alpha_spearman/",dataset,"/"), showWarnings = FALSE)
  
  if (core_only) {
    ggsave(paste0("figures/alpha_spearman/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=13, height=11, units = 'cm', dpi=1000, bg='#ffffff')
  } else {
    ggsave(paste0("figures/alpha_spearman/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_all_tools.png"), width=19.5, height=17.5, units = 'cm', dpi=1000, bg='#ffffff')
  }
}
make_all_grid_alpha_spearman <- function() {
  for (dataset in c(list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")], "all")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/alpha_spearman/",dataset,"/",dataset,"_",level_name,"_NCBI_only.png"))) {
        grid_alpha_mantel(dataset, level, TRUE)
      }
      if (!file.exists(paste0("figures/alpha_spearman/",dataset,"/",dataset,"_",level_name,"_NCBI_only_all_tools.png"))) {
        grid_alpha_mantel(dataset, level, TRUE, FALSE)
      }
      if (!file.exists(paste0("figures/alpha_spearman/",dataset,"/",dataset,"_",level_name,"_all_taxa.png"))) {
        grid_alpha_mantel(dataset, level, FALSE)
      }
    }
  }
}
make_all_grid_alpha_spearman()

PERMANOVA_example <- function(dataset, level, ncbi_only) {
  meta = read.csv(paste0("metadata/", dataset, "/merged_metadata.tsv"), sep="\t")

  data <- preprocess(dataset, level, real = TRUE, by_dataset = TRUE)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  for (i in 1:length(profiles)) {
    tmp <- profiles[[i]]
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    profiles[[i]] <- tmp
  }
  
  sample_ids = read.csv(paste0("real_data_outputs/inputs/", dataset, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)[1,] %>% 
    gsub(pattern="_taxonomic_profile", replacement = "") %>% gsub(pattern="^re", replacement="")
  
  sample_ids <- sort(sample_ids[-1])
  
  meta <- meta[meta$sample_id %in% sample_ids,]
  meta <- meta[order(meta$sample_id),]
  rownames(meta) <- meta$sample_id
  
  place_holders = rownames(profiles[[1]])
  for (i in 1:length(profiles)) {
    rownames(profiles[[i]]) <- mapvalues(rownames(profiles[[i]]), place_holders, sample_ids, warn_missing = F)
    if (ncol(profiles[[i]]) == 1) {
      colname_tmp <- colnames(profiles[[i]])
      rownames_tmp <- rownames(profiles[[i]])
      profiles[[i]] <- data.frame(profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,])
      colnames(profiles[[i]]) <- colname_tmp
      rownames(profiles[[i]]) <- rownames_tmp
      subset_vec <- rowSums(profiles[[i]]) != 0
      profiles[[i]] <- data.frame(profiles[[i]][subset_vec,])
      colnames(profiles[[i]]) <- colname_tmp
      rownames(profiles[[i]]) <- rownames_tmp[subset_vec]
    } else {
      profiles[[i]] <- profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,]
      profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
    }
    if (nrow(profiles[[i]]) > 0) {
      rownames(profiles[[i]]) <- paste0(rownames(profiles[[i]]), "_", i)
    }
  }
  
  dat_taxa_species <- bind_rows(profiles[[1]], profiles[[2]])
  for (i in 3:length(profiles)) {
    dat_taxa_species <- bind_rows(dat_taxa_species, profiles[[i]])
  }
  dat_taxa_species[is.na(dat_taxa_species)] <- 0
  
  # Calculate Bray dissimilarity
  bray <- vegdist(dat_taxa_species, method="bray")

  stat_meta <- data.frame(matrix(nrow = 0, ncol = nrow(profiles[[1]])))
  for (i in 1:length(profiles)) {
    tmp <- meta
    tmp$sample_id <- paste0(tmp$sample_id, "_", i)
    tmp$method <- tool_names[i]
    tmp <- tmp[tmp$sample_id %in% rownames(profiles[[i]]),]
    stat_meta <- rbind(stat_meta, tmp)
  }
  
  meta_wo_na = stat_meta[, which(names(stat_meta) != "sample_id")]
  rownames(meta_wo_na) <- stat_meta$sample_id
  
  meta_wo_na <- meta_wo_na[rownames(meta_wo_na) %in% rownames(dat_taxa_species), ]
  meta_wo_na$method_type <- case_when(meta_wo_na$method=="Centrifuge" | meta_wo_na$method=="Kraken 2 / Bracken 2" ~ "kmer",
                                      grepl("MetaPhlAn|mOTUs|Metaxa", meta_wo_na$method) ~ "marker",
                                      TRUE ~ "assembly")
  
  if (dataset != "human") {
    # All data for methods-based
    adonis_res_rsq = vector()
    adonis_res_pval = vector()
    for (col in c("method", "method_type")) {
      adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na)
      adonis_res_rsq[col] = adonis.univ[1,"R2"]
      adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
    }
    
    # Only single copy of data for everything else
    colnames_tmp <- colnames(meta_wo_na)
    meta_wo_na <- data.frame(meta_wo_na[meta_wo_na$method=="MetaPhlAn 4",-which(colnames(meta_wo_na) %in% c("method", "method_type"))])
    colnames(meta_wo_na) <- colnames_tmp[!(colnames_tmp %in% c("method", "method_type"))]
    new_dat_taxa_species <- profiles[[which(tool_names == "MetaPhlAn 4")]]
    for (col in names(meta_wo_na)){
      bray <- vegdist(new_dat_taxa_species[!is.na(meta_wo_na[,col]),])
      if (ncol(meta_wo_na) == 1) {
        colnames_tmp <- colnames(meta_wo_na)
        in_df <- data.frame(meta_wo_na[!is.na(meta_wo_na[,col]),])
        colnames(in_df) <- colnames_tmp
        adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = in_df, permutations = 9999)
      } else {
        adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999)
      }
      
      adonis_res_rsq[col] = adonis.univ[1,"R2"]
      adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
    }
  } else {
    # All data for methods-based
    adonis_res_rsq = vector()
    adonis_res_pval = vector()
    for (col in c("method", "method_type")) {
      adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999, strata = meta_wo_na[!is.na(meta_wo_na[,col]),]$participant)
      adonis_res_rsq[col] = adonis.univ[1,"R2"]
      adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
    }
    
    # Only single copy of data for everything else
    colnames_tmp <- colnames(meta_wo_na)
    meta_wo_na <- data.frame(meta_wo_na[meta_wo_na$method=="MetaPhlAn 4",-which(colnames(meta_wo_na) %in% c("method", "method_type"))])
    colnames(meta_wo_na) <- colnames_tmp[!(colnames_tmp %in% c("method", "method_type"))]
    new_dat_taxa_species <- profiles[[which(tool_names == "MetaPhlAn 4")]]
    for (col in names(meta_wo_na)[names(meta_wo_na) != "participant"]){
      bray <- vegdist(new_dat_taxa_species[!is.na(meta_wo_na[,col]),])
      
      if (nrow(matrix(bray)) != 0) {
        adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999, strata = meta_wo_na[!is.na(meta_wo_na[,col]),]$participant)
        adonis_res_rsq[col] = adonis.univ[1,"R2"]
        adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
      } else {
        adonis_res_rsq[col] = 0
        adonis_res_pval[col] = 1
      }
      
    }
    col = "participant"
    bray <- vegdist(new_dat_taxa_species[!is.na(meta_wo_na[,col]),])
    
    if (nrow(matrix(bray)) != 0) {
      adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999)
      adonis_res_rsq[col] = adonis.univ[1,"R2"]
      adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
    } else {
      adonis_res_rsq[col] = 0
      adonis_res_pval[col] = 1
    }
  }
  
  # Format adonis output for plotting
  univar_res_tax = data.frame(adonis_res_pval, adonis_res_rsq)
  colnames(univar_res_tax) = c("P-Value", "R2")
  univar_res_tax$`P-Value` = as.numeric(univar_res_tax$`P-Value`)
  univar_res_tax$R2 = as.numeric(univar_res_tax$R2)
  proper_names <- gsub("_", " ", names(adonis_res_rsq))
  proper_names <- sapply(proper_names, simpleCap)
  univar_res_tax$variables <- proper_names
  rownames(univar_res_tax) <- univar_res_tax$variables
  univar_res_tax$stars = cut(univar_res_tax$`P-Value`, c(0, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))
  
  univar_res_tax$`P-Value` = round(univar_res_tax$`P-Value`, 3)
  univar_res_tax$R2 = univar_res_tax$R2 *100
  dodge = position_dodge(width = 0.8)
  
  study_name <- case_when(dataset == "acid_mine" ~ "acid mine runoff", 
                          dataset == "animal_gut" ~ "wild animal gut", 
                          dataset == "forest_soil" ~ "forest soil", 
                          dataset == "gator_soil" ~ "gator nest",
                          dataset == "saltmarsh" ~ "salt marsh", 
                          dataset == "tara_polar" ~ "polar ocean",
                          dataset == "human" ~ "human gut")
  
  # Plot
  ggplot(data = univar_res_tax, aes(reorder(variables, R2), y = R2, label = stars)) + 
    geom_bar(stat = "identity", position = dodge, fill="brown") + 
    geom_text(position = dodge, vjust = 0.8, hjust = -0.1, size = 4) + 
    theme_classic(base_size = 10) + ylab("Univarate R-squared") + 
    coord_flip() + ylim(0, max(univar_res_tax$R2) * 1.2) + xlab("") + 
    labs(fill = "Tool") + 
    ggtitle(paste0("PERMANOVA on ", study_name, " samples")) +
    guides(fill = guide_legend(reverse=T)) + 
    theme(text = element_text(size = 18))
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/PERMANOVA/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/PERMANOVA/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_PERMANOVA.png"), width=18, height=15, units = 'cm', dpi=1000)
}
make_all_PERMANOVA_example <- function() {
  for (dataset in c("acid_mine", "animal_gut", "forest_soil", "gator_soil", "saltmarsh", "tara_polar", "human")) {
    for (level in 2:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/PERMANOVA/", dataset,"/", dataset,"_",level_name,"_", "NCBI_only", "_PERMANOVA.png"))) {
        PERMANOVA_example(dataset, level, TRUE)
      }
      if (!file.exists(paste0("figures/PERMANOVA/", dataset,"/", dataset,"_",level_name,"_", "all_taxa", "_PERMANOVA.png"))) {
        PERMANOVA_example(dataset, level, FALSE)
      }
    }
  }
}
make_all_PERMANOVA_example()

PERMANOVA_comparison <- function(dataset, level, ncbi_only = FALSE, core_only = TRUE) {
  df = data.frame(matrix(nrow = 0, ncol = 6))
  colnames(df) <- c("method", "variable", "R2", "p_val", "dataset")
  if (dataset == "all") {
    for (dataset_tmp in c("acid_mine", "animal_gut", "forest_soil", "gator_soil", "saltmarsh", "tara_polar", "human")) {
      meta = read.csv(paste0("metadata/", dataset_tmp, "/merged_metadata.tsv"), sep="\t")
      
      data <- preprocess(dataset_tmp, level, real = TRUE, by_dataset = TRUE)
      
      if (ncbi_only) {
        profiles <- remove_non_ncbi(data)
      } else {
        profiles <- remove_unknown(data)
      }
      
      profiles = lapply(profiles, renormalize)
      profiles = lapply(profiles, threshold_sample, threshold)
      profiles = lapply(profiles, renormalize)
      
      for (i in 1:length(profiles)) {
        tmp <- profiles[[i]]
        taxa_names <- tmp$TaxIDs
        tmp <- data.frame(t(tmp[,-1]))
        colnames(tmp) <- as.character(taxa_names)
        profiles[[i]] <- tmp
      }
      
      sample_ids = read.csv(paste0("real_data_outputs/inputs/", dataset_tmp, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)[1,] %>% 
        gsub(pattern="_taxonomic_profile", replacement = "") %>% gsub(pattern="^re", replacement="")
      
      sample_ids <- sort(sample_ids[-1])
      
      meta <- meta[meta$sample_id %in% sample_ids,]
      meta <- meta[order(meta$sample_id),]
      rownames(meta) <- meta$sample_id
      
      place_holders = rownames(profiles[[1]])
      for (i in 1:length(profiles)) {
        rownames(profiles[[i]]) <- mapvalues(rownames(profiles[[i]]), place_holders, sample_ids, warn_missing = F)
        if (ncol(profiles[[i]]) == 1) {
          colname_tmp <- colnames(profiles[[i]])
          rownames_tmp <- rownames(profiles[[i]])
          profiles[[i]] <- data.frame(profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,])
          colnames(profiles[[i]]) <- colname_tmp
          rownames(profiles[[i]]) <- rownames_tmp
          subset_vec <- rowSums(profiles[[i]]) != 0
          profiles[[i]] <- data.frame(profiles[[i]][subset_vec,])
          colnames(profiles[[i]]) <- colname_tmp
          rownames(profiles[[i]]) <- rownames_tmp[subset_vec]
        } else {
          profiles[[i]] <- profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,]
          profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
        }
      }
      
      for (i in 1:length(profiles)) {
        meta_wo_na <- meta
        meta_wo_na$sample_id <- NULL
        colnames_tmp <- colnames(meta_wo_na)
        meta_wo_na <- data.frame(meta_wo_na[rownames(meta_wo_na) %in% rownames(profiles[[i]]), ])
        colnames(meta_wo_na) <- colnames_tmp
        
        # Only single copy of data for everything else
        adonis_res_rsq = vector()
        adonis_res_pval = vector()
        if (dataset != "human") {
          for (col in names(meta_wo_na)){
            bray <- vegdist(profiles[[i]][!is.na(meta_wo_na[,col]),])
            if (nrow(matrix(bray)) != 0) {
              data_in <- data.frame(meta_wo_na[!is.na(meta_wo_na[,col]),])
              colnames(data_in) <- colnames(meta_wo_na)
              if (length(unique(data_in[,col])) > 1) {
                adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = data_in, permutations = 9999)
                adonis_res_rsq[col] = adonis.univ[1,"R2"]
                adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
              } else {
                adonis_res_rsq[col] = 0
                adonis_res_pval[col] = 1
              }
            } else {
              adonis_res_rsq[col] = 0
              adonis_res_pval[col] = 1
            }
          }
        } else {
          for (col in names(meta_wo_na)[names(meta_wo_na) != "participant"]){
            bray <- vegdist(profiles[[i]][!is.na(meta_wo_na[,col]),])
            
            if (nrow(matrix(bray)) != 0) {
              adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999, strata = meta_wo_na[!is.na(meta_wo_na[,col]),]$participant)
              adonis_res_rsq[col] = adonis.univ[1,"R2"]
              adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
            } else {
              adonis_res_rsq[col] = 0
              adonis_res_pval[col] = 1
            }
            
          }
          col = "participant"
          bray <- vegdist(profiles[[i]][!is.na(meta_wo_na[,col]),])
          
          if (nrow(matrix(bray)) != 0) {
            adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999)
            adonis_res_rsq[col] = adonis.univ[1,"R2"]
            adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
          } else {
            adonis_res_rsq[col] = 0
            adonis_res_pval[col] = 1
          }
        }
        
        df_addition = data.frame("method"=i, "variable"=names(adonis_res_pval), "R2"=adonis_res_rsq, "p_val"=adonis_res_pval, "dataset" = dataset_tmp)
        df <- rbind(df, df_addition)
      }
    }
  } else {
    meta = read.csv(paste0("metadata/", dataset, "/merged_metadata.tsv"), sep="\t")
    
    data <- preprocess(dataset, level, real = TRUE, by_dataset = TRUE)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    for (i in 1:length(profiles)) {
      tmp <- profiles[[i]]
      taxa_names <- tmp$TaxIDs
      tmp <- data.frame(t(tmp[,-1]))
      colnames(tmp) <- as.character(taxa_names)
      profiles[[i]] <- tmp
    }
    
    sample_ids = read.csv(paste0("real_data_outputs/inputs/", dataset, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)[1,] %>% 
      gsub(pattern="_taxonomic_profile", replacement = "") %>% gsub(pattern="^re", replacement="")
    
    sample_ids <- sort(sample_ids[-1])
    
    meta <- meta[meta$sample_id %in% sample_ids,]
    meta <- meta[order(meta$sample_id),]
    rownames(meta) <- meta$sample_id
    
    place_holders = rownames(profiles[[1]])
    for (i in 1:length(profiles)) {
      rownames(profiles[[i]]) <- mapvalues(rownames(profiles[[i]]), place_holders, sample_ids, warn_missing = F)
      if (ncol(profiles[[i]]) == 1) {
        colname_tmp <- colnames(profiles[[i]])
        rownames_tmp <- rownames(profiles[[i]])
        profiles[[i]] <- data.frame(profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,])
        colnames(profiles[[i]]) <- colname_tmp
        rownames(profiles[[i]]) <- rownames_tmp
        subset_vec <- rowSums(profiles[[i]]) != 0
        profiles[[i]] <- data.frame(profiles[[i]][subset_vec,])
        colnames(profiles[[i]]) <- colname_tmp
        rownames(profiles[[i]]) <- rownames_tmp[subset_vec]
      } else {
        profiles[[i]] <- profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,]
        profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
      }
    }
    
    for (i in 1:length(profiles)) {
      meta_wo_na <- meta
      meta_wo_na$sample_id <- NULL
      colnames_tmp <- colnames(meta_wo_na)
      meta_wo_na <- data.frame(meta_wo_na[rownames(meta_wo_na) %in% rownames(profiles[[i]]), ])
      colnames(meta_wo_na) <- colnames_tmp
      
      # Only single copy of data for everything else
      adonis_res_rsq = vector()
      adonis_res_pval = vector()
      if (dataset != "human") {
        for (col in names(meta_wo_na)){
          bray <- vegdist(profiles[[i]][!is.na(meta_wo_na[,col]),])
          if (nrow(matrix(bray)) != 0) {
            data_in <- data.frame(meta_wo_na[!is.na(meta_wo_na[,col]),])
            colnames(data_in) <- colnames(meta_wo_na)
            if (length(unique(data_in[,col])) > 1) {
              adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = data_in, permutations = 9999)
              adonis_res_rsq[col] = adonis.univ[1,"R2"]
              adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
            } else {
              adonis_res_rsq[col] = 0
              adonis_res_pval[col] = 1
            }
          } else {
            adonis_res_rsq[col] = 0
            adonis_res_pval[col] = 1
          }
        }
      } else {
        for (col in names(meta_wo_na)[names(meta_wo_na) != "participant"]){
          bray <- vegdist(profiles[[i]][!is.na(meta_wo_na[,col]),])
          
          if (nrow(matrix(bray)) != 0) {
            adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999, strata = meta_wo_na[!is.na(meta_wo_na[,col]),]$participant)
            adonis_res_rsq[col] = adonis.univ[1,"R2"]
            adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
          } else {
            adonis_res_rsq[col] = 0
            adonis_res_pval[col] = 1
          }
          
        }
        col = "participant"
        bray <- vegdist(profiles[[i]][!is.na(meta_wo_na[,col]),])
        
        if (nrow(matrix(bray)) != 0) {
          adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999)
          adonis_res_rsq[col] = adonis.univ[1,"R2"]
          adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
        } else {
          adonis_res_rsq[col] = 0
          adonis_res_pval[col] = 1
        }
      }
      
      df_addition = data.frame("method"=i, "variable"=names(adonis_res_pval), "R2"=adonis_res_rsq, "p_val"=adonis_res_pval, "dataset" = dataset)
      df <- rbind(df, df_addition)
    }
  }
  
  df$method <- tool_names[df$method]
  if (core_only) {
    df <- df[df$method %in% tool_core,]
  }
  df$variable <- gsub("_", " ", df$variable)
  df$variable <- sapply(df$variable, simpleCap)
  
  df$dataset <- case_when(df$dataset == "acid_mine" ~ "Acid mine runoff", 
                          df$dataset == "animal_gut" ~ "Wild animal gut", 
                          df$dataset == "forest_soil" ~ "Forest soil", 
                          df$dataset == "gator_soil" ~ "Gator nest",
                          df$dataset == "saltmarsh" ~ "Salt marsh", 
                          df$dataset == "tara_polar" ~ "Polar ocean",
                          df$dataset == "human" ~ "Human gut")
  
  col = scale_color_manual(values = colAdd)
  
  df$`P-value below` <- df$p_val
  df$`P-value below` <- case_when(df$`P-value below` < 0.001 ~ "0.001",
                                  df$`P-value below` < 0.01 ~ "0.01",
                                  df$`P-value below` < 0.05 ~ "0.05",
                                  TRUE ~ "1")
  df$`P-value below` <- factor(df$`P-value below`, c("1", "0.05", "0.01", "0.001"))
  
  if (core_only) {
    colAdd <- colAdd[tool_names %in% tool_core]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  
  # In-text
  
  range_df <- df %>%
    group_by(dataset, variable) %>%
    summarize(range = max(R2) - min(R2))
  mean(range_df[range_df$dataset %in% c("Human gut", "Wild animal gut"),]$range)
  mean(range_df[!range_df$dataset %in% c("Human gut", "Wild animal gut"),]$range)
  mean(range_df[!range_df$dataset %in% c("Human gut", "Wild animal gut", "Salt marsh", "Acid mine runoff"),]$range)
  
  ggplot(df, aes(x = variable, y = R2, color=factor(method, tool_names), size=`P-value below`)) + 
    geom_beeswarm(cex=3, stroke=2) + 
    theme_linedraw() + 
    ylab("PERMANOVA R-squared") + 
    xlab("Variable") + 
    col + 
    scale_size_manual(values=c(0.5, 1.25, 2, 2.75)) + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(1, 'cm'),
          legend.key.width = unit(1, 'cm'),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Method")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 12, hjust = 1)) + 
    coord_flip() + 
    facet_wrap(~factor(dataset, c("Human gut", "Forest soil", "Gator nest", "Acid mine runoff", "Wild animal gut", "Salt marsh", "Polar ocean")), ncol = 4, scales = "free") + 
    theme(legend.box = "horizontal")
  
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/PERMANOVA_comparison/",dataset,"/"), showWarnings = FALSE)
  if (dataset == "all") {
    if (core_only) {
      ggsave(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_PERMANOVA.png"), width=54, height=16, units = 'cm', dpi=1000, bg='#ffffff')
    } else {
      ggsave(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_PERMANOVA_all_tools.png"), width=54, height=16, units = 'cm', dpi=1000, bg='#ffffff')
    }
  } else {
    if (core_only) {
      ggsave(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_PERMANOVA.png"), width=24, height=16, units = 'cm', dpi=1000, bg='#ffffff')
    } else {
      ggsave(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_PERMANOVA_all_tools.png"), width=24, height=16, units = 'cm', dpi=1000, bg='#ffffff')
    }
  }
}
make_all_PERMANOVA_comparison <- function() {
  for (dataset in c("acid_mine", "animal_gut", "forest_soil", "gator_soil", "saltmarsh", "tara_polar", "human", "all")) {
    for (level in 2:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", "NCBI_only", "_PERMANOVA.png"))) {
        PERMANOVA_comparison(dataset, level, TRUE)
      }
      if (!file.exists(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", "all_taxa", "_PERMANOVA_all_tools.png"))) {
        PERMANOVA_comparison(dataset, level, FALSE, FALSE)
      }
      if (!file.exists(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", "all_taxa", "_PERMANOVA.png"))) {
        PERMANOVA_comparison(dataset, level, FALSE)
      }
    }
  }
}
make_all_PERMANOVA_comparison()

Maaslin_comparison <- function(dataset, level, ncbi_only = TRUE, core_only = TRUE) {
  if (dataset == "all") {
    maaslin_out_final <- c()
    melted_cormat <- c()
    for (dataset_tmp in c("acid_mine", "animal_gut", "forest_soil", "gator_soil", "saltmarsh", "tara_polar", "human")) {
      if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names_tmp <- tool_names[-c(3, 11)]
      } else {
        tool_names_tmp <- tool_names
      }
      
      meta = read.csv(paste0("metadata/", dataset_tmp, "/merged_metadata.tsv"), sep="\t")
      
      data <- preprocess(dataset_tmp, level, real = TRUE, by_dataset = TRUE)
      
      if (ncbi_only) {
        profiles <- remove_non_ncbi(data)
      } else {
        profiles <- remove_unknown(data)
      }
      
      profiles = lapply(profiles, renormalize)
      profiles = lapply(profiles, threshold_sample, threshold)
      profiles = lapply(profiles, renormalize)
      
      for (i in 1:length(profiles)) {
        tmp <- profiles[[i]]
        taxa_names <- tmp$TaxIDs
        tmp <- data.frame(t(tmp[,-1]))
        colnames(tmp) <- as.character(taxa_names)
        profiles[[i]] <- tmp
      }
      
      sample_ids = read.csv(paste0("real_data_outputs/inputs/", dataset_tmp, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)[1,] %>% 
        gsub(pattern="_taxonomic_profile", replacement = "") %>% gsub(pattern="^re", replacement="")
      
      sample_ids <- sort(sample_ids[-1])
      
      meta <- meta[meta$sample_id %in% sample_ids,]
      meta <- meta[order(meta$sample_id),]
      rownames(meta) <- meta$sample_id
      
      place_holders = rownames(profiles[[1]])
      for (i in 1:length(profiles)) {
        rownames(profiles[[i]]) <- mapvalues(rownames(profiles[[i]]), place_holders, sample_ids, warn_missing = F)
        if (ncol(profiles[[i]]) == 1) {
          colname_tmp <- colnames(profiles[[i]])
          rownames_tmp <- rownames(profiles[[i]])
          profiles[[i]] <- data.frame(profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,])
          colnames(profiles[[i]]) <- colname_tmp
          rownames(profiles[[i]]) <- rownames_tmp
          subset_vec <- rowSums(profiles[[i]]) != 0
          profiles[[i]] <- data.frame(profiles[[i]][subset_vec,])
          colnames(profiles[[i]]) <- colname_tmp
          rownames(profiles[[i]]) <- rownames_tmp[subset_vec]
        } else {
          profiles[[i]] <- profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,]
          profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
        }
      }
      
      meta$sample_id <- NULL
      maaslin_out <- data.frame(matrix(nrow = 0, ncol = 7))
      reference_list <- list()
      reference_vec <- c()
      for (colname in colnames(meta)) {
        if (sum(!is.na(as.numeric(meta[,colname]))) == 0) {
          # Categorical
          if (length(unique(meta[,colname])) > 2) {
            reference_list[[colname]] <- unique(meta[,colname])
            reference_vec <- c(reference_vec, paste0("(Reference: ", names(which.max(table(meta[,colname]))), ")"))
          } else {
            reference_vec <- c(reference_vec, paste0("(Reference: ", sort(unique(meta[,colname]))[1], ")"))
          }
        } else {
          reference_vec <- c(reference_vec, "")
        }
      }
      
      names(reference_vec) <- colnames(meta)
      
      for (i in 1:length(profiles)) {
        maaslin_tmp <- data.frame(matrix(nrow = 0, ncol = 5))
        if (length(reference_list) > 0) {
          reference_string = ""
          for (j in 1:length(reference_list)) {
            reference_string = paste0(reference_string, ";", names(reference_list[j]), ",", names(which.max(table(meta[,names(reference_list[j])]))))
          }
          reference_string = gsub("^;", "", reference_string)
          
          if (ncol(profiles[[i]]) > 1) {
            if (dataset != "human") {
              maaslin_results <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", reference = reference_string, plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1)$results[,c(1, 2, 3, 4, 6)]
            } else {
              maaslin_results <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", reference = reference_string, plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1, random_effects = "participant")$results[,c(1, 2, 3, 4, 6)]
            }
            maaslin_tmp <- rbind(maaslin_tmp, maaslin_results)
          }
        } else {
          if (ncol(profiles[[i]]) > 1) {
            if (dataset != "human") {
              maaslin_tmp <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1)$results[,c(1, 2, 3, 4, 6)]
            } else {
              maaslin_tmp <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1, random_effects = "participant")$results[,c(1, 2, 3, 4, 6)]
            }        }
        }
        
        if (nrow(maaslin_tmp) > 0) {
          maaslin_tmp <- distinct(maaslin_tmp, feature, metadata, value, .keep_all = TRUE)
          maaslin_tmp$qval <- p.adjust(maaslin_tmp$pval, "BH")
          maaslin_tmp$method = i
          maaslin_tmp <- maaslin_tmp[maaslin_tmp$qval < 0.25,]
          
          maaslin_out <- rbind(maaslin_out, maaslin_tmp)
        }
      }
      
      maaslin_out$value <- paste0(maaslin_out$value, "\n", reference_vec[maaslin_out$metadata])
      
      dataset_tmp <- case_when(dataset_tmp == "acid_mine" ~ "Acid mine runoff", 
                               dataset_tmp == "animal_gut" ~ "Wild animal gut", 
                               dataset_tmp == "forest_soil" ~ "Forest soil", 
                               dataset_tmp == "gator_soil" ~ "Gator nest",
                               dataset_tmp == "saltmarsh" ~ "Salt marsh", 
                               dataset_tmp == "tara_polar" ~ "Polar ocean",
                               dataset_tmp == "human" ~ "Human gut")
      
      maaslin_out$joined_name <- paste0(dataset_tmp, ": ", maaslin_out$feature, ": ", maaslin_out$value)
      
      cormat = matrix(nrow = length(data), ncol = length(data))
      for (i in 1:nrow(cormat)) {
        for (j in 1:ncol(cormat)) {
          cormat[i,j] <- sum(maaslin_out$joined_name[maaslin_out$method == i][!is.na(maaslin_out$joined_name[maaslin_out$method == i])] %in% maaslin_out$joined_name[maaslin_out$method == j]) / 
            length(unique(c(maaslin_out$joined_name[maaslin_out$method == i][!is.na(maaslin_out$joined_name[maaslin_out$method == i])], maaslin_out$joined_name[maaslin_out$method == j][!is.na(maaslin_out$joined_name[maaslin_out$method == j])])))
        }
      }
      if (max(cormat, na.rm = T) > 1) {
        stop()
      }
      
      rownames(cormat) <- tool_names_tmp
      colnames(cormat) <- tool_names_tmp
      
      melted_cormat_tmp <- reshape2::melt(get_upper_tri(cormat), na.rm = T)
      if (dataset_tmp != "Human gut") {
        melted_cormat <- rbind(melted_cormat, melted_cormat_tmp)
      }
      
      maaslin_out$method <- tool_names_tmp[maaslin_out$method]
      maaslin_out$dataset <- dataset_tmp
      maaslin_out_final <- rbind(maaslin_out_final, maaslin_out)
    }
    melted_cormat <- aggregate(value ~ Var1 + Var2, melted_cormat, FUN = mean)
  } else {
    melted_cormat <- c()
    if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
      tool_names_tmp <- tool_names[-c(3, 11)]
    } else {
      tool_names_tmp <- tool_names
    }
    
    meta = read.csv(paste0("metadata/", dataset, "/merged_metadata.tsv"), sep="\t")
    
    data <- preprocess(dataset, level, real = TRUE, by_dataset = TRUE)
    
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    profiles = lapply(profiles, renormalize)
    
    for (i in 1:length(profiles)) {
      tmp <- profiles[[i]]
      taxa_names <- tmp$TaxIDs
      tmp <- data.frame(t(tmp[,-1]))
      colnames(tmp) <- as.character(taxa_names)
      profiles[[i]] <- tmp
    }
    
    sample_ids = read.csv(paste0("real_data_outputs/inputs/", dataset, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)[1,] %>% 
      gsub(pattern="_taxonomic_profile", replacement = "") %>% gsub(pattern="^re", replacement="")
    
    sample_ids <- sort(sample_ids[-1])
    
    meta <- meta[meta$sample_id %in% sample_ids,]
    meta <- meta[order(meta$sample_id),]
    rownames(meta) <- meta$sample_id
    
    place_holders = rownames(profiles[[1]])
    for (i in 1:length(profiles)) {
      rownames(profiles[[i]]) <- mapvalues(rownames(profiles[[i]]), place_holders, sample_ids, warn_missing = F)
      if (ncol(profiles[[i]]) == 1) {
        colname_tmp <- colnames(profiles[[i]])
        rownames_tmp <- rownames(profiles[[i]])
        profiles[[i]] <- data.frame(profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,])
        colnames(profiles[[i]]) <- colname_tmp
        rownames(profiles[[i]]) <- rownames_tmp
        subset_vec <- rowSums(profiles[[i]]) != 0
        profiles[[i]] <- data.frame(profiles[[i]][subset_vec,])
        colnames(profiles[[i]]) <- colname_tmp
        rownames(profiles[[i]]) <- rownames_tmp[subset_vec]
      } else {
        profiles[[i]] <- profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,]
        profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
      }
    }
    
    meta$sample_id <- NULL
    maaslin_out <- data.frame(matrix(nrow = 0, ncol = 7))
    reference_list <- list()
    reference_vec <- c()
    for (colname in colnames(meta)) {
      if (sum(!is.na(as.numeric(meta[,colname]))) == 0) {
        # Categorical
        if (length(unique(meta[,colname])) > 2) {
          reference_list[[colname]] <- unique(meta[,colname])
          reference_vec <- c(reference_vec, paste0("(Reference: ", names(which.max(table(meta[,colname]))), ")"))
        } else {
          reference_vec <- c(reference_vec, paste0("(Reference: ", sort(unique(meta[,colname]))[1], ")"))
        }
      } else {
        reference_vec <- c(reference_vec, "")
      }
    }
    
    names(reference_vec) <- colnames(meta)
    
    for (i in 1:length(profiles)) {
      maaslin_tmp <- data.frame(matrix(nrow = 0, ncol = 5))
      if (length(reference_list) > 0) {
        reference_string = ""
        for (j in 1:length(reference_list)) {
          reference_string = paste0(reference_string, ";", names(reference_list[j]), ",", names(which.max(table(meta[,names(reference_list[j])]))))
        }
        reference_string = gsub("^;", "", reference_string)
        
        if (ncol(profiles[[i]]) > 1) {
          if (dataset != "human") {
            maaslin_results <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", reference = reference_string, plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1)$results[,c(1, 2, 3, 4, 6)]
          } else {
            maaslin_results <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", reference = reference_string, plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1, random_effects = "participant")$results[,c(1, 2, 3, 4, 6)]
          }
          maaslin_tmp <- rbind(maaslin_tmp, maaslin_results)
        }
      } else {
        if (ncol(profiles[[i]]) > 1) {
          if (dataset != "human") {
            maaslin_tmp <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1)$results[,c(1, 2, 3, 4, 6)]
          } else {
            maaslin_tmp <- Maaslin2(profiles[[i]], meta, "scripts/cache/maaslin", plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 1, random_effects = "participant")$results[,c(1, 2, 3, 4, 6)]
          }        }
      }
      
      if (nrow(maaslin_tmp) > 0) {
        maaslin_tmp <- distinct(maaslin_tmp, feature, metadata, value, .keep_all = TRUE)
        maaslin_tmp$qval <- p.adjust(maaslin_tmp$pval, "BH")
        maaslin_tmp$method = i
        maaslin_tmp <- maaslin_tmp[maaslin_tmp$qval < 0.25,]
        
        maaslin_out <- rbind(maaslin_out, maaslin_tmp)
      }
    }
    maaslin_out$value <- paste0(maaslin_out$value, "\n", reference_vec[maaslin_out$metadata])
    
    study_name <- case_when(dataset == "acid_mine" ~ "Acid mine runoff", 
                         dataset == "animal_gut" ~ "Wild animal gut", 
                         dataset == "forest_soil" ~ "Forest soil", 
                         dataset == "gator_soil" ~ "Gator nest",
                         dataset == "saltmarsh" ~ "Salt marsh", 
                         dataset == "tara_polar" ~ "Polar ocean",
                         dataset == "human" ~ "Human gut")
    
    if (nrow(maaslin_out) > 0) {
      maaslin_out$joined_name <- paste0(study_name, ": ", maaslin_out$feature, ": ", maaslin_out$value)
    }
    
    cormat = matrix(nrow = length(data), ncol = length(data))
    for (i in 1:nrow(cormat)) {
      for (j in 1:ncol(cormat)) {
        cormat[i,j] <- sum(maaslin_out$joined_name[maaslin_out$method == i][!is.na(maaslin_out$joined_name[maaslin_out$method == i])] %in% 
                             maaslin_out$joined_name[maaslin_out$method == j][!is.na(maaslin_out$joined_name[maaslin_out$method == j])]) / 
          length(unique(c(maaslin_out$joined_name[maaslin_out$method == i][!is.na(maaslin_out$joined_name[maaslin_out$method == i])], 
                          maaslin_out$joined_name[maaslin_out$method == j][!is.na(maaslin_out$joined_name[maaslin_out$method == j])])))
      }
    }
    
    rownames(cormat) <- tool_names_tmp
    colnames(cormat) <- tool_names_tmp
    
    melted_cormat_tmp <- reshape2::melt(get_upper_tri(cormat), na.rm = T)
    melted_cormat <- rbind(melted_cormat, melted_cormat_tmp)
    
    if (nrow(melted_cormat) > 0) {
      melted_cormat <- aggregate(value ~ Var1 + Var2, melted_cormat, FUN = mean)
    }
    
    if (nrow(maaslin_out) > 0) {
      maaslin_out$method <- tool_names_tmp[maaslin_out$method]
      maaslin_out$dataset <- study_name
    }
    maaslin_out_final <- maaslin_out
  }
  
  maaslin_out_final$joined_name <- gsub("\n$", "", maaslin_out_final$joined_name)
  
  if (core_only) {
    melted_cormat <- melted_cormat[melted_cormat$Var1 %in% c(tool_core) & melted_cormat$Var2 %in% c(tool_core),]
  }
  melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2,]
  
  melted_cormat$value <- round(melted_cormat$value, 2)
  ggheatmap <- ggplot(melted_cormat, aes(factor(Var2, tool_names), factor(Var1, tool_names), fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), space = "Lab", 
                         name="Intersection\n over union") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(Var2, tool_names), factor(Var1, tool_names), label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.5, 0.7),
      legend.direction = "horizontal",
      legend.title=element_text(size=14)) +
    guides(fill = guide_colorbar(barwidth = 8, barheight = 1.5,
                                 title.position = "top", title.hjust = 0.5))
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/Maaslin_comparison/",dataset,"/"), showWarnings = FALSE)
  if (core_only) {
    ggsave(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_Maaslin_grid.png"), width=18, height=15, units = 'cm', dpi=1000, bg='#ffffff')
  } else {
    ggsave(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_Maaslin_grid_all_tools.png"), width=18, height=15, units = 'cm', dpi=1000, bg='#ffffff')
  }
  
  if (core_only) {
    maaslin_out_final <- maaslin_out_final[maaslin_out_final$method %in% tool_core,]
  }
  
  if (dataset == "all") {
    signif_assoc <- table(maaslin_out_final$method, maaslin_out_final$dataset)
    signif_df <- reshape2::melt(signif_assoc)
    signif_df$org_value <- round(signif_df$value, 2)
    signif_df$value <- log(signif_df$value)
    signif_df$value[signif_df$value == -Inf] <- 0
    
    ggheatmap <- ggplot(signif_df, aes(factor(Var1, tool_names), factor(Var2, rev(c("Human gut", "Wild animal gut", "Forest soil", "Gator nest", "Acid mine runoff", "Salt marsh", "Polar ocean"))), fill = value)) +
      geom_tile(color = "white")+
      scale_fill_viridis_c(option = "plasma", begin = 0.3) + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                       size = 12, hjust = 1),
            axis.text.y = element_text(size = 12))+
      coord_fixed()
    
    ggheatmap + 
      geom_text(aes(factor(Var1, tool_names), factor(Var2, rev(c("Human gut", "Wild animal gut", "Forest soil", "Gator nest", "Acid mine runoff", "Salt marsh", "Polar ocean"))), label = org_value), color = "black", size = 4) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
    
    if (core_only) {
      ggsave(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_Maaslin_counts.png"), width=24, height=12, units = 'cm', dpi=1000, bg='#ffffff')
    } else {
      ggsave(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_Maaslin_counts_all_tools.png"), width=24, height=12, units = 'cm', dpi=1000, bg='#ffffff')
    }
  }
  
  tmp_sort_list <- sort(table(maaslin_out_final[maaslin_out_final$dataset != "Human gut",]$joined_name), decreasing=TRUE)
  keep_sort_list <- c()
  for (dataset_item in unique(gsub("\\:.*", "", names(tmp_sort_list)))) {
    keep_sort_list <- c(keep_sort_list, tmp_sort_list[gsub("\\:.*", "", names(tmp_sort_list)) == dataset_item][1:2])
  }
  keep_sort_list <- sort(keep_sort_list, decreasing = T)
  
  overlap_maaslin <- maaslin_out_final[maaslin_out_final$joined_name %in% names(keep_sort_list),]
  
  level_list <- fread(paste0("databases/ncbi_taxdump/names.dmp"), sep="\t", header = F)
  level_list <- level_list[,c(1,3,7)]
  colnames(level_list) <- c("TaxID", "Name", "Status")
  level_list <- level_list[level_list$Status == "scientific name" | (!duplicated(level_list$Status) & !duplicated(level_list$Status, fromLast = T)),]
  tax_list <- level_list$Name
  tax_list <- c(tax_list, overlap_maaslin$feature[is.na(as.numeric(gsub("^X", "", overlap_maaslin$feature)))])
  names(tax_list) <- c(level_list$TaxID, overlap_maaslin$feature[is.na(as.numeric(gsub("^X", "", overlap_maaslin$feature)))])
  rm(level_list)
  
  if (nrow(overlap_maaslin) > 0) {
    overlap_maaslin$joined_name <- paste0(gsub(":.*", "", overlap_maaslin$joined_name), ": ",
                                          tax_list[gsub("^X", "", overlap_maaslin$feature)], ":", 
                                          sub(".*:(?!\n)", "", overlap_maaslin$joined_name, perl = TRUE))
  }
  
  colAdd <- c("#CF9FFF", "#A15BE4", "#00FFFF", "#0061FE", "#1434A4", "#81C784", "#2E7D32", "#FF7300", "#FBB15B", "#AC0911", "#E95420")
  col = scale_color_manual(values = colAdd[tool_names %in% unique(overlap_maaslin$method)])
  
  overlap_maaslin$`Q-value below` <- overlap_maaslin$qval
  overlap_maaslin$`Q-value below` <- case_when(overlap_maaslin$`Q-value below` < 0.001 ~ "0.001",
                                             overlap_maaslin$`Q-value below` < 0.01 ~ "0.01",
                                             overlap_maaslin$`Q-value below` < 0.05 ~ "0.05",
                                             TRUE ~ "0.25")
  overlap_maaslin$`Q-value below` <- factor(overlap_maaslin$`Q-value below`, c("0.25", "0.05", "0.01", "0.001"))
  
  if (nrow(overlap_maaslin) > 0) {
    ggplot(overlap_maaslin, aes(y=coef, x=factor(joined_name, rev(sort(unique(overlap_maaslin$joined_name)))), color = factor(method, tool_names), size = `Q-value below`)) + 
      geom_beeswarm(cex=2.5, stroke=2) + 
      theme_bw() + 
      ylab("Effect size") + 
      xlab("Association") + 
      col + 
      scale_size_manual(values=c(1, 1.75, 2.5, 3.25)) + 
      theme(text=element_text(size=14),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=12),
            legend.position = "bottom",
            legend.direction = "horizontal") + 
      guides(color="none") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                       size = 14, hjust = 1)) + 
      geom_hline(yintercept = 0, linetype="dashed") + 
      coord_flip()
    
    if (core_only) {
      ggsave(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_Maaslin_effects.png"), width=28, height=20, units = 'cm', dpi=1000, bg='#ffffff')
    } else {
      ggsave(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_Maaslin_effects_all_tools.png"), width=28, height=20, units = 'cm', dpi=1000, bg='#ffffff')
    }
  } else {
    file.create(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_Maaslin_effects.png"))
  }
}
make_all_Maaslin_comparison <- function() {
  for (dataset in c("acid_mine", "animal_gut", "forest_soil", "gator_soil", "saltmarsh", "tara_polar", "human", "all")) {
    for (level in 2:7) {
      ncbi_only <- TRUE
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", "NCBI_only", "_Maaslin_effects.png"))) {
        Maaslin_comparison(dataset, level)
      }
      if (!file.exists(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", "NCBI_only", "_Maaslin_effects_all_tools.png"))) {
        Maaslin_comparison(dataset, level, TRUE, FALSE)
      }
      if (!file.exists(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", "all_taxa", "_Maaslin_effects.png"))) {
        Maaslin_comparison(dataset, level, ncbi_only = FALSE)
      }
      if (!file.exists(paste0("figures/Maaslin_comparison/",dataset,"/",dataset,"_",level_name,"_", "all_taxa", "_Maaslin_effects_all_tools.png"))) {
        Maaslin_comparison(dataset, level, ncbi_only = FALSE, core_only = FALSE)
      }
    }
  }
}
make_all_Maaslin_comparison()

# Plot taxa called per method
taxa_num_comparison <- function(dataset, level, ncbi_only) {
  df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(df) <- c("method", "count", "sample_id", "threshold")
  for (threshold_status in c("yes", "no")) {
    if (dataset == "all") {
      
      for (dataset_tmp in list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("human", "overall_dists")]) {
        tool_names_tmp <- tool_names
        if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
          tool_names_tmp <- tool_names[-c(3, 11)]
        }
        set.seed(1)
        
        data <- preprocess(dataset_tmp, level)
        
        if (ncbi_only) {
          profiles <- remove_non_ncbi(data)
        } else {
          profiles <- remove_unknown(data)
        }
        
        if (threshold_status == "yes") {
          profiles = lapply(profiles, renormalize)
          profiles = lapply(profiles, threshold_sample, threshold)
        }
        profiles = lapply(profiles, renormalize)
        
        df_tmp <- data.frame(matrix(nrow = 0, ncol = 3))
        colnames(df_tmp) <- c("method", "count", "sample_id")
        for (i in 1:length(profiles)) {
          if (threshold_status == "yes") {
            count_taxa = sapply(list_taxa_by_sample(profiles[[i]], threshold), length)
          } else {
            count_taxa = sapply(list_taxa_by_sample(profiles[[i]], 0), length)
          }
          df_append <- data.frame("method"=i, 
                                  "count" = count_taxa,
                                  "sample_id" = colnames(profiles[[i]])[-1],
                                  "threshold" = threshold_status)
          df_tmp <- rbind(df_tmp, df_append)
        }
        
        df_tmp$method <- tool_names_tmp[df_tmp$method]
        
        df <- rbind(df, df_tmp)
      }
    } else {
      if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names <- tool_names[-c(9, 11)]
        colAdd <- colAdd[-c(9, 11)]
        col = scale_color_manual(values = colAdd)
        fil = scale_fill_manual(values = colAdd)
      }
      set.seed(1)
      
      data <- preprocess(dataset, level)
      
      if (ncbi_only) {
        profiles <- remove_non_ncbi(data)
      } else {
        profiles <- remove_unknown(data)
      }
      
      if (threshold_status == "yes") {
        profiles = lapply(profiles, renormalize)
        profiles = lapply(profiles, threshold_sample, threshold)
      }
      profiles = lapply(profiles, renormalize)
      
      for (i in 1:length(profiles)) {
        if (threshold_status == "yes") {
          count_taxa = sapply(list_taxa_by_sample(profiles[[i]], threshold), length)
        } else {
          count_taxa = sapply(list_taxa_by_sample(profiles[[i]], 0), length)
        }
        df_append <- data.frame("method"=i, 
                                "count" = count_taxa, 
                                "sample_id" = colnames(profiles[[i]])[-1],
                                "threshold" = threshold_status)
        df <- rbind(df, df_append)
      }
      
      
    }
  }
  
  if (dataset != "all") {
    df$method <- tool_names[df$method]
  }
  
  study = case_when(dataset == "acid_mine" ~ "Acid Mine Drainage",
                    dataset == "animal_gut" ~ "Wild Animal Gut",
                    dataset == "cat_gut" ~ "Cat Gut",
                    dataset == "coastal_sediment" ~ "Coastal Sediment",
                    dataset == "dog_gut" ~ "Dog Gut",
                    dataset == "forest_soil" ~ "Forest Soil",
                    dataset == "gator_soil" ~ "Gator Nest",
                    dataset == "human" ~ "Human Gut",
                    dataset == "saltmarsh" ~ "Salt Marsh",
                    dataset == "tara_polar" ~ "Polar ocean",
                    dataset == "all" ~ "All datasets")
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  study <- paste0(study, ", ", level_name, " level")
  
  df$threshold <- ifelse(df$threshold == "yes", 
                         paste0("Threshold at ", threshold),
                         paste0("No threshold"))
  
  # Plot the overlaps
  ggplot(df, aes(y=count, x=factor(method, c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")), group=sample_id)) + 
    geom_line(aes(alpha=0.5, color=sample_id)) + 
    geom_point(aes(alpha=0.5, color=sample_id)) + 
    theme_linedraw() +
    theme(legend.position="none") + 
    ggtitle(study) +
    xlab("Method") +
    ylab(paste0("Number of ", level_name, " assigned")) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    facet_wrap(~threshold, scale = "free")
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("Figures/taxa_num/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/taxa_num/",dataset,"/",dataset,"_",level_name,"_", ncbi_only,".png"), width=10, height=5)
}
make_all_taxa_num_comparison <- function() {
  for (dataset in c(list.files("real_data_outputs/inputs")[!list.files("real_data_outputs/inputs") %in% c("overall_dists")], "all")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("Figures/taxa_num/",dataset,"/",dataset,"_",level_name,"_", "NCBI_only",".png"))) {
        taxa_num_comparison(dataset, level, TRUE)
      }
      if (!file.exists(paste0("Figures/taxa_num/",dataset,"/",dataset,"_",level_name,"_", "all_taxa",".png"))) {
        taxa_num_comparison(dataset, level, FALSE)
      }
    }
  }
}
make_all_taxa_num_comparison()

