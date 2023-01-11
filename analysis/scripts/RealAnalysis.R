rm(list = ls())
setwd('..')

source("scripts/helpers.R")

preprocess_all_real()

threshold = 0.08
tool_names = c("Centrifuge", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "Kraken 2 / Braken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "mOTUs 3", "Metaxa 2", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
colAdd <- c("#FF96D0", "#FF7300", "#FBB15B", "#A15BE4", "#00FFFF", "#00CCFF", "#0061FE", "#81C784", "#2E7D32", "#AC0911", "#E95420")
col = scale_color_manual(values = colAdd)
fil = scale_fill_manual(values = colAdd)

# Plot PCoA from Bray Curtis dissimilarity on normalized abundances
PCoABray <- function(dataset, level, ncbi_only) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
    tool_names <- tool_names[-c(3, 11)]
    colAdd <- colAdd[-c(3, 11)]
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
  distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
  distances[is.na(distances)] <- 1
  
  pc <- capscale(distances~1, comm = distances, na.action = na.omit)
  cap = data.frame(pc$CA$u)
  cap$Method <- tool_names[gsub(".*_", "", rownames(cap)) %>% as.numeric()]

  cap$Sample <- gsub("_.*", "", rownames(cap))
  s = summary(pc)
  
  # Swap study name for plotting
  study = case_when(dataset == "acid_mine" ~ "Acid Mine Drainage",
                     dataset == "animal_gut" ~ "Wild Animal Gut",
                     dataset == "cat_gut" ~ "Cat Gut",
                     dataset == "coastal_sediment" ~ "Coastal Sediment",
                     dataset == "dog_gut" ~ "Dog Gut",
                     dataset == "forest_soil" ~ "Forest Soil",
                     dataset == "gator_soil" ~ "Gator Nest Soil",
                     dataset == "human" ~ "Human Gut",
                     dataset == "saltmarsh" ~ "Salt Marsh",
                     dataset == "tara_polar" ~ "Tara Polar",)
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  study <- paste0(study, ", ", level_name)
  
  # Plot ordination
  ggplot(cap, aes(MDS1, MDS2)) + 
    geom_point(aes(color=factor(Method, levels=tool_names)), size = 1, alpha = 0.7) + 
    theme_classic() +coord_fixed() + 
    ggtitle(study) +
    labs( color = "Tool") + 
    labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
         y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = '')) +
    col + guides(color=guide_legend(title="Method"))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/PCoABray/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_", ncbi_only,"_PCoA.png"), width=14, height=12, units = 'cm', dpi=1000)
}

run_all_PCoABray <- function() {
  for (dataset in list.files("real_data_outputs/inputs")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdoms",
                              level == 2 ~ "phyla", 
                              level == 3 ~ "classes",
                              level == 4 ~ "orders",
                              level == 5 ~ "families",
                              level == 6 ~ "genera",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_NCBI_only_PCoA.png")) | !file.exists(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_all_taxa_PCoA.png"))) {
        PCoABray(dataset, level, TRUE)
        PCoABray(dataset, level, FALSE)
      }
    }
  }
}













# Return list of taxonomic names overlapping from k methods
get_k_way_overlaps <- function(input_list, k) {
  combns <- combn(length(input_list), k)
  total <- c()
  for (i in 1:ncol(combns)) {
    intersecting <- combns[,i]
    nonintersecting <- setdiff(1:length(input_list), intersecting)
    temp <- input_list[[intersecting[1]]]
    for (j in 2:length(intersecting)) {
      temp <- intersect(temp, input_list[[intersecting[j]]])
    }
    for (j in 1:length(nonintersecting)) {
      temp <- setdiff(temp, input_list[[nonintersecting[j]]])
    }
    total <- c(total, temp)
  }
  
  return(total)
}

renormalize <- function(column, coef) {
  if (is.infinite(coef)) {coef = 0}
  return (as.numeric(column) * coef)
}



# Get dataframe of samples and reads
preprocessReads <- function(dataset) {
  if (dataset=="All") {
    datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar")
    reads <- matrix(nrow = 0, ncol = 2)
    for (name in datasets) {
      reads <- rbind(reads, read.csv(paste(name, "/reads.tsv",sep = ""), sep="\t", header = T))
    }
    reads <- as.data.frame(reads)
    colnames(reads) <- c("ID", "reads")
    return (data.frame("Sample"=reads$ID, "Reads"=reads$reads))
  }
  
  reads <- read.csv(paste(dataset, "/reads.tsv",sep = ""), sep="\t", header = T)
  return (data.frame("Sample"=reads$ID, "Reads"=reads$reads))
}

# Plot the overlap between method assignments with p-values
# For level, kingdom = 1, phylum = 2, class = 3, order = 4, family = 5, genus = 6, species = 7
# Not updated with kraken
multipleFoundPlot <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  taxa_names <- list()
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    for (j in 2:ncol(tmp)) {
      tmp[,j] <- ifelse(tmp[,j]!=0, tmp[,1], NA_character_)
    }
    tmp <- tmp[!tmp$Taxa=="UNCLASSIFIED",-1]
    taxa_names <- append(taxa_names, list(tmp))
  }
  
  # 11 rows because 11 categories of overlap or nonoverlap
  whole_dataset <- data.frame(matrix(nrow = 11, ncol = ncol(taxa_names[[1]])))
  colnames(whole_dataset) <- colnames(taxa_names[[1]])
  
  # Proportion that could be called by all databases
  all_potential <- as.list(vector(length = 11))
  
  for (i in 1:ncol(taxa_names[[1]])) {
    MPA2 <- setdiff(taxa_names[[1]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
                       setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA2 <- MPA2[!is.na(MPA2)]
    MPA3 <- setdiff(taxa_names[[2]][,i], taxa_names[[1]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA3 <- MPA3[!is.na(MPA3)]
    MPA4 <- setdiff(taxa_names[[3]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[1]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA4 <- MPA4[!is.na(MPA4)]
    metaxa <- setdiff(taxa_names[[4]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[1]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    metaxa <- metaxa[!is.na(metaxa)]
    motus <- setdiff(taxa_names[[5]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[1]][,i]) %>% setdiff(taxa_names[[6]][,i])
    motus <- MPA2[!is.na(motus)]
    phylophlan <- setdiff(taxa_names[[6]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[1]][,i])
    phylophlan <- phylophlan[!is.na(phylophlan)]
    
    way2 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                  taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 2)
    way2 <- way2[!is.na(way2)]
    way3 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                  taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 3)
    way3 <- way3[!is.na(way3)]
    way4 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                  taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 4)
    way4 <- way4[!is.na(way4)]
    way5 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                  taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 5)
    way5 <- way5[!is.na(way5)]
    all_overlap <- Reduce(intersect, list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                         taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]))
    all_overlap <- all_overlap[!is.na(all_overlap)]
    
    # Need to return proportion of all database occurrences
    whole_dataset[,i] <- c(length(MPA2), length(MPA3), length(MPA4), length(metaxa), length(motus), length(phylophlan),
                           length(way2), length(way3), length(way4), length(way5), length(all_overlap))
    
    all_potential[[1]] <- c(all_potential[[1]], MPA2)
    all_potential[[2]] <- c(all_potential[[2]], MPA3)
    all_potential[[3]] <- c(all_potential[[3]], MPA4)
    all_potential[[4]] <- c(all_potential[[4]], metaxa)
    all_potential[[5]] <- c(all_potential[[5]], motus)
    all_potential[[6]] <- c(all_potential[[6]], phylophlan)
    all_potential[[7]] <- c(all_potential[[7]], way2)
    all_potential[[8]] <- c(all_potential[[8]], way3)
    all_potential[[9]] <- c(all_potential[[9]], way4)
    all_potential[[10]] <- c(all_potential[[10]], way5)
    all_potential[[11]] <- c(all_potential[[11]], all_overlap)
  }
  
  for (i in 1:length(all_potential)) {
    all_potential[[i]] <- all_potential[[i]][all_potential[[i]]!="FALSE"]
  }
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  all_dbs <- read.csv(paste0("Databases/DB_merged/", level_name, ".csv"), sep = ",")
  tmp <- all_dbs[,1]
  for (i in 2:length(ncol(all_dbs))) {
    tmp <- intersect(tmp, all_dbs[,i])
  }
  
  prop_in_all <- vector(length = length(all_potential))
  for (i in 1:length(all_potential)) {
    prop_in_all[i] <- 1 - length(all_potential[[i]][all_potential[[i]] %in% setdiff(all_potential[[i]], tmp)]) / length(all_potential[[i]])
  }
  prop_in_all <- ifelse(is.na(prop_in_all), 0, prop_in_all)
  
  chyper_db <- read.csv(paste0("Databases/DB_merged/for_chyper.csv"), sep = ",")
  
  for_chyper <- chyper_db[,(level + 1)]
  
  # Update to use pchyper
  ns <- rep(list(for_chyper[-1]), ncol(taxa_names[[1]]))
  s <- rep(list(for_chyper[1]), ncol(taxa_names[[1]]))
  
  tmp <- t(data.frame(lapply(taxa_names, function (method) colSums(!is.na(method)))))
  ms <- as.list(vector(length = ncol(tmp)))
  for (i in 1:ncol(tmp)) {
    ms[[i]] <- unname(tmp[,i])
  }
  
  overlapall <- as.list(whole_dataset[nrow(whole_dataset),])
  
  ps <- mapply(pvalchyper, overlapall, s, ns, ms, "upper", verbose=F)
  ps <- signif(ps, 1)
  
  orders <- (as.numeric(whole_dataset[nrow(whole_dataset),]) + 0.00001)/(colSums(whole_dataset) + 0.00001)
  
  # Make dataframe for plotting
  whole_dataset$labels <- c("Unique to MPA2", "Unique to MPA3", "Unique to MPA4", 
                               "Unique to metaxa", "Unique to mOTUs3", "Unique to Phylophlan",
                               "2-way overlaps", "3-way overlaps", "4-way overlaps", "5-way overlaps",
                               "All overlapping")
  
  df <- melt(whole_dataset, id.vars = "labels")
  
  
  
  # Sort and add p-values
  tmp <- data.frame(names=names(orders))
  tmp$orders <- unname(orders)
  tmp$ps <- ps
  
  df <- merge(df, tmp, by.x = "variable", by.y = "names")
  
  df <- df %>% arrange(-orders)
  df$ps <- ifelse(df$ps == 0, "<1e-16", df$ps)
  
  label <- paste("Proportion of ",level_name," overlapping")
  
  # Change study name for plot
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_fill_manual(values = colAdd[1:11])
  
  # Plot the overlaps
  ggplot(df, aes(fill=factor(labels, levels=c("All overlapping", "5-way overlaps", 
                                              "4-way overlaps", "3-way overlaps", 
                                              "2-way overlaps", "Unique to Phylophlan",
                                              "Unique to mOTUs3", "Unique to metaxa",
                                              "Unique to MPA4", "Unique to MPA3",
                                              "Unique to MPA2")), 
                             y=value, x=reorder(variable, -orders))) + 
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    theme(axis.text.x=element_text(angle = 90),
          axis.ticks.x = element_blank(),
          legend.title = element_blank()) + 
    ggtitle(study) +
    xlab("Sample (labeled by p-value)") +
    ylab(label) + 
    scale_x_discrete(labels=df$ps[df$labels=="Unique to MPA2"]) + 
    col
    dir.create(file.path("Figures/Prevalence/",dataset,"/"), showWarnings = FALSE)
    ggsave(paste0("Figures/Prevalence/",dataset,"/",dataset,"_",level_name,"_Misclassified.png"), width=10, height=5)
    
    df2 <- data.frame(value=prop_in_all)
    df2$labels <- whole_dataset$labels
    
    # Plot how many of each overlap type could be called by all tools
    ggplot(data = df2, aes(x = factor(labels, levels=c("All overlapping", "5-way overlaps", 
                                                       "4-way overlaps", "3-way overlaps", 
                                                       "2-way overlaps", "Unique to Phylophlan",
                                                       "Unique to mOTUs3", "Unique to metaxa",
                                                       "Unique to MPA4", "Unique to MPA3",
                                                       "Unique to MPA2")),
                           y = value, fill=factor(labels, levels=c("All overlapping", "5-way overlaps", 
                                                                   "4-way overlaps", "3-way overlaps", 
                                                                   "2-way overlaps", "Unique to Phylophlan",
                                                                   "Unique to mOTUs3", "Unique to metaxa",
                                                                   "Unique to MPA4", "Unique to MPA3",
                                                                   "Unique to MPA2")))) + 
      ggtitle(study) +
      geom_bar(stat = "identity") + 
      theme_classic(base_size = 10)  + 
      theme(legend.title = element_blank(),
            axis.text.x = element_blank()) + 
      ylab(paste0("Proportion of ",level_name, " that could\nhave been identified by all tools")) + 
      ylim(0,1) + xlab("Classification") +
      col
    
    dir.create(file.path("Figures/Misclassified/",dataset,"/"), showWarnings = FALSE)
    ggsave(paste0("Figures/Misclassified/",dataset,"/",dataset,"_",level_name,"_Misclassified.png"), width=9, height=12, units = 'cm', dpi=1000)
    
}
make_all_prevalence_plots <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      multipleFoundPlot(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
}

# Database effect
databaseEffect <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  all_dbs <- read.csv(paste0("Databases/DB_merged/", level_name, ".csv"), sep = ",")
  tmp <- all_dbs[,1]
  for (i in 2:ncol(all_dbs)) {
    tmp <- intersect(tmp, all_dbs[,i])
  }
  
  meanUnique <- vector()
  for (i in 1:length(data)) {
    renormed <- data[[i]]
    unknownAmts <- as.numeric(renormed[renormed$TaxIDs=="UNCLASSIFIED",-1])
    classified <- (100 - as.numeric(renormed[renormed$TaxIDs=="UNCLASSIFIED",-1])) / 100
    renormed <- renormed[renormed$TaxIDs!="UNCLASSIFIED",]
    renormed[,-1] <- mapply(renormalize, renormed[,-1], 1/classified)
    
    cut_data <- colSums(renormed[renormed$TaxIDs %in% setdiff(renormed$TaxIDs, tmp),][,-1])
    meanUnique <- append(meanUnique, mean(cut_data))
    
    
  }
  
  whole_dataset <- data.frame(Proportion = meanUnique)
  
  # Make dataframe for plotting
  whole_dataset$labels <- c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "metaxa2")
  
  label <- paste("Proportion of ",level_name," overlapping")
  
  colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_fill_manual(values = colAdd[1:11])
  
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  # Plot how many of each overlap type could be called by all tools
  ggplot(data = whole_dataset, aes(x = factor(labels, levels=c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "metaxa2")),
                         y = Proportion, fill=factor(labels, levels=c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "metaxa2")))) + 
    ggtitle(study) +
    geom_bar(stat = "identity") + 
    theme_classic(base_size = 10)  + 
    theme(legend.title = element_blank(),
          axis.text.x = element_blank()) + 
    ylab(paste0("Abundance of ",level_name, " that could not \nhave been identified by all tools")) + 
    ylim(0,100) + xlab("Classification") +
    col
  
  dir.create(file.path("Figures/Misclassified/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/Misclassified/",dataset,"/",dataset,"_",level_name,"_DatabaseEffect.png"), width=9, height=12, units = 'cm', dpi=1000)
  
}
make_all_database_effect <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      databaseEffect(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
}

# Plot an ECDF of unknown percentage by sample
unknownECDF <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  df <- data.frame(matrix(ncol = ncol(data[[1]]) - 1, nrow = 0))
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    df[nrow(df) + 1,] <- as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])
  }
  colnames(df) <- colnames(data[[1]])[-1]
  df <- 100 - df
  df$Methods <- c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "GTDBTk", "metaxa2")
  
  df <- melt(df, id.vars = c("Methods"))
  
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(df$Methods)])
  
  ggplot(df, aes(x=value, color=factor(Methods, levels=c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "GTDBTk", "metaxa2")))) + 
    stat_ecdf(geom = "step", size = 1) +
    theme_classic() + ggtitle("Empirical cumulative density over all samples") +
    ylab("Cumulative density") + xlab(paste0("Percent known at the ", level_name, " level")) +
    col +
    guides(color=guide_legend(title="Method"))
  dir.create(file.path("Figures/unknownECDF/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0(path_to_data, "Figures/unknownECDF/",dataset,"/",dataset,"_",level_name,"_ECDFUnknowns.png"), width=12, height=10, units = 'cm', dpi=1000)
}
make_all_unknownECDF_plots <- function() {
  #datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  datasets <- c("All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      unknownECDF(dataset, level)
    }
  }
}

# Calculate Bray Curtis between methods, abundances normalized
brayCurtisNormalized <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
  }
  
  pairs <- combn(1:length(data), 2)
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    distances[is.na(distances)] <- 1
    output <- vector()
    for (colname in unique(gsub(pattern="_.$", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
  }
  
  colnames(pairwise_diffs) <- c("Method 1", "Method 2", rownames(data[[1]]))
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MPA2",
                                          pairwise_diffs[,1]==2 ~ "MPA3",
                                          pairwise_diffs[,1]==3 ~ "MPA4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "Phylophlan",
                                          pairwise_diffs[,1]==6 ~ "Kraken",
                                          pairwise_diffs[,1]==7 ~ "GTDBTk",
                                          pairwise_diffs[,1]==8 ~ "metaxa2"),
                                "/",
                                case_when(pairwise_diffs[,2]==1 ~ "MPA2",
                                          pairwise_diffs[,2]==2 ~ "MPA3",
                                          pairwise_diffs[,2]==3 ~ "MPA4",
                                          pairwise_diffs[,2]==4 ~ "mOTUs3",
                                          pairwise_diffs[,2]==5 ~ "Phylophlan",
                                          pairwise_diffs[,2]==6 ~ "Kraken",
                                          pairwise_diffs[,2]==7 ~ "GTDBTk",
                                          pairwise_diffs[,2]==8 ~ "metaxa2"))
  
  pairwise_diffs[,1:2] <- NULL
  
  bray_df <- melt(pairwise_diffs, id.vars = "Pair")
  bray_df$Method <- gsub("\\/.*","",bray_df$Pair)
  
  # Reformat Bray Curtis matrix
  
  # Swap study for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  bray_df$dataset <- gsub("_.*", "", bray_df$variable)
  bray_df$assembly <- ifelse(grepl("Phylophlan|GTDBTk", bray_df$Pair), "assembly", "reference")
  
  #colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  colAdd <- c("#38761d")
  colAdd <- colAdd[1:length(unique(gsub("_.*", "", rownames(data[[1]]))))]
  names(colAdd) <- unique(gsub("_.*", "", rownames(data[[1]])))
  col = scale_color_manual(values = colAdd)
  
  # Plot
  ggplot(data = bray_df, aes(x=Pair, y=value)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 1, aes(color=dataset)) + 
    col +
    ggtitle(study) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab(paste0("Bray-Curtis Dissimilarity at the ", level_name, " level (Normalized)")) + 
    ylim(0, 1) + 
    facet_grid(~assembly,scales='free')
  
  dir.create(file.path(path_to_data, "Figures/BrayCurtisNormalizedBoxplots/",dataset,"/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/BrayCurtisNormalizedBoxplots/",dataset,"/",dataset,"_",level_name,"_BrayCurtisNormalized.png"), width=20, height=16, units = 'cm', dpi=1000)
}
make_all_brayCurtisNormalizedBoxplots <- function() {
  #datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  datasets <- c("All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      brayCurtisNormalized(dataset, level)
    }
  }
}

# Calculate Bray Curtis between methods, let unclassified be a category
brayCurtis <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
  }
  
  pairs <- combn(1:length(data), 2)
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    output <- vector()
    for (colname in unique(gsub(pattern="_.*", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
  }
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MPA2",
                                          pairwise_diffs[,1]==2 ~ "MPA3",
                                          pairwise_diffs[,1]==3 ~ "MPA4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "Phylophlan",
                                          pairwise_diffs[,1]==6 ~ "Kraken",
                                          pairwise_diffs[,1]==7 ~ "metaxa2"),
                                "/",
                                case_when(pairwise_diffs[,2]==1 ~ "MPA2",
                                          pairwise_diffs[,2]==2 ~ "MPA3",
                                          pairwise_diffs[,2]==3 ~ "MPA4",
                                          pairwise_diffs[,2]==4 ~ "mOTUs3",
                                          pairwise_diffs[,2]==5 ~ "Phylophlan",
                                          pairwise_diffs[,2]==6 ~ "Kraken",
                                          pairwise_diffs[,2]==7 ~ "metaxa2"))
  
  pairwise_diffs[,1:2] <- NULL
  
  bray_df <- melt(pairwise_diffs, id.vars = "Pair")
  bray_df$Method <- gsub("\\/.*","",bray_df$Pair)
  
  # Reformat Bray Curtis matrix
  
  # Swap study for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  # Plot
  ggplot(data = bray_df, aes(x=Pair, y=value)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.1) + 
    ggtitle(study) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab(paste0("Bray-Curtis Dissimilarity at the ", level_name, " level with UNCLASSIFIED as a group")) + 
    ylim(0, 1)
  
  dir.create(file.path("Figures/BrayCurtisBoxplots/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/BrayCurtisBoxplots/",dataset,"/",dataset,"_",level_name,"_BrayCurtis.png"), width=20, height=16, units = 'cm', dpi=1000)
}
make_all_brayCurtisBoxplots <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      brayCurtis(dataset, level)
    }
  }
}

brayCurtisShared <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  shared_taxa <- data[[1]]$TaxIDs
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
    shared_taxa <- intersect(shared_taxa, colnames(data[[i]]))
  }
  
  for (i in 1:length(data)) {
    row_names <- rownames(data[[i]])
    data[[i]] <- as.data.frame(data[[i]][,shared_taxa])
    rownames(data[[i]]) <- row_names
    colnames(data[[i]]) <- shared_taxa
  }
  
  pairs <- combn(1:length(data), 2)
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    output <- vector()
    for (colname in unique(gsub(pattern="_.*", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
  }
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MPA2",
                                          pairwise_diffs[,1]==2 ~ "MPA3",
                                          pairwise_diffs[,1]==3 ~ "MPA4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "Phylophlan",
                                          pairwise_diffs[,1]==6 ~ "Kraken",
                                          pairwise_diffs[,1]==7 ~ "metaxa2"),
                                "/",
                                case_when(pairwise_diffs[,2]==1 ~ "MPA2",
                                          pairwise_diffs[,2]==2 ~ "MPA3",
                                          pairwise_diffs[,2]==3 ~ "MPA4",
                                          pairwise_diffs[,2]==4 ~ "mOTUs3",
                                          pairwise_diffs[,2]==5 ~ "Phylophlan",
                                          pairwise_diffs[,2]==6 ~ "Kraken",
                                          pairwise_diffs[,2]==7 ~ "metaxa2"))
  
  pairwise_diffs[,1:2] <- NULL
  
  bray_df <- melt(pairwise_diffs, id.vars = "Pair")
  bray_df$Method <- gsub("\\/.*","",bray_df$Pair)
  
  # Reformat Bray Curtis matrix
  
  # Swap study for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  # Plot
  ggplot(data = bray_df, aes(x=Pair, y=value)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.1) + 
    ggtitle(paste0(study, " (", length(shared_taxa), " common taxa)")) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab(paste0("Bray Curtis at the ", level_name, " level (UNCLASSIFIED removed)")) + 
    ylim(0, 1)
  
  dir.create(file.path("Figures/brayCurtisShared/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/brayCurtisShared/",dataset,"/",dataset,"_",level_name,"_brayCurtisShared.png"), width=20, height=16, units = 'cm', dpi=1000)
}
make_all_brayCurtisSharedBoxplots <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      brayCurtisShared(dataset, level)
    }
  }
}

# Plot PCoA from Bray Curtis dissimilarity on normalized abundances
PCoABray <- function(dataset, level) {
  set.seed(1)
  # Read in metadata
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
  }
  
  dat_taxa_species <- bind_rows(data[[1]], data[[2]])
  for (i in 3:length(data)) {
    dat_taxa_species <- bind_rows(dat_taxa_species, data[[i]])
  }
  dat_taxa_species[is.na(dat_taxa_species)] <- 0
  distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
  distances[is.na(distances)] <- 1
  
  pc <- capscale(distances~1, comm = distances, na.action = na.omit)
  cap = data.frame(pc$CA$u)
  cap$Method <- case_when(grepl("_1", rownames(cap)) ~ "MPA2",
                          grepl("_2", rownames(cap)) ~ "MPA3",
                          grepl("_3", rownames(cap)) ~ "MPA4",
                          grepl("_4", rownames(cap)) ~ "mOTUs3",
                          grepl("_5", rownames(cap)) ~ "Phylophlan",
                          grepl("_6", rownames(cap)) ~ "Kraken",
                          grepl("_7", rownames(cap)) ~ "metaxa2")
  cap$Sample <- gsub("_.*", "", rownames(cap))
  s = summary(pc)
  
  # Swap study name for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  study <- paste0(study, ", ", level_name)
  
  colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:11])
  
  # Plot ordination
  ggplot(cap, aes(MDS1, MDS2)) + 
    geom_point(aes(color=factor(Method, levels=c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "metaxa2"))), size = 1, alpha = 0.7) + 
    theme_classic() +coord_fixed() + 
    ggtitle(study) +
    labs( color = "Tool") + 
    labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
         y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = '')) +
    col + guides(color=guide_legend(title="Method"))
  dir.create(file.path("Figures/PCoABray/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/PCoABray/",dataset,"/",dataset,"_",level_name,"_PCoA.png"), width=14, height=12, units = 'cm', dpi=1000)
}

PCoAJaccard <- function(dataset, level) {
  set.seed(1)
  # Read in metadata
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
  }
  
  dat_taxa_species <- bind_rows(data[[1]], data[[2]])
  for (i in 3:length(data)) {
    dat_taxa_species <- bind_rows(dat_taxa_species, data[[i]])
  }
  dat_taxa_species[is.na(dat_taxa_species)] <- 0
  distances <- vegdist(dat_taxa_species, method = "jaccard") %>% as.matrix(labels=TRUE)
  distances[is.na(distances)] <- 1
  
  pc <- capscale(distances~1, comm = distances, na.action = na.omit)
  cap = data.frame(pc$CA$u)
  cap$Method <- case_when(grepl("_1", rownames(cap)) ~ "MPA2",
                          grepl("_2", rownames(cap)) ~ "MPA3",
                          grepl("_3", rownames(cap)) ~ "MPA4",
                          grepl("_4", rownames(cap)) ~ "mOTUs3",
                          grepl("_5", rownames(cap)) ~ "Phylophlan",
                          grepl("_6", rownames(cap)) ~ "Kraken",
                          grepl("_7", rownames(cap)) ~ "metaxa2")
  cap$Sample <- gsub("_.*", "", rownames(cap))
  s = summary(pc)
  
  # Swap study name for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  study <- paste0(study, ", ", level_name)
  
  colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:11])
  
  # Plot ordination
  ggplot(cap, aes(MDS1, MDS2)) + 
    geom_point(aes(color=factor(Method, levels=c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "metaxa2"))), size = 1, alpha = 0.7) + 
    theme_classic() +coord_fixed() + 
    ggtitle(study) +
    labs( color = "Tool") + 
    labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
         y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = '')) +
    col + guides(color=guide_legend(title="Method"))
  dir.create(file.path("Figures/PCoAJaccard/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/PCoAJaccard/",dataset,"/",dataset,"_",level_name,"_PCoA.png"), width=14, height=12, units = 'cm', dpi=1000)
}
make_all_PCoAs <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      PCoABray(dataset, level)
      PCoAJaccard(dataset, level)
    }
  }
}

# Plot overlaps vs read count
overlapVsReads <- function(dataset, level) {
  data <- preprocess(dataset, level)
  reads <- preprocessReads(dataset)
  
  taxa_names <- list()
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    for (j in 2:ncol(tmp)) {
      tmp[,j] <- ifelse(tmp[,j]!=0, tmp[,1], NA_character_)
    }
    tmp <- tmp[!tmp$TaxIDs=="UNCLASSIFIED",-1]
    taxa_names <- append(taxa_names, list(tmp))
  }
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  all_dbs <- read.csv(paste0("Databases/DB_merged/", level_name, ".csv"), sep = ",")
  dbs <- all_dbs[,1]
  for (i in 2:length(ncol(all_dbs))) {
    dbs <- intersect(dbs, all_dbs[,i])
  }
  
  prop_overlap <- data.frame(matrix(nrow = 0, ncol = 2))
  
  for (i in 1:ncol(taxa_names[[1]])) {
    MPA2 <- setdiff(taxa_names[[1]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA2 <- MPA2[!is.na(MPA2)]
    MPA3 <- setdiff(taxa_names[[2]][,i], taxa_names[[1]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA3 <- MPA3[!is.na(MPA3)]
    MPA4 <- setdiff(taxa_names[[3]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[1]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA4 <- MPA4[!is.na(MPA4)]
    metaxa <- setdiff(taxa_names[[4]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[1]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    metaxa <- metaxa[!is.na(metaxa)]
    motus <- setdiff(taxa_names[[5]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[1]][,i]) %>% setdiff(taxa_names[[6]][,i])
    motus <- MPA2[!is.na(motus)]
    phylophlan <- setdiff(taxa_names[[6]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[1]][,i])
    phylophlan <- phylophlan[!is.na(phylophlan)]
    
    way2 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 2)
    way2 <- way2[!is.na(way2)]
    way3 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 3)
    way3 <- way3[!is.na(way3)]
    way4 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 4)
    way4 <- way4[!is.na(way4)]
    way5 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 5)
    way5 <- way5[!is.na(way5)]
    all_overlap <- Reduce(intersect, list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                          taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]))
    all_overlap <- all_overlap[!is.na(all_overlap)]
    
    all_not_overlap <- c(MPA2, MPA3, MPA4, metaxa, motus, phylophlan, way2, way3, way4, way5)
    
    prop_overlap <- rbind(prop_overlap, c(colnames(taxa_names[[1]])[i], 
                                          length(all_overlap)/
                             (length(all_overlap) + length(intersect(dbs, all_not_overlap)))))
  }
  
  colnames(prop_overlap) <- c("Sample", "Proportion")
  prop_overlap$Proportion <- as.numeric(prop_overlap$Proportion)
  
  df <- left_join(prop_overlap, reads, by=c("Sample"))
  
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  ggplot(df, aes(x=Reads, y=Proportion)) + geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() + ggtitle("Reads vs Overlap Proportion") +
    ylab("Proportion of taxa assigned by any tool \n(that were in all the databases)\nthat were assigned by all of the tools") + xlab("Reads")
  dir.create(file.path("Figures/readsVsOverlap/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/readsVsOverlap/",dataset,"/",dataset,"_",level_name,"_readsVsOverlap.png"), width=12, height=10, units = 'cm', dpi=1000)
}
make_all_overlapVsReads <- function() {
  datasets <- c("All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      overlapVsReads(dataset, level)
    }
  }
}

# Plot taxa called per method
taxaCalledSlopePlot <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  taxa_names <- list()
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    for (j in 2:ncol(tmp)) {
      tmp[,j] <- ifelse(tmp[,j]!=0, tmp[,1], NA_character_)
    }
    tmp <- tmp[!tmp$TaxIDs=="UNCLASSIFIED",-1]
    taxa_names <- append(taxa_names, list(tmp))
  }
  
  whole_dataset <- data.frame(matrix(nrow = 7, ncol = ncol(taxa_names[[1]])))
  colnames(whole_dataset) <- colnames(taxa_names[[1]])
  
  
  for (i in 1:nrow(whole_dataset)) {
    whole_dataset[i,] <- colSums(!is.na(taxa_names[[i]]))
  }
  
  whole_dataset$Method <- c("MPA2", "MPA3", "MPA4", "mOTUs3", "Phylophlan", "Kraken", "metaxa2")
  
  df <- melt(whole_dataset, id.vars = "Method")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  # Change study name for plot
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  # Plot the overlaps
  ggplot(df, aes(y=value, x=Method, group=variable)) + 
    geom_line(aes(alpha=0.5, color=variable)) + 
    geom_point(aes(alpha=0.5, color=variable)) + 
    theme_classic() +
    theme(legend.position="none") + 
    ggtitle(study) +
    xlab("Method") +
    ylab(paste0("Number of ", level_name, " assigned"))
  dir.create(file.path("Figures/TaxaCalls/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/TaxaCalls/",dataset,"/",dataset,"_",level_name,"_TaxaCalls.png"), width=10, height=5)
}
make_all_taxaCalled_plots <- function() {
  datasets <- c("All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      taxaCalledSlopePlot(dataset, level)
    }
  }
}

# Plot reads vs bins
ReadsVsBins <- function(dataset) {
  df <- read.csv(paste(dataset, "/reads_vs_bins.tsv",sep = ""), sep="\t", header = T)
  ggplot(df, aes(x=reads, y=bins)) + 
    geom_point() + theme_classic() + xlab("Reads") + 
    ylab("Bins")
  ggsave(paste0("Figures/ReadsVsBins/",dataset,"_ReadsVsBins.png"), width=9, height=8, units = 'cm', dpi=1000)
}
make_all_ReadsVsBins <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  
  for (dataset in datasets) {
    if (dataset == "All") {
      tmp_datasets <- datasets[datasets!=dataset]
      
      df <- read.csv(paste(tmp_datasets[1], "/reads_vs_bins.tsv",sep = ""), sep="\t", header = T)
      for (i in 2:length(tmp_datasets)) {
        df <- rbind(df, 
                    read.csv(paste(tmp_datasets[i], "/reads_vs_bins.tsv",sep = ""), sep="\t", header = T))
      }
      ggplot(df, aes(x=reads, y=bins)) + 
        geom_point() + theme_classic() + xlab("Reads") + 
        geom_smooth(method="lm") +
        ylab("Bins")
      ggsave(paste0("Figures/ReadsVsBins/",dataset,"_ReadsVsBins.png"), width=9, height=8, units = 'cm', dpi=1000)
      
    }
    else {
    }
  }
}

# Calculate Bray Curtis differences against reads
brayNormVsReads <- function(dataset, level) {
  data <- preprocess(dataset, level)
  reads <- preprocessReads(dataset)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
  }
  
  pairs <- combn(1:length(data), 2)
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    output <- vector()
    for (colname in unique(gsub(pattern="_.*", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
  }
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MPA2",
                                          pairwise_diffs[,1]==2 ~ "MPA3",
                                          pairwise_diffs[,1]==3 ~ "MPA4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "Phylophlan",
                                          pairwise_diffs[,1]==6 ~ "Kraken",
                                          pairwise_diffs[,1]==7 ~ "metaxa2"),
                                "/",
                                case_when(pairwise_diffs[,2]==1 ~ "MPA2",
                                          pairwise_diffs[,2]==2 ~ "MPA3",
                                          pairwise_diffs[,2]==3 ~ "MPA4",
                                          pairwise_diffs[,2]==4 ~ "mOTUs3",
                                          pairwise_diffs[,2]==5 ~ "Phylophlan",
                                          pairwise_diffs[,2]==6 ~ "Kraken",
                                          pairwise_diffs[,2]==7 ~ "metaxa2"))
  
  pairwise_diffs[,1:2] <- NULL
  
  colnames(pairwise_diffs) <- c(gsub(pattern="_.*", "", rownames(data[[1]])), "Pair")
  rownames(pairwise_diffs) <- pairwise_diffs$Pair
  pairwise_diffs$Pair <- NULL
  
  pairwise_diffs <- data.frame(t(pairwise_diffs))
  
  pairwise_diffs$Sample <- rownames(pairwise_diffs)
  df <- left_join(pairwise_diffs, reads, by=c("Sample"))
  colnames(df) <- gsub("\\.", "\\/", colnames(df))
  
  bray_df <- melt(df, id.vars = c("Sample", "Reads"))
  colnames(bray_df) <- c("Sample", "Reads", "Pair", "value")

  # Swap study for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  # Plot
  ggplot(data = bray_df, aes(x=Reads, y=value)) +
    geom_point() + geom_smooth(method='lm', formula= y~x, color="black") +
    ggtitle(study) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab(paste0("Bray-Curtis Dissimilarity for ", level_name," level (normalized)")) + 
    ylim(0, 1)
  
  dir.create(file.path("Figures/BrayNormVsReads/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/BrayNormVsReads/",dataset,"/",dataset,"_",level_name,"_BrayCurtisNormalized.png"), width=20, height=12, units = 'cm', dpi=1000)
}
make_all_brayNormVsReads <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      brayNormVsReads(dataset, level)
    }
  }
}

brayVsReads <- function(dataset, level) {
  data <- preprocess(dataset, level)
  reads <- preprocessReads(dataset)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
  }
  
  pairs <- combn(1:length(data), 2)
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    output <- vector()
    for (colname in unique(gsub(pattern="_.*", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
  }
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MPA2",
                                          pairwise_diffs[,1]==2 ~ "MPA3",
                                          pairwise_diffs[,1]==3 ~ "MPA4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "Phylophlan",
                                          pairwise_diffs[,1]==6 ~ "Kraken",
                                          pairwise_diffs[,1]==7 ~ "metaxa2"),
                                "/",
                                case_when(pairwise_diffs[,2]==1 ~ "MPA2",
                                          pairwise_diffs[,2]==2 ~ "MPA3",
                                          pairwise_diffs[,2]==3 ~ "MPA4",
                                          pairwise_diffs[,2]==4 ~ "mOTUs3",
                                          pairwise_diffs[,2]==5 ~ "Phylophlan",
                                          pairwise_diffs[,2]==6 ~ "Kraken",
                                          pairwise_diffs[,2]==7 ~ "metaxa2"))
  
  pairwise_diffs[,1:2] <- NULL
  
  colnames(pairwise_diffs) <- c(gsub(pattern="_.*", "", rownames(data[[1]])), "Pair")
  rownames(pairwise_diffs) <- pairwise_diffs$Pair
  pairwise_diffs$Pair <- NULL
  
  pairwise_diffs <- data.frame(t(pairwise_diffs))
  
  pairwise_diffs$Sample <- rownames(pairwise_diffs)
  df <- left_join(pairwise_diffs, reads, by=c("Sample"))
  colnames(df) <- gsub("\\.", "\\/", colnames(df))
  
  bray_df <- melt(df, id.vars = c("Sample", "Reads"))
  colnames(bray_df) <- c("Sample", "Reads", "Pair", "value")
  
  # Swap study for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  # Plot
  ggplot(data = bray_df, aes(x=Reads, y=value)) +
    geom_point() + geom_smooth(method='lm', formula= y~x, color="black") +
    ggtitle(study) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab(paste0("Bray-Curtis Dissimilarity for ", level_name," level with UNCLASSIFIED")) + 
    ylim(0, 1)
  
  dir.create(file.path("Figures/BrayVsReads/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/BrayVsReads/",dataset,"/",dataset,"_",level_name,"_BrayCurtis.png"), width=20, height=12, units = 'cm', dpi=1000)
}
make_all_brayVsReads <- function() {
  datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      brayVsReads(dataset, level)
    }
  }
}

# Compare PERMANOVA results among tools
PERMANOVAAcrossTools <- function(dataset, level, cutoff) {
  # Import metadata
  meta = read.csv("Metadata/all_metadata.csv")
  
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    taxa_names <- tmp$TaxIDs
    tmp <- data.frame(t(tmp[,-1]))
    colnames(tmp) <- as.character(taxa_names)
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- tmp
  }
  
  meta <- meta[meta$sample_accession_WGS %in% gsub("_.*", "", rownames(data[[1]])),]
  
  all_data <- data[[1]]
  all_data$Sample <- rownames(all_data)
  for (i in 2:length(data)) {
    data[[i]]$Sample <- rownames(data[[i]])
    all_data <- merge(all_data,data[[i]], all = T)
  }
  
  all_data[is.na(all_data)] <- 0
  
  # Keep only high abundance taxa
  all_barcode = all_data > 0.0001
  common_index = apply(all_barcode,2,mean) > cutoff

  # Fill in abundance table and apply abundance threshold
  dat_taxa_species_common = all_data[, common_index]
  rownames(dat_taxa_species_common) <- dat_taxa_species_common$Sample
  dat_taxa_species_common <- dat_taxa_species_common[order(rownames(dat_taxa_species_common)),]
  dat_taxa_species_common$Sample <- NULL
  
  # Calculate Bray dissimilarity
  bray <- vegdist(dat_taxa_species_common, method="bray")
  
  stat_meta <- data.frame(matrix(nrow = 0, ncol = nrow(data[[1]])))
  for (i in 1:length(data)) {
    tmp <- meta
    tmp$sample_accession_WGS <- paste0(tmp$sample_accession_WGS, "_", i)
    tmp$method <- case_when(i==1 ~ "MPA2",
                            i==2 ~ "MPA3",
                            i==3 ~ "MPA4",
                            i==4 ~ "mOTUs3",
                            i==5 ~ "Phylophlan",
                            i==6 ~ "Kraken",
                            i==7 ~ "metaxa2")
    stat_meta <- rbind(stat_meta, tmp)
  }
  
  # Keep only fully known metadata
  na = c("dataset_name", "PMID", "location_accession", "sample_accession", "sample_accession_WGS", "database", "study_accession_db", "location_accession_db", "sample_accession_db", "sample_type", "date", "latitude", "longitude", "water_temp", "salinity", "oxygen", "phosphate", "nitrate.nitrite", "silicon", "sulfate", "redox", "conductivity", "carbon", "prior_precip", "pH", "iron", "freeze", "ID", "version")
  meta_wo_na = stat_meta[, -which(names(stat_meta) %in% na)]
  rownames(meta_wo_na) <- stat_meta$sample_accession_WGS
  
  meta_wo_na <- meta_wo_na[rownames(meta_wo_na) %in% rownames(dat_taxa_species_common), ]
  
  bray[is.na(bray)] <- 1
  
  ### Statistics Beta Diversity Taxonomy  
  #### Univariate on complete cases
  adonis_res_rsq = vector()
  adonis_res_pval = vector()
  for (col in names(meta_wo_na)){
    adonis.univ = adonis(as.formula(paste("bray ~ ", col)), data = meta_wo_na)
    adonis_res_rsq[col] = adonis.univ$aov.tab[1,]$R2
    adonis_res_pval[col] = adonis.univ$aov.tab[1,]$`Pr(>F)`
  }
  
  # Format adonis output for plotting
  univar_res_tax = rbind(adonis_res_pval, adonis_res_rsq)
  univar_res_tax = as.data.frame(t(univar_res_tax))
  names(univar_res_tax) = c("P-Value", "R2")
  univar_res_tax$`P-Value` = as.numeric(univar_res_tax$`P-Value`)
  univar_res_tax$R2 = as.numeric(univar_res_tax$R2)
  univar_res_tax$p_adj = p.adjust(univar_res_tax$`P-Value`, "fdr")
  univar_res_tax = subset(univar_res_tax, row.names(univar_res_tax) != "batch")
  univar_res_tax = subset(univar_res_tax, row.names(univar_res_tax) != "hla_b27")
  univar_res_tax$clinical <- c("Elevation", "Prior high", "Prior low", "Annual precipitation", "Method")
  rownames(univar_res_tax) <- c("Elevation", "Prior high", "Prior low", "Annual precipitation", "Method")
  univar_res_tax$stars = cut(univar_res_tax$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", ".", ""))
  
  univar_res_tax$p_adj = round(univar_res_tax$p_adj, 3)
  univar_res_tax$R2 = univar_res_tax$R2 *100
  dodge = position_dodge(width = 0.8)
  
  # Swap study name for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  # Plot
  ggplot(data = univar_res_tax, aes(reorder(clinical, R2), y = R2, label = stars)) + 
    geom_bar(stat = "identity", position = dodge, fill="brown") + 
    geom_text(position = dodge, vjust = 0.8, hjust = -0.1, size = 4) + 
    theme_classic(base_size = 10) + ylab("Univarate R-squared") + 
    coord_flip() + ylim(0, 25) + xlab("") + 
    labs(fill = "Tool") + 
    ggtitle(study) +
    guides(fill = guide_legend(reverse=T))
  
  dir.create(file.path("Figures/PERMANOVA/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("Figures/PERMANOVA/",dataset,"/",dataset,"_",level_name,"_MethodPERMANOVA.png"), width=9, height=9, units = 'cm', dpi=1000)
}
make_all_PERMANOVA <- function() {
  datasets <- c("Coastal", "Mine", "Gators", "TaraPolar", "Trees")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      tryCatch({
        PERMANOVAAcrossTools(dataset, level, 0.001)
      }, error=function(e){})
      
    }
  }
}

# Intersection over union plots 
intersection_over_union <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    data[[i]] <- data[[i]][data[[i]]$TaxIDs!="UNCLASSIFIED",]
    data[[i]][,-1][data[[i]][,-1]<0.00001]=0
    tmp <- data[[i]]$TaxIDs
    data[[i]] <- data.frame(lapply(data[[i]], function(x){ifelse(x==0,NA, tmp)}))
    data[[i]]$TaxIDs <- NULL
  }
  
  pairs <- combn(1:length(data), 2)
  
  iou <- function(x, y) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    output <- length(intersect(x,y))/length(unique(c(x,y)))
    if (is.na(output)) {output <- 0}
    return(output)
  }
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + ncol(data[[1]])))
  for (i in 1:ncol(pairs)) {
    output <- mapply(iou, data[[pairs[1,i]]], data[[pairs[2,i]]])
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
  }
  
  colnames(pairwise_diffs) <- c("Method 1", "Method 2", colnames(data[[1]]))
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MPA2",
                                          pairwise_diffs[,1]==2 ~ "MPA3",
                                          pairwise_diffs[,1]==3 ~ "MPA4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "Phylophlan",
                                          pairwise_diffs[,1]==6 ~ "Kraken",
                                          pairwise_diffs[,1]==7 ~ "GTDBTk",
                                          pairwise_diffs[,1]==8 ~ "metaxa2"),
                                "/",
                                case_when(pairwise_diffs[,2]==1 ~ "MPA2",
                                          pairwise_diffs[,2]==2 ~ "MPA3",
                                          pairwise_diffs[,2]==3 ~ "MPA4",
                                          pairwise_diffs[,2]==4 ~ "mOTUs3",
                                          pairwise_diffs[,2]==5 ~ "Phylophlan",
                                          pairwise_diffs[,2]==6 ~ "Kraken",
                                          pairwise_diffs[,2]==7 ~ "GTDBTk",
                                          pairwise_diffs[,2]==8 ~ "metaxa2"))
  
  pairwise_diffs[,1:2] <- NULL
  
  iou_df <- melt(pairwise_diffs, id.vars = "Pair")
  iou_df$Method <- gsub("\\/.*","",iou_df$Pair)
  
  # Reformat Bray Curtis matrix
  
  # Swap study for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  iou_df$dataset <- gsub("_.*", "", iou_df$variable)
  
  colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  colAdd <- colAdd[1:length(unique(gsub("_.*", "", colnames(data[[1]]))))]
  names(colAdd) <- unique(gsub("_.*", "", colnames(data[[1]])))
  col = scale_color_manual(values = colAdd)
  
  # Plot
  ggplot(data = iou_df, aes(x=Pair, y=value)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 1, aes(color=dataset)) + 
    col +
    ggtitle(study) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab(paste0("Intersection Over Union at the ", level_name, " level")) + 
    ylim(0, 1)
  
  dir.create(file.path(path_to_data, "Figures/IOUBoxplots/",dataset,"/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/IOUBoxplots/",dataset,"/",dataset,"_",level_name,"_IOU.png"), width=20, height=16, units = 'cm', dpi=1000)
}
make_all_intersection_over_union_Boxplots <- function() {
  #datasets <- c("Coastal", "Gators", "Mine", "Trees", "TaraPolar", "All")
  datasets <- c("All")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      intersection_over_union(dataset, level)
    }
  }
}

# Figures for Herchel Smith
herchelSmith <- function() {
  dataset <- "All"
  level <- 5
  data <- preprocess(dataset, level)
  
  taxa_names <- list()
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    for (j in 2:ncol(tmp)) {
      tmp[,j] <- ifelse(tmp[,j]!=0, tmp[,1], NA_character_)
    }
    tmp <- tmp[!tmp$TaxIDs=="UNCLASSIFIED",-1]
    taxa_names <- append(taxa_names, list(tmp))
  }
  
  # 11 rows because 11 categories of overlap or nonoverlap
  whole_dataset <- data.frame(matrix(nrow = 11, ncol = ncol(taxa_names[[1]])))
  colnames(whole_dataset) <- colnames(taxa_names[[1]])
  
  # Proportion that could be called by all databases
  all_potential <- as.list(vector(length = 11))
  
  for (i in 1:ncol(taxa_names[[1]])) {
    MPA2 <- setdiff(taxa_names[[1]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA2 <- MPA2[!is.na(MPA2)]
    MPA3 <- setdiff(taxa_names[[2]][,i], taxa_names[[1]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA3 <- MPA3[!is.na(MPA3)]
    MPA4 <- setdiff(taxa_names[[3]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[1]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    MPA4 <- MPA4[!is.na(MPA4)]
    metaxa <- setdiff(taxa_names[[4]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[1]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[6]][,i])
    metaxa <- metaxa[!is.na(metaxa)]
    motus <- setdiff(taxa_names[[5]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[1]][,i]) %>% setdiff(taxa_names[[6]][,i])
    motus <- MPA2[!is.na(motus)]
    phylophlan <- setdiff(taxa_names[[6]][,i], taxa_names[[2]][,i]) %>% setdiff(taxa_names[[3]][,i]) %>% 
      setdiff(taxa_names[[4]][,i]) %>% setdiff(taxa_names[[5]][,i]) %>% setdiff(taxa_names[[1]][,i])
    phylophlan <- phylophlan[!is.na(phylophlan)]
    
    way2 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 2)
    way2 <- way2[!is.na(way2)]
    way3 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 3)
    way3 <- way3[!is.na(way3)]
    way4 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 4)
    way4 <- way4[!is.na(way4)]
    way5 <- get_k_way_overlaps(list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                    taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]), 5)
    way5 <- way5[!is.na(way5)]
    all_overlap <- Reduce(intersect, list(taxa_names[[1]][,i], taxa_names[[2]][,i], taxa_names[[3]][,i],
                                          taxa_names[[4]][,i], taxa_names[[5]][,i], taxa_names[[6]][,i]))
    all_overlap <- all_overlap[!is.na(all_overlap)]
    
    # Need to return proportion of all database occurrences
    whole_dataset[,i] <- c(length(MPA2), length(MPA3), length(MPA4), length(metaxa), length(motus), length(phylophlan),
                           length(way2), length(way3), length(way4), length(way5), length(all_overlap))
    
    all_potential[[1]] <- c(all_potential[[1]], MPA2)
    all_potential[[2]] <- c(all_potential[[2]], MPA3)
    all_potential[[3]] <- c(all_potential[[3]], MPA4)
    all_potential[[4]] <- c(all_potential[[4]], metaxa)
    all_potential[[5]] <- c(all_potential[[5]], motus)
    all_potential[[6]] <- c(all_potential[[6]], phylophlan)
    all_potential[[7]] <- c(all_potential[[7]], way2)
    all_potential[[8]] <- c(all_potential[[8]], way3)
    all_potential[[9]] <- c(all_potential[[9]], way4)
    all_potential[[10]] <- c(all_potential[[10]], way5)
    all_potential[[11]] <- c(all_potential[[11]], all_overlap)
  }
  
  for (i in 1:length(all_potential)) {
    all_potential[[i]] <- all_potential[[i]][all_potential[[i]]!="FALSE"]
  }
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  all_dbs <- read.csv(paste0("Databases/DB_merged/", level_name, ".csv"), sep = ",")
  tmp <- all_dbs[,1]
  for (i in 2:length(ncol(all_dbs))) {
    tmp <- intersect(tmp, all_dbs[,i])
  }
  
  prop_in_all <- vector(length = length(all_potential))
  for (i in 1:length(all_potential)) {
    prop_in_all[i] <- 1 - length(all_potential[[i]][all_potential[[i]] %in% setdiff(all_potential[[i]], tmp)]) / length(all_potential[[i]])
  }
  prop_in_all <- ifelse(is.na(prop_in_all), 0, prop_in_all)
  
  chyper_db <- read.csv(paste0("Databases/DB_merged/for_chyper.csv"), sep = ",")
  
  for_chyper <- chyper_db[,(level + 1)]
  
  # Update to use pchyper
  ns <- rep(list(for_chyper[-1]), ncol(taxa_names[[1]]))
  s <- rep(list(for_chyper[1]), ncol(taxa_names[[1]]))
  
  tmp <- t(data.frame(lapply(taxa_names, function (method) colSums(!is.na(method)))))
  ms <- as.list(vector(length = ncol(tmp)))
  for (i in 1:ncol(tmp)) {
    ms[[i]] <- unname(tmp[,i])
  }
  
  overlapall <- as.list(whole_dataset[nrow(whole_dataset),])
  
  ps <- mapply(pvalchyper, overlapall, s, ns, ms, "upper", verbose=F)
  ps <- signif(ps, 1)
  
  orders <- (as.numeric(whole_dataset[7,]) + 0.00001)/(colSums(whole_dataset) + 0.00001)
  
  # Make dataframe for plotting
  whole_dataset$labels <- c("Unique to MPA2", "Unique to MPA3", "Unique to MPA4", 
                            "Unique to metaxa", "Unique to mOTUs3", "Unique to Phylophlan",
                            "2-way overlaps", "3-way overlaps", "4-way overlaps", "5-way overlaps",
                            "All overlapping")
  
  df <- melt(whole_dataset, id.vars = "labels")
  
  
  
  # Sort and add p-values
  tmp <- data.frame(names=names(orders))
  tmp$orders <- unname(orders)
  tmp$ps <- ps
  
  df <- merge(df, tmp, by.x = "variable", by.y = "names")
  
  df <- df %>% arrange(-orders)
  df$ps <- ifelse(df$ps == 0, "<1e-16", df$ps)
  
  label <- paste("Proportion of",level_name)
  
  # Change study name for plot
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  colAdd <- c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_fill_manual(values = colAdd[1:11])
  
  # Plot the overlaps
  ggp1 <- ggplot(df, aes(fill=factor(labels, levels=c("All overlapping", "5-way overlaps", 
                                              "4-way overlaps", "3-way overlaps", 
                                              "2-way overlaps", "Unique to Phylophlan",
                                              "Unique to mOTUs3", "Unique to metaxa",
                                              "Unique to MPA4", "Unique to MPA3",
                                              "Unique to MPA2")), 
                 y=value, x=reorder(variable, -orders))) + 
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10)) + 
    xlab("Sample") +
    ylab(label) + guides(shape = guide_legend(override.aes = list(size = 0.5)),
                         fill = guide_legend(override.aes = list(size = 0.1))) + 
    scale_x_discrete(labels=df$ps[df$labels=="Unique to MPA2"]) + 
    col
  
  data <- preprocess(dataset, level)
  
  df <- data.frame(matrix(ncol = ncol(data[[1]]) - 1, nrow = 0))
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    df[nrow(df) + 1,] <- as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])
  }
  df$Methods <- c("MPA2", "MPA3", "MPA4", "metaxa2", "mOTUs3", "Phylophlan")
  
  df <- melt(df, id.vars = c("Methods"))
  
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  ggp2 <- ggplot(df, aes(x=value, color=Methods)) + stat_ecdf(geom = "step") +
    theme_classic() + ggtitle("Empirical cumulative density over all samples") +
    ylab("Cumulative density") + xlab("Percent unknown")
  
  
  df <- read.csv("HerchelSmith/Datasets.csv")
  colnames(df) <- c("Dataset", "Location", "Number of samples", "Sequencing range", "Mean Sequencing Depth")
  rownames(df) <- NULL
  ggp1 <- tableGrob(df, rows = rep("",5))
  
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    taxa_names <- tmp$TaxIDs
    tmp <- t(tmp[,-1])
    colnames(tmp) <- taxa_names
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    data[[i]] <- data.frame(tmp)
  }
  
  pairs <- combn(1:length(data), 2)
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    output <- vector()
    for (colname in unique(gsub(pattern="_.*", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
  }
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MPA2",
                                          pairwise_diffs[,1]==2 ~ "MPA3",
                                          pairwise_diffs[,1]==3 ~ "MPA4",
                                          pairwise_diffs[,1]==4 ~ "metaxa2",
                                          pairwise_diffs[,1]==5 ~ "mOTUs3",
                                          pairwise_diffs[,1]==6 ~ "Phylophlan"),
                                "/",
                                case_when(pairwise_diffs[,2]==1 ~ "MPA2",
                                          pairwise_diffs[,2]==2 ~ "MPA3",
                                          pairwise_diffs[,2]==3 ~ "MPA4",
                                          pairwise_diffs[,2]==4 ~ "metaxa2",
                                          pairwise_diffs[,2]==5 ~ "mOTUs3",
                                          pairwise_diffs[,2]==6 ~ "Phylophlan"))
  
  pairwise_diffs[,1:2] <- NULL
  
  bray_df <- melt(pairwise_diffs, id.vars = "Pair")
  bray_df$Method <- gsub("\\/.*","",bray_df$Pair)
  
  # Reformat Bray Curtis matrix
  
  # Swap study for plotting
  study <- case_when(dataset=="TaraPolar" ~ "Artic ocean",
                     dataset=="Mine" ~ "Acid mine runoff",
                     dataset=="Trees" ~ "North American forest soils",
                     dataset=="Coastal" ~ "Baltic Coastal Sediment",
                     dataset=="Gators" ~ "Gator nest soil",
                     dataset=="All" ~ "All Datasets")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")
  
  # Plot
  ggp3 <- ggplot(data = bray_df, aes(x=Pair, y=value)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.1) + 
    ggtitle(study) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("Bray-Curtis Dissimilarity") + 
    ylim(0, 1)
  
  figure <- ggarrange(
    ggarrange(ggp3, ggp2, ncol = 2, labels = c("A", "B")), 
    ggp1,
    nrow = 2, 
    heights = c(3,1),
    labels = c("", "C")
  ) 
  figure
  
  dir.create(file.path("Figures/HerchelSmith/"), showWarnings = FALSE)
  ggsave(plot = figure, filename = "Figures/HerchelSmith/HerchelSmith_Figure.png", width=24, height=20, units = 'cm', dpi=1000)
}

runAll <- function() {
  make_all_database_effect()
  print("Done Database Effect")
  make_all_unknownECDF_plots()
  print("Done ECDF")
  make_all_brayCurtisNormalizedBoxplots()
  make_all_brayCurtisBoxplots()
  make_all_brayCurtisSharedBoxplots()
  print("Done BC")
  make_all_PCoAs()
  print("Done PCoA")
  make_all_overlapVsReads()
  print("Done OverlapVsReads")
  make_all_taxaCalled_plots()
  print("Done TaxaCalled")
  make_all_ReadsVsBins()
  print("Done Reads vs bins")
  make_all_brayVsReads()
  make_all_brayNormVsReads()
  print("Done Bray vs reads")
  #make_all_PERMANOVA()
  #print("Done PERMANOVA")
}




