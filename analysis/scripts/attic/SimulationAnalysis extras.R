rm(list = ls())
library(grid)
library(ggplot2)
library(gridExtra)
library(ggplotify)
library(reshape2)
library(tidyr)
library(ggh4x)
library(drc)
library(forcats)

source("analysis/scripts/helpers.R")

preprocess_all_simulated()
preprocess_all_truths()

# Check the F1 for many tools at many abundance thresholds
# Outputs dataframe of (number of tools) x (number of thresholds)
# Using NCBI names only
test_abundance_thresholds <- function(thresholds) {
  level = 7 # Species level
  out_mat = matrix(nrow = length(thresholds), ncol = 11)
  k = 1
  for (threshold in thresholds) {
    f1s_avg <- matrix(nrow = length(list.files("analysis/simulation_outputs")), ncol = 11)
    j = 1
    
    # Compute F1s for each dataset
    for (dataset in list.files("analysis/simulation_outputs")) {
      taxa = remove_non_ncbi(preprocess(dataset, level, FALSE))
      truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
      taxa <- lapply(taxa, renormalize)
      f1s = vector(length = length(taxa))
      # Compute average F1 for each tool
      for (i in 1:length(taxa)) {
        taxa_names <- list_taxa_by_sample(taxa[[i]], threshold)
        truth_names <- list_taxa_by_sample(truth, 0)
        f1s[i] <- mean(calculate_f1(taxa_names, truth_names))
      }
      f1s_avg[j,] <- f1s
      j = j + 1
    }
    # Compute average F1 across the datasets for each tool
    out_mat[k,] <- colMeans(f1s_avg)
    k = k + 1
  }
  rownames(out_mat) <- thresholds
  colnames(out_mat) <- c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", 
                         "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
  return(out_mat)
}

# Save optimal thresholds; use 0.05% as threshold throughout
out_df <- data.frame(test_abundance_thresholds(c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)), check.rows = F)
out_df <- data.frame(cbind("Thresholds" = rownames(out_df), out_df), check.rows = F)
dir.create(paste0("analysis/figures/thresholding/"), showWarnings = FALSE, recursive = T)
write.table(out_df, "analysis/figures/thresholding/thresholds.tsv", sep="\t", row.names = F)

threshold = 0.05
tool_names = c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
tool_core = c("Kraken 2 / Bracken 2", "MetaPhlAn 4", "mOTUs 3", "GTDB-Tk MEGAHIT", "PhyloPhlAn MEGAHIT")
colAdd <- c("#CF9FFF", "#A15BE4", "#00FFFF", "#0061FE", "#1434A4", "#81C784", "#2E7D32", "#FF7300", "#FBB15B", "#AC0911", "#E95420")
col = scale_color_manual(values = colAdd)
fil = scale_fill_manual(values = colAdd)

# Create an example of the phyla present in the simulated datasets
phyla_example <- function() {
  plot_df = data.frame(matrix(nrow = 0, ncol = 3))
  colnames(plot_df) <- c("TaxID", "abundance", "Dataset")
  level = 2 # Phylum level
  
  # Import each simulated dataset core profile
  for (dataset in list.files("analysis/simulation_outputs")) {
    unknown_proportion <- ifelse(dataset=="gut", as.character(5), as.character(75))
    truth = renormalize(preprocess_simulation_truths(dataset, level))
    standards = truth[,grepl(paste0("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_"), colnames(truth)) | colnames(truth) == "TaxID"]
    join_df <- data.frame("TaxID" = standards$TaxID, "abundance" = rowMeans(standards[,colnames(standards) != "TaxID"]), "Dataset" = dataset)
    plot_df <- rbind(plot_df, join_df)
  }
  
  plot_df$Dataset <- case_when(plot_df$Dataset == "soil" ~ "Soil",
                               plot_df$Dataset == "ocean" ~ "Ocean",
                               plot_df$Dataset == "gut" ~ "Animal gut")
  
  tmp <- aggregate(abundance~TaxID, data=plot_df, FUN=sum)
  tmp <- tmp[order(tmp$abundance, decreasing = T),]
  
  nphylum = 15
  plot_df <- plot_df[plot_df$TaxID %in% tmp$TaxID[1:nphylum],]
  
  colAdd <- c('#f58231', '#2f4f4f', '#eee8aa', '#4363d8', '#e6194b', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075')
  colAdd <- c(colAdd[1:(nphylum - 1)], '#808080', '#000000')
  col = scale_color_manual(values = colAdd)
  fil = scale_fill_manual(values = colAdd)
  
  plot_df <- rbind(plot_df, data.frame("TaxID" = "Other", 
                                       "abundance" = 100 - aggregate(abundance~Dataset, data=plot_df, FUN=sum)$abundance,
                                       "Dataset" = aggregate(abundance~Dataset, data=plot_df, FUN=sum)$Dataset)
                   )
  
  plot_df$TaxID[plot_df$TaxID == "UNCLASSIFIED"] <- "Unclassified"
  
  # Convert from NCBI IDs to known names
  level_list <- fread(paste0("analysis/databases/ncbi_taxdump/names.dmp"), sep="\t", header = F)
  level_list <- level_list[,c(1,3,7)]
  colnames(level_list) <- c("TaxID", "Name", "Status")
  level_list <- level_list[level_list$Status == "scientific name" | (!duplicated(level_list$Status) & !duplicated(level_list$Status, fromLast = T)),]
  tax_list <- level_list$Name
  tax_list <- c(tax_list, "Unclassified", "Other")
  names(tax_list) <- c(level_list$TaxID, "Unclassified", "Other")
  rm(level_list)
  
  plot_df$TaxID <- tax_list[plot_df$TaxID]
  
  abun_order = aggregate(abundance~TaxID, data=plot_df, FUN=sum)$TaxID[
    order(aggregate(abundance~TaxID, data=plot_df, FUN=sum)$abundance, decreasing = T)]
  
  abun_order <- c(abun_order[!(abun_order %in% c("Unclassified", "Other"))], "Other", "Unclassified")
  
  p1 <- ggplot(plot_df, aes(x=Dataset, y=abundance, fill=factor(TaxID, abun_order))) + 
    geom_bar(stat="identity") + 
    theme_bw() +
    col + fil + 
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(fill = "Phylum") + 
    ylab("Abundance") + 
    xlab("Simulated dataset")
  
  ggsave(paste0("analysis/figures/overview_figures/sim_phylum.png"), width=15, height=13.5, units = 'cm', dpi=1000)
}
phyla_example()

# Show how much abundance comes from known and unknown taxa
simulation_explanation <- function() {
  unknown_taxa = read.csv('analysis/databases/ncbi_taxdump/unknown_taxa.tsv', header = F)[,1]
  plot_df = data.frame(matrix(nrow = 0, ncol = 5))
  colnames(plot_df) <- c("Dataset", "Level", "Type", "Value", "file_name")
  
  # Read in core files for all datasets
  for (dataset in list.files("analysis/simulation_outputs")) {
    unknown_proportion <- ifelse(dataset=="gut", as.character(5), as.character(75))
    in_files = list.files(paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles"))
    in_files = paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles/", in_files)
    in_files = in_files[grepl(paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles/profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9].txt"), in_files)]
    
    # Read in core profiles
    for (file in in_files) {
      profile = read.csv(file, sep='\t', skip = 4)
      
      # No eukaryotes
      profile <- profile[!grepl("^2759",profile[,3]),c(1,2,5)]
      
      colnames(profile) <- c("TaxID", "rank", "abundance")
      profile$TaxID <- as.integer(profile$TaxID)
      
      for (level in c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) {
        df_append <- data.frame("Dataset" = dataset,
                                "Level" = level, 
                                "Type" = c("NCBI available", "New NCBI (Unknown)", "New SGB (Unknown)"),
                                "Value" = c(sum(profile$abundance[profile$rank == level & !(profile$TaxID %in% unknown_taxa)]),
                                            sum(profile$abundance[profile$rank == level & profile$TaxID %in% unknown_taxa]),
                                            100 - sum(c(profile$abundance[profile$rank == level & !(profile$TaxID %in% unknown_taxa)],
                                                        profile$abundance[profile$rank == level & profile$TaxID %in% unknown_taxa]))
                                            ),
                                "file_name" = file
        )
        
        plot_df <- rbind(plot_df, df_append)
      }
    }
  }
  
  plot_df <- aggregate(Value~Dataset + Level + Type, data=plot_df, FUN=mean)
  
  col <- c("#38A061", "#AAAAAA", "#888888")
  names(col) <- c("NCBI available", "New NCBI (Unknown)", "New SGB (Unknown)")
  fil = scale_fill_manual(values = col)
  
  plot_df$Level <- case_when(plot_df$Level == "superkingdom" ~ "Kingdom",
                             plot_df$Level == "phylum" ~ "Phylum",
                             plot_df$Level == "class" ~ "Class",
                             plot_df$Level == "order" ~ "Order",
                             plot_df$Level == "family" ~ "Family",
                             plot_df$Level == "genus" ~ "Genus",
                             plot_df$Level == "species" ~ "Species"
  )
  
  plot_df$Dataset <- case_when(plot_df$Dataset == "gut" ~ "Animal gut",
                               plot_df$Dataset == "ocean" ~ "Ocean",
                               plot_df$Dataset == "soil" ~ "Soil")
  
  p1 <- ggplot(plot_df, aes(fill=Type, y=Value, x=factor(Level, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))) + 
    geom_bar(stat="identity") + 
    fil + 
    facet_wrap(~Dataset, scale="fixed", ncol=1) + 
    theme_linedraw() + 
    xlab("Taxonomic level") + 
    ylab("Percent abundance in core samples") + 
    theme(text=element_text(size=13),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title = element_blank(),
          legend.direction = "vertical",
          legend.position="bottom",
          strip.background = element_blank(),
          strip.text = element_text(colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA))
  
  dir.create(paste0("analysis/figures/simulation_setup/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/simulation_setup/","simulation_setup.png"), p1, width=6, height=18, units = 'cm', dpi=1000)
}
simulation_explanation()

# Show how much abundance comes from known and unknown taxa for each tool individually
simulation_explanation_by_tool <- function() {
  unknown_taxa = read.csv('analysis/databases/ncbi_taxdump/unknown_taxa.tsv', header = F)[,1]
  plot_df = data.frame(matrix(nrow = 0, ncol = 6))
  colnames(plot_df) <- c("Dataset", "Level", "Type", "Value", "file_name", "Tool")
  for (dataset in list.files("analysis/simulation_outputs")) {
    unknown_proportion <- ifelse(dataset=="gut", as.character(5), as.character(75))
    in_files = list.files(paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles"))
    in_files = paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles/", in_files)
    in_files = in_files[grepl(paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles/profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9].txt"), in_files)]
    
    for (file in in_files) {
      profile = read.csv(file, sep='\t', skip = 4)
      
      # No eukaryotes
      profile <- profile[!grepl("^2759",profile[,3]),c(1,2,5)]
      
      colnames(profile) <- c("TaxID", "rank", "abundance")
      profile$TaxID <- as.integer(profile$TaxID)
      
      for (level in c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) {
        level_to_read <- ifelse(level == "superkingdom", "kingdom", level)
        level_list <- read.csv(paste0("analysis/databases/taxa_lists/", level_to_read, ".tsv"), sep="\t")
        
        for (i in 1:length(tool_names)) {
          database_to_look = case_when(i == 1 ~ "centrifuge",
                                       i == 2 ~ "kraken",
                                       i == 3 ~ "metaphlan2",
                                       i == 4 ~ "metaphlan3",
                                       i == 5 ~ "metaphlan4",
                                       i == 6 ~ "metaxa2",
                                       i == 7 ~ "mOTUs3",
                                       i == 8 ~ "gtdbtk",
                                       i == 9 ~ "gtdbtk",
                                       i == 10 ~ "phylophlan3",
                                       i == 11 ~ "phylophlan3")
          df_append <- data.frame("Dataset" = dataset,
                                  "Level" = level, 
                                  "Type" = c("NCBI available and in database", "NCBI available but not in database", "New NCBI (Unknown)", "New SGB (Unknown)"),
                                  "Value" = c(sum(profile$abundance[profile$rank == level & !(profile$TaxID %in% unknown_taxa) & profile$TaxID %in% level_list[,database_to_look][!is.na(level_list[,database_to_look])]]),
                                              sum(profile$abundance[profile$rank == level & !(profile$TaxID %in% unknown_taxa) & !(profile$TaxID %in% level_list[,database_to_look][!is.na(level_list[,database_to_look])])]),
                                              sum(profile$abundance[profile$rank == level & profile$TaxID %in% unknown_taxa]),
                                              100 - sum(c(profile$abundance[profile$rank == level & !(profile$TaxID %in% unknown_taxa)],
                                                          profile$abundance[profile$rank == level & profile$TaxID %in% unknown_taxa]))
                                  ),
                                  "file_name" = file,
                                  "Tool" = i
          )
          plot_df <- rbind(plot_df, df_append)
        }
          
        
      }
    }
  }
  
  plot_df$Tool <- tool_names[plot_df$Tool]
  plot_df <- aggregate(Value~Dataset + Level + Type + Tool, data=plot_df, FUN=mean)
  
  col <- c("#38A061", "#ADD8E6", "#AAAAAA", "#888888")
  names(col) <- c("NCBI available and in database", "NCBI available but not in database",  "New NCBI (Unknown)", "New SGB (Unknown)")
  fil = scale_fill_manual(values = col)
  
  plot_df$Level <- case_when(plot_df$Level == "superkingdom" ~ "Kingdom",
                             plot_df$Level == "phylum" ~ "Phylum",
                             plot_df$Level == "class" ~ "Class",
                             plot_df$Level == "order" ~ "Order",
                             plot_df$Level == "family" ~ "Family",
                             plot_df$Level == "genus" ~ "Genus",
                             plot_df$Level == "species" ~ "Species"
  )
  
  plot_df$Dataset <- case_when(plot_df$Dataset == "gut" ~ "Animal gut",
                               plot_df$Dataset == "ocean" ~ "Ocean",
                               plot_df$Dataset == "soil" ~ "Soil")
  
  plot_df <- plot_df[plot_df$Tool %in% c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "PhyloPhlAn MEGAHIT"),]
  plot_df$Tool <- case_when(plot_df$Tool == "GTDB-Tk MEGAHIT" ~ "GTDB-Tk 2",
                            plot_df$Tool == "PhyloPhlAn MEGAHIT" ~ "PhyloPhlAn 3",
                            TRUE ~ plot_df$Tool
                            )
  
  p1 <- ggplot(plot_df, aes(fill=Type, y=Value, x=factor(Level, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))) + 
    geom_bar(stat="identity") + 
    fil + 
    facet_grid(Dataset~Tool, scale="fixed") + 
    theme_linedraw() + 
    xlab("Taxonomic level") + 
    ylab("Percent abundance in core samples") + 
    theme(text=element_text(size=16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.title = element_blank(),
          legend.position = "bottom")
  
  dir.create(paste0("analysis/figures/simulation_setup/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/simulation_setup/","simulation_setup_by_tool.png"), p1, width=44, height=18, units = 'cm', dpi=1000)
}
simulation_explanation_by_tool()

# Create a summary figure for precision, recall, and F1 for the core samples
prf1_summary <- function(ncbi_only, core_five = FALSE) {
  plot_df = data.frame(matrix(nrow = 0, ncol = 5))
  colnames(plot_df) <- c("Tool", "Dataset", "Level", "Metric", "Value")
  # Phylum, family, species
  for (level in c(2, 5, 7)) {
    for (dataset in list.files("analysis/simulation_outputs")) {
      if (ncbi_only) {
        profiles <- remove_non_ncbi(preprocess(dataset, level, FALSE))
      } else {
        profiles <- remove_unknown(preprocess(dataset, level, FALSE))
      }
      truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
      profiles = lapply(profiles, renormalize)
      profiles = lapply(profiles, threshold_sample, threshold)
      profiles = lapply(profiles, renormalize)
      
      prf1_table <- matrix(nrow = 3, ncol = length(profiles))
      for (i in 1:length(profiles)) {
        taxa_names <- list_taxa_by_sample(profiles[[i]], threshold)
        truth_names <- list_taxa_by_sample(truth, 0)
        
        # Use core samples
        sample_names = names(taxa_names)[2:length(taxa_names)]
        if (dataset %in% c("soil", "ocean")) {
          standards = grepl("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.75\\.mut_rate\\.0\\.0_", sample_names)
        } else {
          standards = grepl("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.5\\.mut_rate\\.0\\.0_", sample_names)
        }
        
        df_addition <- data.frame(cbind(tool_names[i], dataset, level, "F1", calc_precision_recall_f1(taxa_names, truth_names)[3,standards]))
        colnames(df_addition) <- c("Tool", "Dataset", "Level", "Metric", "Value")
        plot_df <- rbind(plot_df, df_addition)
        
        df_addition <- data.frame(cbind(tool_names[i], dataset, level, "BC", 1 - calc_bc_sim(profiles[[i]], truth)[standards]))
        colnames(df_addition) <- c("Tool", "Dataset", "Level", "Metric", "Value")
        plot_df <- rbind(plot_df, df_addition)
      }
    }
  }
  plot_df$Dataset <- case_when(plot_df$Dataset == "soil" ~ "Soil",
                               plot_df$Dataset == "ocean" ~ "Ocean",
                               plot_df$Dataset == "gut" ~ "Animal gut")
  plot_df$Level <- case_when(plot_df$Level == "2" ~ "Phylum",
                             plot_df$Level == "5" ~ "Family",
                             plot_df$Level == "7" ~ "Species")
  plot_df$Level <- factor(plot_df$Level, levels = c("Phylum", "Family", "Species"))
  
  plot_df$Value <- as.numeric(plot_df$Value)

  plot_df2 <- plot_df %>%
    group_by(Dataset, Level, Tool, Metric) %>%
    summarise_at(vars(Value), list(mean=mean, sd=sd)) %>% 
    as.data.frame()
  
  plot_df <- left_join(plot_df, plot_df2)
  plot_df$Tool <- factor(plot_df$Tool, c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes"))
  
  if (core_five) {
    plot_df <- plot_df[plot_df$Tool %in% tool_core,]
    colAdd <- colAdd[tool_names %in% tool_core]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  
  plot_df$Metric[plot_df$Metric == "BC"] <- "1 - Bray Curtis dissimilarity"
  
  p1 <- ggplot(plot_df, aes(x=Tool, y=Value, fill=Tool)) + 
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0), size=3, alpha=0.7) +
    facet_nested(Dataset ~ factor(Metric, c("F1", "1 - Bray Curtis dissimilarity")) + Level) +
    theme_linedraw() + 
    ylim(0, 1) + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 12, hjust = 1),
          legend.position = "none",
          strip.text = element_text(colour = 'black'),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    col + fil

  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  core_five <- ifelse(core_five, "core_tools", "all_tools")
  
  ggsave(paste0("analysis/figures/simulation_core/f1_and_bc_", ncbi_only, "_", core_five, ".png"), p1, width = 10, height = 6)
}

if (!file.exists(paste0("analysis/figures/simulation_core/f1_and_bc_", "NCBI_only", "_core_tools.png"))) {
  prf1_summary(TRUE)
  prf1_summary(TRUE, TRUE)
}
if (!file.exists(paste0("analysis/figures/simulation_core/f1_and_bc_", "all_taxa", "_all_tools.png"))) {
  prf1_summary(FALSE)
  prf1_summary(FALSE, TRUE)
}

# Show how BC dissimilarity from true profile changes with sample features
bray_vs_parameter_individual <- function(dataset, level, ncbi_only=TRUE) {
  data <- preprocess(dataset, level)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  # Check over all parameters
  bray_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(bray_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  for (i in 1:length(profiles)) {
    variable_names = gsub("_sample_.*", "", colnames(profiles[[i]])[-1])
    df_append <- data.frame(tool_names[i],
                             gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                             gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                             gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                             gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                             gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                             gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                             1 - calc_bc_sim(profiles[[i]], truth)
                             )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
    bray_df <- data.frame(rbind(bray_df, df_append))
  }
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  colAdd <- colAdd[tool_names %in% bray_df$Method]
  col = scale_color_manual(values = colAdd)
  
  bray_df$Method <- factor(bray_df$Method, c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes"))
  
  # One plot for each parameter
  p1 <- ggplot(bray_df[bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=n, y=value, color=Method)) + geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
    xlab("Species count") + 
    col + 
    scale_x_continuous(trans='log10', breaks = c(75, 150, 300, 600)) + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(paste0("analysis/figures/BrayVsParameter/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/BrayVsParameter/",dataset,"/n/",dataset,"_",level_name,"_",ncbi_only, ".png"), p1, width=20, height=20, units = 'cm', dpi=1000)
  
  p2 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=sample_size, y=value, color=Method)) + 
    geom_jitter(width = 0.05) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
    xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10', breaks = c(0.05, 0.5, 1.5, 7.5, 30)) + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(paste0("analysis/figures/BrayVsParameter/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/BrayVsParameter/",dataset,"/sample_size/",dataset,"_",level_name,"_",ncbi_only, ".png"), p2, width=20, height=20, units = 'cm', dpi=1000)
  
  if (level == 7) {
    p3 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
                 aes(x=100 * unknown, y=value, color=Method)) + 
      geom_jitter(width = 1) +
      geom_smooth(method = 'glm', formula = 'y~log(100-x+0.001)', size=1, alpha=0.15) +
      theme_classic() + 
      ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
      xlab("Unknown species (%)") + 
      col + 
      theme(text=element_text(size=12),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10)) + 
      guides(color=guide_legend(title="Tool"))
  } else {
    p3 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
                 aes(x=100 * unknown, y=value, color=Method)) + 
      geom_jitter(width = 1) +
      geom_smooth(method = "glm", 
                  method.args = list(family = "binomial"), 
                  se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
      theme_classic() + 
      ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
      xlab("Unknown species (%)") + 
      col + 
      theme(text=element_text(size=12),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10)) + 
      guides(color=guide_legend(title="Tool"))
  }
  
  dir.create(paste0("analysis/figures/BrayVsParameter/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/BrayVsParameter/",dataset,"/unknown_prop/",dataset,"_",level_name,"_",ncbi_only, ".png"), p3, width=20, height=20, units = 'cm', dpi=1000)
  
  p4 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"],], 
               aes(x=mut_rate, y=value, color=Method)) + 
    geom_jitter(width = 0.001) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
    xlab("Mutation rate") + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(paste0("analysis/figures/BrayVsParameter/",dataset,"/mut_rate/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/BrayVsParameter/",dataset,"/mut_rate/",dataset,"_",level_name,"_",ncbi_only, ".png"), p4, width=20, height=20, units = 'cm', dpi=1000)
  
  p5 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=k, y=value, color=Method)) + 
    geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
    xlab("Genome size distribution\n(lower k is more skewed)") + 
    scale_x_continuous(trans='log10', breaks = c(40, 100, 250)) + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(paste0("analysis/figures/BrayVsParameter/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/BrayVsParameter/",dataset,"/k/",dataset,"_",level_name,"_",ncbi_only, ".png"), p5, width=20, height=20, units = 'cm', dpi=1000)
  
  p6 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=sigma, y=value, color=Method)) + 
    geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
    xlab("Abundance skew") +
    scale_x_continuous(trans='log10', breaks = c(0.5, 1, 2, 4)) + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  dir.create(paste0("analysis/figures/BrayVsParameter/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/BrayVsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_",ncbi_only, ".png"), p6, width=20, height=20, units = 'cm', dpi=1000)
}
make_all_bray_vs_parameter_individual <- function() {
  for (dataset in list.files("analysis/simulation_outputs")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("analysis/figures/BrayVsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_", "NCBI_only", ".png"))) {
        bray_vs_parameter_individual(dataset, level)
      }
      if (!file.exists(paste0("analysis/figures/BrayVsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_", "all_taxa", ".png"))) {
        bray_vs_parameter_individual(dataset, level, FALSE)
      }
    }
  }
}
make_all_bray_vs_parameter_individual()

# Show how F1, precision, and recall from true profile changes with sample features
f1_vs_parameter_individual <- function(dataset, level, ncbi_only=TRUE, metric="F1") {
  if (!(metric %in% c("F1", "Precision", "Recall"))) {
    stop("Metric not allowed")
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
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  f1_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(f1_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  for (i in 1:length(profiles)) {
    taxa_names <- list_taxa_by_sample(profiles[[i]], threshold)
    truth_names <- list_taxa_by_sample(truth, 0)
    
    variable_names = gsub("_sample_.*", "", colnames(profiles[[i]])[-1])
    
    extract_index = case_when(metric == "F1" ~ 3,
                              metric == "Precision" ~ 1,
                              metric == "Recall" ~ 2)
    
    df_append <- data.frame(tool_names[i],
                            gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                            calc_precision_recall_f1(taxa_names, truth_names)[extract_index,]
    )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
    f1_df <- data.frame(rbind(f1_df, df_append))
  }
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  colAdd <- colAdd[tool_names %in% f1_df$Method]
  col = scale_color_manual(values = colAdd)
  
  f1_df$Method <- factor(f1_df$Method, c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes"))
  
  p1 <- ggplot(f1_df[f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],], 
               aes(x=n, y=value, color=Method)) + geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0(metric, " at the \n", level_name, " level")) + 
    xlab("Species count") + 
    col + 
    scale_x_continuous(trans='log10', breaks = c(75, 150, 300, 600)) + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/n/",dataset,"_",level_name,"_",ncbi_only, ".png"), p1, width=20, height=20, units = 'cm', dpi=1000)
  
  p2 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],], 
               aes(x=sample_size, y=value, color=Method)) + 
    geom_jitter(width = 0.05) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0(metric, " at the \n", level_name, " level")) + 
    xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10', breaks = c(0.05, 0.5, 1.5, 7.5, 30)) + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/sample_size/",dataset,"_",level_name,"_",ncbi_only, ".png"), p2, width=20, height=20, units = 'cm', dpi=1000)
  
  if (level == 7) {
    p3 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],], 
                 aes(x=100 * unknown, y=value, color=Method)) + 
      geom_jitter(width = 1) +
      geom_smooth(method = 'glm', formula = 'y~log(100-x+0.001)', size=1, alpha=0.15, se = FALSE) +
      theme_classic() + 
      ylab(paste0(metric, " at the \n", level_name, " level")) + 
      xlab("Unknown species (%)") + 
      col + 
      theme(text=element_text(size=12),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10)) + 
      guides(color=guide_legend(title="Tool"))
  } else {
    p3 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],], 
                 aes(x=100 * unknown, y=value, color=Method)) + 
      geom_jitter(width = 1) +
      geom_smooth(method = "glm", 
                  method.args = list(family = "binomial"), 
                  se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
      theme_classic() + 
      ylab(paste0(metric, " at the \n", level_name, " level")) + 
      xlab("Unknown species (%)") + 
      col + 
      theme(text=element_text(size=12),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10)) + 
      guides(color=guide_legend(title="Tool"))
  }
  
  dir.create(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/unknown_prop/",dataset,"_",level_name,"_",ncbi_only, ".png"), p3, width=20, height=20, units = 'cm', dpi=1000)
  
  p4 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"],], 
               aes(x=mut_rate, y=value, color=Method)) + 
    geom_jitter(width = 0.001) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0(metric, " at the \n", level_name, " level")) + 
    xlab("Mutation rate") + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/mut_rate/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/mut_rate/",dataset,"_",level_name,"_",ncbi_only, ".png"), p4, width=20, height=20, units = 'cm', dpi=1000)
  
  p5 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],], 
               aes(x=k, y=value, color=Method)) + 
    geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0(metric, " at the \n", level_name, " level")) + 
    xlab("Genome size distribution\n(lower k is more skewed)") + 
    scale_x_continuous(trans='log10', breaks = c(40, 100, 250)) + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/k/",dataset,"_",level_name,"_",ncbi_only, ".png"), p5, width=20, height=20, units = 'cm', dpi=1000)
  
  p6 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$mut_rate==default_params["mut_rate"],], 
               aes(x=sigma, y=value, color=Method)) + 
    geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE, linewidth=1, alpha=0.15, formula='y~x') +
    theme_classic() + 
    ylab(paste0(metric, " at the \n", level_name, " level")) + 
    xlab("Abundance skew") +
    scale_x_continuous(trans='log10', breaks = c(0.5, 1, 2, 4)) + 
    col + 
    theme(text=element_text(size=12),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10)) + 
    guides(color=guide_legend(title="Tool"))
  dir.create(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figures/", metric, "VsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_",ncbi_only, ".png"), p6, width=20, height=20, units = 'cm', dpi=1000)
}
make_all_f1_vs_parameter_individual <- function() {
  for (dataset in list.files("analysis/simulation_outputs")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("analysis/figures/F1VsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_", "NCBI_only", ".png"))) {
        f1_vs_parameter_individual(dataset, level, TRUE, "Precision")
        f1_vs_parameter_individual(dataset, level, TRUE, "Recall")
        f1_vs_parameter_individual(dataset, level)
      }
      if (!file.exists(paste0("analysis/figures/F1VsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_", "all_taxa", ".png"))) {
        f1_vs_parameter_individual(dataset, level, FALSE, "Precision")
        f1_vs_parameter_individual(dataset, level, FALSE, "Recall")
        f1_vs_parameter_individual(dataset, level, FALSE)
      }
      
      
    }
  }
}
make_all_f1_vs_parameter_individual()


bc_by_rank <- function(dataset, unknown_proportion, use_title, ncbi_only=TRUE) {
  data_tmp <- preprocess(dataset, 1)
  
  output_df <- data.frame(matrix(ncol = 4, nrow = 7 * (length(data_tmp))))
  colnames(output_df) <- c("level", "method", "bc", "sd")
  
  for (level in 1:7) {
    data <- preprocess(dataset, level)
    
    # Normalize
    if (ncbi_only) {
      profiles <- remove_non_ncbi(data)
    } else {
      profiles <- remove_unknown(data)
    }
    
    profiles = lapply(profiles, renormalize)
    profiles = lapply(profiles, threshold_sample, threshold)
    data = lapply(profiles, renormalize)
    
    truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
    
    # data[[length(data)]][data[[length(data)]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[length(data)]]))] <- 0
    for (i in 1:length(data)) {
      #data[[i]][data[[i]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[i]]))] <- 0
      
      vals = 1 - calc_bc_sim(data[[i]], truth)[grepl(paste0("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up.", unknown_proportion, "[\\.0]*\\.mut_rate\\.0\\.0_"),colnames(data[[i]])[-1])]
      
      output_df[((level - 1) * length(data) + (i - 1) + 1):((level - 1) * length(data) + i ), 1:4] <- 
        cbind(level, i, mean(vals), sd(vals))
    }
    
  }
  
  output_df$Method <- tool_names[output_df$method]
  
  output_df$level <- case_when(output_df$level == 1 ~ "Kingdom",
                               output_df$level == 2 ~ "Phylum", 
                               output_df$level == 3 ~ "Class",
                               output_df$level == 4 ~ "Order",
                               output_df$level == 5 ~ "Family",
                               output_df$level == 6 ~ "Genus",
                               output_df$level == 7 ~ "Species")
  
  output_df$Method <- factor(output_df$Method, c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes"))
  
  if (use_title) {
    p <- ggplot(output_df, aes(y=bc, x=factor(level, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), 
                               group=Method)) + 
      geom_line(aes(color=Method), size=3) + 
      geom_errorbar(aes(ymin=bc-sd, ymax=bc+sd), width=2,
                    position=position_dodge(0.05), size=3) +
      theme_classic() +
      xlab(paste0("Unknown Proportion: ", unknown_proportion)) +
      ylab(paste0("1 - Bray Curtis dissimilarity")) + 
      col + 
      guides(color=guide_legend(title="Method")) + 
      theme(text=element_text(size=80),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1))
  } else {
    p <- ggplot(output_df, aes(y=bc, x=factor(level, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), 
                               group=Method)) + 
      geom_line(aes(color=Method), size=3) + 
      geom_errorbar(aes(ymin=bc-sd, ymax=bc+sd), width=2,
                    position=position_dodge(0.05), size=3) +
      theme_classic() +
      xlab(paste0("Unknown Proportion: ", unknown_proportion)) +
      ylab(paste0("1 - Bray Curtis dissimilarity")) + 
      col + 
      theme(legend.position = "none", text=element_text(size=80),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1))
  }
  return(p)
}
make_bc_by_rank <- function() {
  for (ncbi_only in c(TRUE, FALSE)) {
    for (dataset in list.files("simulation_outputs/inputs")) {
      ncbi_only_tmp <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
      if (!file.exists(paste0("figures/BC_By_Rank/",dataset,"/",dataset,"_",ncbi_only,".png"))) {
        dir.create(file.path("figures/BC_By_Rank/",dataset,"/"), showWarnings = FALSE, recursive = T)
        png(file=paste0("figures/BC_By_Rank/",dataset,"/",dataset,"_",ncbi_only_tmp,".png"), width = 6400, height = 2200)
        grid.arrange(as.grob(bc_by_rank(dataset, 0, F, ncbi_only)), 
                     as.grob(bc_by_rank(dataset, 0.25, F, ncbi_only)), 
                     as.grob(bc_by_rank(dataset, 0.5, F, ncbi_only)),
                     as.grob(bc_by_rank(dataset, 0.75, F, ncbi_only)),
                     as.grob(bc_by_rank(dataset, 1, T, ncbi_only)), nrow=1, ncol=5,
                     top = textGrob("Bray Curtis Dissimilarity vs True Profile",gp=gpar(fontsize=100,font=1)), widths=c(2,2,2,2,2.8))
        dev.off()
      }
    }
  }
}
make_bc_by_rank()

genome_size_densities <- function(dataset) {
  unknown_proportion <- ifelse(dataset=="gut", as.character(5), as.character(75))
  lengths_df = read.csv(paste0("simulation_outputs/inputs/", dataset, "/genome_sizes/genome_lengths.tsv"), sep = "\t")
  lengths_df$name <- case_when(!grepl("sgb_", lengths_df$name) ~ gsub("_genomic", "", lengths_df$name),
                               grepl("sgb_", lengths_df$name) ~ lengths_df$genome_ID)
  
  in_files = list.files(paste0("simulation_outputs/inputs/", dataset, "/true_profile/profiles"))
  in_files = paste0("simulation_outputs/inputs/", dataset, "/true_profile/profiles/", in_files)
  in_files = in_files[grepl(paste0("simulation_outputs/inputs/", dataset, "/true_profile/profiles/profile\\.n\\.300\\.size\\.7\\.5\\.k\\.[0-9]+\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9].txt"), in_files)]
  
  lengths_map = lengths_df$length
  names(lengths_map) = lengths_df$name
  
  df_out <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(df_out) <- c("k", "length", "sample")
  for (file in in_files) {
    profile = read.csv(file, sep='\t', skip = 4)
    genomes = profile$X_CAMI_genomeID[profile$X_CAMI_genomeID != ""]
    genomes <- case_when(grepl("BIN_UNKNOWN", genomes) ~ genomes,
                         grepl("UNKNOWN_NEW", genomes) ~ gsub("UNKNOWN_NEW_[A-Z_a-z]+\\.", "", genomes),
                         TRUE ~ gsub("[A-Z_a-z]+\\.", "",  genomes))
    df_append <- data.frame("k"=gsub(paste0("simulation_outputs/inputs/", dataset, "/true_profile/profiles/profile\\.n\\.300\\.size\\.7\\.5\\.k\\."), "", file) %>% gsub(pattern="\\..*", replacement="") %>% as.numeric(),
                            "length" = unname(lengths_map[genomes]),
                            "sample" = gsub(".*sample_", "", file) %>% gsub(pattern="\\..*", replacement="") %>% as.numeric())
    
    df_out <- rbind(df_out, df_append)
  }
  
  df_out$Distribution <- case_when(df_out$k == 40 ~ "Widest spread (k = 40)", 
                                   df_out$k == 100 ~ "Moderate spread (k = 100)",
                                   df_out$k == 250 ~ "Narrow spread (k = 250)")
  
  df_out$Distribution <- factor(df_out$Distribution, levels = c("Widest spread (k = 40)",
                                                                "Moderate spread (k = 100)",
                                                                "Narrow spread (k = 250)"))
  
  p1 <- ggplot(df_out[df_out$sample == 0,], aes(length/10^6, fill = Distribution, colour = Distribution)) +
    geom_density(alpha = 0.1, ) + 
    scale_x_continuous(trans="log10") + 
    theme_classic() + 
    ylab("Density") + 
    xlab("Genome size (Mb)") + 
    ggtitle("Genome size distribution")
  
  dir.create(paste0("figures/genome_skew/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("figures/genome_skew/",dataset,".png"), p1, width=15, height=12, units = 'cm', dpi=1000)
}
make_genome_size_densities <- function() {
  for (dataset in list.files("simulation_outputs/inputs")) {
    genome_size_densities(dataset)
  }
}
make_genome_size_densities()

unknown_evaluation <- function(dataset, level, core_only = TRUE) {
  data <- preprocess(dataset, level, FALSE, TRUE)
  data = lapply(data, renormalize)
  data = lapply(data, threshold_sample, threshold)
  data = lapply(data, renormalize)
  data = remove_non_ncbi(data)
  for (i in 1:length(data)) {
    data[[i]][nrow(data[[i]]) + 1, ] <- c("UNCLASSIFIED", 100 - colSums(data[[i]][,-1]))
  }
  
  truth = renormalize(preprocess_simulation_truths(dataset, level))
  truth <- truth[,grepl("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.[0-9.]+\\.mut_rate\\.0\\.0_sample_|TaxID", colnames(truth))]
  
  for (i in 1:length(data)) {
    data[[i]] <- data[[i]][,grepl("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.[0-9.]+\\.mut_rate\\.0\\.0_sample_|TaxID", colnames(data[[i]]))]
  }
  
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    df_append = data.frame("method" = i, 
                           "assigned" = c(as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1]), rep(NA, length(truth[truth$TaxID=="UNCLASSIFIED",][-1]) - ncol(tmp) + 1)),
                           "truth" = unname(as.numeric(truth[truth$TaxID=="UNCLASSIFIED",][-1])))
    df <- rbind(df, df_append)
  }
  
  df$method <- tool_names[df$method]
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  colAdd <- colAdd[tool_names %in% df$method]
  col = scale_color_manual(values = colAdd)
  
  for (tool_name in tool_names) {
    tmp <- df[df$method == tool_name,]
    m <- lm(assigned ~ truth, tmp)
    df$method[df$method == tool_name] <- paste0(tool_name, "\n(", "y = ", 
                                                format(abs(round(unname(coef(m)[1]), digits = 2)), nsmall = 2), 
                                                ifelse(round(unname(coef(m)[2]), digits = 2) > 0, " + ", " - "), 
                                                format(abs(round(unname(coef(m)[2]), digits = 2)), nsmall = 2), "x)\n")
  }
  
  c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
  df$method <- factor(df$method, c(unique(df$method)[grepl("Centrifuge", unique(df$method))],
                                   unique(df$method)[grepl("Kraken 2 \\/ Bracken 2", unique(df$method))],
                                   unique(df$method)[grepl("MetaPhlAn 2", unique(df$method))],
                                   unique(df$method)[grepl("MetaPhlAn 3", unique(df$method))],
                                   unique(df$method)[grepl("MetaPhlAn 4", unique(df$method))],
                                   unique(df$method)[grepl("Metaxa 2", unique(df$method))],
                                   unique(df$method)[grepl("mOTUs 3", unique(df$method))],
                                   unique(df$method)[grepl("GTDB-Tk MEGAHIT", unique(df$method))],
                                   unique(df$method)[grepl("GTDB-Tk metaSPAdes", unique(df$method))],
                                   unique(df$method)[grepl("PhyloPhlAn MEGAHIT", unique(df$method))],
                                   unique(df$method)[grepl("PhyloPhlAn metaSPAdes", unique(df$method))]))
  
  if (core_only) {
    df <- df[gsub("\n\\(.*", "", df$method) %in% tool_core,]
    colAdd <- colAdd[tool_names %in% tool_core]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  
  p1 <- ggplot(df, aes(x=truth, y=assigned, color=method)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y ~ x', size=1, alpha = 0.25) +
    theme_bw() + 
    ylab(paste0("Estimated unknown\nat the ", level_name, " level")) + 
    xlab(paste0("True unknown\nat the ", level_name, " level")) + 
    col + 
    ylim(0, 100) + 
    geom_abline() +
    theme(text=element_text(size=20),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(1, 'cm'),
          legend.key.width = unit(1, 'cm'),
          legend.position = "bottom",
          ) + 
    guides(color=guide_legend(title="Method", ncol=1, title.position = "top"))

  dir.create(file.path("figures/unknown_evaluation/",dataset,"/"), showWarnings = FALSE)
  if (core_only) {
    ggsave(paste0("figures/unknown_evaluation/",dataset,"/",dataset,"_",level_name,".png"), width=14, height=26, units = 'cm', dpi=1000)
  } else {
    ggsave(paste0("figures/unknown_evaluation/",dataset,"/",dataset,"_",level_name,"_all_tools.png"), width=20, height=24, units = 'cm', dpi=1000)
  }
}
make_all_unknown_evaluation <- function() {
  for (dataset in list.files("simulation_outputs/inputs")) {
    for (level in 1:7) {
      unknown_evaluation(dataset, level)
      unknown_evaluation(dataset, level, FALSE)
    }
  }
}
make_all_unknown_evaluation()

grid_analysis <- function(dataset, level, ncbi_only = TRUE) {
  data <- preprocess(dataset, level, FALSE, TRUE)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  # unknown_proportion <- ifelse(dataset=="gut", as.character(5), as.character(75))
  # for (i in 1:length(profiles)) {
  #   profiles[[i]] = profiles[[i]][,grepl(paste0("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.[0-9]+\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9]"), colnames(profiles[[i]]))]
  # }
  # truth = truth[,grepl(paste0("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.[0-9]+\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9]"), colnames(truth))]
  
  
  cormat = matrix(nrow = 2, ncol = length(data))
  signif_mat = matrix(nrow = 2, ncol = length(data))
  
  for (i in 1:length(profiles)) {
    bc_out <- run_bc_mantel(profiles[[i]], truth)
    alpha_out <- run_alpha_cor(profiles[[i]], truth)
    cormat[1, i] <- bc_out[1]
    signif_mat[1, i] <- bc_out[2]
    cormat[2, i] <- alpha_out[1]
    signif_mat[2, i] <- alpha_out[2]
  }
  
  colnames(cormat) <- tool_names
  rownames(cormat) <- c("Beta diversity", "Alpha diversity")
  melted_cormat <- melt(cormat, na.rm = T)
  
  colnames(signif_mat) <- tool_names
  rownames(signif_mat) <- c("Beta diversity", "Alpha diversity")
  melted_signif <- melt(signif_mat, na.rm = T)
  melted_cormat$`P-value below` <- melted_signif$value
  melted_cormat$`P-value below` <- case_when(melted_cormat$`P-value below` < 0.001 ~ "0.001",
                                             melted_cormat$`P-value below` < 0.01 ~ "0.01",
                                             melted_cormat$`P-value below` < 0.05 ~ "0.05",
                                             TRUE ~ "1")
  melted_cormat$`P-value below` <- factor(melted_cormat$`P-value below`, c("1", "0.05", "0.01", "0.001"))
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  ggplot(melted_cormat, aes(y=value, x=Var1, color = factor(Var2, tool_names), size=factor(`P-value below`, c("1", "0.05", "0.01", "0.001")))) +  
    geom_beeswarm(cex=7, stroke=2) + 
    theme_linedraw() + 
    xlab("Metric") + 
    col + 
    ylim(0, 1) + 
    theme(text=element_text(size=20),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=18), #change legend title font size
          legend.text = element_text(size=14),
          legend.position = "bottom",
          legend.direction = "vertical") + 
    guides(color = "none") + 
    scale_size_manual(values = c(1.5,3,4.5,6), drop = FALSE) +
    labs(size="P-value below") + 
    labs(color="Method") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 20, hjust = 1)) + 
    scale_y_continuous(
      "Mantel r (Beta diversity)", 
      sec.axis = sec_axis(~ ., name = "Spearman correlation (Alpha diversity)"),
      limits = c(0,1)
    )
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/simulation_structure_grid/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/simulation_structure_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=10, height=36, units = 'cm', dpi=1000, bg='#ffffff')
  
  melted_cormat <- melted_cormat[melted_cormat$Var2 %in% tool_core,]
  colAdd <- colAdd[tool_names %in% tool_core]
  col = scale_color_manual(values = colAdd)
  fil = scale_fill_manual(values = colAdd)
  
  ggplot(melted_cormat, aes(y=value, x=Var1, color = factor(Var2, tool_names), size=factor(`P-value below`, c("1", "0.05", "0.01", "0.001")))) + 
    geom_beeswarm(cex=11, stroke=2) + 
    theme_linedraw() + 
    xlab("Metric") + 
    col + 
    theme(text=element_text(size=20),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=18), #change legend title font size
          legend.text = element_text(size=14),
          legend.position = "bottom",
          legend.direction = "vertical") + 
    labs(color="Method") + 
    guides(size = "none",
           color = guide_legend(override.aes = list(size=6))) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 20, hjust = 1)) + 
    scale_size_manual(values = c(1.5,3,4.5,6), drop = FALSE) +
    scale_y_continuous(
      "Mantel r (Beta diversity)", 
      sec.axis = sec_axis(~ ., name = "Spearman correlation (Alpha diversity)"),
      limits = c(0,1)
    )
  
  ggsave(paste0("figures/simulation_structure_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_core.png"), width=10, height=36, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_grid_analysis <- function() {
  for (dataset in list.files("simulation_outputs/inputs")) {
    for (level in c(7)) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/simulation_structure_grid/",dataset,"/",dataset,"_",level_name,"_all_taxa.png"))) {
        grid_analysis(dataset, level, FALSE)
      }
    }
  }
}
make_all_grid_analysis()

grid_analysis_horizontal <- function(dataset, level, ncbi_only = TRUE) {
  data <- preprocess(dataset, level, FALSE, TRUE)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  cormat = matrix(nrow = 2, ncol = length(data))
  signif_mat = matrix(nrow = 2, ncol = length(data))
  
  for (i in 1:length(profiles)) {
    bc_out <- run_bc_mantel(profiles[[i]], truth)
    alpha_out <- run_alpha_cor(profiles[[i]], truth)
    cormat[1, i] <- bc_out[1]
    signif_mat[1, i] <- bc_out[2]
    cormat[2, i] <- alpha_out[1]
    signif_mat[2, i] <- alpha_out[2]
  }
  
  colnames(cormat) <- tool_names
  rownames(cormat) <- c("Beta diversity", "Alpha diversity")
  melted_cormat <- melt(cormat, na.rm = T)
  
  colnames(signif_mat) <- tool_names
  rownames(signif_mat) <- c("Beta diversity", "Alpha diversity")
  melted_signif <- melt(signif_mat, na.rm = T)
  melted_cormat$`P-value below` <- melted_signif$value
  melted_cormat$`P-value below` <- case_when(melted_cormat$`P-value below` < 0.001 ~ "0.001",
                                             melted_cormat$`P-value below` < 0.01 ~ "0.01",
                                             melted_cormat$`P-value below` < 0.05 ~ "0.05",
                                             TRUE ~ "1")
  melted_cormat$`P-value below` <- factor(melted_cormat$`P-value below`, c("1", "0.05", "0.01", "0.001"))
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  ggplot(melted_cormat, aes(y=value, x=Var1, color = factor(Var2, tool_names), size=factor(`P-value below`, c("1", "0.05", "0.01", "0.001")))) +  
    geom_jitter(width = 0.4) +
    theme_linedraw() + 
    xlab("Metric") + 
    col + 
    ylim(0, 1) + 
    theme(text=element_text(size=20),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=18), #change legend title font size
          legend.text = element_text(size=14),
          legend.position = "right",
          legend.direction = "vertical") + 
    scale_size_manual(values = c(1.5,3,4.5,6), drop = FALSE) +
    labs(size="P-value below") + 
    labs(color="Method") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 20, hjust = 1)) + 
    scale_y_continuous(
      "Mantel r (Beta diversity)", 
      sec.axis = sec_axis(~ ., name = "Spearman correlation (Alpha diversity)"),
      limits = c(0,1)
    ) + 
    coord_flip()
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/simulation_structure_grid/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/simulation_structure_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_horizontal.png"), width=36, height=12, units = 'cm', dpi=1000, bg='#ffffff')
  
  melted_cormat <- melted_cormat[melted_cormat$Var2 %in% tool_core,]
  colAdd <- colAdd[tool_names %in% tool_core]
  col = scale_color_manual(values = colAdd)
  fil = scale_fill_manual(values = colAdd)
  
  ggplot(melted_cormat, aes(y=value, x=Var1, color = factor(Var2, tool_names), size=factor(`P-value below`, c("1", "0.05", "0.01", "0.001")))) + 
    geom_jitter(width = 0.2) +
    theme_linedraw() + 
    xlab("Metric") + 
    col + 
    theme(text=element_text(size=20),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=18), #change legend title font size
          legend.text = element_text(size=14),
          legend.position = "right",
          legend.direction = "vertical") + 
    labs(color="Method") + 
    guides(size = "none",
           color = guide_legend(override.aes = list(size=6))) + 
    scale_size_manual(values = c(1.5,3,4.5,6), drop = FALSE) +
    scale_y_continuous(
      "Mantel r (Beta diversity)", 
      sec.axis = sec_axis(~ ., name = "Spearman correlation (Alpha diversity)"),
      limits = c(0,1)
    ) +
    coord_flip()
  
  ggsave(paste0("figures/simulation_structure_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_core_horizontal.png"), width=50, height=10, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_grid_analysis <- function() {
  for (dataset in list.files("simulation_outputs/inputs")) {
    for (level in c(7)) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      if (!file.exists(paste0("figures/simulation_structure_grid/",dataset,"/",dataset,"_",level_name,"_all_taxa_horizontal.png"))) {
        grid_analysis_horizontal(dataset, level, FALSE)
      }
    }
  }
}
make_all_grid_analysis()

f1_vs_parameter_merged <- function(dataset, level, ncbi_only=TRUE, metric="F1") {
  if (!(metric %in% c("F1", "Precision", "Recall"))) {
    stop("Metric not allowed")
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
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  f1_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(f1_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  for (i in 1:length(profiles)) {
    taxa_names <- list_taxa_by_sample(profiles[[i]], threshold)
    truth_names <- list_taxa_by_sample(truth, 0)
    
    variable_names = gsub("_sample_.*", "", colnames(profiles[[i]])[-1])
    
    extract_index = case_when(metric == "F1" ~ 3,
                              metric == "Precision" ~ 1,
                              metric == "Recall" ~ 2)
    
    df_append <- data.frame(tool_names[i],
                            gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                            calc_precision_recall_f1(taxa_names, truth_names)[extract_index,]
    )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
    f1_df <- data.frame(rbind(f1_df, df_append))
  }
  
  f1_df$Method <- factor(f1_df$Method, rev(c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")))
  
  f1_df <- f1_df[f1_df$Method %in% tool_core,]
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  colAdd <- colAdd[tool_names %in% tool_core]
  tool_names <- tool_names[tool_names %in% tool_core]
  col = scale_color_manual(values = colAdd)
  fil = scale_fill_manual(values = colAdd)
  
  p1 <- ggplot(f1_df[f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],],
               aes(x=n, y=value, color=factor(Method, tool_core))) + geom_jitter(width = 0.01) +
    geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
    theme_bw() +
    ylab(paste0(metric, " at the ", level_name, " level")) +
    xlab("Species count") +
    col +
    scale_x_continuous(trans='log10', breaks = c(75, 150, 300, 600)) +
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") +
    guides(color=guide_legend(title="Tool"))
  
  p2 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],],
               aes(x=sample_size, y=value, color=factor(Method, tool_core))) +
    geom_jitter(width = 0.05) +
    geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
    theme_bw() +
    ylab("") +
    xlab("Sample size (GB)") +
    scale_x_continuous(trans='log10', breaks = c(0.05, 0.5, 1.5, 7.5, 30)) +
    col +
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") +
    guides(color=guide_legend(title="Tool"))
  
  if (level == 7) {
    p3 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"] & f1_df$unknown != 1,],
                 aes(x=100 * unknown, y=value, color=factor(Method, tool_core))) +
      geom_jitter(width = 1) +
      geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
      theme_bw() +
      ylab("") +
      xlab("Unknown species (%)") +
      col +
      theme(text=element_text(size=16),
            plot.margin = unit(c(1,0,1,0), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10),
            legend.position = "none") +
      guides(color=guide_legend(title="Tool"))
  } else {
    p3 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],],
                 aes(x=100 * unknown, y=value, color=factor(Method, tool_core))) +
      geom_jitter(width = 1) +
      geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
      theme_bw() +
      ylab("") +
      xlab("Unknown species (%)") +
      col +
      theme(text=element_text(size=16),
            plot.margin = unit(c(1,0,1,0), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10),
            legend.position = "none") +
      guides(color=guide_legend(title="Tool"))
  }
  
  p4 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"],],
               aes(x=mut_rate, y=value, color=factor(Method, tool_core))) +
    geom_jitter(width = 0.001) +
    geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
    theme_bw() +
    ylab("") +
    xlab("Mutation rate") +
    col +
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") +
    guides(color=guide_legend(title="Tool"))
  
  g <- arrangeGrob(
    grobs = list(p1, p2, p3, p4),
    widths = c(1.1, 1, 1, 1), nrow = 1)
  
  dir.create(file.path("figures/F1VsParameter/",dataset,"/merged/"), showWarnings = FALSE)
  ggsave(paste0("figures/F1VsParameter/",dataset,"/merged/",dataset,"_",level_name,"_", ncbi_only, ".png"), g, width=44, height=10, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_f1_vs_parameter_merged <- function() {
  dataset <- "soil"
  for (level in c(2, 7)) {
    print(paste0("Processing ", dataset, " at level ", level))
    level_name <- case_when(level == 1 ~ "kingdom",
                            level == 2 ~ "phylum", 
                            level == 3 ~ "class",
                            level == 4 ~ "order",
                            level == 5 ~ "family",
                            level == 6 ~ "genus",
                            level == 7 ~ "species")
    f1_vs_parameter_merged(dataset, level)
  }
}
make_all_f1_vs_parameter_merged()

f1_vs_parameter_merged_heatmap <- function(dataset, level, ncbi_only=TRUE, metric="F1") {
  if (!(metric %in% c("F1", "Precision", "Recall"))) {
    stop("Metric not allowed")
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
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  f1_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(f1_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  for (i in 1:length(profiles)) {
    taxa_names <- list_taxa_by_sample(profiles[[i]], threshold)
    truth_names <- list_taxa_by_sample(truth, 0)
    
    variable_names = gsub("_sample_.*", "", colnames(profiles[[i]])[-1])
    
    extract_index = case_when(metric == "F1" ~ 3,
                              metric == "Precision" ~ 1,
                              metric == "Recall" ~ 2)
    
    df_append <- data.frame(tool_names[i],
                            gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                            calc_precision_recall_f1(taxa_names, truth_names)[extract_index,]
    )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
    f1_df <- data.frame(rbind(f1_df, df_append))
  }
  
  f1_df$Method <- factor(f1_df$Method, rev(c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")))
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  colAdd <- colAdd[tool_names %in% f1_df$Method]
  col = scale_color_manual(values = colAdd)
  
  p1 <- f1_df[f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],] %>%
    complete(n,Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(n,Method) %>%
    ggplot(aes(x=factor(n), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=metric) +
    theme_classic() + 
    xlab("Species count") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") + 
    guides(color=guide_legend(title="Tool"))
  
  p2 <- f1_df[f1_df$n==default_params["n"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],] %>%
    complete(sample_size, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(sample_size, Method) %>%
    ggplot(aes(x=factor(sample_size), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=metric, na.value = "gray") +
    theme_classic() + 
    xlab("Sample size (GB)") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()
    ) + 
    guides(color=guide_legend(title="Tool"))
  
  p3 <- f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],] %>%
    complete(unknown, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(unknown, Method) %>%
    ggplot(aes(x=factor(100 * unknown), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=metric, na.value = "gray") +
    theme_classic() + 
    xlab("Unknown species (%)") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p4 <- f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"],] %>%
    complete(mut_rate, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(mut_rate, Method) %>%
    ggplot(aes(x=factor(mut_rate), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0(metric, " at the\n", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Mutation rate") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p5 <- f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$mut_rate==default_params["mut_rate"] & f1_df$sigma==default_params["sigma"],] %>%
    complete(k, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(k, Method) %>%
    ggplot(aes(x=factor(k), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0(metric, " at the\n", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Genome skew") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p6 <- f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$mut_rate==default_params["mut_rate"] & f1_df$k==default_params["k"],] %>%
    complete(sigma, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(sigma, Method) %>%
    ggplot(aes(x=factor(sigma), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0(metric, " at the\n", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Abundance skew") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  g <- arrangeGrob(
    grobs = list(p1, p2, p3, p4, p5, p6),
    widths = c(1.9, 1.2, 1.2, 1, 0.75, 1.55), nrow = 1)
  
  dir.create(file.path("figures/", metric, "VsParameter/",dataset,"/merged/"), showWarnings = FALSE)
  ggsave(paste0("figures/", metric, "VsParameter/",dataset,"/merged/",dataset,"_",level_name,"_", ncbi_only, "_heatmap.png"), g, width=50, height=14, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_f1_vs_parameter_merged_heatmap <- function(metric = "F1") {
  for (dataset in c("soil", "ocean", "gut")) {
    for (level in c(2, 7)) {
      print(paste0("Processing ", dataset, " at level ", level))
      level_name <- case_when(level == 1 ~ "kingdom",
                              level == 2 ~ "phylum", 
                              level == 3 ~ "class",
                              level == 4 ~ "order",
                              level == 5 ~ "family",
                              level == 6 ~ "genus",
                              level == 7 ~ "species")
      f1_vs_parameter_merged_heatmap(dataset, level, TRUE, metric)
    }
  }
}
make_all_f1_vs_parameter_merged_heatmap()
make_all_f1_vs_parameter_merged_heatmap("Precision")
make_all_f1_vs_parameter_merged_heatmap("Recall")

bray_vs_parameter_merged <- function(dataset, level, ncbi_only=TRUE) {
  data <- preprocess(dataset, level)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  bray_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(bray_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  for (i in 1:length(profiles)) {
    variable_names = gsub("_sample_.*", "", colnames(profiles[[i]])[-1])
    df_append <- data.frame(tool_names[i],
                            gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                            1 - calc_bc_sim(profiles[[i]], truth)
    )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
    bray_df <- data.frame(rbind(bray_df, df_append))
  }
  
  bray_df$Method <- factor(bray_df$Method, rev(c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")))
  
  bray_df <- bray_df[bray_df$Method %in% tool_core,]
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  colAdd <- colAdd[tool_names %in% tool_core]
  tool_names <- tool_names[tool_names %in% tool_core]
  col = scale_color_manual(values = colAdd)
  fil = scale_fill_manual(values = colAdd)
  
  p1 <- ggplot(bray_df[bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],],
               aes(x=n, y=value, color=factor(Method, tool_core))) + geom_jitter(width = 0.01) +
    geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
    theme_bw() +
    ylab(paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level")) +
    xlab("Species count") +
    col +
    scale_x_continuous(trans='log10', breaks = c(75, 150, 300, 600)) +
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") +
    guides(color=guide_legend(title="Tool"))
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  
  p2 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],],
               aes(x=sample_size, y=value, color=factor(Method, tool_core))) +
    geom_jitter(width = 0.05) +
    geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
    theme_bw() +
    ylab("") +
    xlab("Sample size (GB)") +
    scale_x_continuous(trans='log10', breaks = c(0.05, 0.5, 1.5, 7.5, 30)) +
    col +
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") +
    guides(color=guide_legend(title="Tool"))
  
  if (level == 7) {
    p3 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"] & bray_df$unknown != 1,],
                 aes(x=100 * unknown, y=value, color=factor(Method, tool_core))) +
      geom_jitter(width = 1) +
      geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
      theme_bw() +
      ylab("") +
      xlab("Unknown species (%)") +
      scale_x_continuous(breaks=c(0,25,50,75)) + 
      col +
      theme(text=element_text(size=16),
            plot.margin = unit(c(1,0,1,0), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10),
            legend.position = "none") +
      guides(color=guide_legend(title="Tool"))
  } else {
    p3 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],],
                 aes(x=100 * unknown, y=value, color=factor(Method, tool_core))) +
      geom_jitter(width = 1) +
      geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
      theme_bw() +
      ylab("") +
      scale_x_continuous(breaks=c(0,25,50,75,100)) + 
      xlab("Unknown species (%)") +
      col +
      theme(text=element_text(size=16),
            plot.margin = unit(c(1,0,1,0), "cm"),
            legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=10),
            legend.position = "none") +
      guides(color=guide_legend(title="Tool"))
  }
  
  p4 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"],],
               aes(x=mut_rate, y=value, color=factor(Method, tool_core))) +
    geom_jitter(width = 0.001) +
    geom_smooth(method = drm, method.args = list(fct = L.4()), se = F) +
    theme_bw() +
    ylab("") +
    xlab("Mutation rate") +
    col +
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") +
    guides(color=guide_legend(title="Tool"))
  
  g <- arrangeGrob(
    grobs = list(p1, p2, p3, p4),
    widths = c(1.1, 1, 1, 1), nrow = 1)
  
  dir.create(file.path("figures/BrayVsParameter/",dataset,"/merged/"), showWarnings = FALSE)
  ggsave(paste0("figures/BrayVsParameter/",dataset,"/merged/",dataset,"_",level_name,"_", ncbi_only, ".png"), g, width=44, height=10, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_bray_vs_parameter_merged <- function() {
  dataset = "soil"
  for (level in c(2,7)) {
    bray_vs_parameter_merged(dataset, level)
  }
}
make_all_bray_vs_parameter_merged()

bray_vs_parameter_merged_heatmap <- function(dataset, level, ncbi_only=TRUE) {
  data <- preprocess(dataset, level)
  
  if (ncbi_only) {
    profiles <- remove_non_ncbi(data)
  } else {
    profiles <- remove_unknown(data)
  }
  
  profiles = lapply(profiles, renormalize)
  profiles = lapply(profiles, threshold_sample, threshold)
  profiles = lapply(profiles, renormalize)
  
  truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
  
  bray_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(bray_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  for (i in 1:length(profiles)) {
    variable_names = gsub("_sample_.*", "", colnames(profiles[[i]])[-1])
    df_append <- data.frame(tool_names[i],
                            gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                            1 - calc_bc_sim(profiles[[i]], truth)
    )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
    bray_df <- data.frame(rbind(bray_df, df_append))
  }
  
  bray_df$Method <- factor(bray_df$Method, rev(c("Centrifuge", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")))
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  colAdd <- colAdd[tool_names %in% bray_df$Method]
  col = scale_color_manual(values = colAdd)
  
  p1 <- bray_df[bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],] %>%
    complete(n, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(n, Method) %>%
    ggplot(aes(x=factor(n), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0("1 - Bray Curtis Dissimilarity\nat the ", level_name, " level")) +
    theme_classic() + 
    xlab("Species count") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none") + 
    guides(color=guide_legend(title="Tool"))
  
  p2 <- bray_df[bray_df$n==default_params["n"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],] %>%
    complete(sample_size, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(sample_size, Method) %>%
    ggplot(aes(x=factor(sample_size), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0("1 - Bray Curtis Dissimilarity\nat the ", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Sample size (GB)") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()
    ) + 
    guides(color=guide_legend(title="Tool"))
  
  p3 <- bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],] %>%
    complete(unknown, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(unknown, Method) %>%
    ggplot(aes(x=factor(100 * unknown), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0("1 - Bray Curtis Dissimilarity\nat the ", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Unknown species (%)") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p4 <- bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"],] %>%
    complete(mut_rate,Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(mut_rate, Method) %>%
    ggplot(aes(x=factor(mut_rate), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0("1 - Bray Curtis Dissimilarity\nat the ", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Mutation rate") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p5 <- bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$mut_rate==default_params["mut_rate"] & bray_df$sigma==default_params["sigma"],] %>%
    complete(k, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(k, Method) %>%
    ggplot(aes(x=factor(k), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0("1 - Bray Curtis Dissimilarity\nat the ", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Genome skew") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p6 <- bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$mut_rate==default_params["mut_rate"] & bray_df$k==default_params["k"],] %>%
    complete(sigma, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(sigma, Method) %>%
    ggplot(aes(x=factor(sigma), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), 
                         name=paste0("1 - Bray Curtis dissimilarity\nat the ", level_name, " level"), na.value = "gray") +
    theme_classic() + 
    xlab("Abundance skew") + 
    theme(text=element_text(size=16),
          plot.margin = unit(c(1,0,1,0), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  g <- arrangeGrob(
    grobs = list(p1, p2, p3, p4, p5, p6),
    widths = c(1.9, 1.2, 1.2, 1, 0.75, 2), nrow = 1)
  
  dir.create(file.path("figures/BrayVsParameter/",dataset,"/merged/"), showWarnings = FALSE)
  ggsave(paste0("figures/BrayVsParameter/",dataset,"/merged/",dataset,"_",level_name,"_", ncbi_only, "_heatmap.png"), g, width=50, height=14, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_bray_vs_parameter_merged_heatmap <- function() {
  for (dataset in c("soil", "ocean", "gut")) {
    for (level in c(2,7)) {
      bray_vs_parameter_merged_heatmap(dataset, level)
    }
  }
}
make_all_bray_vs_parameter_merged_heatmap()

assess_reconstructed <- function(dataset, core_only = TRUE) {
  plot_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(plot_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  for (file in c(paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT PhyloPhlAn 3/genome_matching.tsv"),
                 paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT GTDB-Tk/genome_matching.tsv"),
                 paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes PhyloPhlAn 3/genome_matching.tsv"),
                 paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes GTDB-Tk/genome_matching.tsv"))) {
    genomes_matching <- read.csv(file, sep = "\t")
    
    if (grepl("GTDB-Tk", file)) {
      genomes_matching$reconstructed_tax <- gsub(pattern=" ",replacement="_", genomes_matching$reconstructed_tax) %>% tolower()
      genomes_matching$reconstructed_tax[genomes_matching$reconstructed_tax == "unknown"] <- "UNCLASSIFIED"
      genomes_matching$reconstructed_tax <- gsub("\\|p__\\|.*|\\|c__\\|.*|\\|o__\\|.*|\\|f__\\|.*|\\|g__\\|.*|\\|s__$|", "", genomes_matching$reconstructed_tax)
      
      map <- read.csv(paste0("databases/standardized_databases/GTDBTk.tsv"), sep = "\t")
      map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
      genomes_matching$reconstructed_tax <- mapvalues(genomes_matching$reconstructed_tax, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
    } else {
      genomes_matching$reconstructed_tax <- gsub(pattern=" ",replacement="_", genomes_matching$reconstructed_tax) %>% tolower()
      genomes_matching$reconstructed_tax[grepl("^sgb", genomes_matching$reconstructed_tax)] <- paste0("s__", genomes_matching$reconstructed_tax[grepl("^sgb", genomes_matching$reconstructed_tax)], "_MEGAHIT")
      
      map <- read.csv(paste0("databases/standardized_databases/PhyloPhlAn3.tsv"), sep = "\t")
      map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
      genomes_matching$reconstructed_tax <- mapvalues(genomes_matching$reconstructed_tax, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
    }
    
    keep_only_ncbi_name <- function(x) {
      if (is.na(x)) {
        return (NA)
      }
      tmp <- unlist(str_split(x, "\\|"))
      return (paste0(tmp[!is.na(as.numeric(tmp))], collapse = "|"))
    }
    
    genomes_matching$reconstructed_tax <- sapply(genomes_matching$reconstructed_tax, keep_only_ncbi_name)
    
    reconstructed <- function(x) {
      if (is.na(x[5])) {
        return (FALSE)
      }
      tmp2 <- unlist(str_split(x[5], "\\|"))
      tmp2 <- tmp2[!is.na(as.numeric(tmp2))]
      return (grepl(paste0("^", paste0(tmp2, collapse = "\\|")), x[4]))
    }
    
    genomes_matching$reconstructed <- apply(genomes_matching, 1, reconstructed)
    
    genomes_matching$unknown <- grepl("UNKNOWN", genomes_matching$real_org_id)
    genomes_matching$abundance_reconstructed[is.na(genomes_matching$abundance_reconstructed)] <- 0

    for_norm <- aggregate(abundance_real~sample, data=genomes_matching[genomes_matching$reconstructed & genomes_matching$unknown,], FUN=sum, drop=FALSE)
    genomes_matching <- aggregate(abundance_real~sample, data=genomes_matching[genomes_matching$reconstructed & genomes_matching$unknown,], FUN=sum, drop=FALSE)
    
    # 0.05 and 0.5 had nothing
    genomes_matching <- rbind(genomes_matching, data.frame("sample" = paste0("profile.n.300.size.", rep(c("0.05", "0.5"), each=5), ".k.100.sigma.1.0.up.", default_params["unknown"], ".mut_rate.0.0_sample_", rep(0:4, 2)),
                                                           "abundance_real" = 0))
    
    
    variable_names = gsub("_sample_.*", "", genomes_matching$sample)
    
    cur_tool_name <- case_when(file == paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT PhyloPhlAn 3/genome_matching.tsv") ~ "PhyloPhlAn MEGAHIT",
                               file == paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT GTDB-Tk/genome_matching.tsv") ~ "GTDB-Tk MEGAHIT",
                               file == paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes PhyloPhlAn 3/genome_matching.tsv") ~ "PhyloPhlAn metaSPAdes",
                               file == paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes GTDB-Tk/genome_matching.tsv") ~ "GTDB-Tk metaSPAdes")
    
    df_append <- data.frame(cur_tool_name,
                            gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                            genomes_matching$abundance_real
    )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
    plot_df <- data.frame(rbind(plot_df, df_append))
  }
  
  if (core_only) {
    plot_df <- plot_df[plot_df$Method %in% c("PhyloPhlAn MEGAHIT", "GTDB-Tk MEGAHIT"),]
  }
  
  level_name <- "species"
  
  max_val_col <- max(aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, plot_df, FUN = mean)$value)
  
  p1 <- plot_df[plot_df$sample_size==default_params["sample_size"] & plot_df$unknown==default_params["unknown"] & plot_df$k==default_params["k"] & plot_df$sigma==default_params["sigma"] & plot_df$mut_rate==default_params["mut_rate"],] %>%
    complete(n,Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    ggplot(aes(x=fct_rev(factor(n)), y=Method, fill=value)) +
    geom_tile(color="gray", linewidth=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max_val_col/2,
                         limit = c(0,max_val_col), 
                         name="Recovered abundance") +
    theme_classic() + 
    xlab("Species count") + 
    coord_flip() +
    theme(text=element_text(size=16),
          plot.margin = unit(c(0,1,0,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
    
  
  p2 <- plot_df[plot_df$n==default_params["n"] & plot_df$unknown==default_params["unknown"] & plot_df$k==default_params["k"] & plot_df$sigma==default_params["sigma"] & plot_df$mut_rate==default_params["mut_rate"],] %>%
    complete(sample_size, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    complete(sample_size, Method) %>%
    ggplot(aes(x=fct_rev(factor(sample_size)), y=Method, fill=value)) +
    geom_tile(color="gray", linewidth=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max_val_col/2,
                         limit = c(0,max_val_col)) +
    theme_classic() + 
    xlab("Sample size (GB)") + 
    coord_flip() +
    theme(text=element_text(size=16),
          plot.margin = unit(c(0,1,0,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p3 <- plot_df[plot_df$n==default_params["n"] & plot_df$sample_size==default_params["sample_size"] & plot_df$k==default_params["k"] & plot_df$sigma==default_params["sigma"] & plot_df$mut_rate==default_params["mut_rate"],] %>%
    complete(unknown, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    ggplot(aes(x=fct_rev(factor(100 * unknown)), y=Method, fill=value)) +
    geom_tile(color="gray", linewidth=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max_val_col/2,
                         limit = c(0,max_val_col)) +
    theme_classic() + 
    xlab("Unknown species (%)") + 
    coord_flip() +
    theme(text=element_text(size=16),
          plot.margin = unit(c(0,1,0,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p4 <- plot_df[plot_df$n==default_params["n"] & plot_df$sample_size==default_params["sample_size"] & plot_df$unknown==default_params["unknown"] & plot_df$k==default_params["k"] & plot_df$sigma==default_params["sigma"],] %>%
    complete(mut_rate, Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    ggplot(aes(x=fct_rev(factor(mut_rate)), y=Method, fill=value)) +
    geom_tile(color="gray", linewidth=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max_val_col/2,
                         limit = c(0,max_val_col), 
                         name=paste0("Percent of\nthe unknown\nabundance\nreconstructed")) +
    theme_classic() + 
    xlab("Mutation rate") + 
    coord_flip() +
    theme(text=element_text(size=16),
          plot.margin = unit(c(0,1,0,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.line=element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    guides(fill = guide_legend(title.position = "left"))
  
  p5 <- plot_df[plot_df$n==default_params["n"] & plot_df$sample_size==default_params["sample_size"] & plot_df$unknown==default_params["unknown"] & plot_df$mut_rate==default_params["mut_rate"] & plot_df$sigma==default_params["sigma"],] %>%
    complete(k,Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>%
    ggplot(aes(x=fct_rev(factor(k)), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max_val_col/2,
                         limit = c(0,max_val_col), 
                         name=paste0("Percent of\nthe unknown\nabundance\nreconstructed")) +
    theme_classic() + 
    xlab("Genome skew") + 
    coord_flip() +
    theme(text=element_text(size=16),
          plot.margin = unit(c(0,1,0,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  p6 <- plot_df[plot_df$n==default_params["n"] & plot_df$sample_size==default_params["sample_size"] & plot_df$unknown==default_params["unknown"] & plot_df$mut_rate==default_params["mut_rate"] & plot_df$k==default_params["k"],] %>%
    complete(sigma,Method) %>% aggregate(value ~ n + sample_size + unknown + k + sigma + mut_rate + Method, FUN = mean) %>% 
    ggplot(aes(x=fct_rev(factor(sigma)), y=Method, fill=value)) +
    geom_tile(color="gray", size=1) + 
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = max_val_col/2,
                         limit = c(0,max_val_col), 
                         name=paste0("Percent of\nthe unknown\nabundance\nreconstructed")) +
    theme_classic() + 
    xlab("Abundance skew") + 
    coord_flip() +
    theme(text=element_text(size=16),
          plot.margin = unit(c(0,1,0,1), "cm"),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line=element_blank()) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(file.path("figures/ReconstructedVsParameter/",dataset,"/merged/"), showWarnings = FALSE)
  if (core_only) {
    g <- arrangeGrob(
      grobs = list(p1, p2, p3, p4),
      heights = c(1, 1.25, 1, 3.1), ncol = 1)
    ggsave(paste0("figures/ReconstructedVsParameter/",dataset,"/merged/",dataset,".png"), g, width=9, height=36, units = 'cm', dpi=1000, bg='#ffffff')
  } else {
    g <- arrangeGrob(
      grobs = list(p1, p2, p3, p6, p5, p4),
      heights = c(1, 1.25, 1, 1, 0.75, 3.1), ncol = 1)
    ggsave(paste0("figures/ReconstructedVsParameter/",dataset,"/merged/",dataset,"_all_tools.png"), g, width=9, height=50, units = 'cm', dpi=1000, bg='#ffffff')
  }
}
make_all_assess_reconstructed <- function() {
  for (dataset in c("soil", "ocean", "gut")) {
    assess_reconstructed(dataset)
    assess_reconstructed(dataset, FALSE)
  }
}
make_all_assess_reconstructed()

find_reconstruction_limit <- function(dataset, core_only = TRUE) {
  plot_df <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(plot_df) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "value")
  
  if (dataset == "gut") {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.5, "mut_rate"=0)
  } else {
    default_params <- c("n"=300, "sample_size"=7.5, "k"=100, "sigma"=1, "unknown"=0.75, "mut_rate"=0)
  }
  
  for (file in c(paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT PhyloPhlAn 3/genome_matching.tsv"),
                 paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT GTDB-Tk/genome_matching.tsv"),
                 paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes PhyloPhlAn 3/genome_matching.tsv"),
                 paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes GTDB-Tk/genome_matching.tsv"))) {
    genomes_matching <- read.csv(file, sep = "\t")
    genomes_matching$SGB <- ifelse(grepl("^GCF_", genomes_matching$SGB), paste0(genomes_matching$SGB, "_genomic"), genomes_matching$SGB)  
    genome_lengths <- read.csv(paste0("simulation_outputs/inputs/", dataset, "/genome_sizes/genome_lengths.tsv"), sep="\t")
    genomes_matching <- left_join(genomes_matching, genome_lengths, by=c("SGB" = "name"))
    
    if (grepl("GTDB-Tk", file)) {
      genomes_matching$reconstructed_tax <- gsub(pattern=" ",replacement="_", genomes_matching$reconstructed_tax) %>% tolower()
      genomes_matching$reconstructed_tax[genomes_matching$reconstructed_tax == "unknown"] <- "UNCLASSIFIED"
      genomes_matching$reconstructed_tax <- gsub("\\|p__\\|.*|\\|c__\\|.*|\\|o__\\|.*|\\|f__\\|.*|\\|g__\\|.*|\\|s__$|", "", genomes_matching$reconstructed_tax)
      
      map <- read.csv(paste0("databases/standardized_databases/GTDBTk.tsv"), sep = "\t")
      map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
      genomes_matching$reconstructed_tax <- mapvalues(genomes_matching$reconstructed_tax, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
    } else {
      genomes_matching$reconstructed_tax <- gsub(pattern=" ",replacement="_", genomes_matching$reconstructed_tax) %>% tolower()
      genomes_matching$reconstructed_tax[grepl("^sgb", genomes_matching$reconstructed_tax)] <- paste0("s__", genomes_matching$reconstructed_tax[grepl("^sgb", genomes_matching$reconstructed_tax)], "_MEGAHIT")
      
      map <- read.csv(paste0("databases/standardized_databases/PhyloPhlAn3.tsv"), sep = "\t")
      map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
      genomes_matching$reconstructed_tax <- mapvalues(genomes_matching$reconstructed_tax, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
    }
    
    keep_only_ncbi_name <- function(x) {
      if (is.na(x)) {
        return (NA)
      }
      tmp <- unlist(str_split(x, "\\|"))
      return (paste0(tmp[!is.na(as.numeric(tmp))], collapse = "|"))
    }
    
    genomes_matching$reconstructed_tax <- sapply(genomes_matching$reconstructed_tax, keep_only_ncbi_name)
    
    reconstructed <- function(x) {
      if (is.na(x[5])) {
        return (FALSE)
      }
      tmp2 <- unlist(str_split(x[5], "\\|"))
      tmp2 <- tmp2[!is.na(as.numeric(tmp2))]
      return (grepl(paste0("^", paste0(tmp2, collapse = "\\|")), x[4]))
    }
    
    genomes_matching$reconstructed <- apply(genomes_matching, 1, reconstructed)
    
    genomes_matching$unknown <- grepl("UNKNOWN", genomes_matching$real_org_id)
    genomes_matching$abundance_reconstructed[is.na(genomes_matching$abundance_reconstructed)] <- 0
    
    # for_norm <- aggregate(abundance_real~sample, data=genomes_matching[genomes_matching$reconstructed & genomes_matching$unknown,], FUN=sum, drop=FALSE)
    # genomes_matching <- aggregate(abundance_real~sample, data=, FUN=sum, drop=FALSE)
    
    genomes_matching <- genomes_matching[genomes_matching$unknown,]
    
    # 0.05 and 0.5 had nothing
    # genomes_matching <- rbind(genomes_matching, data.frame("sample" = paste0("profile.n.300.size.", rep(c("0.05", "0.5"), each=5), ".k.100.sigma.1.0.up.", default_params["unknown"], ".mut_rate.0.0_sample_", rep(0:4, 2)),
    #                                                        "abundance_real" = 0))
    
    
    variable_names = gsub("_sample_.*", "", genomes_matching$sample)
    
    cur_tool_name <- case_when(file == paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT PhyloPhlAn 3/genome_matching.tsv") ~ "PhyloPhlAn MEGAHIT",
                               file == paste0("simulation_outputs/inputs/", dataset, "/MEGAHIT GTDB-Tk/genome_matching.tsv") ~ "GTDB-Tk MEGAHIT",
                               file == paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes PhyloPhlAn 3/genome_matching.tsv") ~ "PhyloPhlAn metaSPAdes",
                               file == paste0("simulation_outputs/inputs/", dataset, "/metaSPAdes GTDB-Tk/genome_matching.tsv") ~ "GTDB-Tk metaSPAdes")
    
    df_append <- data.frame(cur_tool_name,
                            gsub(".*\\.n\\.", "", variable_names) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.size\\.", "", variable_names) %>% gsub(pattern="\\.k.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.k\\.", "", variable_names) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.sigma\\.", "", variable_names) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.up\\.", "", variable_names) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric(),
                            gsub(".*\\.mut_rate\\.", "", variable_names) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric(),
                            genomes_matching$abundance_real,
                            genomes_matching$reconstructed,
                            genomes_matching$length
    )
    colnames(df_append) <- c("Method", "n", "sample_size", "k", "sigma", "unknown", "mut_rate", "abundance_real", "reconstructed", "length")
    
    df_append$coverage <- df_append$sample_size * 1000000000 * df_append$abundance_real / 100 / df_append$length
    
    print(ggplot(df_append[df_append$coverage < 100,], aes(x = reconstructed, y = coverage)) + 
      geom_boxplot())
  }
}








