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
library(ROCR)
library(data.table)
library(ggbeeswarm)
library(patchwork)

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
dir.create(paste0("analysis/figure_aitchison/thresholding2/"), showWarnings = FALSE, recursive = T)
write.table(out_df, "analysis/figure_aitchison/thresholding2/thresholds.tsv", sep="\t", row.names = F)

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
  level_list <- fread(paste0("analysis/databases/ncbi_taxdump/names.dmp"), sep="\t", header = F, quote="")
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
    xlab("Simulated dataset") + 
    ggtitle("Average abundance\nin core samples")
  
  ggsave(paste0("analysis/figure_aitchison/overview_figure2/sim_phylum.png"), width=15, height=13.5, units = 'cm', dpi=1000)
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
  
  dir.create(paste0("analysis/figure_aitchison/simulation_setup/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figure_aitchison/simulation_setup/","simulation_setup.png"), p1, width=6, height=18, units = 'cm', dpi=1000)
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
  
  # In-text
  plot_df <- aggregate(Value~Dataset + Level + Type + Tool, data=plot_df, FUN=mean)
  print(max(plot_df[plot_df$Level != "superkingdom" & plot_df$Type == "NCBI available but not in database",]$Value))
  print(mean(plot_df[plot_df$Level != "superkingdom" & plot_df$Type == "NCBI available but not in database",]$Value))
  
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
  
  dir.create(paste0("analysis/figure_aitchison/simulation_setup/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("analysis/figure_aitchison/simulation_setup/","simulation_setup_by_tool.png"), p1, width=44, height=18, units = 'cm', dpi=1000)
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
      ## grab the truth simulated profile.
      truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
      ## grab the profiles from each tool, and re-normalize them
      profiles = lapply(profiles, renormalize)
      ## apply a threshold if there is one, threshold is a global variable
      profiles = lapply(profiles, threshold_sample, threshold)
      ## re-nomarlize after applying the threshold
      profiles = lapply(profiles, renormalize)
      
      ## table that collects the performance data
      prf1_table <- matrix(nrow = 3, ncol = length(profiles))
      for (i in 1:length(profiles)) {
        ##grab taxa names from the profile of interest
        taxa_names <- list_taxa_by_sample(profiles[[i]], threshold)
        ##grab the taxa names of the truth
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
        
        df_addition <- data.frame(cbind(tool_names[i], dataset, level, "BC", calc_bc_sim(profiles[[i]], truth, dist_type = "robust.aitchison")[standards]))
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
  
  # In-text
  print(aggregate(Value ~ Dataset + Level + Metric, plot_df[plot_df$Level == 'Phylum',], mean))
  plot_df$Tool_type <- ifelse(grepl("GTDB|PhyloPhlAn", plot_df$Tool), "Assembly", "Reference")
  print(aggregate(Value ~ Level + Metric + Tool_type, plot_df[plot_df$Level == 'Family',], mean))
  print(aggregate(Value ~ Level + Metric + Tool, plot_df[plot_df$Level == 'Species',], mean))
  print(aggregate(Value ~ Level + Metric + Tool, plot_df[plot_df$Level == 'Species' & plot_df$Tool_type == 'Assembly',], mean))

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
  
  plot_df$Metric[plot_df$Metric == "BC"] <- "Robust Aitchison"
  
  
  
  p1 <- plot_df %>% filter(Metric=="F1") %>% 
    ggplot(aes(x=Tool, y=Value, fill=Tool)) + 
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0), size=3, alpha=0.7) +
    facet_nested(Dataset ~ Level) +
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
    col + fil +
    ggtitle("F1")
  
  
  p2 <- plot_df %>% filter(Metric=="Robust Aitchison") %>% 
    ggplot(aes(x=Tool, y=Value, fill=Tool)) + 
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0), size=3, alpha=0.7) +
    facet_nested(Dataset ~ Level) +
    theme_linedraw() + 
    #ylim(0, 1) + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 12, hjust = 1),
          legend.position = "none",
          strip.text = element_text(colour = 'black'),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    col + fil +
    ggtitle("Robust Aitchison")
  
  
  merged <- p1 + p2

  #rename variable for output file name
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  core_five <- ifelse(core_five, "core_tools", "all_tools")
  
  ggsave(paste0("analysis/figure_aitchison/simulation_core/f1_and_jc_", ncbi_only, "_", core_five, ".png"), merged, width = 10, height = 6)
}

if (!file.exists(paste0("analysis/figure_aitchison/simulation_core/f1_and_jc_", "NCBI_only", "_core_tools.png"))) {
  prf1_summary(TRUE)
  prf1_summary(TRUE, TRUE)
}
if (!file.exists(paste0("analysis/figure_aitchison/simulation_core/f1_and_jc_", "all_taxa", "_all_tools.png"))) {
  prf1_summary(FALSE)
  prf1_summary(FALSE, TRUE)
}

