rm(list = ls())
setwd('..')
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggplotify)

source("scripts/helpers.R")

preprocess_all_real()
preprocess_all_simulated()
preprocess_all_truths()

threshold = 0.05
tool_names = c("Centrifuge", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "Kraken 2 / Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "mOTUs 3", "Metaxa 2", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
colAdd <- c("#FF96D0", "#FF7300", "#FBB15B", "#A15BE4", "#00FFFF", "#0061FE", "#142755", "#81C784", "#2E7D32", "#AC0911", "#E95420")
col = scale_color_manual(values = colAdd)
fil = scale_fill_manual(values = colAdd)

simulation_explanation <- function() {
  unknown_taxa = read.csv('databases/ncbi_taxdump/unknown_taxa.tsv', header = F)[,1]
  plot_df = data.frame(matrix(nrow = 0, ncol = 5))
  colnames(plot_df) <- c("Dataset", "Level", "Type", "Value", "file_name")
  for (dataset in list.files("simulation_outputs/inputs")) {
    unknown_proportion <- ifelse(dataset=="gut", as.character(5), as.character(75))
    in_files = list.files(paste0("simulation_outputs/inputs/", dataset, "/true_profile/profiles"))
    in_files = paste0("simulation_outputs/inputs/", dataset, "/true_profile/profiles/", in_files)
    in_files = in_files[grepl(paste0("simulation_outputs/inputs/", dataset, "/true_profile/profiles/profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9].txt"), in_files)]
    
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
  
  p1 <- ggplot(plot_df, aes(fill=Type, y=Value, x=factor(Level, c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")))) + 
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
          legend.position="bottom")
  
  dir.create(paste0("figures/simulation_setup/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0("figures/simulation_setup/","simulation_setup.png"), p1, width=6, height=18, units = 'cm', dpi=1000)
}
simulation_explanation()

prf1_summary <- function(ncbi_only) {
  plot_df = data.frame(matrix(nrow = 0, ncol = 5))
  colnames(plot_df) <- c("Tool", "Dataset", "Level", "Metric", "Value")
  for (level in c(2, 5, 7)) {
    for (dataset in list.files("simulation_outputs/inputs")) {
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
  
  plot_df_bc = plot_df[plot_df$Metric == "BC",]
  colnames(plot_df_bc)[5] <- "Bray Curtis"
  plot_df_bc$`Bray Curtis` <- as.numeric(plot_df_bc$`Bray Curtis`)
  
  plot_df_f1 = plot_df[plot_df$Metric == "F1",]
  colnames(plot_df_f1)[5] <- "F1"
  plot_df_f1$F1 <- as.numeric(plot_df_f1$F1)
  
  p1 <- ggplot(plot_df_f1, aes(x=Dataset, y=F1, fill=Tool)) + 
    geom_errorbar(aes(ymin = mean - sd,
                      ymax = mean + sd,
                      color=Tool), width=1.5,
                  position=position_dodge(0.75),
                  size=1) +
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0), size=1) +
    facet_wrap(~Level, scale="fixed") + 
    theme_linedraw() + 
    theme(text = element_text(size = 15)) + 
    col + fil
  
  p2 <- ggplot(plot_df_bc, aes(x=Dataset, y=`Bray Curtis`, fill=Tool)) + 
    geom_errorbar(aes(ymin = mean - sd,
                      ymax = mean + sd,
                      color=Tool), width=1.5,
                  position=position_dodge(0.75),
                  size=1) +
    ylab("1 - Bray Curtis dissimilarity") +
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0), size=1) +
    facet_wrap(~Level, scale="fixed") + 
    theme_linedraw() + 
    theme(text = element_text(size = 15)) + 
    col + fil
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  
  ggsave(paste0("figures/simulation_core/f1_", ncbi_only, ".png"), p1, width = 12, height = 4)
  ggsave(paste0("figures/simulation_core/bc_", ncbi_only, ".png"), p2, width = 12, height = 4)
}
prf1_summary(TRUE)

grid_analysis <- function(dataset, level, ncbi_only = TRUE) {
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
  
  unknown_proportion <- ifelse(dataset=="gut", as.character(5), as.character(75))
  for (i in 1:length(profiles)) {
    profiles[[i]] = profiles[[i]][,grepl(paste0("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.[0-9]+\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9]"), colnames(profiles[[i]]))]
  }
  truth = truth[,grepl(paste0("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.[0-9]+\\.sigma\\.1\\.0\\.up\\.0\\.", unknown_proportion, "\\.mut_rate\\.0\\.0_sample_[0-9]"), colnames(truth))]
  
  
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
  rownames(cormat) <- c("Beta Diversity (Mantel r)", "Alpha Diversity (Spearman correlation)")
  melted_cormat <- melt(cormat, na.rm = T)
  
  colnames(signif_mat) <- tool_names
  rownames(signif_mat) <- c("Beta Diversity (Mantel r)", "Alpha Diversity (Spearman correlation)")
  melted_signif <- melt(signif_mat, na.rm = T)
  melted_signif$value <- case_when(melted_signif$value < 0.001 ~ "***",
                                   melted_signif$value < 0.01 ~ "**",
                                   melted_signif$value < 0.05 ~ "*",
                                   TRUE ~ "")
  
  melted_cormat$stars <- melted_signif$value
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  melted_cormat$value <- round(melted_cormat$value, 2)
  ggheatmap <- ggplot(melted_cormat, aes(factor(Var1, c("Beta Diversity (Mantel r)", "Alpha Diversity (Spearman correlation)")), factor(Var2, rev(tool_names)), fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                         limit = c(-1,1), space = "Lab", 
                         name="Spearman correlation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 20, hjust = 1),
          axis.text.y = element_text(angle = 45, vjust = 0, 
                                     size = 20, hjust = 1),
          text = element_text(size=20))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(Var1, c("Beta Diversity (Mantel r)", "Alpha Diversity (Spearman correlation)")), factor(Var2, rev(tool_names)), label = paste0(value, "\n", stars)), color = "black", size = 8) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      text = element_text(size=16))
  #legend.justification = c(1, 0),
  #legend.position = c(0.5, 0.7),
  #legend.direction = "horizontal",
  #legend.title=element_text(size=14)) +
  #guides(fill = guide_colorbar(barwidth = 8, barheight = 1.5,
  #                             title.position = "top", title.hjust = 0.5))
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/simulation_analysis_grid/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/simulation_analysis_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=14, height=42, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_grid_analysis <- function() {
  dataset = "gut"
  level = 7
  grid_analysis(dataset, level, TRUE)
}

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
  
  p1 <- ggplot(f1_df[f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],], 
               aes(x=n, y=value, color=Method)) + geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE) +
    theme_classic() + 
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
               aes(x=sample_size, y=value, color=Method)) + 
    geom_jitter(width = 0.05) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE) +
    theme_classic() + 
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
          legend.position = "none",
          ) + 
    guides(color=guide_legend(title="Tool"))
  
  p3 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$k==default_params["k"] & f1_df$sigma==default_params["sigma"] & f1_df$mut_rate==default_params["mut_rate"],], 
               aes(x=100 * unknown, y=value, color=Method)) + 
    geom_jitter(width = 1) +
    geom_smooth(method = 'glm', formula = 'y~log(100-x+0.001)', size=1, alpha=0.15) +
    theme_classic() + 
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
  
  p6 <- ggplot(f1_df[f1_df$n==default_params["n"] & f1_df$sample_size==default_params["sample_size"] & f1_df$unknown==default_params["unknown"] & f1_df$k==default_params["k"] & f1_df$mut_rate==default_params["mut_rate"],], 
               aes(x=sigma, y=value, color=Method)) + 
    geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE) +
    theme_classic() + 
    ylab("") + 
    xlab("Abundance skew") +
    scale_x_continuous(trans='log10', breaks = c(0.5, 1, 2, 4)) + 
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
    grobs = list(p1, p2, p3, p6),
    widths = c(1, 1, 1, 1), nrow = 1)
  
  dir.create(file.path("figures/F1VsParameter/",dataset,"/merged/"), showWarnings = FALSE)
  ggsave(paste0("figures/F1VsParameter/",dataset,"/merged/",dataset,"_",level_name,"_", ncbi_only, ".png"), g, width=48, height=14, units = 'cm', dpi=1000, bg='#ffffff')
  
}
make_all_f1_vs_parameter_individual <- function() {
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
      f1_vs_parameter_individual(dataset, level)
      f1_vs_parameter_individual(dataset, level, FALSE)
  }
}

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
  
  p1 <- ggplot(bray_df[bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=n, y=value, color=Method)) + geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE) +
    theme_classic() + 
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
  
  p2 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=sample_size, y=value, color=Method)) + 
    geom_jitter(width = 0.05) +
    geom_smooth(method = "loess", span=1.5) +
    theme_classic() + 
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
  
  p3 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$k==default_params["k"] & bray_df$sigma==default_params["sigma"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=100 * unknown, y=value, color=Method)) + 
    geom_jitter(width = 1) +
    geom_smooth(method = 'glm', formula = 'y~log(100-x+0.001)', size=1, alpha=0.15) +
    theme_classic() + 
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
  
  p6 <- ggplot(bray_df[bray_df$n==default_params["n"] & bray_df$sample_size==default_params["sample_size"] & bray_df$unknown==default_params["unknown"] & bray_df$k==default_params["k"] & bray_df$mut_rate==default_params["mut_rate"],], 
               aes(x=sigma, y=value, color=Method)) + 
    geom_jitter(width = 0.01) +
    geom_smooth(method = "glm", 
                method.args = list(family = "binomial"), 
                se = FALSE) +
    theme_classic() + 
    ylab("") + 
    xlab("Abundance skew") +
    scale_x_continuous(trans='log10', breaks = c(0.5, 1, 2, 4)) + 
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
    grobs = list(p1, p2, p3, p6),
    widths = c(1, 1, 1, 1), nrow = 1)
  
  dir.create(file.path("figures/BrayVsParameter/",dataset,"/merged/"), showWarnings = FALSE)
  ggsave(paste0("figures/BrayVsParameter/",dataset,"/merged/",dataset,"_",level_name,"_", ncbi_only, ".png"), g, width=48, height=14, units = 'cm', dpi=1000, bg='#ffffff')
  
}
make_all_bray_vs_parameter_individual <- function() {
  dataset = "soil"
  for (level in c(2,7)) {
    bray_vs_parameter_individual(dataset, level)
  }
}

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
  dat_taxa_species <- dat_taxa_species[rowSums(dat_taxa_species) != 0,]
  distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
  
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
                    dataset == "gator_soil" ~ "Gator Nest Soil",
                    dataset == "human" ~ "Human Gut",
                    dataset == "saltmarsh" ~ "Salt Marsh",
                    dataset == "tara_polar" ~ "Tara Polar",)
  
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
  
  if (dataset != "human") {
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
      theme(legend.position = "none")
    
    ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
    dir.create(file.path("figures/PCoABray/",dataset,"/"), showWarnings = FALSE)
    ggsave(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_", ncbi_only,"_PCoA.png"), width=12, height=12, units = 'cm', dpi=1000)
  } else {
    ggplot(cap, aes(MDS1, MDS2)) + 
      geom_point(aes(color=factor(Method, levels=tool_names)), size = 2, alpha = 0.7) + 
      theme_classic() +
      ggtitle(study) +
      labs( color = "Tool") + 
      labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
           y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = '')) +
      col + guides(color=guide_legend(title="Method")) +
      geom_line(subset(cap,!is.na(Sample)), mapping=aes(MDS1, MDS2, linetype=Sample), alpha=0.5, linewidth=0.5) + 
      guides(linetype = "none")

    ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
    dir.create(file.path("figures/PCoABray/",dataset,"/"), showWarnings = FALSE)
    ggsave(paste0("figures/PCoABray/",dataset,"/",dataset,"_",level_name, "_", ncbi_only,"_PCoA.png"), width=18, height=12, units = 'cm', dpi=1000)
  }
  
}

run_all_PCoABray <- function() {
  for (dataset in c("forest_soil", "cat_gut", "human")) {
    level = 6
    PCoABray(dataset, level, TRUE)
  }
}

grid_iou <- function(dataset, level, ncbi_only) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(3, 11)]
    colAdd <- colAdd[-c(3, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    melted_cormat = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(melted_cormat) <- c("Var1", "Var2")
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    for (dataset_tmp in list.files("real_data_outputs/inputs")[list.files("real_data_outputs/inputs")!="human"]) {
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
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), space = "Lab", 
                         name="Intersection\n over union") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 16, hjust = 1),
          axis.text.y = element_text(size = 16))+
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
      legend.title=element_text(size=16)) +
    guides(fill = guide_colorbar(barwidth = 8, barheight = 1.5,
                                 title.position = "top", title.hjust = 0.5))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/iou/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/iou/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=18, height=18, units = 'cm', dpi=1000, bg='#ffffff')
  
}
grid_iou("all", 6, TRUE)

grid_bc <- function(dataset, level, ncbi_only) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(3, 11)]
    colAdd <- colAdd[-c(3, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    melted_cormat = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(melted_cormat) <- c("Var1", "Var2")
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    for (dataset_tmp in list.files("real_data_outputs/inputs")[list.files("real_data_outputs/inputs")!="human"]) {
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
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.5,
                         limit = c(0,1), space = "Lab", 
                         name="1 - Bray Curtis\ndissimilarity") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 16, hjust = 1),
          axis.text.y = element_text(size = 16))+
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
      legend.title=element_text(size=16)) +
    guides(fill = guide_colorbar(barwidth = 8, barheight = 1.5,
                                 title.position = "top", title.hjust = 0.5))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/bray_grid/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/bray_grid/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=18, height=18, units = 'cm', dpi=1000, bg='#ffffff')
  
}
grid_bc("all", 6, TRUE)

grid_bc_mantel <- function(dataset, level, ncbi_only) {
  if (level == 7) {
    tool_names <- c(tool_names, "Simka")
  }
  
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(3, 11)]
    colAdd <- colAdd[-c(3, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    melted_cormat = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(melted_cormat) <- c("Var1", "Var2")
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    for (dataset_tmp in list.files("real_data_outputs/inputs")[list.files("real_data_outputs/inputs")!="human"]) {
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
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                         limit = c(-1,1), space = "Lab", 
                         name="Mantel r") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     size = 16, hjust = 1),
          axis.text.y = element_text(size = 16))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(Var2, tool_names), factor(Var1, tool_names), label = paste0(value, "\n", stars)), color = "black", size = 4) +
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
      legend.title=element_text(size=16)) +
    guides(fill = guide_colorbar(barwidth = 8, barheight = 1.5,
                                 title.position = "top", title.hjust = 0.5))
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  dir.create(file.path("figures/beta_mantel/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/beta_mantel/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, ".png"), width=18, height=18, units = 'cm', dpi=1000, bg='#ffffff')
  
}
grid_bc_mantel("all", 6, TRUE)

unknownECDF <- function(dataset, level) {
  if (!(dataset %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh", "all"))) {
    tool_names <- tool_names[-c(3, 11)]
    colAdd <- colAdd[-c(3, 11)]
    col = scale_color_manual(values = colAdd)
    fil = scale_fill_manual(values = colAdd)
  }
  set.seed(1)
  
  if (dataset == "all") {
    df <- data.frame(matrix(ncol = 3, nrow = 0))
    for (dataset_tmp in list.files("real_data_outputs/inputs")[list.files("real_data_outputs/inputs")!="human"]) {
      if (!(dataset_tmp %in% c("acid_mine", "animal_gut", "cat_gut", "coastal_sediment", "dog_gut", "forest_soil", "human", "saltmarsh"))) {
        tool_names_tmp <- tool_names[-c(3, 11)]
      } else {
        tool_names_tmp <- tool_names
      }
      
      data <- preprocess(dataset_tmp, level)
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
    data <- preprocess(dataset, level)
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
                    dataset == "gator_soil" ~ "gator nest soil",
                    dataset == "human" ~ "human gut",
                    dataset == "saltmarsh" ~ "salt marsh",
                    dataset == "tara_polar" ~ "tara polar",
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
  
  ggplot(df, aes(x=value, color=Methods)) + 
    stat_ecdf(geom = "step", size = 1) +
    theme_classic() + 
    ylab("Empirical cumulative density") + xlab(paste0("Percent unknown at the ", level_name, " level")) +
    col +
    guides(color=guide_legend(title="Method"))
  dir.create(file.path("figures/unknownECDF/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/unknownECDF/",dataset,"/",dataset,"_",level_name,"_ECDFUnknowns.png"), width=14.5, height=10, units = 'cm', dpi=1000)
}
make_all_unknownECDF_plots <- function() {
  dataset = "all"
  for (level in c(2, 5, 7)) {
    unknownECDF(dataset, level)
  }
}
make_all_unknownECDF_plots()

PERMANOVA_comparison <- function(level, ncbi_only) {
  meta = read.csv("metadata/forest_soil/merged_metadata.tsv", sep="\t")
  dataset <- "forest_soil"
  
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
  
  place_holders = rownames(profiles[[1]])
  for (i in 1:length(profiles)) {
    rownames(profiles[[i]]) <- mapvalues(rownames(profiles[[i]]), place_holders, sample_ids, warn_missing = F)
    profiles[[i]] <- profiles[[i]][rownames(profiles[[i]]) %in% meta$sample_id,]
    profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
  }
  
  df = data.frame(matrix(nrow = 0, ncol = 4))
  colnames(df) <- c("method", "variable", "R2", "p_val")
  for (i in 1:length(profiles)) {
    meta_wo_na <- meta
    rownames(meta_wo_na) <- meta$sample_id
    meta_wo_na$sample_id <- NULL
    meta_wo_na$pH <- as.numeric(meta_wo_na$pH)
    meta_wo_na <- meta_wo_na[rownames(meta_wo_na) %in% rownames(profiles[[i]]), ]
    
    # Only single copy of data for everything else
    adonis_res_rsq = vector()
    adonis_res_pval = vector()
    for (col in names(meta_wo_na)){
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
    df_addition = data.frame("method"=i, "variable"=names(adonis_res_pval), "R2"=adonis_res_rsq, "p_val"=adonis_res_pval)
    df <- rbind(df, df_addition)
  }
  
  df$method <- tool_names[df$method]
  df$variable <- case_when(df$variable == "soil_layer" ~ "Soil layer",
                           df$variable == "location" ~ "Location",
                           df$variable == "pH" ~ "pH")
  df$stars <- case_when(df$p_val < 0.001 ~ "***",
                        df$p_val < 0.01 ~ "**",
                        df$p_val < 0.05 ~ "*",
                        TRUE ~ "")
  
  ggheatmap <- ggplot(df, aes(factor(method, tool_names), variable, fill = R2))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", mid = "#FFBBBB", midpoint = 0.2,
                         limit = c(0,0.4), space = "Lab", 
                         name="PERMANOVA\nUnivariate R-squared") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(factor(method, tool_names), variable, label = paste0(round(R2,2), "\n", stars)), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
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
  dir.create(file.path("figures/PERMANOVA_comparison/",dataset,"/"), showWarnings = FALSE)
  ggsave(paste0("figures/PERMANOVA_comparison/",dataset,"/",dataset,"_",level_name,"_", ncbi_only, "_PERMANOVA.png"), width=24, height=9, units = 'cm', dpi=1000, bg='#ffffff')
}
make_all_PERMANOVA_comparison <- function() {
  level = 7
  PERMANOVA_comparison(level, FALSE)
}
make_all_PERMANOVA_comparison()

PERMANOVA_example <- function(level, ncbi_only) {
  meta = read.csv("metadata/forest_soil/merged_metadata.tsv", sep="\t")
  dataset <- "forest_soil"
  
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
    rownames(tmp) <- paste0(rownames(tmp), "_", i)
    profiles[[i]] <- tmp
  }
  
  sample_ids = read.csv(paste0("real_data_outputs/inputs/", dataset, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)[1,] %>% 
    gsub(pattern="_taxonomic_profile", replacement = "") %>% gsub(pattern="^re", replacement="")
  
  sample_ids <- sort(sample_ids[-1])
  
  meta <- meta[meta$sample_id %in% sample_ids,]
  meta <- meta[order(meta$sample_id),]
  
  place_holders = gsub("_[0-9]+$", "", rownames(profiles[[1]]))
  for (i in 1:length(profiles)) {
    rownames(profiles[[i]]) <- paste0(mapvalues(gsub("_[0-9]+$", "", rownames(profiles[[i]])), place_holders, sample_ids, warn_missing = F), "_", i)
    profiles[[i]] <- profiles[[i]][gsub("_.*", "", rownames(profiles[[i]])) %in% meta$sample_id,]
    profiles[[i]] <- profiles[[i]][rowSums(profiles[[i]]) != 0,]
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
  
  meta_wo_na = stat_meta[, which(names(stat_meta) %in% c("pH", "location", "soil_layer", "method"))]
  meta_wo_na$pH <- as.numeric(meta_wo_na$pH)
  rownames(meta_wo_na) <- stat_meta$sample_id
  
  meta_wo_na <- meta_wo_na[rownames(meta_wo_na) %in% rownames(dat_taxa_species), ]
  meta_wo_na$method_type <- case_when(meta_wo_na$method=="Centrifuge" | meta_wo_na$method=="Kraken 2 / Bracken 2" ~ "kmer",
                                      grepl("MetaPhlAn|mOTUs|Metaxa", meta_wo_na$method) ~ "marker",
                                      TRUE ~ "assembly")
  
  # All data for methods-based
  adonis_res_rsq = vector()
  adonis_res_pval = vector()
  for (col in c("method", "method_type")) {
    adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na)
    adonis_res_rsq[col] = adonis.univ[1,"R2"]
    adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
  }
  
  # Only single copy of data for everything else
  meta_wo_na <- meta_wo_na[meta_wo_na$method=="MetaPhlAn 4",-which(colnames(meta_wo_na) %in% c("method", "method_type"))]
  new_dat_taxa_species <- profiles[[7]]
  for (col in names(meta_wo_na)){
    
    bray <- vegdist(new_dat_taxa_species[!is.na(meta_wo_na[,col]),])
    adonis.univ = adonis2(as.formula(paste("bray ~ ", col)), data = meta_wo_na[!is.na(meta_wo_na[,col]),], permutations = 9999)
    adonis_res_rsq[col] = adonis.univ[1,"R2"]
    adonis_res_pval[col] = adonis.univ[1,"Pr(>F)"]
  }
  
  # Format adonis output for plotting
  univar_res_tax = data.frame(adonis_res_pval, adonis_res_rsq)
  colnames(univar_res_tax) = c("P-Value", "R2")
  univar_res_tax$`P-Value` = as.numeric(univar_res_tax$`P-Value`)
  univar_res_tax$R2 = as.numeric(univar_res_tax$R2)
  univar_res_tax$variables <- c("Method", "Method type", "Soil layer", "Location", "pH")
  rownames(univar_res_tax) <- univar_res_tax$variables
  univar_res_tax$stars = cut(univar_res_tax$`P-Value`, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", ".", ""))
  
  univar_res_tax$`P-Value` = round(univar_res_tax$`P-Value`, 3)
  univar_res_tax$R2 = univar_res_tax$R2 *100
  dodge = position_dodge(width = 0.8)
  
  # Plot
  ggplot(data = univar_res_tax, aes(reorder(variables, R2), y = R2, label = stars)) + 
    geom_bar(stat = "identity", position = dodge, fill="brown") + 
    geom_text(position = dodge, vjust = 0.8, hjust = -0.1, size = 4) + 
    theme_classic(base_size = 10) + ylab("PERMANOVA Univarate R-squared") + 
    coord_flip() + ylim(0, 50) + xlab("") + 
    labs(fill = "Tool") + 
    guides(fill = guide_legend(reverse=T)) + 
    theme(text = element_text(size = 22))
  
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
  level <- 7
  PERMANOVA_example(level, FALSE)
}
make_all_PERMANOVA_example()


