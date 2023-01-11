rm(list = ls())
setwd('..')

source("scripts/helpers.R")

preprocess_all_simulated()
preprocess_all_truths()

# Outputs dataframe of (number of tools) x (number of thresholds)
# Using NCBI names only
test_abundance_thresholds <- function(thresholds) {
  level = 7
  out_mat = matrix(nrow = length(thresholds), ncol = 11)
  k = 1
  for (threshold in thresholds) {
    f1s_avg <- matrix(nrow = length(list.files("simulation_outputs/inputs")), ncol = 11)
    j = 1
    for (dataset in list.files("simulation_outputs/inputs")) {
      taxa = remove_non_ncbi(preprocess(dataset, level, FALSE))
      truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset, level)))
      taxa <- lapply(taxa, renormalize)
      f1s = vector(length = length(taxa))
      for (i in 1:length(taxa)) {
        taxa_names <- list_taxa_by_sample(taxa[[i]], threshold)
        truth_names <- list_taxa_by_sample(truth, 0)
        f1s[i] <- mean(calculate_f1(taxa_names, truth_names))
      }
      f1s_avg[j,] <- f1s
      j = j + 1
    }
    out_mat[k,] <- colMeans(f1s_avg)
    k = k + 1
  }
  rownames(out_mat) <- thresholds
  colnames(out_mat) <- c("Centrifuge", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "Kraken 2 / Braken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "mOTUs 3", "Metaxa 2", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
  return(out_mat)
}

# Save optimal thresholds; use 0.03 as threshold throughout
out_df <- data.frame(test_abundance_thresholds(c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1, 1)), check.rows = F)
out_df <- data.frame(cbind("Thresholds" = rownames(out_df), out_df), check.rows = F)
write.table(out_df, "figures/thresholding/thresholds.tsv", sep="\t", row.names = F)

threshold = 0.08
tool_names = c("Centrifuge", "GTDB-Tk MEGAHIT", "GTDB-Tk metaSPAdes", "Kraken 2 / Braken 2", "MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "mOTUs 3", "Metaxa 2", "PhyloPhlAn MEGAHIT", "PhyloPhlAn metaSPAdes")
colAdd <- c("#FF96D0", "#FF7300", "#FBB15B", "#A15BE4", "#00FFFF", "#00CCFF", "#0061FE", "#81C784", "#2E7D32", "#AC0911", "#E95420")
col = scale_color_manual(values = colAdd)
fil = scale_fill_manual(values = colAdd)

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
        
        sample_names = names(truth)[2:length(truth)]
        if (dataset %in% c("soil", "ocean")) {
          standards = grepl("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.75\\.mut_rate\\.0\\.0_", sample_names)
        } else {
          standards = grepl("profile\\.n\\.300\\.size\\.7\\.5\\.k\\.100\\.sigma\\.1\\.0\\.up\\.0\\.5\\.mut_rate\\.0\\.0_", sample_names)
        }
        
        df_addition <- data.frame(cbind(tool_names[i], dataset, level, "F1", calc_precision_recall_f1(taxa_names, truth_names)[3,standards]))
        colnames(df_addition) <- c("Tool", "Dataset", "Level", "Metric", "Value")
        plot_df <- rbind(plot_df, df_addition)
        
        df_addition <- data.frame(cbind(tool_names[i], dataset, level, "BC", calc_bc_sim(profiles[[i]], truth)[standards]))
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
  plot_df$Level <- factor(plot_df$Level, levels = c("Phylum", "Genus", "Species"))
  
  plot_df$Value <- as.numeric(plot_df$Value)

  plot_df2 <- plot_df %>%
    group_by(Dataset, Level, Tool, Metric) %>%
    summarise_at(vars(Value), list(mean=mean, sd=sd)) %>% 
    as.data.frame()
  
  plot_df <- left_join(plot_df, plot_df2)
  
  plot_df_bc = plot_df[plot_df$Metric == "BC",]
  colnames(plot_df_bc)[5] <- "Bray-Curtis"
  plot_df_bc$`Bray-Curtis` <- as.numeric(plot_df_bc$`Bray-Curtis`)
  
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
  
  p2 <- ggplot(plot_df_bc, aes(x=Dataset, y=`Bray-Curtis`, fill=Tool)) + 
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
  
  ncbi_only <- ifelse(ncbi_only, "NCBI_only", "all_taxa")
  
  ggsave(paste0("figures/simulation_core/f1_", ncbi_only, ".png"), p1, width = 12, height = 6)
  ggsave(paste0("figures/simulation_core/bc_", ncbi_only, ".png"), p2, width = 12, height = 6)
}

prf1_summary(TRUE)
prf1_summary(FALSE)

precision_recall_f1 <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    # Fix one taxa slipping through
    data[[i]][data[[i]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[i]]))] <- 0
    data[[i]] <- data[[i]][data[[i]]$TaxIDs!="UNCLASSIFIED",]
    data[[i]][,-1][data[[i]][,-1]<0.0001]=0
    tmp <- data[[i]]$TaxIDs
    data[[i]] <- data.frame(lapply(data[[i]], function(x){ifelse(x==0,NA, tmp)}))
    data[[i]]$TaxIDs <- NULL
  }
  
  pairs <- rbind(1:(length(data) - 1), length(data))
  
  precision_recall_f1 <- function(x, y) {
    x_set = x[!is.na(x)]
    y_set = y[!is.na(y)]
    
    precision <- length(intersect(x_set,y_set))/length(x_set)
    recall <- length(intersect(x_set,y_set))/length(y_set)
    f1 <- 2 * precision * recall / (precision + recall)
    
    return(c(precision, recall, f1))
  }
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(data[[1]][,-1]) * ncol(pairs), ncol = 5))
  for (i in 1:ncol(pairs)) {
    output <- t(mapply(precision_recall_f1, data[[pairs[1,i]]][,-1], data[[pairs[2,i]]][,-1]))
    pairwise_diffs[((i-1) * ncol(data[[1]][,-1]) + 1):(i * ncol(data[[1]][,-1])),1] <- i
    pairwise_diffs[((i-1) * ncol(data[[1]][,-1]) + 1):(i * ncol(data[[1]][,-1])),2] <- rownames(output)
    pairwise_diffs[((i-1) * ncol(data[[1]][,-1]) + 1):(i * ncol(data[[1]][,-1])),3:5] <- output
  }
  
  colnames(pairwise_diffs) <- c("Method", "sample", "precision", "recall", "f1")
  
  pairwise_diffs$Method <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MetaPhlAn2",
                                          pairwise_diffs[,1]==2 ~ "MetaPhlAn3",
                                          pairwise_diffs[,1]==3 ~ "MetaPhlAn4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "PhyloPhlAn3",
                                          pairwise_diffs[,1]==6 ~ "Kraken2",
                                          pairwise_diffs[,1]==7 ~ "GTDBTk",
                                          pairwise_diffs[,1]==8 ~ "Metaxa2"))
  
  pairwise_diffs$n = gsub(".*\\.n\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric()
  pairwise_diffs$sample_size = gsub(".*\\.size\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  pairwise_diffs$k = gsub(".*\\.k\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  pairwise_diffs$sigma = gsub(".*\\.sigma\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric()
  pairwise_diffs$unknown_prop = gsub(".*\\.up\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric()
  pairwise_diffs$mut_rate = gsub(".*\\.mut_rate\\.", "", pairwise_diffs$sample) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()
  
  # Reformat Bray Curtis matrix
  
  # Swap study for plotting
  study <- case_when(dataset=="Soil" ~ "Soil Simulations")
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  colAdd <- c("#00ffff", "#00ccff", "#0066cc", "#81C784", "#2E7D32", "#7B1FA2", "#FF8F00", "#FBC02D", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(unique(pairwise_diffs$Method))])
  
  # Recall
  p1 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=n, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab("Recall") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  #dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/n/",dataset,"_",level_name,"_RecallVsSpeciesNumber.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p2 <- ggplot(pairwise_diffs[pairwise_diffs$n==80 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=sample_size, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    col + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    scale_x_continuous(trans='log10') + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/sample_size/",dataset,"_",level_name,"_RecallVsSampleSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p3 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=unknown_prop, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    col + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  
  p4 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], 
               aes(x=mut_rate, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)  
  # ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/unknown_prop/",dataset,"_",level_name,"_RecallVsUnknown.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n == 80 & pairwise_diffs$sigma==1,], aes(x=k, y=recall, color=Method)) + geom_point() +
  #   geom_smooth(method = "lm", formula = 'y ~ x') +
  #   theme_classic() + ggtitle("Recall") +
  #   ylab("Recall") + xlab("Genome size distribution (lower k is more skewed)")
  # dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/k/",dataset,"_",level_name,"_RecallVsGenomeSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$n==80,], aes(x=sigma, y=recall, color=Method)) + geom_point() +
  #   geom_smooth(method = "lm", formula = 'y ~ x') +
  #   theme_classic() + ggtitle("Recall") +
  #   ylab("Recall") + xlab("Abundance skew (larger sigma is more skewed)")
  # dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/sigma/",dataset,"_",level_name,"_RecallVsSigma.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # Precision
  p5 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=n, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab("Precision") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/n/",dataset,"_",level_name,"_PrecisionVsSpeciesNumber.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p6 <- ggplot(pairwise_diffs[pairwise_diffs$n==80 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=sample_size, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    col + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    scale_x_continuous(trans='log10') + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/sample_size/",dataset,"_",level_name,"_PrecisionVsSampleSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p7 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=unknown_prop, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    col + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  
  p8 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], 
               aes(x=mut_rate, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)  
  # ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/unknown_prop/",dataset,"_",level_name,"_PrecisionVsUnknown.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n == 80 & pairwise_diffs$sigma==1,], aes(x=k, y=precision, color=Method)) + geom_point() +
  #   geom_smooth(method = "lm", formula = 'y ~ x') +
  #   theme_classic() + ggtitle("Precision") +
  #   ylab("Precision") + xlab("Genome size distribution (lower k is more skewed)")
  # dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/k/",dataset,"_",level_name,"_PrecisionVsGenomeSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$n==80,], aes(x=sigma, y=precision, color=Method)) + geom_point() +
  #   geom_smooth(method = "lm", formula = 'y ~ x') +
  #   theme_classic() + ggtitle("Precision") +
  #   ylab("Precision") + xlab("Abundance skew (larger sigma is more skewed)")
  # dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/sigma/",dataset,"_",level_name,"_PrecisionVsSigma.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # F1
  p9 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=n, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab("F1") + 
    xlab("Number of species") +
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/n/",dataset,"_",level_name,"_F1VsSpeciesNumber.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p10 <- ggplot(pairwise_diffs[pairwise_diffs$n==80 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=sample_size, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    xlab("Sample size (GB)") + 
    col + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.y=element_blank()) + 
    scale_x_continuous(trans='log10') + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/sample_size/",dataset,"_",level_name,"_F1VsSampleSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p11 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=unknown_prop, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    xlab("Abundance unknown") +
    col + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"), 
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  
  p12 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], 
               aes(x=mut_rate, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    xlab("Insertion/deletion rate") +
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=100),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.y=element_blank()) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  # dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)  
  # ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/unknown_prop/",dataset,"_",level_name,"_F1VsUnknown.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n == 80 & pairwise_diffs$sigma==1,], aes(x=k, y=f1, color=Method)) + geom_point() +
  #   geom_smooth(method = "lm", formula = 'y ~ x') +
  #   theme_classic() + ggtitle("F1") +
  #   ylab("F1") + xlab("Genome size distribution (lower k is more skewed)")
  # dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/k/",dataset,"_",level_name,"_F1VsGenomeSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$n==80,], aes(x=sigma, y=f1, color=Method)) + geom_point() +
  #   geom_smooth(method = "lm", formula = 'y ~ x') +
  #   theme_classic() + ggtitle("F1") +
  #   ylab("F1") + xlab("Abundance skew (larger sigma is more skewed)")
  # dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  # ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/sigma/",dataset,"_",level_name,"_F1VsSigma.png"), width=12, height=10, units = 'cm', dpi=1000)

  ggexport(ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow=3, ncol=4,common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/",dataset,"_",level_name,"_All_Combined_Precision_Recall_F1_Vs_Parameters.png"), width = 5200, height = 3600)

}
make_all_precision_recall_f1 <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      precision_recall_f1(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
}

precision_recall_f1_individual <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    # Fix one taxa slipping through
    data[[i]][data[[i]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[i]]))] <- 0
    data[[i]] <- data[[i]][data[[i]]$TaxIDs!="UNCLASSIFIED",]
    data[[i]][,-1][data[[i]][,-1]<0.0001]=0
    tmp <- data[[i]]$TaxIDs
    data[[i]] <- data.frame(lapply(data[[i]], function(x){ifelse(x==0,NA, tmp)}))
    data[[i]]$TaxIDs <- NULL
  }
  
  pairs <- rbind(1:(length(data) - 1), length(data))
  
  precision_recall_f1 <- function(x, y) {
    x_set = x[!is.na(x)]
    y_set = y[!is.na(y)]
    
    precision <- length(intersect(x_set,y_set))/length(x_set)
    recall <- length(intersect(x_set,y_set))/length(y_set)
    f1 <- 2 * precision * recall / (precision + recall)
    
    return(c(precision, recall, f1))
  }
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(data[[1]][,-1]) * ncol(pairs), ncol = 5))
  for (i in 1:ncol(pairs)) {
    output <- t(mapply(precision_recall_f1, data[[pairs[1,i]]][,-1], data[[pairs[2,i]]][,-1]))
    pairwise_diffs[((i-1) * ncol(data[[1]][,-1]) + 1):(i * ncol(data[[1]][,-1])),1] <- i
    pairwise_diffs[((i-1) * ncol(data[[1]][,-1]) + 1):(i * ncol(data[[1]][,-1])),2] <- rownames(output)
    pairwise_diffs[((i-1) * ncol(data[[1]][,-1]) + 1):(i * ncol(data[[1]][,-1])),3:5] <- output
  }
  
  colnames(pairwise_diffs) <- c("Method", "sample", "precision", "recall", "f1")
  
  pairwise_diffs$Method <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MetaPhlAn2",
                                            pairwise_diffs[,1]==2 ~ "MetaPhlAn3",
                                            pairwise_diffs[,1]==3 ~ "MetaPhlAn4",
                                            pairwise_diffs[,1]==4 ~ "mOTUs3",
                                            pairwise_diffs[,1]==5 ~ "PhyloPhlAn3",
                                            pairwise_diffs[,1]==6 ~ "Kraken2",
                                            pairwise_diffs[,1]==7 ~ "GTDBTk",
                                            pairwise_diffs[,1]==8 ~ "Metaxa2"))
  
  pairwise_diffs$n = gsub(".*\\.n\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric()
  pairwise_diffs$sample_size = gsub(".*\\.size\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  pairwise_diffs$k = gsub(".*\\.k\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  pairwise_diffs$sigma = gsub(".*\\.sigma\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric()
  pairwise_diffs$unknown_prop = gsub(".*\\.up\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric()
  pairwise_diffs$mut_rate = gsub(".*\\.mut_rate\\.", "", pairwise_diffs$sample) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()
  
  # Reformat Bray Curtis matrix
  
  # Swap study for plotting
  study <- case_when(dataset=="Soil" ~ "Soil Simulations")
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  colAdd <- c("#00ffff", "#00ccff", "#0066cc", "#81C784", "#2E7D32", "#7B1FA2", "#FF8F00", "#FBC02D", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(unique(pairwise_diffs$Method))])
  
  # Recall
  p1 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=n, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Recall at the\n", level_name, " level")) + 
    xlab("Number of species") +
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/n/",dataset,"_",level_name,"_RecallVsSpeciesNumber.png"), p1, width=20, height=20, units = 'cm', dpi=1000)
  
  p2 <- ggplot(pairwise_diffs[pairwise_diffs$n==80 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=sample_size, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Recall at the\n", level_name, " level")) + 
    xlab("Sample size (GB)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    scale_x_continuous(trans='log10') + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/sample_size/",dataset,"_",level_name,"_RecallVsSampleSize.png"), p2, width=20, height=20, units = 'cm', dpi=1000)
  
  p3 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
               aes(x=unknown_prop, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Recall at the\n", level_name, " level")) + 
    xlab("Proportion unknown") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/unknown_prop/",dataset,"_",level_name,"_RecallVsUnknown.png"), p3, width=20, height=20, units = 'cm', dpi=1000)
  
  p4 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], 
               aes(x=mut_rate, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Recall at the\n", level_name, " level")) + 
    xlab("Mutation rate") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/mut_rate/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/mut_rate/",dataset,"_",level_name,"_RecallVsMutRate.png"), p4, width=20, height=20, units = 'cm', dpi=1000)
  
  p5 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n == 80 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=k, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Recall at the\n", level_name, " level")) + 
    xlab("Genome size distribution\n(lower k is more skewed)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/k/",dataset,"_",level_name,"_RecallVsGenomeSize.png"), p5, width=20, height=20, units = 'cm', dpi=1000)
  
  p6 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$n==80 & pairwise_diffs$mut_rate==0,], 
         aes(x=sigma, y=recall, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Recall at the\n", level_name, " level")) + 
    xlab("Abundance skew\n(larger sigma is more skewed)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Recall/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/sigma/",dataset,"_",level_name,"_RecallVsSigma.png"), p6, width=20, height=20, units = 'cm', dpi=1000)
  
  # Precision
  p7 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=n, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Precision at the\n", level_name, " level")) + 
    xlab("Number of species") +
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/n/",dataset,"_",level_name,"_PrecisionVsSpeciesNumber.png"), p7, width=20, height=20, units = 'cm', dpi=1000)
  
  p8 <- ggplot(pairwise_diffs[pairwise_diffs$n==80 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=sample_size, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Precision at the\n", level_name, " level")) + 
    xlab("Sample size (GB)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    scale_x_continuous(trans='log10') + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/sample_size/",dataset,"_",level_name,"_PrecisionVsSampleSize.png"), p8, width=20, height=20, units = 'cm', dpi=1000)
  
  p9 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=unknown_prop, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Precision at the\n", level_name, " level")) + 
    xlab("Proportion unknown") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/unknown_prop/",dataset,"_",level_name,"_PrecisionVsUnknown.png"), p9, width=20, height=20, units = 'cm', dpi=1000)
  
  p10 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], 
         aes(x=mut_rate, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Precision at the\n", level_name, " level")) + 
    xlab("Mutation rate") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/mut_rate/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/mut_rate/",dataset,"_",level_name,"_PrecisionVsMutRate.png"), p10, width=20, height=20, units = 'cm', dpi=1000)
  
  p11 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n == 80 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=k, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Precision at the\n", level_name, " level")) + 
    xlab("Genome size distribution\n(lower k is more skewed)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/k/",dataset,"_",level_name,"_PrecisionVsGenomeSize.png"), p11, width=20, height=20, units = 'cm', dpi=1000)
  
  p12 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$n==80 & pairwise_diffs$mut_rate==0,], 
         aes(x=sigma, y=precision, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("Precision at the\n", level_name, " level")) + 
    xlab("Abundance skew\n(larger sigma is more skewed)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/Precision/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/sigma/",dataset,"_",level_name,"_PrecisionVsSigma.png"), p12, width=20, height=20, units = 'cm', dpi=1000)
  
  # F1
  p13 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=n, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("F1 at the\n", level_name, " level")) + 
    xlab("Number of species") +
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/n/",dataset,"_",level_name,"_F1VsSpeciesNumber.png"), p13, width=20, height=20, units = 'cm', dpi=1000)
  
  p14 <- ggplot(pairwise_diffs[pairwise_diffs$n==80 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=sample_size, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("F1 at the\n", level_name, " level")) + 
    xlab("Sample size (GB)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    scale_x_continuous(trans='log10') + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/sample_size/",dataset,"_",level_name,"_F1VsSampleSize.png"), p14, width=20, height=20, units = 'cm', dpi=1000)
  
  p15 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=unknown_prop, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("F1 at the\n", level_name, " level")) + 
    xlab("Proportion unknown") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/unknown_prop/",dataset,"_",level_name,"_F1VsUnknown.png"), p15, width=20, height=20, units = 'cm', dpi=1000)
  
  p16 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], 
         aes(x=mut_rate, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("F1 at the\n", level_name, " level")) + 
    xlab("Mutation rate") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/mut_rate/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/mut_rate/",dataset,"_",level_name,"_F1VsMutRate.png"), p16, width=20, height=20, units = 'cm', dpi=1000)
  
  p17 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n == 80 & pairwise_diffs$sigma==1 & pairwise_diffs$mut_rate==0,], 
         aes(x=k, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("F1 at the\n", level_name, " level")) + 
    xlab("Genome size distribution\n(lower k is more skewed)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/k/",dataset,"_",level_name,"_F1VsGenomeSize.png"), p17, width=20, height=20, units = 'cm', dpi=1000)
  
  p18 <- ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$n==80 & pairwise_diffs$mut_rate==0,], 
         aes(x=sigma, y=f1, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = 'y ~ x', size=3) +
    theme_classic() + 
    ylab(paste0("F1 at the\n", level_name, " level")) + 
    xlab("Abundance skew\n(larger sigma is more skewed)") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method")) + 
    scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1.1))
  dir.create(file.path(path_to_data, "Figures/Simulations/F1/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/sigma/",dataset,"_",level_name,"_F1VsSigma.png"), p18, width=20, height=20, units = 'cm', dpi=1000)
  
  ggexport(ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, nrow=3, ncol=6,common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/F1/",dataset,"/",dataset,"_",level_name,"_All_Combined_Precision_Recall_F1_Vs_Parameters.png"), width = 8000, height = 4000)
  ggexport(ggarrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3,common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/Recall/",dataset,"/",dataset,"_",level_name,"_All_Combined_Recall_Vs_Parameters.png"), width = 4500, height = 3000)
  ggexport(ggarrange(p7, p8, p9, p10, p11, p12, nrow=2, ncol=3,common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/",dataset,"_",level_name,"_All_Combined_Precision_Vs_Parameters.png"), width = 4500, height = 3000)
  return(list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12))
}
make_all_precision_recall_f1_individual <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      precision_recall_f1_individual(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
  for (dataset in datasets) {
    phylum_list <- precision_recall_f1_individual(dataset, 2)
    species_list <- precision_recall_f1_individual(dataset, 7)
    plot_list <- append(phylum_list[1:6], species_list[1:6])
    plot_list <- append(plot_list, phylum_list[7:12])
    plot_list <- append(plot_list, species_list[7:12])
    ggexport(ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[6]], 
                       plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[12]], 
                       plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]], plot_list[[18]], 
                       plot_list[[19]], plot_list[[20]], plot_list[[21]], plot_list[[22]], plot_list[[24]], 
                       nrow=4, ncol=5, common.legend = TRUE, legend="right"), 
             filename = paste0(path_to_data, "Figures/Simulations/Precision/",dataset,"/",dataset,"_Herchel_Smith_All_Combined_Precision_Vs_Parameters.png"), width = 8000, height = 6000)
  }
}

absolute_error <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    data[[i]] <- tmp
  }
  
  pairs <- rbind(1:(length(databases) - 1), length(databases))
  
  absolute_error <- function(x, y, x_names, y_names) {
    x <- data.frame("taxa"=x_names, "x"=x)
    y <- data.frame("taxa"=y_names, "y"=y)
    
    merged <- full_join(x,y, by="taxa")
    merged[is.na(merged)] <- 0
    
    return(sum(abs(merged$x - merged$y)))
  }
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(data[[1]][,-1]) * ncol(pairs), ncol = 3))
  for (i in 1:ncol(pairs)) {
    output <- mapply(absolute_error, data[[pairs[1,i]]][,-1], data[[pairs[2,i]]][,-1], data.frame(rep.col(data[[pairs[1,i]]]$TaxIDs, dim(data[[pairs[1,i]]][,-1])[2])), data.frame(rep.col(data[[pairs[2,i]]]$TaxIDs, dim(data[[pairs[2,i]]][,-1])[2])))
    n_entries <- dim(data[[pairs[1,i]]][,-1])[2]
    pairwise_diffs[((i-1) * n_entries + 1):(i * n_entries),1] <- i
    pairwise_diffs[((i-1) * n_entries + 1):(i * n_entries),2] <- colnames(data[[pairs[1,i]]][,-1])
    pairwise_diffs[((i-1) * n_entries + 1):(i * n_entries),3] <- output
  }
  
  colnames(pairwise_diffs) <- c("Method", "sample", "error")
  
  pairwise_diffs$Method <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MetaPhlAn2",
                                            pairwise_diffs[,1]==2 ~ "MetaPhlAn3",
                                            pairwise_diffs[,1]==3 ~ "MetaPhlAn4",
                                            pairwise_diffs[,1]==4 ~ "mOTUs3",
                                            pairwise_diffs[,1]==5 ~ "PhyloPhlAn3",
                                            pairwise_diffs[,1]==6 ~ "Kraken2",
                                            pairwise_diffs[,1]==7 ~ "GTDBTk",
                                            pairwise_diffs[,1]==8 ~ "Metaxa2"))
  
  pairwise_diffs$n = gsub(".*\\.n\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.sample_size.*", replacement = "") %>% as.numeric()
  pairwise_diffs$sample_size = gsub(".*\\.sample_size\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  pairwise_diffs$k = gsub(".*\\.k\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  pairwise_diffs$sigma = gsub(".*\\.sigma\\.", "", pairwise_diffs$sample) %>% gsub(pattern="\\.unknown_prop.*", replacement = "") %>% as.numeric()
  pairwise_diffs$unknown_prop = gsub(".*\\.unknown_prop\\.", "", pairwise_diffs$sample) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()
  
  # Swap study for plotting
  study <- case_when(dataset=="Soil" ~ "Soil Simulations")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")  
  
  ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], aes(x=n, y=error, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Absolute Error") +
    ylab("Absolute Error") + xlab("Number of species")
  dir.create(file.path(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/n/",dataset,"_",level_name,"_AbsoluteErrorVsSpeciesNumber.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(pairwise_diffs[pairwise_diffs$n==80 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], aes(x=sample_size, y=error, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Absolute Error") +
    ylab("Absolute Error") + xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10')
  dir.create(file.path(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/sample_size/",dataset,"_",level_name,"_AbsoluteErrorVsSampleSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$n==80 & pairwise_diffs$k == 100 & pairwise_diffs$sigma==1,], aes(x=unknown_prop, y=error, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Absolute Error") +
    ylab("Absolute Error") + xlab("Unknown proportion")
  dir.create(file.path(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)  
  ggsave(paste0(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/unknown_prop/",dataset,"_",level_name,"_AbsoluteErrorVsUnknown.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$n == 80 & pairwise_diffs$sigma==1,], aes(x=k, y=error, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Absolute Error") +
    ylab("Absolute Error") + xlab("Genome size distribution (lower k is more skewed)")
  dir.create(file.path(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/k/",dataset,"_",level_name,"_AbsoluteErrorVsGenomeSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(pairwise_diffs[pairwise_diffs$sample_size==7.5 & pairwise_diffs$unknown_prop==0.4 & pairwise_diffs$k == 100 & pairwise_diffs$n==80,], aes(x=sigma, y=error, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Absolute Error") +
    ylab("Absolute Error") + xlab("Abundance skew (larger sigma is more skewed)")
  dir.create(file.path(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/Absolute Error/",dataset,"/sigma/",dataset,"_",level_name,"_AbsoluteErrorVsSigma.png"), width=12, height=10, units = 'cm', dpi=1000)
}
make_all_absolute_error <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      absolute_error(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
}

bray_vs_parameter <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    data[[i]][data[[i]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[i]]))] <- 0
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
  
  pairs <- as.matrix(rbind(1:(length(data) - 1), length(data)))
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    output <- vector()
    for (colname in unique(gsub(pattern="_[0-9]$", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
    colnames(pairwise_diffs) <- c("Method 1", "Method 2", unique(gsub(pattern="_[0-9]$", "", colnames(distances))))
  }
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MetaPhlAn2",
                                          pairwise_diffs[,1]==2 ~ "MetaPhlAn3",
                                          pairwise_diffs[,1]==3 ~ "MetaPhlAn4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "PhyloPhlAn3",
                                          pairwise_diffs[,1]==6 ~ "Kraken2",
                                          pairwise_diffs[,1]==7 ~ "GTDBTk",
                                          pairwise_diffs[,1]==8 ~ "Metaxa2"),
                                "/",
                                case_when(pairwise_diffs[,2]==9 ~ "Correct"))
  
  pairwise_diffs[,1:2] <- NULL
  
  bray_df <- melt(pairwise_diffs, id.vars = "Pair")
  bray_df$Method <- gsub("\\/.*","",bray_df$Pair)
  
  bray_df$variable <- gsub("_sample_.*", "", bray_df$variable)
  
  bray_df$n = gsub(".*\\.n\\.", "", bray_df$variable) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric()
  bray_df$sample_size = gsub(".*\\.size\\.", "", bray_df$variable) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  bray_df$k = gsub(".*\\.k\\.", "", bray_df$variable) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  bray_df$sigma = gsub(".*\\.sigma\\.", "", bray_df$variable) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric()
  bray_df$unknown_prop = gsub(".*\\.up\\.", "", bray_df$variable) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric()
  bray_df$mut_rate = gsub(".*\\.mut_rate\\.", "", bray_df$variable) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()

  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  colAdd <- c("#00ffff", "#00ccff", "#0066cc", "#81C784", "#2E7D32", "#7B1FA2", "#FF8F00", "#FBC02D", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(unique(bray_df$Method))])
  
  p1 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$unknown_prop==0.4 & bray_df$k == 100 & bray_df$sigma==1 & bray_df$mut_rate==0,], 
               aes(x=n, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Bray Curtis Dissimilarity") + xlab("Number of species") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  
  #dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/n/",dataset,"_",level_name,"_BrayVsSpeciesNumber.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p2 <- ggplot(bray_df[bray_df$n==80 & bray_df$unknown_prop==0.4 & bray_df$k == 100 & bray_df$sigma==1 & bray_df$mut_rate==0,], 
               aes(x=sample_size, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Bray Curtis Dissimilarity") + xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10') + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  
  #dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sample_size/",dataset,"_",level_name,"_BrayVsSampleSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p3 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$n==80 & bray_df$k == 100 & bray_df$sigma==1 & bray_df$mut_rate==0,], 
               aes(x=unknown_prop, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Bray Curtis Dissimilarity") + xlab("Proportion of unknown species") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
    
  
  #dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/unknown_prop/",dataset,"_",level_name,"_BrayVsUnknown.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p4 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$unknown_prop==0.4 & bray_df$n == 80 & bray_df$k == 100 & bray_df$sigma==1,], 
               aes(x=mut_rate, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Bray Curtis Dissimilarity") + xlab("Random insertion/deletion rate") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  
  #dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/k/",dataset,"_",level_name,"_BrayVsGenomeSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # p5 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$unknown_prop==0.4 & bray_df$k == 100 & bray_df$n==80 & bray_df$mut_rate==0,], 
  #              aes(x=sigma, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
  #   geom_point(size=5) +
  #   geom_smooth(method = "lm", formula = "y~x", size=3) +
  #   theme_classic() + 
  #   ylab("Bray Curtis Dissimilarity") + xlab("Abundance skew (larger sigma is more skewed)") + 
  #   col + 
  #   theme(text=element_text(size=80)) + 
  #   guides(color=guide_legend(title="Method"))
  #dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_BrayVsSigma.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggexport(ggarrange(p1, p2, p3, p4, nrow=1, ncol=4,common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/",dataset,"_",level_name,"_All_Combined_BC_Vs_Parameters.png"), width = 4200, height = 1800)
  

}
make_all_bray_vs_parameter <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      bray_vs_parameter(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
}

bray_vs_parameter_individual <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  for (i in 1:length(data)) {
    data[[i]][data[[i]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[i]]))] <- 0
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
  
  pairs <- as.matrix(rbind(1:(length(data) - 1), length(data)))
  
  pairwise_diffs <- data.frame(matrix(nrow = ncol(pairs), ncol = 2 + nrow(data[[1]])))
  for (i in 1:ncol(pairs)) {
    dat_taxa_species <- bind_rows(data[[pairs[1,i]]], data[[pairs[2,i]]])
    dat_taxa_species[is.na(dat_taxa_species)] <- 0
    distances <- vegdist(dat_taxa_species, method="bray") %>% as.matrix(labels=TRUE)
    output <- vector()
    for (colname in unique(gsub(pattern="_[0-9]$", "", colnames(distances)))) {
      same_ids <- colnames(distances)[grepl(colname, colnames(distances))]
      output <- append(output, distances[same_ids[1], same_ids[2]])
    }
    pairwise_diffs[i,] <- c("Method 1"=pairs[1,i], "Method 2"=pairs[2,i], output)
    colnames(pairwise_diffs) <- c("Method 1", "Method 2", unique(gsub(pattern="_[0-9]$", "", colnames(distances))))
  }
  
  pairwise_diffs$Pair <- paste0(case_when(pairwise_diffs[,1]==1 ~ "MetaPhlAn 2",
                                          pairwise_diffs[,1]==2 ~ "MetaPhlAn 3",
                                          pairwise_diffs[,1]==3 ~ "MetaPhlAn 4",
                                          pairwise_diffs[,1]==4 ~ "mOTUs3",
                                          pairwise_diffs[,1]==5 ~ "PhyloPhlAn 3",
                                          pairwise_diffs[,1]==6 ~ "Kraken 2",
                                          pairwise_diffs[,1]==7 ~ "GTDB-Tk 2",
                                          pairwise_diffs[,1]==8 ~ "Metaxa 2"),
                                "/",
                                case_when(pairwise_diffs[,2]==9 ~ "Correct"))
  
  pairwise_diffs[,1:2] <- NULL
  
  bray_df <- melt(pairwise_diffs, id.vars = "Pair")
  bray_df$Method <- gsub("\\/.*","",bray_df$Pair)
  
  bray_df$variable <- gsub("_sample_.*", "", bray_df$variable)
  
  bray_df$n = gsub(".*\\.n\\.", "", bray_df$variable) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric()
  bray_df$sample_size = gsub(".*\\.size\\.", "", bray_df$variable) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  bray_df$k = gsub(".*\\.k\\.", "", bray_df$variable) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  bray_df$sigma = gsub(".*\\.sigma\\.", "", bray_df$variable) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric()
  bray_df$unknown_prop = gsub(".*\\.up\\.", "", bray_df$variable) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric()
  bray_df$mut_rate = gsub(".*\\.mut_rate\\.", "", bray_df$variable) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")
  
  colAdd <- c("#00ffff", "#00ccff", "#0066cc", "#81C784", "#2E7D32", "#7B1FA2", "#FF8F00", "#FBC02D", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(unique(bray_df$Method))])
  
  p1 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$unknown_prop==0.4 & bray_df$k == 100 & bray_df$sigma==1 & bray_df$mut_rate==0,], 
               aes(x=n, y=value, color=factor(Method, levels=c("MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs3", "Kraken 2", "GTDB-Tk 2", "PhyloPhlAn 3")))) + geom_point() +
    geom_point(size=7) +
    geom_smooth(method = "lm", formula = "y~x", size=5) +
    theme_classic() + 
    ylab(paste0("Bray Curtis dissimilarity\nat the ", level_name, " level")) + 
    xlab("Species count") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=75),
          plot.margin = unit(c(2,1,1,1), "cm"),
          legend.key.size = unit(3, 'cm'), #change legend key size
          legend.key.height = unit(3, 'cm'), #change legend key height
          legend.key.width = unit(3, 'cm'), #change legend key width
          legend.title = element_text(size=120), #change legend title font size
          legend.text = element_text(size=80)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/n/",dataset,"_",level_name,"_BrayVsSpeciesNumber.png"), p1, width=20, height=20, units = 'cm', dpi=1000)
  
  p2 <- ggplot(bray_df[bray_df$n==80 & bray_df$unknown_prop==0.4 & bray_df$k == 100 & bray_df$sigma==1 & bray_df$mut_rate==0,], 
               aes(x=sample_size, y=value, color=factor(Method, levels=c("MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs3", "Kraken 2", "GTDB-Tk 2", "PhyloPhlAn 3")))) + geom_point() +
    geom_point(size=7) +
    geom_smooth(method = "lm", formula = "y~x", size=5) +
    theme_classic() + 
    xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10') + 
    col + 
    theme(text=element_text(size=75),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.y = element_blank(),
          legend.key.size = unit(3, 'cm'), #change legend key size
          legend.key.height = unit(3, 'cm'), #change legend key height
          legend.key.width = unit(3, 'cm'), #change legend key width
          legend.title = element_text(size=120), #change legend title font size
          legend.text = element_text(size=80)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sample_size/",dataset,"_",level_name,"_BrayVsSampleSize.png"), p2, width=20, height=20, units = 'cm', dpi=1000)
  
  p3 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$n==80 & bray_df$k == 100 & bray_df$sigma==1 & bray_df$mut_rate==0,], 
               aes(x=100 * unknown_prop, y=value, color=factor(Method, levels=c("MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs3", "Kraken 2", "GTDB-Tk 2", "PhyloPhlAn 3")))) + geom_point() +
    geom_point(size=7) +
    geom_smooth(method = "lm", formula = "y~x", size=5) +
    theme_classic() + 
    xlab("Unknown species (%)") + 
    col + 
    theme(text=element_text(size=75),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.y = element_blank(),
          legend.key.size = unit(3, 'cm'), #change legend key size
          legend.key.height = unit(3, 'cm'), #change legend key height
          legend.key.width = unit(3, 'cm'), #change legend key width
          legend.title = element_text(size=120), #change legend title font size
          legend.text = element_text(size=80)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/unknown_prop/",dataset,"_",level_name,"_BrayVsUnknown.png"), p3, width=20, height=20, units = 'cm', dpi=1000)
  
  p4 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$unknown_prop==0.4 & bray_df$n == 80 & bray_df$k == 100 & bray_df$sigma==1,], 
               aes(x=mut_rate, y=value, color=factor(Method, levels=c("MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs3", "Kraken 2", "GTDB-Tk 2", "PhyloPhlAn 3")))) + geom_point() +
    geom_point(size=7) +
    geom_smooth(method = "lm", formula = "y~x", size=5) +
    theme_classic() + 
    xlab("Mutation rate") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=75),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.y = element_blank(),
          legend.key.size = unit(3, 'cm'), #change legend key size
          legend.key.height = unit(3, 'cm'), #change legend key height
          legend.key.width = unit(3, 'cm'), #change legend key width
          legend.title = element_text(size=120), #change legend title font size
          legend.text = element_text(size=80)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/mut_rate/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/mut_rate/",dataset,"_",level_name,"_BrayVsMutRate.png"), p4, width=20, height=20, units = 'cm', dpi=1000)
  
  p5 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$unknown_prop==0.4 & bray_df$n == 80 & bray_df$sigma==1 & bray_df$mut_rate==0,], 
               aes(x=k, y=value, color=factor(Method, levels=c("MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs3", "Kraken 2", "GTDB-Tk 2", "PhyloPhlAn 3")))) + geom_point() +
    geom_point(size=7) +
    geom_smooth(method = "lm", formula = "y~x", size=5) +
    theme_classic() + 
    xlab("Genome size distribution\n(lower k is more skewed)") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=75),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.y = element_blank(),
          legend.key.size = unit(3, 'cm'), #change legend key size
          legend.key.height = unit(3, 'cm'), #change legend key height
          legend.key.width = unit(3, 'cm'), #change legend key width
          legend.title = element_text(size=120), #change legend title font size
          legend.text = element_text(size=80)) + 
    guides(color=guide_legend(title="Tool"))
  
  dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/k/",dataset,"_",level_name,"_BrayVsGenomeSize.png"), p5, width=20, height=20, units = 'cm', dpi=1000)
  
  p6 <- ggplot(bray_df[bray_df$sample_size==7.5 & bray_df$unknown_prop==0.4 & bray_df$k == 100 & bray_df$n==80 & bray_df$mut_rate==0,],
               aes(x=sigma, y=value, color=factor(Method, levels=c("MetaPhlAn 2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs3", "Kraken 2", "GTDB-Tk 2", "PhyloPhlAn 3")))) +
    geom_point(size=7) +
    geom_smooth(method = "lm", formula = "y~x", size=5) +
    theme_classic() +
    xlab("Abundance skew") +
    col +
    theme(text=element_text(size=75),
          plot.margin = unit(c(2,1,1,1), "cm"),
          axis.title.y = element_blank(),
          legend.key.size = unit(3, 'cm'), #change legend key size
          legend.key.height = unit(3, 'cm'), #change legend key height
          legend.key.width = unit(3, 'cm'), #change legend key width
          legend.title = element_text(size=120), #change legend title font size
          legend.text = element_text(size=80)) +
    guides(color=guide_legend(title="Tool"))
  dir.create(file.path(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/sigma/",dataset,"_",level_name,"_BrayVsSigma.png"), p6, width=20, height=20, units = 'cm', dpi=1000)
  
  ggexport(ggarrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3,common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/",dataset,"_",level_name,"_All_Combined_BC_Vs_Parameters.png"), width = 4500, height = 3000)
  return(list(p1, p2, p3, p4, p5, p6))
}
make_all_bray_vs_parameter_individual <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      bray_vs_parameter_individual(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
  for (dataset in datasets) {
    phylum_list <- bray_vs_parameter_individual(dataset, 2)
    species_list <- bray_vs_parameter_individual(dataset, 7)
    plot_list <- append(phylum_list[1:6], species_list[1:6])
    ggexport(ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[6]], 
                       plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[12]],
                       nrow=2, ncol=5, legend="right", legend.grob = get_legend(plot_list[[1]]), widths = c(1.2, 1, 1, 1, 1)), 
             filename = paste0(path_to_data, "Figures/Simulations/BrayVsParameter/",dataset,"/",dataset,"_Herchel_Smith_All_Combined_BC_Vs_Parameters.png"), width = 5000, height = 2000)
  }
}

unknown_absolute_error <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  output_matrix <- matrix(nrow = length(data) - 1, ncol = ncol(data[[1]]) - 1)
  for (i in 1:(length(data) - 1)) {
    output_matrix[i,] <- abs(unlist(data[[i]][data[[i]]$TaxIDs=="UNCLASSIFIED",-1]) - unlist(data[[length(data)]][data[[length(data)]]$TaxIDs=="UNCLASSIFIED",-1]))
  }
  
  output_df <- data.frame(t(output_matrix))
  colnames(output_df) <- c("MPA2", "MPA3", "MPA4", "Kraken", "metaxa2")
  output_df$sample <- colnames(data[[1]][-1])
  
  output_df$n = gsub(".*\\.n\\.", "", output_df$sample) %>% gsub(pattern="\\.sample_size.*", replacement = "") %>% as.numeric()
  output_df$sample_size = gsub(".*\\.sample_size\\.", "", output_df$sample) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  output_df$k = gsub(".*\\.k\\.", "", output_df$sample) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  output_df$sigma = gsub(".*\\.sigma\\.", "", output_df$sample) %>% gsub(pattern="\\.unknown_prop.*", replacement = "") %>% as.numeric()
  output_df$unknown_prop = gsub(".*\\.unknown_prop\\.", "", output_df$sample) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()
  
  melted_df <- melt(output_df, id.vars=c("sample", "n", "sample_size", "k", "sigma", "unknown_prop"), variable.name = "Method")
  
  # Swap study for plotting
  study <- case_when(dataset=="Soil" ~ "Soil Simulations")
  
  level_name <- case_when(level == 1 ~ "kingdoms",
                          level == 2 ~ "phyla", 
                          level == 3 ~ "classes",
                          level == 4 ~ "orders",
                          level == 5 ~ "families",
                          level == 6 ~ "genera",
                          level == 7 ~ "species")  
  
  ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$sigma==1,], aes(x=n, y=value, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Correct vs Estimated Unknown") +
    ylab("Absolute Error") + xlab("Number of species")
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/n/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSpeciesNumber.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(melted_df[melted_df$n==80 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$sigma==1,], aes(x=sample_size, y=value, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Correct vs Estimated Unknown") +
    ylab("Absolute Error") + xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10')
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sample_size/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSampleSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$n==80 & melted_df$k == 100 & melted_df$sigma==1,], aes(x=unknown_prop, y=value, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Correct vs Estimated Unknown") +
    ylab("Absolute Error") + xlab("Unknown proportion")
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)  
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/unknown_prop/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsUnknown.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$n == 80 & melted_df$sigma==1,], aes(x=k, y=value, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Correct vs Estimated Unknown") +
    ylab("Absolute Error") + xlab("Genome size distribution (lower k is more skewed)")
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/k/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsGenomeSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$n==80,], aes(x=sigma, y=value, color=Method)) + geom_point() +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    theme_classic() + ggtitle("Correct vs Estimated Unknown") +
    ylab("Absolute Error") + xlab("Abundance skew (larger sigma is more skewed)")
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sigma/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSigma.png"), width=12, height=10, units = 'cm', dpi=1000)
}
make_all_unknown_absolute_error <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      unknown_absolute_error(dataset, level)
      print(paste0(dataset, "/", level))
    }
  }
}

unknown_error <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  output_matrix <- matrix(nrow = length(data) - 1, ncol = ncol(data[[1]]) - 1)
  for (i in 1:(length(data) - 1)) {
    output_matrix[i,] <- unlist(data[[i]][data[[i]]$TaxIDs=="UNCLASSIFIED",-1]) - unlist(data[[length(data)]][data[[length(data)]]$TaxIDs=="UNCLASSIFIED",-1])
  }
  
  output_df <- data.frame(t(output_matrix))
  colnames(output_df) <- c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "mOTUs3", "PhyloPhlAn3", "Kraken2", "GTDBTk", "Metaxa2")
  output_df$sample <- colnames(data[[1]][-1])
  
  output_df$n = gsub(".*\\.n\\.", "", output_df$sample) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric()
  output_df$sample_size = gsub(".*\\.size\\.", "", output_df$sample) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  output_df$k = gsub(".*\\.k\\.", "", output_df$sample) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  output_df$sigma = gsub(".*\\.sigma\\.", "", output_df$sample) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric()
  output_df$unknown_prop = gsub(".*\\.up\\.", "", output_df$sample) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric()
  output_df$mut_rate = gsub(".*\\.mut_rate\\.", "", output_df$sample) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()
  
  melted_df <- melt(output_df, id.vars=c("sample", "n", "sample_size", "k", "sigma", "unknown_prop", "mut_rate"), variable.name = "Method")
  
  # Swap study for plotting
  study <- case_when(dataset=="Soil" ~ "Soil Simulations")
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")  
  
  colAdd <- c("#00ffff", "#00ccff", "#0066cc", "#81C784", "#2E7D32", "#7B1FA2", "#FF8F00", "#FBC02D", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(unique(melted_df$Method))])
  
  p1 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$sigma==1 & melted_df$mut_rate==0,], 
         aes(x=n, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Number of species") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  #dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/n/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSpeciesNumber.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p2 <- ggplot(melted_df[melted_df$n==80 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$sigma==1 & melted_df$mut_rate==0,], 
               aes(x=sample_size, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10') + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  #dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sample_size/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSampleSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p3 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$n==80 & melted_df$k == 100 & melted_df$sigma==1 & melted_df$mut_rate==0,], 
               aes(x=unknown_prop, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Species proportion unknown") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  #dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)  
  #ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/unknown_prop/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsUnknown.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  p4 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$n == 80 & melted_df$k == 100 & melted_df$sigma==1,], 
               aes(x=mut_rate, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Random insertion/deletion rate") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  
  #dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/k/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsGenomeSize.png"), width=12, height=10, units = 'cm', dpi=1000)
  
  # p5 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$n==80,], 
  #        aes(x=sigma, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
  #   geom_point(size=5) +
  #   geom_smooth(method = "lm", formula = 'y ~ x', size=3)+
  #   theme_classic() + 
  #   ylab("Estimated Unknown - True Unknown") + xlab("Abundance skew (larger sigma is more skewed)")
  #dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  #ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sigma/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSigma.png"), width=12, height=10, units = 'cm', dpi=1000)

  ggexport(ggarrange(p1, p2, p3, p4, nrow=1, ncol=4,common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSigma.png"), width = 5200, height = 1800)
  
}
make_all_unknown_error <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      unknown_error(dataset, level)
    }
  }
}

unknown_error_individual <- function(dataset, level) {
  data <- preprocess(dataset, level)
  
  output_matrix <- matrix(nrow = length(data) - 1, ncol = ncol(data[[1]]) - 1)
  for (i in 1:(length(data) - 1)) {
    output_matrix[i,] <- unlist(data[[i]][data[[i]]$TaxIDs=="UNCLASSIFIED",-1]) - unlist(data[[length(data)]][data[[length(data)]]$TaxIDs=="UNCLASSIFIED",-1])
  }
  
  output_df <- data.frame(t(output_matrix))
  colnames(output_df) <- c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "mOTUs3", "PhyloPhlAn3", "Kraken2", "GTDBTk", "Metaxa2")
  output_df$sample <- colnames(data[[1]][-1])
  
  output_df$n = gsub(".*\\.n\\.", "", output_df$sample) %>% gsub(pattern="\\.size.*", replacement = "") %>% as.numeric()
  output_df$sample_size = gsub(".*\\.size\\.", "", output_df$sample) %>% gsub(pattern="\\.gc.*", replacement = "") %>% as.numeric()
  output_df$k = gsub(".*\\.k\\.", "", output_df$sample) %>% gsub(pattern="\\.sigma.*", replacement = "") %>% as.numeric()
  output_df$sigma = gsub(".*\\.sigma\\.", "", output_df$sample) %>% gsub(pattern="\\.up.*", replacement = "") %>% as.numeric()
  output_df$unknown_prop = gsub(".*\\.up\\.", "", output_df$sample) %>% gsub(pattern="\\.mut_rate.*", replacement = "") %>% as.numeric()
  output_df$mut_rate = gsub(".*\\.mut_rate\\.", "", output_df$sample) %>% gsub(pattern="_sample_.*", replacement = "") %>% as.numeric()
  
  melted_df <- melt(output_df, id.vars=c("sample", "n", "sample_size", "k", "sigma", "unknown_prop", "mut_rate"), variable.name = "Method")
  
  # Swap study for plotting
  study <- case_when(dataset=="Soil" ~ "Soil Simulations")
  
  level_name <- case_when(level == 1 ~ "kingdom",
                          level == 2 ~ "phylum", 
                          level == 3 ~ "class",
                          level == 4 ~ "order",
                          level == 5 ~ "family",
                          level == 6 ~ "genus",
                          level == 7 ~ "species")  
  
  colAdd <- c("#00ffff", "#00ccff", "#0066cc", "#81C784", "#2E7D32", "#7B1FA2", "#FF8F00", "#FBC02D", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(unique(melted_df$Method))])
  
  p1 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$sigma==1 & melted_df$mut_rate==0,], 
               aes(x=n, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=5) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Number of species") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/n/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/n/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSpeciesNumber.png"), p1, width=20, height=20, units = 'cm', dpi=1000)
  
  p2 <- ggplot(melted_df[melted_df$n==80 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$sigma==1 & melted_df$mut_rate==0,], 
               aes(x=sample_size, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Sample size (GB)") + 
    scale_x_continuous(trans='log10') + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sample_size/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sample_size/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSampleSize.png"), p2, width=20, height=20, units = 'cm', dpi=1000)
  
  p3 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$n==80 & melted_df$k == 100 & melted_df$sigma==1 & melted_df$mut_rate==0,], 
               aes(x=unknown_prop, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Species proportion unknown") + 
    col + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/unknown_prop/"), showWarnings = FALSE, recursive = T)  
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/unknown_prop/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsUnknown.png"), p3, width=20, height=20, units = 'cm', dpi=1000)
  
  p4 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$n == 80 & melted_df$k == 100 & melted_df$sigma==1,], 
               aes(x=mut_rate, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Random insertion/deletion rate") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/mut_rate/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/mut_rate/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsMutRate.png"), width=20, height=20, units = 'cm', dpi=1000)
  
  p5 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$n == 80 & melted_df$sigma==1 & melted_df$mut_rate==0,], 
               aes(x=k, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + geom_point() +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Genome size distribution\n(lower k is more skewed)") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/k/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/k/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsGenomeSize.png"), width=20, height=20, units = 'cm', dpi=1000)
  
  p6 <- ggplot(melted_df[melted_df$sample_size==7.5 & melted_df$unknown_prop==0.4 & melted_df$k == 100 & melted_df$n==80 & melted_df$mut_rate==0,],
         aes(x=sigma, y=value, color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) +
    geom_point(size=5) +
    geom_smooth(method = "lm", formula = "y~x", size=3) +
    theme_classic() + 
    ylab("Estimated Unknown - True Unknown") + xlab("Abundance skew (larger sigma is more skewed)") + 
    col + 
    scale_x_continuous(trans='log10') + 
    theme(text=element_text(size=80),
          plot.margin = unit(c(2,1,1,1), "cm")) + 
    guides(color=guide_legend(title="Method"))
  dir.create(file.path(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sigma/"), showWarnings = FALSE, recursive = T)
  ggsave(paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/sigma/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSigma.png"), width=20, height=20, units = 'cm', dpi=1000)
  
  ggexport(ggarrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3, common.legend = TRUE, legend="right"), 
           filename = paste0(path_to_data, "Figures/Simulations/CorrectEstimatedUnknown/",dataset,"/",dataset,"_",level_name,"_CorrectEstimatedUnknownVsSigma.png"), width = 4500, height = 3000)
  
  return(list(p1, p2, p3, p4, p5, p6))
}
make_all_unknown_error_individual <- function() {
  datasets <- c("Soil")
  levels <- 1:7
  for (dataset in datasets) {
    for (level in levels) {
      unknown_error_individual(dataset, level)
    }
  }
}

bc_by_rank <- function(dataset, sample_settings, use_title) {
  
  data_tmp <- preprocess(dataset, 1)
  
  output_df <- data.frame(matrix(ncol = 4, nrow = 7 * (length(data_tmp) - 1)))
  colnames(output_df) <- c("level", "method", "bc", "sd")
  
  for (level in 1:7) {
    data <- preprocess(dataset, level)
    data[[length(data)]][data[[length(data)]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[length(data)]]))] <- 0
    for (i in 1:(length(data) - 1)) {
      data[[i]][data[[i]]$TaxIDs=="254246",grepl("up\\.1\\.0", colnames(data[[i]]))] <- 0
      
      # Normalize
      tmp <- data[[i]]
      classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
      tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
      tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
      data[[i]] <- tmp
      
      extract_df = data[[i]][,grepl(paste0(sample_settings, "|TaxIDs"), colnames(data[[i]]))]
      colnames(extract_df)[-1] <- paste0(colnames(extract_df)[-1], "_extracted")
      correct_df = data[[length(data)]][,grepl(paste0(sample_settings, "|TaxIDs"), colnames(data[[length(data)]]))]
      colnames(correct_df)[-1] <- paste0(colnames(correct_df)[-1], "_correct")
      merged <- merge(extract_df, correct_df, by="TaxIDs", all = T)
      merged[is.na(merged)] <- 0
      TaxIDs <- merged$TaxIDs
      merged <- t(merged[,-1])
      colnames(merged) <- TaxIDs
      
      tmp <- vegdist(merged, method="bray") %>% as.matrix(labels=TRUE)
      vals <- diag(tmp[grepl("extracted", rownames(tmp)),grepl("correct", rownames(tmp))])
      output_df[((level - 1) * (length(data) - 1) + (i - 1) + 1):((level - 1) * (length(data) - 1) + i ), 1:4] <- 
        cbind(level, i, mean(vals), sd(vals))
    }
    
  }
  
  output_df$Method <- case_when(output_df$method==1 ~ "MetaPhlAn2",
                                output_df$method==2 ~ "MetaPhlAn3",
                                output_df$method==3 ~ "MetaPhlAn4",
                                output_df$method==4 ~ "mOTUs3",
                                output_df$method==5 ~ "PhyloPhlAn3",
                                output_df$method==6 ~ "Kraken2",
                                output_df$method==7 ~ "GTDBTk",
                                output_df$method==8 ~ "Metaxa2")
  
  output_df$level <- case_when(output_df$level == 1 ~ "Kingdom",
                               output_df$level == 2 ~ "Phylum", 
                               output_df$level == 3 ~ "Class",
                               output_df$level == 4 ~ "Order",
                               output_df$level == 5 ~ "Family",
                               output_df$level == 6 ~ "Genus",
                               output_df$level == 7 ~ "Species")
  
  colAdd <- c("#00ffff", "#00ccff", "#0066cc", "#81C784", "#2E7D32", "#7B1FA2", "#FF8F00", "#FBC02D", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142")
  col = scale_color_manual(values = colAdd[1:length(unique(output_df$Method))])
  
  if (use_title) {
    p <- ggplot(output_df, aes(y=bc, x=factor(level, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), 
                               group=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
      geom_line(aes(color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3"))), size=3) + 
      geom_errorbar(aes(ymin=bc-sd, ymax=bc+sd), width=2,
                    position=position_dodge(0.05), size=3) +
      theme_classic() +
      xlab(paste0("Unknown Proportion: ", gsub(".*up\\.", "", sample_settings))) +
      ylab(paste0("Bray Curtis dissimilarity")) + 
      col + 
      guides(color=guide_legend(title="Method")) + 
      theme(text=element_text(size=80),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1))
  } else {
    p <- ggplot(output_df, aes(y=bc, x=factor(level, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), 
                               group=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3")))) + 
      geom_line(aes(color=factor(Method, levels=c("MetaPhlAn2", "MetaPhlAn3", "MetaPhlAn4", "Metaxa2", "mOTUs3", "Kraken2", "GTDBTk", "PhyloPhlAn3"))), size=3) + 
      geom_errorbar(aes(ymin=bc-sd, ymax=bc+sd), width=2,
                    position=position_dodge(0.05), size=3) +
      theme_classic() +
      xlab(paste0("Unknown Proportion: ", gsub(".*up\\.", "", sample_settings))) +
      ylab(paste0("Bray Curtis dissimilarity")) + 
      col + 
      theme(legend.position = "none", text=element_text(size=80),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      scale_y_continuous(breaks=seq(0,1,by=.2), limits = c(0,1))
  }
  return(p)
}
make_bc_by_rank <- function() {
  dataset <- "Soil"
  dir.create(file.path(path_to_data, "Figures/Simulations/BC_By_Rank/",dataset,"/"), showWarnings = FALSE, recursive = T)
  png(file=paste0(path_to_data, "Figures/Simulations/BC_By_Rank/",dataset,"/",dataset,"_All_Combined_ECDFUnknowns.png"), width = 4800, height = 2000)
  grid.arrange(as.grob(bc_by_rank("Soil", "profile.n.80.size.7.5.gc.unrestricted.k.100.sigma.1.0.up.0.0", F)), 
               as.grob(bc_by_rank("Soil", "profile.n.80.size.7.5.gc.unrestricted.k.100.sigma.1.0.up.0.2", F)), 
               #as.grob(bc_by_rank("Soil", "profile.n.80.size.7.5.gc.unrestricted.k.100.sigma.1.0.up.0.4", F)),
               #as.grob(bc_by_rank("Soil", "profile.n.80.size.7.5.gc.unrestricted.k.100.sigma.1.0.up.0.6", F)),
               as.grob(bc_by_rank("Soil", "profile.n.80.size.7.5.gc.unrestricted.k.100.sigma.1.0.up.0.8", F)),
               as.grob(bc_by_rank("Soil", "profile.n.80.size.7.5.gc.unrestricted.k.100.sigma.1.0.up.1.0", T)), nrow=1, ncol=4,
               top = textGrob("Bray Curtis Dissimilarity vs True Profile",gp=gpar(fontsize=100,font=1)), widths=c(2,2,2,2.8))
  dev.off()
}

# unknown_saturation <- function(dataset, sample_settings) {
#   
#   data_tmp <- preprocess(dataset, 1)
#   
#   output_df <- data.frame(matrix(ncol = 4, nrow = 7 * length(data_tmp)))
#   colnames(output_df) <- c("level", "method", "bc", "sd")
#   
#   for (level in 1:7) {
#     data <- preprocess(dataset, level)
#     for (i in 1:length(data)) {
#       extract_df = data[[i]][,grepl(paste0(sample_settings, "|TaxIDs"), colnames(data[[i]]))]
#       extract_df$TaxIDs[extract_df$TaxIDs=="UNCLASSIFIED"] <- "UNCLASSIFIED2"
#       colnames(extract_df)[-1] <- paste0(colnames(extract_df)[-1], "_extracted")
#       correct_df = data[[length(data)]][,grepl(paste0(sample_settings, "|TaxIDs"), colnames(data[[length(data)]]))]
#       colnames(correct_df)[-1] <- paste0(colnames(correct_df)[-1], "_correct")
#       merged <- merge(extract_df, correct_df, by="TaxIDs", all = T)
#       merged[is.na(merged)] <- 0
#       TaxIDs <- merged$TaxIDs
#       merged <- t(merged[,-1])
#       colnames(merged) <- TaxIDs
#       
#       tmp <- vegdist(merged, method="bray") %>% as.matrix(labels=TRUE)
#       vals <- diag(tmp[grepl("extracted", rownames(tmp)),grepl("correct", rownames(tmp))])
#       output_df[((level - 1) * length(data) + (i - 1) + 1):((level - 1) * length(data) + i ), 1:4] <- 
#         cbind(level, i, mean(vals), sd(vals))
#     }
#     
#   }
#   
#   output_df$method <- case_when(output_df$method==1 ~ "MPA2",
#                                 output_df$method==2 ~ "MPA3",
#                                 output_df$method==3 ~ "MPA4",
#                                 #output_df$method==4 ~ "mOTUs3",
#                                 #output_df$method==5 ~ "Phylophlan",
#                                 output_df$method==4 ~ "Kraken",
#                                 #output_df$method==7 ~ "GTDBTk",
#                                 output_df$method==5 ~ "metaxa2",
#                                 output_df$method==6 ~ "Best")
#   
#   output_df$level <- case_when(output_df$level == 1 ~ "kingdoms",
#                                output_df$level == 2 ~ "phyla", 
#                                output_df$level == 3 ~ "classes",
#                                output_df$level == 4 ~ "orders",
#                                output_df$level == 5 ~ "families",
#                                output_df$level == 6 ~ "genera",
#                                output_df$level == 7 ~ "species")  
#   
#   # Swap study for plotting
#   study <- case_when(dataset=="Soil" ~ "Soil Simulations")
#   
#   
#   
#   ggplot(output_df, aes(y=bc, x=factor(level, c("kingdoms", "phyla", "classes", "orders", "families", "genera", "species")), group=method)) + 
#     geom_line(aes(color=method), size=1) + 
#     geom_errorbar(aes(ymin=bc-sd, ymax=bc+sd), width=0.5,
#                   position=position_dodge(0.05)) +
#     theme_classic() +
#     ggtitle(study) +
#     xlab("Level") +
#     ylab(paste0("Bray Curtis dissimilarity"))
#   
#   dir.create(file.path(path_to_data, "Figures/Simulations/UnknownSaturation/",dataset,"/"), showWarnings = FALSE, recursive = T)
#   ggsave(paste0(path_to_data, "Figures/Simulations/UnknownSaturation/",dataset,"/",sample_settings,"_UnknownSaturation.png"), width=12, height=10, units = 'cm', dpi=1000)
#   
# }
# make_all_unknown_saturation <- function() {
#   datasets <- c("Soil")
#   for (dataset in datasets) {
#     sample_settings <- unique(gsub("_sample_.*$", "", colnames(preprocess(dataset, 1)[[1]][-1])))
#     for (sample_setting in sample_settings) {
#       unknown_saturation(dataset, sample_setting)
#     }
#   }
# }