### For every NCBI-renormalized sample across all tools, simulation environments,
### and real datasets, compute what fraction of abundance comes from taxa present
### in the phylogenetic tree built from the unique TaxIDs.

source("analysis/scripts/helpers.R")

library(dplyr)
library(ape)

threshold <- .1

phylo_tree <- ape::read.tree("analysis/databases/Phylogenetic_trees/phylogenetic_tree.nwk")
tree_tips <- phylo_tree$tip.label

cat("Tree has", length(tree_tips), "tips\n\n")

label_converts_soil <- c('2599805'=931866, '1404649'=2823807, '2735433'=3014751, '2829818'=2952571, '2840472'=2840469, 
                         '2884447'=2792859, '80870'=80867, '1532'=33035, '2778071'=29523, '285567'=1942, '290385'=104623, '37480'=67362,
                         '423539'=3074428, '50340'=53407, '84292'=162393, '88075'=88074, '96101'=46163, '1497613'=1505087, 
                         '1915400'=67332, '335659'=1404864, '90270'=2754056, '80870'=80867, '1532'=33035, 
                         '285567'=1942, '290385'=104623, '37480'=67362, '423539'=3074428, '50340'=53407,
                         '84292'=162393, '88075'=88074, '96101'=46163, '335659'=1404864, '90270'=2754056)

label_converts_ocean <- c('2841037'=3028070, '29570'=77097, '1404649'=2823807, '2884447'=2792859, '158080'=141390, '50340'=53407, '81037'=77608, 
                          '1304902'=2731756, '1453999'=2954383, '1457154'=2954382, '29570'=77097, '158080'=141390,
                          '50340'=53407, '81037'=77608, '1960941'=2942470)

label_converts_gut <- c('2841037'=3028070, '1404649'=2823807, '2884447'=2792859, '1532'=33035, '1841867'=2897707, 
                        '1118061'=2585118, '1453999'=2954383, '1457154'=2954382, '1532'=33035, '1898205'=2485925,
                        '1960941'=2942470)

tool_names <- c("Centrifuge", "Kraken 2/Bracken 2", "MetaPhlAn 2", "MetaPhlAn 3", 
                "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "GTDB-Tk (MEGAHIT)", 
                "GTDB-Tk (metaSPAdes)", "PhyloPhlAn 3 (MEGAHIT)", "PhyloPhlAn 3 (metaSPAdes)")
tool_names_no_metaspades <- tool_names[1:9]
tool_names_no_metaspades[8] <- "GTDB-Tk (MEGAHIT)"

process_profiles <- function(profiles, label_converts = NULL) {
  profiles <- remove_non_ncbi(profiles)
  profiles <- lapply(profiles, renormalize)
  profiles <- lapply(profiles, threshold_sample, threshold)
  profiles <- lapply(profiles, function(x) {
    if (length(which(rowSums(x[,-1]) == 0)) != 0)
      return(x[-which(rowSums(x[,-1]) == 0),])
    else
      return(x)
  })
  profiles <- lapply(profiles, renormalize)
  if (!is.null(label_converts)) {
    profiles <- lapply(profiles, function(x) update_labels(x, label_converts))
  }
  return(profiles)
}

calc_tree_coverage <- function(profile, tree_tips) {
  in_tree <- profile$TaxIDs %in% tree_tips
  sample_cols <- setdiff(colnames(profile), "TaxIDs")
  coverages <- numeric(length(sample_cols))
  names(coverages) <- sample_cols
  for (j in seq_along(sample_cols)) {
    col <- sample_cols[j]
    total <- sum(as.numeric(profile[[col]]))
    if (total > 0) {
      in_tree_total <- sum(as.numeric(profile[[col]][in_tree]))
      coverages[j] <- in_tree_total / total
    } else {
      coverages[j] <- NA
    }
  }
  return(coverages)
}

all_results <- data.frame(environment = character(), type = character(),
                          tool = character(), sample = character(),
                          abundance_in_tree = numeric(),
                          stringsAsFactors = FALSE)

sim_envs <- list(
  list(name="soil", label_converts=label_converts_soil),
  list(name="Ocean", label_converts=label_converts_ocean),
  list(name="gut", label_converts=label_converts_gut)
)

for (env in sim_envs) {
  cat("Processing simulation environment:", env$name, "\n")
  
  profiles <- process_profiles(preprocess(env$name, 7, FALSE), env$label_converts)
  
  truth <- renormalize(remove_unclassified(preprocess_simulation_truths(env$name, 7)))
  truth <- update_labels(truth, env$label_converts)
  colnames(truth)[1] <- "TaxIDs"
  
  n_tools <- length(profiles)
  tnames <- if (n_tools == 11) tool_names else tool_names_no_metaspades[1:n_tools]
  
  for (i in seq_along(profiles)) {
    covs <- calc_tree_coverage(profiles[[i]], tree_tips)
    covs <- covs[!is.na(covs)]
    if (length(covs) > 0) {
      all_results <- rbind(all_results, data.frame(
        environment = env$name, type = "simulation",
        tool = tnames[i], sample = names(covs),
        abundance_in_tree = covs, stringsAsFactors = FALSE))
    }
  }
  
  truth_covs <- calc_tree_coverage(truth, tree_tips)
  truth_covs <- truth_covs[!is.na(truth_covs)]
  if (length(truth_covs) > 0) {
    all_results <- rbind(all_results, data.frame(
      environment = env$name, type = "simulation_truth",
      tool = "Truth", sample = names(truth_covs),
      abundance_in_tree = truth_covs, stringsAsFactors = FALSE))
  }
}

real_datasets <- list.files("analysis/real_data_outputs")
real_datasets <- real_datasets[real_datasets != "overall_dists"]

for (dataset in real_datasets) {
  cat("Processing real dataset:", dataset, "\n")
  
  profiles <- process_profiles(preprocess(dataset, 7, TRUE))
  
  n_tools <- length(profiles)
  tnames <- if (n_tools == 11) tool_names else tool_names_no_metaspades[1:n_tools]
  
  for (i in seq_along(profiles)) {
    covs <- calc_tree_coverage(profiles[[i]], tree_tips)
    covs <- covs[!is.na(covs)]
    if (length(covs) > 0) {
      all_results <- rbind(all_results, data.frame(
        environment = dataset, type = "real",
        tool = tnames[i], sample = names(covs),
        abundance_in_tree = covs, stringsAsFactors = FALSE))
    }
  }
}

cat("\n========== RESULTS ==========\n\n")

cat(sprintf("Total tool*sample observations: %d\n", nrow(all_results)))
cat(sprintf("Overall: median=%.1f%%, mean=%.1f%%, min=%.1f%%, max=%.1f%%\n",
            median(all_results$abundance_in_tree) * 100,
            mean(all_results$abundance_in_tree) * 100,
            min(all_results$abundance_in_tree) * 100,
            max(all_results$abundance_in_tree) * 100))

cat(sprintf("  5th percentile=%.1f%%, 25th=%.1f%%, 75th=%.1f%%, 95th=%.1f%%\n",
            quantile(all_results$abundance_in_tree, 0.05) * 100,
            quantile(all_results$abundance_in_tree, 0.25) * 100,
            quantile(all_results$abundance_in_tree, 0.75) * 100,
            quantile(all_results$abundance_in_tree, 0.95) * 100))

cat("\n--- By type ---\n")
for (tp in unique(all_results$type)) {
  sub <- all_results[all_results$type == tp,]
  cat(sprintf("  %s: n=%d, median=%.1f%%, mean=%.1f%%, [%.1f%%–%.1f%%]\n",
              tp, nrow(sub),
              median(sub$abundance_in_tree) * 100,
              mean(sub$abundance_in_tree) * 100,
              min(sub$abundance_in_tree) * 100,
              max(sub$abundance_in_tree) * 100))
}

cat("\n--- By environment ---\n")
for (env_name in unique(all_results$environment)) {
  sub <- all_results[all_results$environment == env_name,]
  cat(sprintf("  %s: n=%d, median=%.1f%%, mean=%.1f%%, [%.1f%%–%.1f%%]\n",
              env_name, nrow(sub),
              median(sub$abundance_in_tree) * 100,
              mean(sub$abundance_in_tree) * 100,
              min(sub$abundance_in_tree) * 100,
              max(sub$abundance_in_tree) * 100))
}

cat("\n--- By tool (across all environments) ---\n")
for (tl in unique(all_results$tool)) {
  sub <- all_results[all_results$tool == tl,]
  cat(sprintf("  %s: n=%d, median=%.1f%%, mean=%.1f%%, [%.1f%%–%.1f%%]\n",
              tl, nrow(sub),
              median(sub$abundance_in_tree) * 100,
              mean(sub$abundance_in_tree) * 100,
              min(sub$abundance_in_tree) * 100,
              max(sub$abundance_in_tree) * 100))
}
