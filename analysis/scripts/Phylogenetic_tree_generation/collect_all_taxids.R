### Collect all unique NCBI species-level TaxIDs across simulation profiles,
### simulation truths, and real data profiles for phylogenetic tree construction.

source("analysis/scripts/helpers.R")

library(dplyr)

threshold <- .1

# Label remappings from each environment's taxonomic distance analysis
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

all_taxids <- c()

# Collect TaxIDs from simulation profiles and truths
sim_envs <- list(
  list(name="soil", label_converts=label_converts_soil),
  list(name="Ocean", label_converts=label_converts_ocean),
  list(name="gut", label_converts=label_converts_gut)
)

for (env in sim_envs) {
  cat("Processing simulation environment:", env$name, "\n")
  
  profiles <- remove_non_ncbi(preprocess(env$name, 7, FALSE))
  truth <- renormalize(remove_unclassified(preprocess_simulation_truths(env$name, 7)))
  
  profiles <- lapply(profiles, renormalize)
  profiles <- lapply(profiles, threshold_sample, threshold)
  profiles <- lapply(profiles, function(x) {
    if (length(which(rowSums(x[,-1]) == 0)) != 0)
      return(x[-which(rowSums(x[,-1]) == 0),])
    else
      return(x)
  })
  profiles <- lapply(profiles, renormalize)
  
  truth <- update_labels(truth, env$label_converts)
  profiles <- lapply(profiles, function(x) update_labels(x, env$label_converts))
  
  for (p in profiles) {
    all_taxids <- c(all_taxids, p$TaxIDs)
  }
  
  colnames(truth)[1] <- "TaxIDs"
  all_taxids <- c(all_taxids, truth$TaxIDs)
  
  # Also include the remapped target IDs
  all_taxids <- c(all_taxids, as.character(env$label_converts))
}

# Collect TaxIDs from real data profiles
real_datasets <- list.files("analysis/real_data_outputs")
real_datasets <- real_datasets[real_datasets != "overall_dists"]

for (dataset in real_datasets) {
  cat("Processing real dataset:", dataset, "\n")
  
  profiles <- remove_non_ncbi(preprocess(dataset, 7, TRUE))
  profiles <- lapply(profiles, renormalize)
  profiles <- lapply(profiles, threshold_sample, threshold)
  profiles <- lapply(profiles, function(x) {
    if (length(which(rowSums(x[,-1]) == 0)) != 0)
      return(x[-which(rowSums(x[,-1]) == 0),])
    else
      return(x)
  })
  profiles <- lapply(profiles, renormalize)
  
  for (p in profiles) {
    all_taxids <- c(all_taxids, p$TaxIDs)
  }
}

# Deduplicate and keep only valid numeric NCBI TaxIDs
all_taxids <- as.character(all_taxids)
all_taxids <- unique(all_taxids)
all_taxids <- all_taxids[!is.na(suppressWarnings(as.numeric(all_taxids)))]
all_taxids <- all_taxids[all_taxids != "UNCLASSIFIED"]

cat("Total unique NCBI species TaxIDs:", length(all_taxids), "\n")

dir.create("analysis/databases/Phylogenetic_trees", showWarnings = FALSE, recursive = TRUE)
writeLines(all_taxids, "analysis/databases/Phylogenetic_trees/all_species_taxids.txt")
cat("Written to analysis/databases/Phylogenetic_trees/all_species_taxids.txt\n")
