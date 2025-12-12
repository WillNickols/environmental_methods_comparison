##Script for generating taxonomic trees


library(taxonomizr)
#load in database
taxonomizr::prepareDatabase(getAccessions=FALSE)
source("analysis/scripts/helpers.R")

#make tree for gut data
dataset_gut <- "gut"
#set level to species
level <- 7
#set abundance threshold to 0.1%
threshold <- .1

#generate the profiles and only keep NCBI taxa
profiles <- remove_non_ncbi(preprocess(dataset_gut, level, FALSE))
#grab the truth samples and remove the unclassified portion
truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset_gut, level)))
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


## get all IDs on the profiles
all_ids_gut <- unique(get_all_present_IDs(profiles, truth))

#get taxonomic names
full_taxonomy_gut <- getTaxonomy(all_ids_gut, desiredTaxa = c("domain", "phylum", "class", "order", "family", "genus",
                                                      "species"))
full_taxonomy_gut <- cbind(full_taxonomy_gut, rownames(full_taxonomy_gut))
colnames(full_taxonomy_gut)[8] <- "ID"

#remove viruses
virus_rows <- grep("virus", ignore.case = T, full_taxonomy_gut[,7])
full_taxonomy_gut <- full_taxonomy_gut[-virus_rows,]

#remove phages
phage_rows <- grep("phage", ignore.case = T, full_taxonomy_gut[,7])
full_taxonomy_gut <- full_taxonomy_gut[-phage_rows,]

##start by cleaning up those that are all NA by attempting to remap them..
rows_all_NA <- which(apply(is.na(full_taxonomy_gut[,-8]), 1, all))
full_taxonomy_gut <- full_taxonomy_gut[-rows_all_NA,]

##we can just remap these to their newer IDs?
#262071 is a virus.
new_ids_to_search <- c('2841037'=3028070, '1404649'=2823807, '2884447'=2792859, '1532'=33035, '1841867'=2897707, 
                       '1118061'=2585118, '1453999'=2954383, '1457154'=2954382, '1532'=33035, '1898205'=2485925,
                       '1960941'=2942470)
new_ids_to_search <- getTaxonomy(new_ids_to_search, desiredTaxa = c("domain", "phylum", "class", "order", "family", "genus",
                                                                    "species"))

new_ids_to_search <- cbind(new_ids_to_search, rownames(new_ids_to_search))
colnames(new_ids_to_search)[8] <- "ID"
full_taxonomy_gut <- rbind(full_taxonomy_gut, new_ids_to_search)

#grab rows that need label fixing
rows_with_any_NA <- full_taxonomy_gut[apply(is.na(full_taxonomy_gut), 1, any),]

any_NA_index <- which(apply(is.na(full_taxonomy_gut), 1, any))

full_taxonomy_gut <- full_taxonomy_gut[-any_NA_index,]

fixed_taxonomy_gut <- propogate_labels(rows_with_any_NA)

full_taxonomy_gut <- rbind(full_taxonomy_gut, fixed_taxonomy_gut)


which(duplicated(full_taxonomy_gut))
full_taxonomy_gut <- full_taxonomy_gut[-which(duplicated(full_taxonomy_gut)),]

full_taxonomy_tree <- makeNewick(full_taxonomy_gut, excludeTerminalNAs = T, quote = "'")
writeLines(text = full_taxonomy_tree, "analysis/databases/Taxonomic_trees/gut_simulation_tree.nwk")

