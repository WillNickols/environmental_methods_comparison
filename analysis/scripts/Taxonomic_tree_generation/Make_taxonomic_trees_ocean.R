##Script for generating taxonomic trees


library(taxonomizr)
#load in database
taxonomizr::prepareDatabase(getAccessions=FALSE)
source("analysis/scripts/helpers.R")

#make tree for gut data
dataset_ocean <- "ocean"
#set level to species
level <- 7
#set abundance threshold to 0.1%
threshold <- .1

#generate the profiles and only keep NCBI taxa
profiles <- remove_non_ncbi(preprocess(dataset_ocean, level, FALSE))
#grab the truth samples and remove the unclassified portion
truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset_ocean, level)))
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
all_ids_ocean <- unique(get_all_present_IDs(profiles, truth))

#get taxonomic names
full_taxonomy_ocean <- getTaxonomy(all_ids_ocean, desiredTaxa = c("domain", "phylum", "class", "order", "family", "genus",
                                                      "species"))
full_taxonomy_ocean <- cbind(full_taxonomy_ocean, rownames(full_taxonomy_ocean))
colnames(full_taxonomy_ocean)[8] <- "ID"

#remove viruses
virus_rows <- grep("virus", ignore.case = T, full_taxonomy_ocean[,7])
full_taxonomy_ocean <- full_taxonomy_ocean[-virus_rows,]

#remove phages
phage_rows <- grep("phage", ignore.case = T, full_taxonomy_ocean[,7])
full_taxonomy_ocean <- full_taxonomy_ocean[-phage_rows,]

##start by cleaning up those that are all NA by attempting to remap them..
rows_all_NA <- which(apply(is.na(full_taxonomy_ocean[,-8]), 1, all))
full_taxonomy_ocean <- full_taxonomy_ocean[-rows_all_NA,]

##we can just remap these to their newer IDs?
#262071 is a virus.
new_ids_to_search <- c('2841037'=3028070, '29570'=77097, '1404649'=2823807, '2884447'=2792859, '158080'=141390, '50340'=53407, '81037'=77608, 
                       '1304902'=2731756, '1453999'=2954383, '1457154'=2954382, '29570'=77097, '158080'=141390,
                       '50340'=53407, '81037'=77608, '1960941'=2942470)

new_ids_to_search <- getTaxonomy(new_ids_to_search, desiredTaxa = c("domain", "phylum", "class", "order", "family", "genus",
                                                                    "species"))

new_ids_to_search <- cbind(new_ids_to_search, rownames(new_ids_to_search))
colnames(new_ids_to_search)[8] <- "ID"
full_taxonomy_ocean <- rbind(full_taxonomy_ocean, new_ids_to_search)

#grab rows that need label fixing
rows_with_any_NA <- full_taxonomy_ocean[apply(is.na(full_taxonomy_ocean), 1, any),]

any_NA_index <- which(apply(is.na(full_taxonomy_ocean), 1, any))

full_taxonomy_ocean <- full_taxonomy_ocean[-any_NA_index,]

fixed_taxonomy_ocean <- propogate_labels(rows_with_any_NA)

full_taxonomy_ocean <- rbind(full_taxonomy_ocean, fixed_taxonomy_ocean)


which(duplicated(full_taxonomy_ocean))
full_taxonomy_ocean <- full_taxonomy_ocean[-which(duplicated(full_taxonomy_ocean)),]

full_taxonomy_ocean_df <- data.frame(full_taxonomy_ocean)
full_taxonomy_ocean_df <- cbind(root="root", full_taxonomy_ocean_df)

full_taxonomy_tree <- makeNewick(full_taxonomy_ocean_df, excludeTerminalNAs = T, quote = "'")
writeLines(text = full_taxonomy_tree, "analysis/databases/Taxonomic_trees/ocean_simulation_tree.nwk")
