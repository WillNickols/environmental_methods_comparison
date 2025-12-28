##Script for generating taxonomic trees


library(taxonomizr)
#load in database
taxonomizr::prepareDatabase(getAccessions=FALSE)
#setwd("~/Dropbox_Harvard/hutlab/Jacob/Repos/Enviro_revs/environmental_methods_comparison/")
source("analysis/scripts/helpers.R")

#make tree for gut data
dataset_soil <- "soil"
#set level to species
level <- 7
#set abundance threshold to 0.1%
threshold <- .1

#generate the profiles and only keep NCBI taxa
profiles <- remove_non_ncbi(preprocess(dataset_soil, level, FALSE))
#grab the truth samples and remove the unclassified portion
truth = renormalize(remove_unclassified(preprocess_simulation_truths(dataset_soil, level)))
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
all_ids_soil <- unique(get_all_present_IDs(profiles, truth))

#get taxonomic names
full_taxonomy_soil <- getTaxonomy(all_ids_soil, desiredTaxa = c("domain", "phylum", "class", "order", "family", "genus",
                                                      "species"))
full_taxonomy_soil <- cbind(full_taxonomy_soil, rownames(full_taxonomy_soil))
colnames(full_taxonomy_soil)[8] <- "ID"

#remove viruses
virus_rows <- grep("virus", ignore.case = T, full_taxonomy_soil[,7])
full_taxonomy_soil <- full_taxonomy_soil[-virus_rows,]

#remove phages
phage_rows <- grep("phage", ignore.case = T, full_taxonomy_soil[,7])
full_taxonomy_soil <- full_taxonomy_soil[-phage_rows,]

##start by cleaning up those that are all NA by attempting to remap them..
rows_all_NA <- which(apply(is.na(full_taxonomy_soil[,-8]), 1, all))
full_taxonomy_soil <- full_taxonomy_soil[-rows_all_NA,]

##we can just remap these to their newer IDs?
#262071 is a virus.
new_ids_to_search <- c('2599805'=931866, '1404649'=2823807, '2735433'=3014751, '2829818'=2952571, '2840472'=2840469, 
                       '2884447'=2792859, '80870'=80867, '1532'=33035, '2778071'=29523, '285567'=1942, '290385'=104623, '37480'=67362,
                       '423539'=3074428, '50340'=53407, '84292'=162393, '88075'=88074, '96101'=46163, '1497613'=1505087, 
                       '1915400'=67332, '335659'=1404864, '90270'=2754056, '80870'=80867, '1532'=33035, 
                       '285567'=1942, '290385'=104623, '37480'=67362, '423539'=3074428, '50340'=53407,
                       '84292'=162393, '88075'=88074, '96101'=46163, '335659'=1404864, '90270'=2754056)

new_ids_to_search <- getTaxonomy(new_ids_to_search, desiredTaxa = c("domain", "phylum", "class", "order", "family", "genus",
                                                                    "species"))

new_ids_to_search <- cbind(new_ids_to_search, rownames(new_ids_to_search))
colnames(new_ids_to_search)[8] <- "ID"
full_taxonomy_soil <- rbind(full_taxonomy_soil, new_ids_to_search)

#grab rows that need label fixing
rows_with_any_NA <- full_taxonomy_soil[apply(is.na(full_taxonomy_soil), 1, any),]

any_NA_index <- which(apply(is.na(full_taxonomy_soil), 1, any))

full_taxonomy_soil <- full_taxonomy_soil[-any_NA_index,]

fixed_taxonomy_soil <- propogate_labels(rows_with_any_NA)

full_taxonomy_soil <- rbind(full_taxonomy_soil, fixed_taxonomy_soil)


which(duplicated(full_taxonomy_soil))
full_taxonomy_soil <- full_taxonomy_soil[-which(duplicated(full_taxonomy_soil)),]

full_taxonomy_soil_df <- data.frame(full_taxonomy_soil)
full_taxonomy_soil_df <- cbind(root="root", full_taxonomy_soil_df)

##add 
full_taxonomy_tree <- makeNewick(full_taxonomy_soil_df, excludeTerminalNAs = T, quote = "'")
writeLines(text = full_taxonomy_tree, "analysis/databases/Taxonomic_trees/soil_simulation_tree.nwk")
