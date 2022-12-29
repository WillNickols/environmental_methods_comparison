#!/usr/bin/env Rscript

# Need to deal with strains vs species

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   get_GC_abundances.R [--GC_target <GC_target> --num_species <num_species> --genome_to_id <genome_to_id.tsv> --profile <profile.tsv> --GC_and_lengths <GC_and_lengths.tsv> -o <output>]

Options:
   --GC_target GC_target
   --num_species num_species
   --genome_to_id genome_to_id.tsv
   --profile profile.tsv
   --GC_and_lengths GC_and_lengths.tsv
   -o output

' -> doc

opts <- docopt(doc)

genome_to_id = fread(gsub("/$", "", opts$genome_to_id), header = FALSE)

colnames(genome_to_id) <- c("ID", "filepath")
genome_to_id$whole_ID <- genome_to_id$ID
genome_to_id$ID <- gsub("\\..*", "", genome_to_id$ID)

profile = fread(gsub("/$", "", opts$profile), header = TRUE) %>%
  mutate(genome_ID = gsub(" ", "_", genome_ID))

gc_content = fread(gsub("/$", "", opts$GC_and_lengths), header = TRUE)

merged_genomes = genome_to_id %>%
  left_join(profile, by = c("ID"="genome_ID")) %>% 
  mutate(filepath = coalesce(filepath.x, filepath.y)) %>% 
  select(-filepath.x, -filepath.y) %>%
  select(whole_ID, ID, abundance, filepath) %>%
  left_join(gc_content, by=c("filepath"="filename")) %>%
  mutate(gc_contents = gc_contents/100,
         abundance = ifelse(grepl("UNKNOWN", ID), NA, abundance/unname(table(ID)[ID])))

known_abundance = sum(merged_genomes$abundance[!grepl("UNKNOWN", merged_genomes$ID)], na.rm = T)
unknown_abundance = 1 - known_abundance
known_gc_content = sum(merged_genomes$gc_contents[!grepl("UNKNOWN", merged_genomes$ID)] * 
                         merged_genomes$abundance[!grepl("UNKNOWN", merged_genomes$ID)], na.rm = T)

global_target = as.numeric(opts$GC_target)

target = (global_target - known_gc_content) / unknown_abundance

unknown_subset <- merged_genomes[grepl("UNKNOWN", merged_genomes$ID),]
unknown_subset <- unknown_subset[order(unknown_subset$gc_contents),]

unknown_species = as.numeric(opts$num_species) - sum(!is.na(merged_genomes$abundance[!grepl("UNKNOWN", merged_genomes$ID)]))

print("Dropping unknown genomes to meet genome number requirement")

counter = 0
if (length(unique(unknown_subset$ID)) > unknown_species) {
  while (length(unique(unknown_subset$ID)) > unknown_species) {
    counter = counter + 1
    if (counter > 1000) {
      stop("Overran in part 1")
    }
    to_drop = which.max(abs(unknown_subset$gc_contents - target))
    
    if (max(unknown_subset$gc_contents[-to_drop]) >= target & min(unknown_subset$gc_contents[-to_drop]) <= target) {
      unknown_subset = unknown_subset[-to_drop,]
    } else {
      if (unknown_subset$gc_contents[to_drop] >= target) {
        second_drop = which(unknown_subset$gc_contents==
                              unknown_subset$gc_contents[unknown_subset$gc_contents < target][which.min(unknown_subset$gc_contents[unknown_subset$gc_contents < target] - target)])
        if (max(unknown_subset$gc_contents[-second_drop]) >= target & min(unknown_subset$gc_contents[-second_drop]) <= target) {
          unknown_subset = unknown_subset[-second_drop,]
        }
      } else {
        second_drop = which(unknown_subset$gc_contents==
                              unknown_subset$gc_contents[unknown_subset$gc_contents >= target][which.min(unknown_subset$gc_contents[unknown_subset$gc_contents >= target] - target)])
        if (max(unknown_subset$gc_contents[-second_drop]) >= target & min(unknown_subset$gc_contents[-second_drop]) <= target) {
          unknown_subset = unknown_subset[-second_drop,]
        }
      }
    }
  }
} else if (length(unique(unknown_subset$ID)) < unknown_species) {
  stop(paste0("Invalid number of genomes: ", unknown_species, " needed, but ", length(unique(unknown_subset$ID)), " provided."))
}

if (target > max(unknown_subset$gc_contents) | target < min(unknown_subset$gc_contents)) {
  stop("Fitting impossible")
}

change = 0.01

fit_gc <- function(unknown_subset, target) {
  print("Starting fitting")
  
  weights = rep(1/length(unknown_subset$gc_contents), length(unknown_subset$gc_contents))
  if (sum(unknown_subset$gc_contents*weights) < target) {
    counter = 0
    while(sum(unknown_subset$gc_contents*weights) < target) {
      counter = counter + 1
      weights = weights + (unknown_subset$gc_contents - 
                             unknown_subset$gc_contents[which(unknown_subset$gc_contents>=target)[1]]) * change - 
        unknown_subset$gc_contents[which(unknown_subset$gc_contents>=target)[1]] * change
      weights[weights<0] = 0
      weights = weights/sum(weights)
      if (counter > 10000000) {
        stop("Overran in part 2")
      }
    }
    return(weights)
    
  } else if (sum(unknown_subset$gc_contents*weights) > target) {
    counter = 0
    while(sum(unknown_subset$gc_contents*weights) > target) {
      counter = counter + 1
      weights = weights - (unknown_subset$gc_contents*weights - 
                             unknown_subset$gc_contents * 
                             weights[which(unknown_subset$gc_contents*weights<=target)[sum(unknown_subset$gc_contents*weights<=target)]]) * 
        change - 
        unknown_subset$gc_contents*weights[which(unknown_subset$gc_contents*weights<=target)[sum(unknown_subset$gc_contents*weights<=target)]] * 
        change
      weights[weights<0] = 0.0001
      weights = weights/sum(weights)
      if (counter > 10000000) {
        stop("Overran in part 3")
      }
    }
    
    return(weights)
    
  } else {
    return(weights)
  }
}

weights <- fit_gc(unknown_subset, target)

unknown_subset$abundance <- weights * unknown_abundance

merged_genomes = merged_genomes %>%
  left_join(unknown_subset, by=c("ID", "filepath", "gc_contents", "length", "whole_ID")) %>% 
  mutate(abundance = coalesce(abundance.x, abundance.y),
         abundance = ifelse(is.na(abundance), 0, abundance)) %>% 
  select(whole_ID, abundance)

write.table(merged_genomes, opts$o, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
