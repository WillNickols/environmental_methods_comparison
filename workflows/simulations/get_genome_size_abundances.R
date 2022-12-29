#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   get_genome_size_abundances.R [--k <k> --sigma <sigma> --num_species <num_species> --genome_to_id <genome_to_id.tsv> --profile <profile.tsv> --GC_and_lengths <GC_and_lengths.tsv> -o <output>]

Options:
   --k k
   --sigma sigma
   --num_species num_species
   --genome_to_id genome_to_id.tsv
   --profile profile.tsv
   --GC_and_lengths GC_and_lengths.tsv
   -o output

' -> doc

opts <- docopt(doc)

k <- as.numeric(opts$k)
n <- as.numeric(opts$num_species)

genome_to_id = fread(gsub("/$", "", opts$genome_to_id), header = FALSE)

colnames(genome_to_id) <- c("ID", "filepath")
genome_to_id$whole_ID <- genome_to_id$ID
genome_to_id$ID <- gsub("\\..*", "", genome_to_id$ID) %>% gsub(pattern = " ", replacement = "_")

profile = fread(gsub("/$", "", opts$profile), header = TRUE) %>%
  mutate(genome_ID = gsub(" ", "_", genome_ID))

lengths_df = fread(gsub("/$", "", opts$GC_and_lengths), header = TRUE)

true_unknown = sum(profile$abundance[grepl("UNKNOWN", profile$genome_ID)])

merged_genomes_full = genome_to_id %>%
  left_join(profile, by = c("ID"="genome_ID")) %>% 
  mutate(filepath = coalesce(filepath.x, filepath.y)) %>% 
  select(-filepath.x, -filepath.y) %>%
  select(whole_ID, ID, abundance, filepath) %>%
  left_join(lengths_df, by=c("filepath"="filename")) %>%
  mutate(abundance = ifelse(grepl("UNKNOWN", ID), NA, abundance))

merged_genomes_full$abundance = merged_genomes_full$abundance * (1 - true_unknown) / sum(merged_genomes_full$abundance, na.rm=T)

# Keep only the highest quality SGBs
unknown_needed = 480 - sum(grepl("NEW_", merged_genomes_full$ID))
pot_drops = which(grepl("UNKNOWN", merged_genomes_full$ID))
if (length(pot_drops) >= unknown_needed) {
  pot_drops = c()
} else {
  pot_drops = pot_drops[-c(1:unknown_needed)]
}
merged_genomes_full = merged_genomes_full[!(1:nrow(merged_genomes_full) %in% pot_drops),]

known_species = sum(!grepl("UNKNOWN", merged_genomes_full$ID))
unknown_species = sum(grepl("UNKNOWN", merged_genomes_full$ID))

if (known_species + unknown_species < n) {
  stop((paste0("Not enough species: ", n, " needed, but ", known_species + unknown_species, " provided.")))
}

seq_to_fit = vector(length = n)
seq_to_fit[1] = 1/k
for (j in 2:n) {
  seq_to_fit[j] <- seq_to_fit[j - 1] * (k-1)/k
}

seq_to_fit <- seq_to_fit/mean(seq_to_fit) * median(merged_genomes_full$length)

genome_lengths = merged_genomes_full[!grepl("UNKNOWN", merged_genomes_full$ID),c("ID", "length")]
genome_lengths = genome_lengths[order(genome_lengths$length),]
genome_unknown_lengths = merged_genomes_full[grepl("UNKNOWN", merged_genomes_full$ID),c("ID", "length")]
genome_unknown_lengths = genome_unknown_lengths[order(genome_unknown_lengths$length),]

sequence_fit_mat = cbind(seq_to_fit, rep(NA, length(seq_to_fit)), rep(NA, length(seq_to_fit)))

# Fill with known genome lengths
if (nrow(genome_lengths) > 0) {
  for (i in 1:nrow(genome_lengths)) {
    dist_order = order(abs(as.numeric(sequence_fit_mat[,1]) - as.numeric(genome_lengths[i,"length"])))
    j=1
    while (!is.na(sequence_fit_mat[dist_order[j],2])) {
      j <- j+1
    }
    if (j < length(dist_order)) {
      sequence_fit_mat[dist_order[j],2:3] <- unlist(genome_lengths[i,c("length", "ID")])
    }
  }
}

unknown_length_mat <- cbind(genome_unknown_lengths$length, rep(NA, length(genome_unknown_lengths$length)), genome_unknown_lengths$ID)

for (i in 1:nrow(sequence_fit_mat)) {
  if (is.na(sequence_fit_mat[i,2])) {
    dist_order = order(abs(as.numeric(unknown_length_mat[,1]) - as.numeric(sequence_fit_mat[i,1])))
    j=1
    while (!is.na(unknown_length_mat[dist_order[j],2])) {
      j <- j+1
    }
    sequence_fit_mat[i,2:3] <- unknown_length_mat[dist_order[j],c(1,3)]
    unknown_length_mat[dist_order[j],2] <- 1
  }
}

sequence_fit <- data.frame(sequence_fit_mat)
colnames(sequence_fit) <- c("sequence", "length", "ID")

merged_genomes_full = sequence_fit %>%
  left_join(merged_genomes_full, by=c("ID")) %>% 
  select(-length.x, -length.y, -sequence)

unknown_abundances <- rlnorm(sum(grepl("UNKNOWN", merged_genomes_full$ID)), -3, as.numeric(opts$sigma))
known_abundance <- sum(merged_genomes_full$abundance[!grepl("UNKNOWN", merged_genomes_full$ID)], na.rm = T)
known_abundance <- ifelse(abs(known_abundance - 1) < 0.00001, 1, known_abundance)
unknown_abundances <- unknown_abundances / (sum(unknown_abundances)/(1-known_abundance))

if (known_abundance < 1) {
  abundances_to_merge <- data.frame(ID=merged_genomes_full$ID[grepl("UNKNOWN", merged_genomes_full$ID)], 
                                    abundance=unknown_abundances)
  
  merged_genomes_full <- full_join(abundances_to_merge, merged_genomes_full, by=c("ID")) %>%
    mutate(abundance = coalesce(abundance.x, abundance.y)) %>% 
    select(whole_ID, abundance)
} else {
  merged_genomes_full <- merged_genomes_full %>%
    select(whole_ID, abundance)
}

write.table(merged_genomes_full, opts$o, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
