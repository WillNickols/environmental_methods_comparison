#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   identify_high_quality_unknowns.R [--ext <extension> --tax <tax> --qa <combined> -o <output>]

Options:
   --ext extension
   --tax phylophlan tax
   --qa checkm qa combined
   -o output

' -> doc 

opts <- docopt(doc)

# Merge in bin taxonomies
map <- fread(gsub("/$", "", opts$tax), skip = 3)
map <- map[,-5]
colnames(map) <- c("ID", "Species", "Genus", "Family")
map$Taxonomy <- case_when(sapply(str_split(map$Species, ":"), `[`, 4) %>% as.numeric <= 0.05 ~ sapply(str_split(map$Species, ":"), `[`, 3),
                          sapply(str_split(map$Genus, ":"), `[`, 4) %>% as.numeric <= 0.15 ~ sapply(str_split(map$Genus, ":"), `[`, 3),
                          sapply(str_split(map$Family, ":"), `[`, 4) %>% as.numeric <= 0.3 ~ sapply(str_split(map$Family, ":"), `[`, 3),
                          TRUE ~ "known_kingdom")
map$Taxonomy <- gsub("(.*):(.*):(.*):(.*)", "\\3", map$Taxonomy)
map <- data.frame(ID=map$ID, Taxonomy=map$Taxonomy, NCBI_ID=sapply(str_split(map$Species, ":"), `[`, 3) %>% gsub(pattern = "\\|.*", replacement = ""))

# Keep only high quality
qa <- fread(gsub("/$", "", opts$qa))
map <- left_join(map, qa, by=c("ID"="bin_id"))

# Reformat
map$Taxonomy <- gsub(pattern=" \\(root\\)", replacement="", map$Taxonomy)
map <- mutate(map, Taxonomy = case_when(
  grepl("unbinned", ID) ~ "UNKNOWN",
  grepl("low_quality", quality) ~ "UNKNOWN",
  grepl("medium_quality", quality) ~ "UNKNOWN",
  is.na(quality) ~ "UNKNOWN",
  TRUE ~ Taxonomy),
  NCBI_ID = case_when(
    grepl("k__Bacteria", NCBI_ID) ~ "2",
    grepl("k__Eukaryota", NCBI_ID) ~ "2759",
    grepl("k__Viruses", NCBI_ID) ~ "10239",
    grepl("k__Archaea", NCBI_ID) ~ "2157",
    TRUE ~ toString(NCBI_ID))
  )

map <- map[map$Taxonomy == "known_kingdom",]

if(nrow(map) > 0) {
  map$filepath <- paste0(opts$ext, map$ID, ".fa")
  
  write.table(map, opts$o, quote=FALSE, sep="\t", row.names=FALSE)
} else {
  print("No high quality unknowns")
}

#Rscript identify_high_quality_unknowns.R --ext mags/checkm/ --tax mags/phylophlan/combined.tsv --qa mags/merged/checkm_qa_combined.tsv -o to_move.tsv
