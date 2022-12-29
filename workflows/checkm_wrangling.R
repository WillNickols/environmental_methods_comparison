#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
rm(list = ls())

'Usage:
   checkm_wrangling.R [-i <input> --n50 <n50_table> -o <output>]

Options:
   -i directory containing CheckM QA tables
   --n50 N50 table
   -o output

' -> doc 

opts <- docopt(doc)

in_dir <- gsub("/$", "", opts$i)
files <- list.files(path=in_dir, full.names=TRUE)

checkm <- c()
for (file in files) {
  
  data <- fread(file) %>%
    rename(bin_id=1, completeness=12, contamination=13, strain_heterogeneity=14) %>%
    select(1, 12, 13, 14) %>%
    filter(!grepl("\\.bin\\.unbinned", bin_id)) %>%
    mutate(quality = case_when(
      completeness >= 90 & contamination <= 5 ~ "high_quality",
      (completeness >= 90 & contamination <= 10 & contamination >= 5) | (completeness >= 50 & completeness <= 90 & contamination <= 10) ~ "medium_quality",
      completeness < 50 | contamination > 10 ~ "low_quality"
      )) %>%
    select(bin_id, quality, completeness, contamination, strain_heterogeneity)
  
  checkm <- rbind(checkm, data)
  
}

n50 <- fread(opts$n50) %>%
  rename(bin_id=1, n50=2) %>%
  mutate(bin_id = gsub("\\.fa", "", gsub(".*\\/", "", bin_id)))

combo <- merge(checkm, n50, by="bin_id")

write.table(combo, opts$o, quote=FALSE, sep="\t", row.names=FALSE)

#