#!/usr/bin/env Rscript

library(docopt)
library(tidyverse)
library(data.table)

'Usage:
   merge_sgb_abundances.R [--in_dir <in_dir> --out_dir <out_dir>]

Options:
   --in_dir Input directory
   --out_dir Output directory
   
' -> doc 

opts <- docopt(doc)

# input
in_dir <- gsub("/$", "", opts$in_dir)

# output
out_dir <- paste0(gsub("/$", "", opts$out_dir), "/")

# loop through the inputs
extensions = c("sgb_abundance.tsv", "n_mapped_reads.txt")

for (ext in extensions) {
  
  files <- list.files(
    in_dir,
    pattern = ext,
    full.names = TRUE
  )
  
  merged <- c()
  
  for (file in files) {
    
    if (ext == "sgb_abundance.tsv") {
      
      data <- fread(file)
      
      merged <- rbind(merged, data)
      
    } else {
      
      data <- fread(file)
      
      sample_id <- gsub(".*\\/", "", gsub(paste0("\\.", ext), "", file))
      
      if (identical(colnames(data), c("total_reads", "mapped_reads",  "percent"))) {
        
        data <- data %>%
          mutate(sample = sample_id) %>%
          select(sample, percent, total_reads)
        
        merged <- rbind(merged, data)
        
      }
    }
  }
  
  write.table(merged, paste0(out_dir, "merged_", ext), quote = FALSE, sep = "\t", row.names = FALSE)
  
}