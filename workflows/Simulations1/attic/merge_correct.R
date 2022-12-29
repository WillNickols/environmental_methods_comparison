#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   merge_correct.R [--tax_dir <tax_dir>]

Options:
   --tax_dir tax_dir

' -> doc

opts <- docopt(doc)

file_list <- list.files(opts$tax_dir)

abundances <- fread(paste0(opts$tax_dir, file_list[1]))[,c(3,5)]
colnames(abundances) <- c("Taxa", gsub("\\.txt", "", file_list[1]))
for (file in file_list[-1]) {
  tmp_df <- fread(paste0(opts$tax_dir, file))[,c(3,5)]
  colnames(tmp_df) <- c("Taxa", gsub("\\.txt", "", file))
  abundances <- full_join(abundances, tmp_df, by="Taxa")
}

abundances[is.na(abundances)] <- 0

dir.create(paste0(opts$tax_dir, "merged/"))
write.table(abundances, paste0(opts$tax_dir, "merged/merged.tsv"), sep = "\t", row.names = F, col.names = T)
