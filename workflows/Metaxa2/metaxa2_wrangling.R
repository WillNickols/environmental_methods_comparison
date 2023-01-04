#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(readr)
library(stringr)
rm(list = ls())

'Usage:
   checkm_wrangling.R [-i <input> -o <output>]

Options:
   -i merged.tsv file
   -o output

' -> doc 

opts <- docopt(doc)

in_file <- gsub("/$", "", opts$i)

merged <- read_delim(file = in_file, delim = '\t', col_names = TRUE)

colnames(merged) <- gsub(pattern = "^(?:[^.]*/){1}", "", colnames(merged))

# Get sample names and remove "Taxa"
samples <- unique(gsub(pattern = "[\\.]*[0-9]*$", "", colnames(merged))[-1])

output <- data.frame(cbind(merged$Taxa))

for (sample in samples) {
  output <- cbind(output, rowSums(merged[,colnames(merged)[grepl(sample, colnames(merged))]]))
}

colnames(output) <- c("Taxa", samples)

rownames(output) <- merged$Taxa
output$Taxa <- merged$Taxa

output$taxGroup <- str_count(output$Taxa, ";")
output$Taxa <- gsub(" ", "_", output$Taxa)
output$Taxa <- case_when(output$taxGroup == 0 ~ paste0("k__", output$Taxa),
                         output$taxGroup == 1 ~ paste0("k__", gsub("^(.*?(;.*?){0});", "\\1|p__", output$Taxa)),
                         output$taxGroup == 2 ~ paste0("k__", gsub("^(.*?(;.*?){0});(.*?(;.*?){0});", 
                                                                   "\\1|p__\\3|c__", output$Taxa)),
                         output$taxGroup == 3 ~ paste0("k__", gsub("^(.*?(;.*?){0});(.*?(;.*?){0});(.*?(;.*?){0});", 
                                                                   "\\1|p__\\3|c__\\5|o__", output$Taxa)),
                         output$taxGroup == 4 ~ paste0("k__", gsub("^(.*?(;.*?){0});(.*?(;.*?){0});(.*?(;.*?){0});(.*?(;.*?){0});", 
                                                                   "\\1|p__\\3|c__\\5|o__\\7|f__", output$Taxa)),
                         output$taxGroup == 5 ~ paste0("k__", gsub("^(.*?(;.*?){0});(.*?(;.*?){0});(.*?(;.*?){0});(.*?(;.*?){0});(.*?(;.*?){0});", 
                                                                   "\\1|p__\\3|c__\\5|o__\\7|f__\\9|g__", output$Taxa)),
                         output$taxGroup == 6 ~ paste0("k__", gsub("^(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});", 
                                                                   "\\1|p__\\2|c__\\3|o__\\4|f__\\5|g__\\6|s__", output$Taxa)),
                         output$taxGroup == 7 ~ paste0("k__", gsub("^(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});", 
                                                                   "\\1|p__\\2|c__\\3|o__\\4|f__\\5|g__\\6|s__\\7|t__", output$Taxa)),
                         output$taxGroup == 8 ~ paste0("k__", gsub("^(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});(.*?(?:;.*?){0});", 
                                                                   "\\1|p__\\2|c__\\3|o__\\4|f__\\5|g__\\6|s__\\7|t__", output$Taxa)))


colsums <- colSums(output[output$taxGroup==0,-1])

for (i in 2:length(output)) {
  output[,i] <- output[,i]/colsums[i-1] * 100
}

output$taxGroup <- NULL

colnames(output) <- c("clade_name", paste0(colnames(output), "_taxonomic_profile")[-1])

write.table(output, opts$o, sep="\t", row.names=FALSE, quote = FALSE)

#
