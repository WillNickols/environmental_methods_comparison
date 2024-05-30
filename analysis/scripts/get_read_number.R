#! usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   get_read_number.R [-i <input> -o <output>]

Options:
   -i abundance_by_X directory
   -o output.tsv

' -> doc 

opts <- docopt(doc)

# Read folder of abundances
in_dir <- gsub("/$", "", opts$i)
files <- list.files(path=in_dir, full.names=TRUE)

read_nums <- files[grep("read_num",files)]

mapped_props <- c()
for (read_num in read_nums) {
  text <- readChar(read_num, file.info(read_num)$size)
  temp <- strsplit(text, "\n")[[1]] %>% as.numeric()
  sample <- gsub(".*\\/", "", read_num) %>% gsub(pattern="\\.mapped_read_num.*", replacement="")
  mapped_props <- rbind(mapped_props, c(sample, temp[2]))
}

colnames(mapped_props) <- c("ID", "num")
mapped_props <- as.data.frame(mapped_props)
mapped_props$num <- as.numeric(as.character(mapped_props$num))

write.table(mapped_props, opts$o, quote=FALSE, sep="\t", row.names=FALSE)
