#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
rm(list = ls())

'Usage:
   kraken_tax_wrangling.R [-i <input> --tax <tax> --qa <combined> -o <output>]

Options:
   -i kraken directory
   -o output

' -> doc 

opts <- docopt(doc)

# Read folder of abundances
in_dir <- gsub("/$", "", opts$i)
files <- list.files(path=in_dir, full.names=TRUE)

read_nums <- files[grep("read_num",files)]

# Read in total reads
total_reads <- c()
for (read_num in read_nums) {
  sample <- gsub(".*kraken\\/", "", read_num) %>% gsub(pattern="\\.mapped_read_num.*", replacement="")
  total_reads <- rbind(total_reads, c(sample, readChar(read_num, file.info(read_num)$size) %>% as.numeric()))
}

total_reads <- data.frame(total_reads)
colnames(total_reads) <- c("ID", "reads")

# Read in merged kraken tsv
kraken <- read.csv(files[grep("kraken_merged_tmp.tsv",files)], sep="\t", header = T)
colnames(kraken) <- gsub("_report_kraken_bracken_species.txt", "", colnames(kraken))
colnames(kraken)[1] <- "Taxa"
kraken$Taxa <- tolower(kraken$Taxa) %>% gsub(pattern=" ",replacement="_")

for(i in 2:ncol(kraken)) {
  reads <- sum(kraken[which(kraken$Taxa=="k__bacteria"),i], 
                     kraken[which(kraken$Taxa=="k__archaea"),i],
                     kraken[which(kraken$Taxa=="k__viruses"),i],
                     kraken[which(kraken$Taxa=="k__eukaryota"),i])
  total_reads[total_reads$ID==colnames(kraken)[i],3] <- reads
  kraken[,i] <- 100 * kraken[,i]/reads
}

colnames(total_reads) <- c("ID", "total", "reads")
total_reads$total <- as.numeric(levels(total_reads$total))[total_reads$total]

kraken[nrow(kraken) + 1,1] <- c("UNCLASSIFIED")
kraken[nrow(kraken), 2:ncol(kraken)] <- rep(0, ncol(kraken) - 1)
kraken <- kraken[order(kraken$Taxa),]

for (ID in total_reads$ID) {
  prop <- total_reads$reads[total_reads$ID==ID]/total_reads$total[total_reads$ID==ID]
  kraken[,ID] <- kraken[,ID] * prop
  kraken[kraken$Taxa=="UNCLASSIFIED", ID] <- kraken[kraken$Taxa=="UNCLASSIFIED", ID] + 
    100 * (1-prop)
}

write.table(kraken, opts$o, quote=FALSE, sep="\t", row.names=FALSE)

#
