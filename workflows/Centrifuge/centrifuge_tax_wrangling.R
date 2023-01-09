#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
rm(list = ls())

'Usage:
   centrifuge_tax_wrangling.R [-i <input> --tax <tax> --qa <combined> -o <output>]

Options:
   -i centrifuge directory
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
  sample <- gsub(".*centrifuge\\/", "", read_num) %>% gsub(pattern="\\.mapped_read_num.*", replacement="")
  total_reads <- rbind(total_reads, c(sample, readChar(read_num, file.info(read_num)$size) %>% as.numeric()))
}

total_reads <- data.frame(total_reads)
colnames(total_reads) <- c("ID", "reads")

# Read in merged centrifuge tsv
centrifuge <- read.csv(files[grep("centrifuge_merged_tmp.tsv",files)], sep="\t", header = T)
colnames(centrifuge) <- gsub("_report_centrifuge_species.txt", "", colnames(centrifuge))
colnames(centrifuge)[1] <- "Taxa"
centrifuge$Taxa <- tolower(centrifuge$Taxa) %>% gsub(pattern=" ",replacement="_")

for(i in 2:ncol(centrifuge)) {
  reads <- sum(centrifuge[which(centrifuge$Taxa=="k__bacteria"),i], 
                     centrifuge[which(centrifuge$Taxa=="k__archaea"),i],
                     centrifuge[which(centrifuge$Taxa=="k__viruses"),i],
                     centrifuge[which(centrifuge$Taxa=="k__eukaryota"),i])
  total_reads[total_reads$ID==colnames(centrifuge)[i],3] <- reads
  centrifuge[,i] <- 100 * centrifuge[,i]/reads
}

colnames(total_reads) <- c("ID", "total", "reads")
total_reads$total <- as.numeric(levels(total_reads$total))[total_reads$total]

centrifuge[nrow(centrifuge) + 1,1] <- c("UNCLASSIFIED")
centrifuge[nrow(centrifuge), 2:ncol(centrifuge)] <- rep(0, ncol(centrifuge) - 1)
centrifuge <- centrifuge[order(centrifuge$Taxa),]

for (ID in total_reads$ID) {
  prop <- total_reads$reads[total_reads$ID==ID]/total_reads$total[total_reads$ID==ID]
  centrifuge[,ID] <- centrifuge[,ID] * prop
  centrifuge[centrifuge$Taxa=="UNCLASSIFIED", ID] <- centrifuge[centrifuge$Taxa=="UNCLASSIFIED", ID] + 
    100 * (1-prop)
}

write.table(centrifuge, opts$o, quote=FALSE, sep="\t", row.names=FALSE)

#
