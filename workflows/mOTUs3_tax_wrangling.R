#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
rm(list = ls())

'Usage:
   mOTUs3_tax_wrangling.R [-i <input> -o <output>]

Options:
   -i merged sample profiles
   -o output

' -> doc 

opts <- docopt(doc)

motus <- read.csv(gsub("/$", "", opts$i), sep="\t", header = T)

motus[,1] <- gsub(" \\[ref.*", "", motus[,1])

write.table(motus, opts$o, quote=FALSE, sep="\t", row.names=FALSE)

#
