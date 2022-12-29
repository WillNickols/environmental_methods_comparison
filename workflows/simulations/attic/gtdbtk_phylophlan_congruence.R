#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   gtdbtk_phylophlan_congruence.R [--phylophlan_relab <phylophlan_relab> --arc <gtdbtk.ar53.summary.tsv> --bac <gtdbtk.bac120.summary.tsv> --checkm <checkm_qa_and_n50.tsv> --phylophlan_map <phylophlan_renamed.tsv> --out <output.tsv>]

Options:
   --phylophlan_relab phylophlan_relab.tsv
   --arc gtdbtk.ar53.summary.tsv
   --bac gtdbtk.bac120.summary.tsv
   --checkm checkm_qa_and_n50.tsv
   --phylophlan_map phylophlan_renamed.tsv
   -o output

' -> doc

opts <- docopt(doc)

phylophlan_relab <- fread(gsub("/$", "", opts$phylophlan_relab), header = TRUE)[,c(2, 7)]
colnames(phylophlan_relab) <- c("ID", "phylophlan")

arc <- fread(gsub("/$", "", opts$arc))[,1:2]
bac <- fread(gsub("/$", "", opts$bac))[,1:2]
gtdbtk <- full_join(arc,bac)
colnames(gtdbtk) <- c("ID", "GTDBTK")
gtdbtk$GTDBTK <- gsub(";", "\\|", gtdbtk$GTDBTK) %>% gsub(pattern = "d__", replacement = "k__")

checkm <- fread(gsub("/$", "", opts$checkm), header = TRUE)[,c(1, 6)]

all_phylo <- left_join(gtdbtk, phylophlan_relab, by=c("ID")) %>% left_join(checkm, by=c("ID"="bin_id"))

subset_phylo <- all_phylo[all_phylo$keep=="keep",c(1,2,3)]

subset_phylo$phylophlan <- gsub("Candidatus_", "", subset_phylo$phylophlan)
subset_phylo$GTDBTK <- gsub("Candidatus_", "", subset_phylo$GTDBTK) %>% gsub(pattern="Unclassified", replacement = "UNKNOWN")

write.csv(subset_phylo, "output.csv")

get_overlap <- function(first, second) {
  if (str_count(first, "\\|") > 0 & str_count(second, "\\|") > 0) {
    for (i in 1:str_count(first, "\\|")) {
      if (! grepl(gsub("\\|", ";", substr(first, 1, gregexpr(pattern = '\\|',first)[[1]][i])), gsub("\\|", ";", second))) {
        i = i - 1
        break
      }
    }
    return(substr(first, 1, gregexpr(pattern = '\\|',first)[[1]][i] - 1))
  } else {
    if (first=="UNKNOWN" | second=="UNKNOWN") {
      return("UNKNOWN")
    }
  }
  
}

subset_phylo$agree <- unname(mapply(get_overlap, subset_phylo$GTDBTK, subset_phylo$phylophlan))

phylophlan_map <- fread(gsub("/$", "", opts$phylophlan_map), header = TRUE)
subset_phylo <- left_join(subset_phylo, phylophlan_map, by=c("agree"="Taxa"))
subset_phylo[is.na(subset_phylo$renamed)] <- "UNKNOWN"

write.table(subset_phylo, opts$o, sep="\t", row.names = F)
