remove(list=ls())

library(dplyr)
library(ggplot2)
library(UpSetR)
library(gridExtra)
library(ggplotify)

setwd('..')
make_upset_plot <- function(tax_level, only_NCBI = FALSE) {
  df = read.csv(paste0('taxa_lists/',tax_level,'.tsv'), sep='\t')
  to_plot = as.list(df[, 1:(ncol(df) - 1)])
  to_plot = lapply(to_plot, function(group) return(group[!is.na(group)]) )
  if (only_NCBI) {
    to_plot = lapply(to_plot, function(group) return(group[!is.na(as.numeric(group))]) )
  }
  
  to_plot = lapply(to_plot, function(group) return(group[!grepl("unclassified|noname|incertae|/", group)]) )
  
  names(to_plot) <- c("Centrifuge", "GTDB-Tk 2", "Kraken 2 / Bracken 2", "MetaPhlAn2", "MetaPhlAn 3", "MetaPhlAn 4", "Metaxa 2", "mOTUs 3", "PhyloPhlAn 3")
  level <- case_when(tax_level == "kingdom" ~ "Kingdom",
                     tax_level == "phylum" ~ "Phylum",
                     tax_level == "class" ~ "Class",
                     tax_level == "order" ~ "Order",
                     tax_level == "family" ~ "Family",
                     tax_level == "genus" ~ "Genus",
                     tax_level == "species" ~ "Species")
  return (upset(fromList(to_plot), 
        order.by = "freq", 
        nsets = length(to_plot), 
        text.scale = 5,
        show.numbers = T,
        keep.order = T, 
        nintersects = 10,
        number.angles = 30,
        mainbar.y.label = paste0(level, " Intersection Size"),
       ))
}

phylaPlot = make_upset_plot('phylum')
familiesPlot = make_upset_plot('family')
generaPlot = make_upset_plot('genus')
speciesPlot = make_upset_plot('species')

png(file=paste0("figures/","All_Taxa_Database_Overlap.png"), width = 3600, height = 1200)
grid.arrange(as.grob(phylaPlot), as.grob(familiesPlot), as.grob(generaPlot), as.grob(speciesPlot), nrow=1, ncol=4)
dev.off()

phylaPlot = make_upset_plot('phylum', TRUE)
familiesPlot = make_upset_plot('family', TRUE)
generaPlot = make_upset_plot('genus', TRUE)
speciesPlot = make_upset_plot('species', TRUE)

png(file=paste0("figures/","NCBI_Taxa_Database_Overlap.png"), width = 3600, height = 1200)
grid.arrange(as.grob(phylaPlot), as.grob(familiesPlot), as.grob(generaPlot), as.grob(speciesPlot), nrow=1, ncol=4)
dev.off()
