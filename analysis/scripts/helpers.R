rm(list = ls())

library(plyr)
library(dplyr)
library(stringr)
library(vegan)

# Preprocessing to get lists of taxa
preprocess <- function(dataset, level, real=TRUE, by_dataset=FALSE) {
  # Real or simulated
  if (real) {
    path_prefix = "analysis/real_data_outputs/"
  } else {
    path_prefix = "analysis/simulation_outputs/"
  }
  
  # Check if dataset is cached
  if (by_dataset) {
    if(file.exists(paste0("analysis/scripts/cache/", dataset, level, "_by_dataset.rds"))) {
      return (readRDS(paste0("analysis/scripts/cache/", dataset, level, "_by_dataset.rds")))
    }
  } else {
    if(file.exists(paste0("analysis/scripts/cache/", dataset, level, ".rds"))) {
      return (readRDS(paste0("analysis/scripts/cache/", dataset, level, ".rds")))
    }
  }
  
  # Read Centrifuge and standardize naming
  centrifuge <- read.csv(paste0(path_prefix, dataset, "/Centrifuge/centrifuge_merged.tsv"), sep="\t", header = T)
  centrifuge$Taxa <- gsub(pattern=" ",replacement="_", centrifuge$Taxa) %>% tolower() %>% gsub(pattern="candidatus_", replacement="")
  
  last_tax <- str_split(centrifuge$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(centrifuge$Taxa, "\\|") + 1)) %>% tolower()
  
  double_species_list = centrifuge$Taxa[str_count(centrifuge$Taxa, "s__") > 1]
  remove_string = paste0(unique(str_remove(str_extract(double_species_list,"s__.*\\|"), "\\|$")), collapse='|')
  remove_string = gsub("\\.", "\\\\.", remove_string) %>% gsub(pattern="\\(", replacement="\\\\(") %>% gsub(pattern="\\)", replacement="\\\\)")
  centrifuge$Taxa[str_count(centrifuge$Taxa, "s__") > 1] = gsub("s__.*\\|", "", double_species_list)
  
  centrifuge <- centrifuge[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | centrifuge$Taxa=="UNCLASSIFIED") &
                     !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', centrifuge$Taxa) &
                       !grepl(remove_string, last_tax),]
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/centrifuge.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  centrifuge$TaxIDs <- mapvalues(centrifuge$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  
  centrifuge <- centrifuge[order(centrifuge$Taxa),]

  # Read GTDB-Tk for MEGAHIT
  if (by_dataset) {
    gtdbtk <- read.csv(paste0(path_prefix, dataset, "/MEGAHIT GTDB-Tk/final_profile_by_dataset.tsv"), sep="\t", header = T)
  } else {
    gtdbtk <- read.csv(paste0(path_prefix, dataset, "/MEGAHIT GTDB-Tk/final_profile_by_sample.tsv"), sep="\t", header = T)
  }
  
  gtdbtk$Taxa <- gsub(pattern=" ",replacement="_", gtdbtk$Taxonomy) %>% tolower()
  gtdbtk$Taxonomy <- NULL
  gtdbtk$Taxa[gtdbtk$Taxa == "unknown"] <- "UNCLASSIFIED"
  
  # Remove duplicate columns (bug)
  gtdbtk <- gtdbtk[,!grepl("\\.1$|\\.0$", colnames(gtdbtk))]
  
  # GTDB-Tk never actually leaves sgbs
  gtdbtk[gtdbtk$Taxa=="UNCLASSIFIED", -ncol(gtdbtk)] <- gtdbtk[gtdbtk$Taxa=="UNCLASSIFIED", -ncol(gtdbtk)] + 
    colSums(gtdbtk[grepl("sgb_", gtdbtk$Taxa),-ncol(gtdbtk)])
  gtdbtk <- gtdbtk[!grepl("sgb_", gtdbtk$Taxa),]
  
  last_tax <- str_split(gtdbtk$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(gtdbtk$Taxa, "\\|") + 1)) %>% tolower()
  
  gtdbtk <- gtdbtk[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | gtdbtk$Taxa=="UNCLASSIFIED") &
                     !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', gtdbtk$Taxa),]
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/GTDBTk.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  gtdbtk$TaxIDs <- mapvalues(gtdbtk$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  
  gtdbtk_megahit <- gtdbtk[order(gtdbtk$Taxa),]
  
  # Read GTDB-Tk for metaSPAdes
  run_metaspades = file.exists(paste0(path_prefix, dataset, "/metaSPAdes GTDB-Tk/final_profile_by_sample.tsv"))
  if (run_metaspades) {
    if (by_dataset) {
      gtdbtk <- read.csv(paste0(path_prefix, dataset, "/metaSPAdes GTDB-Tk/final_profile_by_dataset.tsv"), sep="\t", header = T)
    } else {
      gtdbtk <- read.csv(paste0(path_prefix, dataset, "/metaSPAdes GTDB-Tk/final_profile_by_sample.tsv"), sep="\t", header = T)
    }
    
    gtdbtk$Taxa <- gsub(pattern=" ",replacement="_", gtdbtk$Taxonomy) %>% tolower()
    gtdbtk$Taxonomy <- NULL
    gtdbtk$Taxa[gtdbtk$Taxa == "unknown"] <- "UNCLASSIFIED"
    
    # Remove duplicate columns (bug)
    gtdbtk <- gtdbtk[,!grepl("\\.1$|\\.0$", colnames(gtdbtk))]
    
    gtdbtk[gtdbtk$Taxa=="UNCLASSIFIED", -ncol(gtdbtk)] <- gtdbtk[gtdbtk$Taxa=="UNCLASSIFIED", -ncol(gtdbtk)] + 
      colSums(gtdbtk[grepl("sgb_", gtdbtk$Taxa),-ncol(gtdbtk)])
    gtdbtk <- gtdbtk[!grepl("sgb_", gtdbtk$Taxa),]
    
    last_tax <- str_split(gtdbtk$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(gtdbtk$Taxa, "\\|") + 1)) %>% tolower()
    
    gtdbtk <- gtdbtk[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | gtdbtk$Taxa=="UNCLASSIFIED") &
                       !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', gtdbtk$Taxa),]
    
    map <- read.csv(paste0("analysis/databases/standardized_databases/GTDBTk.tsv"), sep = "\t")
    map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
    gtdbtk$TaxIDs <- mapvalues(gtdbtk$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
    
    gtdbtk_metaspades <- gtdbtk[order(gtdbtk$Taxa),]
  } else {
    print("Skipping metaSPAdes")
  }
  
  # Read Kraken and standardize naming
  kraken <- read.csv(paste0(path_prefix, dataset, "/Kraken 2 Bracken 2/kraken_merged.tsv"), sep="\t", header = T)
  kraken$Taxa <- gsub(pattern=" ",replacement="_", kraken$Taxa) %>% tolower() %>% gsub(pattern="candidatus_", replacement="")
  
  last_tax <- str_split(kraken$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(kraken$Taxa, "\\|") + 1)) %>% tolower()
  
  kraken <- kraken[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | kraken$Taxa=="UNCLASSIFIED") &
                 !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', kraken$Taxa),]
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/KrakenBracken.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  kraken$TaxIDs <- mapvalues(kraken$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  
  kraken <- kraken[order(kraken$Taxa),]
  
  # Read MPA2 and standardize naming
  MPA2 <- read.csv(paste0(path_prefix, dataset, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, skip = 2)
  colnames(MPA2) <- read.csv(paste0(path_prefix, dataset, "/MetaPhlAn 2/metaphlan2_taxonomic_profiles.tsv"), sep="\t", header = F, nrows = 1)
  colnames(MPA2) <- gsub("_taxonomic_profile", "", colnames(MPA2))
  colnames(MPA2)[1] <- "Taxa"
  MPA2$Taxa <- gsub(pattern=" ",replacement="_", MPA2$Taxa) %>% tolower()
  if ("unclassified" %in% MPA2$Taxa) {
    MPA2$Taxa[MPA2$Taxa=="unclassified"] <- "UNCLASSIFIED"
  } else {
    MPA2[nrow(MPA2) + 1,] <- rep(0, ncol(MPA2))
    MPA2[nrow(MPA2),which(colnames(MPA2)=="Taxa")] <- "UNCLASSIFIED"
  }
  
  last_tax <- str_split(MPA2$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(MPA2$Taxa, "\\|") + 1)) %>% tolower()
  
  MPA2 <- MPA2[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | MPA2$Taxa=="UNCLASSIFIED") &
                 !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', MPA2$Taxa),]

  map <- read.csv(paste0("analysis/databases/standardized_databases/MetaPhlAn2.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  MPA2$TaxIDs <- mapvalues(MPA2$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  
  MPA2 <- MPA2[order(MPA2$Taxa),]
  
  # Read MPA3 and standardize naming
  MPA3 <- read.csv(paste0(path_prefix, dataset, "/MetaPhlAn 3/metaphlan3_taxonomic_profiles.tsv"), sep="\t", header = T, skip = 1)
  MPA3 <- MPA3[,-2]
  colnames(MPA3) <- gsub("_taxonomic_profile", "", colnames(MPA3))
  colnames(MPA3)[1] <- "Taxa"
  MPA3$Taxa <- gsub(pattern=" ",replacement="_", MPA3$Taxa) %>% tolower()
  MPA3$Taxa[MPA3$Taxa=="unknown"] <- "UNCLASSIFIED"
  
  last_tax <- str_split(MPA3$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(MPA3$Taxa, "\\|") + 1)) %>% tolower()
  MPA3 <- MPA3[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | MPA3$Taxa=="UNCLASSIFIED") &
                 !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', MPA3$Taxa),]
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/MetaPhlAn3.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  MPA3$TaxIDs <- mapvalues(MPA3$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  MPA3 <- MPA3[order(MPA3$Taxa),]
  
  # Read MPA4 and standardize naming
  MPA4 <- read.csv(paste0(path_prefix, dataset, "/MetaPhlAn 4/metaphlan4_taxonomic_profiles.tsv"), sep="\t", header = T)
  colnames(MPA4) <- gsub("_taxonomic_profile", "", colnames(MPA4))
  colnames(MPA4)[1] <- "Taxa"
  MPA4$NCBI_tax_id <- NULL
  MPA4$Taxa <- gsub(pattern=" ",replacement="_", MPA4$Taxa) %>% tolower()
  temp <- MPA4[MPA4$Taxa=="UNCLASSIFIED",]
  MPA4 <- MPA4[!MPA4$Taxa=="UNCLASSIFIED",]
  MPA4[nrow(MPA4) + 1,] <- c(0, colSums(temp[,-c(which(colnames(temp)=="Taxa"), which(colnames(temp)=="TaxIDs"))]))
  MPA4[nrow(MPA4),c(which(colnames(temp)=="Taxa"), which(colnames(temp)=="TaxIDs"))] <- "UNCLASSIFIED"
  
  last_tax <- str_split(MPA4$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(MPA4$Taxa, "\\|") + 1)) %>% tolower()
  MPA4 <- MPA4[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | MPA4$Taxa=="UNCLASSIFIED") &
                 !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', MPA4$Taxa),]
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/MetaPhlAn4.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  MPA4$TaxIDs <- mapvalues(MPA4$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  MPA4 <- MPA4[order(MPA4$Taxa),]
  
  # Read metaxa and standardize naming
  metaxa <- read.csv(paste0(path_prefix, dataset, "/Metaxa 2/merged_reorganized.tsv"), sep="\t", header = T)
  colnames(metaxa) <- gsub("_taxonomic_profile", "", colnames(metaxa))
  colnames(metaxa)[1] <- "Taxa"
  metaxa$Taxa <- gsub(pattern=" ",replacement="_", metaxa$Taxa) %>% tolower()
  metaxa = metaxa[!grepl(pattern="__$", metaxa$Taxa),]
  
  last_tax <- str_split(metaxa$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(metaxa$Taxa, "\\|") + 1)) %>% tolower()
  
  # Keep only bacteria, archaea, and viruses
  metaxa <- metaxa[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | metaxa$Taxa=="UNCLASSIFIED") &
                 !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', metaxa$Taxa) ,]
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/Metaxa2.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  metaxa$TaxIDs <- mapvalues(metaxa$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  metaxa <- metaxa[order(metaxa$Taxa),]
  
  # Read mOTUs3 and standardize naming
  motus <- read.csv(paste0(path_prefix, dataset, "/mOTUs3/mOTUs3_merged.tsv"), sep="\t", header = T)
  colnames(motus)[1] <- "Taxa"
  motus$Taxa <- gsub(pattern=" ",replacement="_", motus$Taxa) %>% tolower() %>% gsub(pattern="_\\[.*\\]",replacement="")
  motus[,-1] <- lapply(motus[,-1], function (column) as.numeric(column) * 100)
  
  # Insert unclassified
  motus[nrow(motus) + 1,] <- c(0, 
                               100 - colSums(motus[motus$Taxa=="k__bacteria"|
                                                     motus$Taxa=="k__archaea"|
                                                     motus$Taxa=="k__viruses",-1]))
  motus[nrow(motus),which(colnames(motus)=="Taxa")] <- "UNCLASSIFIED"
  
  last_tax <- str_split(motus$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(motus$Taxa, "\\|") + 1)) %>% tolower()
  motus <- motus[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | motus$Taxa=="UNCLASSIFIED") &
                     !grepl("noname|incertae|/|t__|sp\\.$", last_tax) & grepl('k__bacteria|k__archaea|k__viruses', motus$Taxa) ,]
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/mOTUs3.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  motus$TaxIDs <- mapvalues(motus$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  
  motus <- motus[order(motus$Taxa),]
  
  last_tax <- str_split(motus$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(motus$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  motus <- motus[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | motus$Taxa=="UNCLASSIFIED") &
                   !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", motus$Taxa),]
  
  # Read PhyloPhlAN and standardize naming
  if (by_dataset) {
    phylophlan <- read.csv(paste0(path_prefix, dataset, "/MEGAHIT PhyloPhlAn 3/final_profile_by_dataset.tsv"), sep="\t", header = T)
  } else {
    phylophlan <- read.csv(paste0(path_prefix, dataset, "/MEGAHIT PhyloPhlAn 3/final_profile_by_sample.tsv"), sep="\t", header = T)
  }
  
  colnames(phylophlan)[1] <- "Taxa"
  phylophlan$Taxa <- gsub(pattern=" ",replacement="_", phylophlan$Taxa) %>% tolower()
  phylophlan <- phylophlan[order(phylophlan$Taxa),]
  phylophlan$Taxa <- c(phylophlan$Taxa[-nrow(phylophlan)], "UNCLASSIFIED")
  
  phylophlan <- phylophlan[,!grepl("\\.1$|\\.0$", colnames(phylophlan))]
  
  last_tax <- str_split(phylophlan$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(phylophlan$Taxa, "\\|") + 1)) %>% tolower()
  phylophlan <- phylophlan[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | phylophlan$Taxa=="UNCLASSIFIED") &
                             !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses|^sgb', phylophlan$Taxa),]
  
  phylophlan$Taxa[grepl("^sgb", phylophlan$Taxa)] <- paste0("s__", phylophlan$Taxa[grepl("^sgb", phylophlan$Taxa)], "_MEGAHIT")
  
  map <- read.csv(paste0("analysis/databases/standardized_databases/PhyloPhlAn3.tsv"), sep = "\t")
  map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
  phylophlan$TaxIDs <- mapvalues(phylophlan$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
  phylophlan_megahit <- phylophlan[order(phylophlan$Taxa),]
  
  if (run_metaspades) {
    # Read phylophlan metaSPAdes
    if (by_dataset) {
      phylophlan <- read.csv(paste0(path_prefix, dataset, "/metaSPAdes PhyloPhlAn 3/final_profile_by_dataset.tsv"), sep="\t", header = T)
    } else {
      phylophlan <- read.csv(paste0(path_prefix, dataset, "/metaSPAdes PhyloPhlAn 3/final_profile_by_sample.tsv"), sep="\t", header = T)
    }
    
    colnames(phylophlan)[1] <- "Taxa"
    phylophlan$Taxa <- gsub(pattern=" ",replacement="_", phylophlan$Taxa) %>% tolower()
    phylophlan <- phylophlan[order(phylophlan$Taxa),]
    phylophlan$Taxa <- c(phylophlan$Taxa[-nrow(phylophlan)], "UNCLASSIFIED")
    
    phylophlan <- phylophlan[,!grepl("\\.1$|\\.0$", colnames(phylophlan))]
    
    last_tax <- str_split(phylophlan$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(phylophlan$Taxa, "\\|") + 1)) %>% tolower()
    phylophlan <- phylophlan[!is.na(last_tax) & last_tax != "" & (!grepl("unclassified", last_tax) | phylophlan$Taxa=="UNCLASSIFIED") &
                               !grepl("noname|incertae|/|t__", last_tax) & grepl('k__bacteria|k__archaea|k__viruses|^sgb', phylophlan$Taxa),]
    
    phylophlan$Taxa[grepl("^sgb", phylophlan$Taxa)] <- paste0("s__", phylophlan$Taxa[grepl("^sgb", phylophlan$Taxa)], "_metaSPAdes")
    
    map <- read.csv(paste0("analysis/databases/standardized_databases/PhyloPhlAn3.tsv"), sep = "\t")
    map$Taxa <- gsub(pattern=" ",replacement="_", map$Taxa) %>% tolower()
    phylophlan$TaxIDs <- mapvalues(phylophlan$Taxa, c(map$Taxa, "UNCLASSIFIED"), c(map$TaxID, "UNCLASSIFIED"), warn_missing = FALSE)
    phylophlan_metaspades <- phylophlan[order(phylophlan$Taxa),]
  }
  
  if (run_metaspades) {
    databases <- list(subset(centrifuge, select=c("Taxa",sort(colnames(centrifuge))[sort(colnames(centrifuge))!="Taxa"])),
                      subset(kraken, select=c("Taxa",sort(colnames(kraken))[sort(colnames(kraken))!="Taxa"])),
                      subset(MPA2, select=c("Taxa",sort(colnames(MPA2))[sort(colnames(MPA2))!="Taxa"])), 
                      subset(MPA3, select=c("Taxa",sort(colnames(MPA3))[sort(colnames(MPA3))!="Taxa"])), 
                      subset(MPA4, select=c("Taxa",sort(colnames(MPA4))[sort(colnames(MPA4))!="Taxa"])), 
                      subset(metaxa, select=c("Taxa",sort(colnames(metaxa))[sort(colnames(metaxa))!="Taxa"])),
                      subset(motus, select=c("Taxa",sort(colnames(motus))[sort(colnames(motus))!="Taxa"])), 
                      subset(gtdbtk_megahit, select=c("Taxa",sort(colnames(gtdbtk_megahit))[sort(colnames(gtdbtk_megahit))!="Taxa"])),
                      subset(gtdbtk_metaspades, select=c("Taxa",sort(colnames(gtdbtk_metaspades))[sort(colnames(gtdbtk_metaspades))!="Taxa"])),
                      subset(phylophlan_megahit, select=c("Taxa",sort(colnames(phylophlan_megahit))[sort(colnames(phylophlan_megahit))!="Taxa"])),
                      subset(phylophlan_metaspades, select=c("Taxa",sort(colnames(phylophlan_metaspades))[sort(colnames(phylophlan_metaspades))!="Taxa"]))
    )
  } else {
    databases <- list(subset(centrifuge, select=c("Taxa",sort(colnames(centrifuge))[sort(colnames(centrifuge))!="Taxa"])),
                      subset(kraken, select=c("Taxa",sort(colnames(kraken))[sort(colnames(kraken))!="Taxa"])),
                      subset(MPA2, select=c("Taxa",sort(colnames(MPA2))[sort(colnames(MPA2))!="Taxa"])), 
                      subset(MPA3, select=c("Taxa",sort(colnames(MPA3))[sort(colnames(MPA3))!="Taxa"])), 
                      subset(MPA4, select=c("Taxa",sort(colnames(MPA4))[sort(colnames(MPA4))!="Taxa"])), 
                      subset(metaxa, select=c("Taxa",sort(colnames(metaxa))[sort(colnames(metaxa))!="Taxa"])),
                      subset(motus, select=c("Taxa",sort(colnames(motus))[sort(colnames(motus))!="Taxa"])), 
                      subset(gtdbtk_megahit, select=c("Taxa",sort(colnames(gtdbtk_megahit))[sort(colnames(gtdbtk_megahit))!="Taxa"])),
                      subset(phylophlan_megahit, select=c("Taxa",sort(colnames(phylophlan_megahit))[sort(colnames(phylophlan_megahit))!="Taxa"]))
    )
  }
  
  
  if (length(unique(lapply(databases, colnames))) > 2) { # Only metaSPAdes samples can be different
    stop("Different column names between tools")
  }
  
  if (real) {
    dataset_name = case_when(dataset == "acid_mine" ~ "Acid Mine Drainage",
                             dataset == "animal_gut" ~ "Wild Animal Gut",
                             dataset == "cat_gut" ~ "Cat Gut",
                             dataset == "coastal_sediment" ~ "Coastal Sediment",
                             dataset == "dog_gut" ~ "Dog Gut",
                             dataset == "forest_soil" ~ "Forest Soil",
                             dataset == "gator_soil" ~ "Gator Nest Soil",
                             dataset == "human" ~ "Human Gut",
                             dataset == "saltmarsh" ~ "Salt Marsh",
                             dataset == "tara_polar" ~ "Tara Polar",)
    
    
    from_names <- colnames(databases[[1]])[!(colnames(databases[[1]]) %in% c("Taxa", "TaxIDs"))]
    replace_col_names <- function(x) {
      replace_names <- paste0(dataset_name, "_", 1:length(colnames(databases[[1]])[!grepl("Taxa|TaxIDs", colnames(databases[[1]]))]))
      colnames(x)[!grepl("Taxa|TaxIDs", colnames(x))] <- mapvalues(colnames(x)[!grepl("Taxa|TaxIDs", colnames(x))], from_names, replace_names, warn_missing = F)
      return (x)
    }
    
    databases <- lapply(databases, replace_col_names)
  }
  
  prefix_vec = c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Process results at the desired level
  
  # Kingdom level
  if (level == 1) {
    kingdoms <- vector("list", length(databases))
    for (i in 1:length(databases)) {
      kingdoms[[i]] <- databases[[i]][grepl(paste0(prefix_vec[1], collapse = "|"), databases[[i]]$Taxa) & !grepl(paste0(prefix_vec[-1], collapse = "|"), databases[[i]]$Taxa),]
      
      kingdoms[[i]] <- kingdoms[[i]][str_count(kingdoms[[i]]$Taxa, "k__") == 1, ]
      
      kingdoms[[i]][nrow(kingdoms[[i]]) + 1,] <- c(0, 100 - colSums(kingdoms[[i]][,-c(which(colnames(kingdoms[[i]])=="Taxa"), which(colnames(kingdoms[[i]])=="TaxIDs"))]), 0)
      kingdoms[[i]][nrow(kingdoms[[i]]),c(which(colnames(kingdoms[[i]])=="Taxa"), which(colnames(kingdoms[[i]])=="TaxIDs"))] <- "UNCLASSIFIED"
      kingdoms[[i]]$Taxa <- NULL
      kingdoms[[i]] <- kingdoms[[i]][, c("TaxIDs", colnames(kingdoms[[i]])[colnames(kingdoms[[i]])!="TaxIDs"] )]
      if (nrow(kingdoms[[i]]) > 1) {
        kingdoms[[i]] <- aggregate(.~TaxIDs,data=kingdoms[[i]],FUN=sum)
      }
    }
    if (by_dataset) {
      saveRDS(kingdoms, paste0("analysis/scripts/cache/", dataset, level, "_by_dataset.rds"))
    } else {
      saveRDS(kingdoms, paste0("analysis/scripts/cache/", dataset, level, ".rds"))
    }
    
    return (kingdoms)
  }
  
  if (level %in% 2:6) {
    taxa <- vector("list", length(databases))
    for (i in 1:length(databases)) {
      taxa[[i]] <- databases[[i]][grepl(paste0(prefix_vec[level], collapse = "|"), databases[[i]]$Taxa) & !grepl(paste0(prefix_vec[-(1:level)], collapse = "|"), databases[[i]]$Taxa),]
      
      taxa[[i]]$TaxIDs <- str_split(taxa[[i]]$TaxIDs, "\\|") %>% mapply(FUN = "[[", as.list(str_count(taxa[[i]]$TaxIDs, "\\|") + 1))
      
      taxa[[i]][nrow(taxa[[i]]) + 1,] <- c(0, 100 - colSums(taxa[[i]][,-c(which(colnames(taxa[[i]])=="Taxa"), which(colnames(taxa[[i]])=="TaxIDs"))]), 0)
      taxa[[i]][nrow(taxa[[i]]),c(which(colnames(taxa[[i]])=="Taxa"), which(colnames(taxa[[i]])=="TaxIDs"))] <- "UNCLASSIFIED"
      taxa[[i]]$Taxa <- NULL
      taxa[[i]] <- taxa[[i]][, c("TaxIDs", colnames(taxa[[i]])[colnames(taxa[[i]])!="TaxIDs"] )]
      if (nrow(taxa[[i]]) > 1) {
        taxa[[i]] <- aggregate(.~TaxIDs,data=taxa[[i]],FUN=sum)
      }
    }
    if (by_dataset) {
      saveRDS(taxa, paste0("analysis/scripts/cache/", dataset, level, "_by_dataset.rds"))
    } else {
      saveRDS(taxa, paste0("analysis/scripts/cache/", dataset, level, ".rds"))
    }
    return (taxa)
  }
  
  if (level == 7) {
    taxa <- vector("list", length(databases))
    for (i in 1:length(databases)) {
      taxa[[i]] <- databases[[i]][grepl(paste0(prefix_vec[level], collapse = "|"), databases[[i]]$Taxa),]
      
      taxa[[i]]$TaxIDs <- str_split(taxa[[i]]$TaxIDs, "\\|") %>% mapply(FUN = "[[", as.list(str_count(taxa[[i]]$TaxIDs, "\\|") + 1))
      
      taxa[[i]][nrow(taxa[[i]]) + 1,] <- c(0, 100 - colSums(taxa[[i]][,-c(which(colnames(taxa[[i]])=="Taxa"), which(colnames(taxa[[i]])=="TaxIDs"))]), 0)
      taxa[[i]][nrow(taxa[[i]]),c(which(colnames(taxa[[i]])=="Taxa"), which(colnames(taxa[[i]])=="TaxIDs"))] <- "UNCLASSIFIED"
      taxa[[i]]$Taxa <- NULL
      taxa[[i]] <- taxa[[i]][, c("TaxIDs", colnames(taxa[[i]])[colnames(taxa[[i]])!="TaxIDs"] )]
      if (nrow(taxa[[i]]) > 1) {
        taxa[[i]] <- aggregate(.~TaxIDs,data=taxa[[i]],FUN=sum)
      }
    }
    if (by_dataset) {
      saveRDS(taxa, paste0("analysis/scripts/cache/", dataset, level, "_by_dataset.rds"))
    } else {
      saveRDS(taxa, paste0("analysis/scripts/cache/", dataset, level, ".rds"))
    }
    return (taxa)
  }
}

# Preprocess all real datasets
preprocess_all_real <- function() {
  for (dataset in list.files("analysis/real_data_outputs")) {
    for (level in 1:7) {
      if (dataset != "overall_dists") {
        print(paste0("Processing ", dataset, " at level ", level))
        preprocess(dataset, level)
        preprocess(dataset, level, real=TRUE, by_dataset=TRUE)
        }
      }
  }
}

# Preprocess all simulated datasets
preprocess_all_simulated <- function() {
  for (dataset in list.files("analysis/simulation_outputs")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      preprocess(dataset, level, FALSE)
      preprocess(dataset, level, FALSE, by_dataset = TRUE)
    }
  }
}

# Remove non-NCBI taxa (those without numeric identifiers)
remove_non_ncbi <- function(taxa) {
  for (i in 1:length(taxa)) {
    tmp <- taxa[[i]]
    tmp <- tmp[!is.na(as.numeric(tmp$TaxID)),]
    taxa[[i]] <- tmp
  }
  return(taxa)
}

# Remove the UNCLASSIFIED portion of taxa
remove_unknown <- function(taxa) {
  for (i in 1:length(taxa)) {
    tmp <- taxa[[i]]
    tmp <- tmp[tmp$TaxID != "UNCLASSIFIED",]
    taxa[[i]] <- tmp
  }
  return(taxa)
}

# Preprocess the true simulated profiles
preprocess_simulation_truths <- function(dataset, level) {
  if (file.exists(paste0("analysis/simulation_outputs/", dataset, "/true_profile/", dataset, level, ".tsv"))) {
    return(read.csv(paste0("analysis/simulation_outputs/", dataset, "/true_profile/", dataset, level, ".tsv"), sep='\t'))
  }
  
  in_files = list.files(paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles"))
  in_files = paste0("analysis/simulation_outputs/", dataset, "/true_profile/profiles/", in_files)
  unknown_taxa = read.csv('analysis/databases/ncbi_taxdump/unknown_taxa.tsv', header = F)[,1]
  full_profile <- data.frame("TaxID" = "2")
  for (file in in_files) {
    profile = read.csv(file, sep='\t', skip = 4)
    colnames(profile) <- c("single_id", "rank", "TaxID", "Taxa", (str_split(file, "/")[[1]][6] %>% str_split("\\.(?!.*\\.)"))[[1]][1], "CAMI_genome", "CAMI_OTU")
    profile <- profile[profile$rank != "strain",]
    profile <- profile[c(3, 5)]
    full_profile <- full_join(full_profile, profile, by='TaxID')
  }
  full_profile[is.na(full_profile)] <- 0
  
  # No Eukaryotes
  full_profile <- full_profile[!grepl("^2759",full_profile$TaxID),]
  
  full_profile <- full_profile[str_count(full_profile$TaxID, "\\|") == level - 1,]
  full_profile$TaxID <- gsub(".*\\|","",full_profile$TaxID)
  
  full_profile <- full_profile[!(full_profile$TaxID %in% unknown_taxa),]
  
  full_profile[nrow(full_profile) + 1,] <- c(0, 100 - colSums(full_profile[,-1]))
  full_profile$TaxID[nrow(full_profile)] <- "UNCLASSIFIED"
  
  write.table(full_profile, paste0("analysis/simulation_outputs/", dataset, "/true_profile/", dataset, level, ".tsv"), sep='\t', row.names = F)
  
  return(full_profile)
}

# Proprocess true simulated profiles for all datasets
preprocess_all_truths <- function() {
  for (dataset in list.files("analysis/simulation_outputs")) {
    for (level in 1:7) {
      print(paste0("Processing ", dataset, " at level ", level))
      preprocess_simulation_truths(dataset, level)
    }
  }
}

# Remove the UNCLASSIFIED portion of taxa for a single dataset
remove_unclassified <- function(taxa) {
  return(taxa[taxa[,1]!="UNCLASSIFIED",])
}

# Total sum scale normalize the taxa in units of percent
renormalize <- function(taxa) {
  if (nrow(taxa) > 0) {
    taxa[,-1] <- t(t(taxa[,-1])/colSums(taxa[,-1])) * 100
    taxa[,-1][is.na(taxa[,-1])] <- 0
  }
  return(taxa)
}

# Return a list of the taxa in each sample passing an abundance threshold
list_taxa_by_sample <- function(profile, threshold) {
  output <- list()
  for (i in 2:ncol(profile)) {
    name_vec = profile$TaxID
    abundance_vec = profile[,i]
    name_vec = ifelse(abundance_vec<=threshold, NA, name_vec)
    name_vec = name_vec[!is.na(name_vec) & name_vec != "UNCLASSIFIED"]
    output = append(output, list(name_vec))
  }
  names(output) <- colnames(profile)[2:length(profile)]
  return(output)
}

# Calculate the F1 given two lists of vectors of taxa
calculate_f1 <- function(names_1, names_2) {
  # Set intersection
  names_1 <- names_1[names(names_2)]
  names_2 <- names_2[names(names_1)]
  
  f1_out = vector(length = length(names_1))
  for (i in 1:length(names_1)) {
    prec <- sum(names_1[[i]] %in% names_2[[i]]) / length(names_1[[i]])
    recall <- sum(names_2[[i]] %in% names_1[[i]]) / length(names_2[[i]])
    f1_out[i] <- 2/(1/prec + 1/recall)
  }
  f1_out[is.na(f1_out)] <- 0
  return(f1_out)
}

# Calculate the F1, precision, and recall given two lists of vectors of taxa (names_2 is truth)
calc_precision_recall_f1 <- function(names_1, names_2) {
  names_1 <- names_1[names(names_1) %in% names(names_2)]
  names_2 <- names_2[names(names_2) %in% names(names_1)]
  
  f1 = vector(length = length(names_1))
  prec = vector(length = length(names_1))
  recall = vector(length = length(names_1))
  for (i in 1:length(names_1)) {
    prec[i] <- sum(names_1[[i]] %in% names_2[[i]]) / length(names_1[[i]])
    recall[i] <- sum(names_2[[i]] %in% names_1[[i]]) / length(names_2[[i]])
    f1[i] <- 2/(1/prec[i] + 1/recall[i])
  }
  prec[is.na(prec)] <- 0
  recall[is.na(recall)] <- 0
  f1[is.na(f1)] <- 0
  return(rbind(prec, recall, f1))
}

# Set all taxa below a threshold to 0
threshold_sample <- function(profile, threshold) {
  profile[,-1][profile[,-1] < threshold] <- 0
  return(profile)
}

# Calculate Bray-Curtis dissimilarity between corresponding samples in estimated and true profiles
calc_bc_sim <- function(guess, correct, dist_type = "bray") {
  colnames(correct)[colnames(correct) == "TaxID"] = "TaxIDs"
  correct <- correct[,colnames(guess)]
  
  colnames(guess)[-1] <- paste0("guess_", colnames(guess)[-1])
  colnames(correct)[-1] <- paste0("correct_", colnames(correct)[-1])
  
  joined_df = full_join(guess, correct, by="TaxIDs")
  joined_df[is.na(joined_df)] <- 0
  joined_df = joined_df[,-1]
  joined_df[joined_df < 0] <- 0
  
  dist_mat <- vegdist(t(joined_df), dist_type) %>% as.matrix(labels=TRUE)
  
  output <- vector()
  for (colname in unique(gsub(pattern="^guess_|^correct_", "", colnames(dist_mat)))) {
    same_ids <- colnames(dist_mat)[grepl(paste0(colname,"$"), colnames(dist_mat))]
    output <- append(output, dist_mat[same_ids[1], same_ids[2]])
  }
  
  return(output)
}

# Calculate intersection over union for two lists of vectors of taxa
calc_iou <- function(names_1, names_2) {
  names_1 <- names_1[names(names_1) %in% names(names_2)]
  names_2 <- names_2[names(names_2) %in% names(names_1)]
  
  iou = vector(length = length(names_1))
  for (i in 1:length(names_1)) {
    iou[i] <- sum(names_1[[i]][!is.na(names_1[[i]])] %in% names_2[[i]]) / length(unique(c(names_1[[i]][!is.na(names_1[[i]])], names_2[[i]][!is.na(names_2[[i]])])))
  }
  return(iou)
}

# Calculate Bray-Curtis dissimilarity between corresponding samples in two profiles
calc_bc <- function(profile_1, profile_2, dist_type = "bray") {
  profile_1 <- profile_1[colnames(profile_1) %in% colnames(profile_2)]
  profile_2 <- profile_2[colnames(profile_2) %in% colnames(profile_1)]
  
  colnames(profile_1)[-1] <- paste0("prof1_", colnames(profile_1)[-1])
  colnames(profile_2)[-1] <- paste0("prof2_", colnames(profile_2)[-1])
  
  profile_1$TaxIDs <- as.character(profile_1$TaxIDs)
  profile_2$TaxIDs <- as.character(profile_2$TaxIDs)
  joined_df = full_join(profile_1, profile_2, by="TaxIDs")
  joined_df[is.na(joined_df)] <- 0
  joined_df = joined_df[,-1]
  joined_df[joined_df < 0] <- 0
  
  dist_mat <- vegdist(t(joined_df), dist_type) %>% as.matrix(labels=TRUE)
  
  output <- vector()
  for (colname in unique(gsub(pattern="^prof1_|^prof2_", "", colnames(dist_mat)))) {
    same_ids <- colnames(dist_mat)[grepl(paste0(colname,"$"), colnames(dist_mat))]
    output <- append(output, dist_mat[same_ids[1], same_ids[2]])
  }
  
  return(output)
}

# Set lower triangle of correlation matrix to NA
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Run a mantel test on two profiles and return the statistic and significance
run_bc_mantel <- function(profile_1, profile_2, dist_type = "bray") {
  profile_1 <- profile_1[,c(TRUE, colSums(profile_1[,-1]) != 0)]
  profile_2 <- profile_2[,c(TRUE, colSums(profile_2[,-1]) != 0)]
  
  if (is.null(nrow(profile_1)) | is.null(nrow(profile_2))) {
    return(c(NA, NA))
  }
  if (nrow(profile_1) == 0 | nrow(profile_2) == 0) {
    return(c(NA,NA))
  }
  
  profile_1 <- profile_1[colnames(profile_1) %in% colnames(profile_2)]
  profile_2 <- profile_2[colnames(profile_2) %in% colnames(profile_1)]
  
  if (ncol(profile_1) > 2) {
    dist_1 <- vegdist(t(profile_1[,-1]), dist_type) %>% as.matrix(labels=TRUE)
    dist_2 <- vegdist(t(profile_2[,-1]), dist_type) %>% as.matrix(labels=TRUE)
    
    mantel_out <- mantel(dist_1, dist_2, permutations = 9999)
    
    return(c(mantel_out$statistic, mantel_out$signif))
  } else {
    return(c(NA,NA))
  }
}

# Get the Spearman correlation between alpha diversities on two profiles and return the correlation and significance
run_alpha_cor <- function(profile_1, profile_2, alpha_type = "invsimpson") {
  profile_1 <- profile_1[colnames(profile_1) %in% colnames(profile_2)]
  profile_2 <- profile_2[colnames(profile_2) %in% colnames(profile_1)]
  
  profile_1$TaxIDs <- NULL
  profile_2$TaxIDs <- NULL
  
  if (nrow(profile_1) > 1) {
    if (alpha_type != "count") {
      diversity_1 <- diversity(t(profile_1), index = alpha_type)
    } else {
      diversity_1 <- apply(profile_1, 2, function(x) sum(x!=0))
    }
  } else {
    return(NA) # If only one taxa, diversity will be all 0, so StDev is 0 and correlation is NA
  }
  
  if (nrow(profile_2) > 1) {
    if (alpha_type != "count") {
      diversity_2 <- diversity(t(profile_2), index = alpha_type)
    } else {
      diversity_2 <- apply(profile_2, 2, function(x) sum(x!=0))
    }
  } else {
    return(NA)
  }
  
  sampled_cor <- vector(length=9999) 
  for (i in 1:9999) {
    diversity_1_tmp <- sample(diversity_1, length(diversity_1), replace = F)
    sampled_cor[i] <- cor(diversity_1_tmp, diversity_2, use = "complete.obs", method = "spearman")
  }

  true_cor <- cor(diversity_1, diversity_2, use = "complete.obs", method = "spearman")
  sampled_cor <- c(sampled_cor, true_cor)
  return(c(true_cor, 1 - mean(abs(true_cor) > abs(sampled_cor))))
}

# Nicely capitalize column names
simpleCap <- function(x) {
  if (x == "pH" | x == "ph") {
    return ("pH")
  }
  if (x == "BMI") {
    return ("BMI")
  }
  paste(toupper(substring(x, 1,1)), substring(x, 2),
        sep="", collapse=" ")
}
