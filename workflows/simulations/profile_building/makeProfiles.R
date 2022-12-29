remove(list=ls())

library(stringr)
library(dplyr)
library(plyr)

path_to_data = "C:/Users/willn/Dropbox (Harvard University)/hutlab/Will/Environmental/Analysis/"

preprocess <- function(dataset, level) {
  if(file.exists(paste0(path_to_data, "cache/", dataset, level, ".rds"))) {
    return (readRDS(paste0(path_to_data, "cache/", dataset, level, ".rds")))
  }
  
  # Read MPA2
  MPA2 <- read.csv(paste0(path_to_data, dataset, "/MPA_2/metaphlan2_taxonomic_profiles.tsv",sep = ""), sep="\t", header = F, skip = 2)
  colnames(MPA2) <- read.csv(paste(path_to_data, dataset, "/MPA_2/metaphlan2_taxonomic_profiles.tsv",sep = ""), sep="\t", header = F, nrows = 1)
  colnames(MPA2) <- gsub("_taxonomic_profile", "", colnames(MPA2))
  colnames(MPA2)[1] <- "Taxa"
  MPA2$Taxa <- gsub(pattern=" ",replacement="_", MPA2$Taxa)
  if ("unclassified" %in% MPA2$Taxa) {
    MPA2$Taxa[MPA2$Taxa=="unclassified"] <- "UNCLASSIFIED"
  } else {
    MPA2[nrow(MPA2) + 1,] <- rep(0, ncol(MPA2))
    MPA2[nrow(MPA2),which(colnames(MPA2)=="Taxa")] <- "UNCLASSIFIED"
  }
  
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/MPA2_renamed.tsv"), sep = "\t")
  MPA2$TaxIDs <- mapvalues(MPA2$Taxa, gsub(pattern=" ",replacement="_", c(map$Taxa, "UNCLASSIFIED")), c(map$renamed, "UNCLASSIFIED"), warn_missing = FALSE)
  
  MPA2 <- MPA2[order(MPA2$Taxa),]
  
  last_tax <- str_split(MPA2$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(MPA2$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  MPA2 <- MPA2[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | MPA2$Taxa=="UNCLASSIFIED") &
                 !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", MPA2$Taxa),]
  
  # Read MPA3
  MPA3 <- read.csv(paste0(path_to_data, dataset, "/MPA_3.0.14_Unknown/metaphlan3_taxonomic_profiles.tsv",sep = ""), sep="\t", header = T, skip = 1)
  MPA3 <- MPA3[,-2]
  colnames(MPA3) <- gsub("_taxonomic_profile", "", colnames(MPA3))
  colnames(MPA3)[1] <- "Taxa"
  MPA3$Taxa <- gsub(pattern=" ",replacement="_", MPA3$Taxa)
  MPA3$Taxa[MPA3$Taxa=="UNKNOWN"] <- "UNCLASSIFIED"
  
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/MPA3_renamed.tsv"), sep = "\t")
  MPA3$TaxIDs <- mapvalues(MPA3$Taxa, gsub(pattern=" ",replacement="_", c(map$Taxa, "UNCLASSIFIED")), c(map$renamed, "UNCLASSIFIED"), warn_missing = FALSE)
  
  MPA3 <- MPA3[order(MPA3$Taxa),]
  
  temp <- MPA3[MPA3$Taxa=="UNCLASSIFIED",]
  MPA3 <- MPA3[!MPA3$Taxa=="UNCLASSIFIED",]
  MPA3[nrow(MPA3) + 1,] <- c(0, colSums(temp[,-c(which(colnames(temp)=="Taxa"), which(colnames(temp)=="TaxIDs"))]), 0)
  MPA3[nrow(MPA3),c(which(colnames(temp)=="Taxa"), which(colnames(temp)=="TaxIDs"))] <- "UNCLASSIFIED"
  
  last_tax <- str_split(MPA3$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(MPA3$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  MPA3 <- MPA3[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | MPA3$Taxa=="UNCLASSIFIED") &
                 !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", MPA3$Taxa),]
  
  # Read MPA4
  MPA4 <- read.csv(paste0(path_to_data, dataset, "/MPA_4/metaphlan4_taxonomic_profiles.tsv",sep = ""), sep="\t", header = T)
  colnames(MPA4) <- gsub("_taxonomic_profile", "", colnames(MPA4))
  colnames(MPA4)[1] <- "Taxa"
  MPA4$Taxa <- gsub(pattern=" ",replacement="_", MPA4$Taxa)
  MPA4$NCBI_tax_id <- NULL
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/MPA4_renamed.tsv"), sep = "\t")
  MPA4$TaxIDs <- mapvalues(MPA4$Taxa, gsub(pattern=" ",replacement="_", c(map$Taxa, "UNCLASSIFIED")), c(map$renamed, "UNCLASSIFIED"), warn_missing = FALSE)
  MPA4 <- MPA4[order(MPA4$Taxa),]
  
  # Deal with MPA4 duplicate unclassified issue
  temp <- MPA4[MPA4$Taxa=="UNCLASSIFIED",]
  MPA4 <- MPA4[!MPA4$Taxa=="UNCLASSIFIED",]
  MPA4[nrow(MPA4) + 1,] <- c(0, colSums(temp[,-c(which(colnames(temp)=="Taxa"), which(colnames(temp)=="TaxIDs"))]), 0)
  MPA4[nrow(MPA4),c(which(colnames(temp)=="Taxa"), which(colnames(temp)=="TaxIDs"))] <- "UNCLASSIFIED"
  
  last_tax <- str_split(MPA4$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(MPA4$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  MPA4 <- MPA4[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | MPA4$Taxa=="UNCLASSIFIED") &
                 !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", MPA4$Taxa),]
  
  # Read metaxa
  metaxa <- read.csv(paste0(path_to_data, dataset, "/metaxa2/merged_reorganized.tsv",sep = ""), sep="\t", header = T)
  colnames(metaxa) <- gsub("_taxonomic_profile", "", colnames(metaxa))
  colnames(metaxa)[1] <- "Taxa"
  metaxa$Taxa <- gsub(pattern=" ",replacement="_", metaxa$Taxa)
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/metaxa2_renamed.tsv"), sep = "\t")
  metaxa$TaxIDs <- mapvalues(metaxa$Taxa, gsub(pattern=" ",replacement="_", map$Taxa),
                             map$renamed %>% trimws("right", whitespace = "\\|"), warn_missing = FALSE)
  metaxa$Taxa <- mapvalues(metaxa$Taxa, gsub(pattern=" ",replacement="_", map$Taxa),
                           map$Taxa_correct %>% trimws("right", whitespace = "\\|"), warn_missing = FALSE)
  
  metaxa <- metaxa[order(metaxa$Taxa),]
  
  # Fix unknown estimate
  tmp <- metaxa[grepl("k__Unclassified", metaxa$Taxa),][1,][-c(1, length(metaxa))] + metaxa[grepl("Unknown", metaxa$Taxa),][1,][-c(1, length(metaxa))]
  metaxa <- metaxa[!grepl("k__Unclassified", metaxa$Taxa) & !grepl("Unknown", metaxa$Taxa),]
  metaxa[nrow(metaxa) + 1,] <- c(0, tmp, 0)
  metaxa[nrow(metaxa),c(which(colnames(metaxa)=="Taxa"), which(colnames(metaxa)=="TaxIDs"))] <- "UNCLASSIFIED"
  
  # Remove unknown taxa
  last_tax <- str_split(metaxa$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(metaxa$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  metaxa <- metaxa[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | metaxa$Taxa=="UNCLASSIFIED") &
                     !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", metaxa$Taxa),]
  
  # Read mOTUs3  
  motus <- read.csv(paste0(path_to_data, dataset, "/mOTUs3/mOTUs3_merged.tsv",sep = ""), sep="\t", header = T)
  colnames(motus)[1] <- "Taxa"
  motus$Taxa <- gsub(pattern=" ",replacement="_", motus$Taxa) %>% 
    gsub(pattern="_\\[.*\\]",replacement="")
  
  motus[,-1] <- lapply(motus[,-1], function (column) as.numeric(column) * 100)
  
  # Insert unclassified
  motus[nrow(motus) + 1,] <- c(0, 
                               100 - colSums(motus[motus$Taxa=="k__Bacteria"|
                                                     motus$Taxa=="k__Archaea"|
                                                     motus$Taxa=="k__Eukaryota",-1]))
  motus[nrow(motus),which(colnames(motus)=="Taxa")] <- "UNCLASSIFIED"
  
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/mOTUs3_renamed.tsv"), sep = "\t")
  motus$TaxIDs <- mapvalues(motus$Taxa, gsub(pattern=" ",replacement="_", c(map$Taxa, "UNCLASSIFIED")), c(map$renamed, "UNCLASSIFIED"), warn_missing = FALSE)
  
  motus <- motus[order(motus$Taxa),]
  
  last_tax <- str_split(motus$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(motus$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  motus <- motus[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | motus$Taxa=="UNCLASSIFIED") &
                   !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", motus$Taxa),]
  
  # Read phylophlan
  phylophlan <- read.csv(paste0(path_to_data, dataset, "/phylophlan/final_profile_by_sample.tsv",sep = ""), sep="\t", header = T)
  colnames(phylophlan)[1] <- "Taxa"
  phylophlan$Taxa <- gsub(pattern=" ",replacement="_", phylophlan$Taxa)
  phylophlan <- phylophlan[order(phylophlan$Taxa),]
  phylophlan$Taxa <- c(phylophlan$Taxa[-nrow(phylophlan)], "UNCLASSIFIED")
  
  # Remove weird duplicate columns
  phylophlan <- phylophlan[,!grepl("\\.1$", colnames(phylophlan))]
  
  # Make SGBs unknown
  phylophlan[phylophlan$Taxa=="UNCLASSIFIED", -1] <- phylophlan[phylophlan$Taxa=="UNCLASSIFIED", -1] + 
    colSums(phylophlan[grepl("sgb_", phylophlan$Taxa),-1])
  phylophlan <- phylophlan[!grepl("sgb_", phylophlan$Taxa),]
  
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/phylophlan_renamed.tsv"), sep = "\t")
  phylophlan$TaxIDs <- mapvalues(phylophlan$Taxa, gsub(pattern=" ",replacement="_", c(map$Taxa, "UNCLASSIFIED")), c(map$renamed, "UNCLASSIFIED"), warn_missing = FALSE)
  
  last_tax <- str_split(phylophlan$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(phylophlan$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  phylophlan <- phylophlan[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | phylophlan$Taxa=="UNCLASSIFIED") &
                             !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", phylophlan$Taxa),]
  
  # Read Kraken
  kraken <- read.csv(paste0(path_to_data, dataset, "/kraken/kraken_merged.tsv",sep = ""), sep="\t", header = T)
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/kraken_renamed.tsv"), sep = "\t")
  kraken$TaxIDs <- mapvalues(kraken$Taxa, gsub(pattern=" ",replacement="_", c(map$Taxa, "UNCLASSIFIED") %>% tolower()), c(map$renamed, "UNCLASSIFIED"), warn_missing = FALSE)
  kraken <- kraken[order(kraken$Taxa),]
  
  last_tax <- str_split(kraken$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(kraken$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  kraken <- kraken[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | kraken$Taxa=="UNCLASSIFIED") &
                     !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", kraken$Taxa),]
  
  # GTDB-Tk
  gtdbtk <- read.csv(paste0(path_to_data, dataset, "/gtdbtk/final_profile_by_sample.tsv"), sep="\t", header = T)
  gtdbtk$Taxa <- gsub(pattern=" ",replacement="_", gtdbtk$Taxonomy) %>% tolower()
  gtdbtk$Taxonomy <- NULL
  gtdbtk$Taxa[gtdbtk$Taxa == "unknown"] <- "UNCLASSIFIED"
  
  # Remove weird duplicate columns
  gtdbtk <- gtdbtk[,!grepl("\\.1$", colnames(gtdbtk))]
  
  gtdbtk[gtdbtk$Taxa=="UNCLASSIFIED", -ncol(gtdbtk)] <- gtdbtk[gtdbtk$Taxa=="UNCLASSIFIED", -ncol(gtdbtk)] + 
    colSums(gtdbtk[grepl("sgb_", gtdbtk$Taxa),-ncol(gtdbtk)])
  gtdbtk <- gtdbtk[!grepl("sgb_", gtdbtk$Taxa),]
  
  map <- read.csv(paste0(path_to_data, "Databases/DB_with_IDs_Propagated/gtdbtk_renamed.tsv"), sep = "\t")
  gtdbtk$TaxIDs <- mapvalues(gtdbtk$Taxa, gsub(pattern=" ",replacement="_", c(map$Taxa, "UNCLASSIFIED") %>% tolower()), c(map$renamed, "UNCLASSIFIED"), warn_missing = FALSE)
  
  gtdbtk <- gtdbtk[order(gtdbtk$Taxa),]
  last_tax <- str_split(gtdbtk$Taxa, "\\|") %>% mapply(FUN = "[[", as.list(str_count(gtdbtk$Taxa, "\\|") + 1)) %>% 
    gsub(pattern = "\\w__", replacement = "") %>% tolower()
  gtdbtk <- gtdbtk[!is.na(last_tax) & last_tax != "na" & last_tax != "" & (!grepl("unclassified", last_tax) | gtdbtk$Taxa=="UNCLASSIFIED") &
                     !grepl("noname", last_tax) & !grepl("/", last_tax) & !grepl("t__", gtdbtk$Taxa),]
  
  
  databases <- list(subset(MPA2, select=c("Taxa",sort(colnames(MPA2))[sort(colnames(MPA2))!="Taxa"])), 
                    subset(MPA3, select=c("Taxa",sort(colnames(MPA3))[sort(colnames(MPA3))!="Taxa"])), 
                    subset(MPA4, select=c("Taxa",sort(colnames(MPA4))[sort(colnames(MPA4))!="Taxa"])), 
                    subset(motus, select=c("Taxa",sort(colnames(motus))[sort(colnames(motus))!="Taxa"])), 
                    subset(phylophlan, select=c("Taxa",sort(colnames(phylophlan))[sort(colnames(phylophlan))!="Taxa"])),
                    subset(kraken, select=c("Taxa",sort(colnames(kraken))[sort(colnames(kraken))!="Taxa"])),
                    subset(gtdbtk, select=c("Taxa",sort(colnames(gtdbtk))[sort(colnames(gtdbtk))!="Taxa"])),
                    subset(metaxa, select=c("Taxa",sort(colnames(metaxa))[sort(colnames(metaxa))!="Taxa"])))
  
  if (length(unique(lapply(databases, colnames))) > 1) {
    stop("Different column names between tools")
  }
  
  replace_col_names <- function(x) {
    colnames(x)[!grepl("Taxa|TaxIDs", colnames(x))] <- paste0(dataset, "_", 1:length(colnames(x)[!grepl("Taxa|TaxIDs", colnames(x))]))
    return (x)
  }
  
  databases <- lapply(databases, replace_col_names)
  
  if (level == 7) {
    species <- vector("list", length(databases))
    for (i in 1:length(databases)) {
      str_counts <- str_split(databases[[i]]$Taxa[grepl("s__", databases[[i]]$Taxa)], "s__", 2) %>% sapply("[[", 1) %>%
        str_count("\\|")
      str_ends <- str_split(databases[[i]]$Taxa[grepl("s__", databases[[i]]$Taxa)], "s__", 2) %>% sapply("[[", 2) %>%
        str_count("\\|")
      tmp <- (str_split(databases[[i]]$TaxIDs[grepl("s__", databases[[i]]$Taxa)], "\\|", n = str_counts + 1) %>% 
                mapply(FUN = "[[", as.list(str_counts + 1)) %>%
                str_split("\\|") %>% sapply("[[", 1))
      species[[i]] <- databases[[i]][grepl("s__", databases[[i]]$Taxa),][intersect(seq_along(tmp)[!duplicated(tmp)], which(str_ends==0)),]
      species[[i]]$Taxa <- gsub(".*\\|", "", species[[i]]$Taxa)
      species[[i]]$TaxIDs <- gsub(".*\\|", "", species[[i]]$TaxIDs)
      species[[i]] <- species[[i]][!is.na(species[[i]]$TaxIDs), ]
      species[[i]] <- species[[i]][!grepl("\\|na", species[[i]]$TaxIDs %>% tolower()), ]
      species[[i]] <- species[[i]][!grepl("unclassified", species[[i]]$TaxIDs %>% tolower()), ]
      species[[i]] <- species[[i]][!grepl("noname", species[[i]]$TaxIDs %>% tolower()), ]
      species[[i]] <- species[[i]][!grepl("incertae", species[[i]]$TaxIDs %>% tolower()), ]
      species[[i]] <- species[[i]][!grepl("/", species[[i]]$TaxIDs %>% tolower()), ]
      
      species[[i]][nrow(species[[i]]) + 1,] <- c(0, 100 - colSums(species[[i]][,-c(which(colnames(species[[i]])=="Taxa"), which(colnames(species[[i]])=="TaxIDs"))]), 0)
      species[[i]][nrow(species[[i]]),c(which(colnames(species[[i]])=="Taxa"), which(colnames(species[[i]])=="TaxIDs"))] <- "UNCLASSIFIED"
      species[[i]]$Taxa <- NULL
      species[[i]] <- species[[i]][, c("TaxIDs", colnames(species[[i]])[colnames(species[[i]])!="TaxIDs"] )]
    }
    saveRDS(species, paste0(path_to_data, "cache/", dataset, level, ".rds"))
    return (species)
  }
}

renormalize <- function(column, coef) {
  if (is.infinite(coef)) {coef = 0}
  return (as.numeric(column) * coef)
}

make_soil_profile <- function() {
  data <- preprocess("Trees", 7)
  
  # Normalize and remove columns with 0 abundance
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    tmp <- tmp[!is.na(as.numeric(tmp$TaxIDs)),]
    tmp[nrow(tmp) + 1,] <- c(0, 100 - colSums(tmp[,-which(colnames(tmp)=="TaxIDs")]))
    tmp[nrow(tmp),which(colnames(tmp)=="TaxIDs")] <- "UNCLASSIFIED"
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    tmp[,-1] <- tmp[,-1][,colSums(tmp[,-1], na.rm = T)==100]
    data[[i]] <- tmp
  }
  
  dfs <- list()
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    tmp$Taxa <- NULL
    dfs[[i]] <- data.frame(cbind(means=rowMeans(tmp[,colnames(tmp)!= "TaxIDs"]),
                                 TaxIDs=tmp$TaxIDs))
    colnames(dfs[[i]]) <- c(paste0("means.", i), "TaxIDs")
  }
  
  combined <- dfs[[1]]
  for (i in 2:length(data)) {
    combined <- merge(combined, dfs[[i]], by=("TaxIDs"), all=T)
  }
  
  combined[is.na(combined)] <- 0
  combined[,colnames(combined)!= "TaxIDs"] <- lapply(combined[,colnames(combined)!= "TaxIDs"], as.numeric)
  
  out_df <- data.frame(cbind(TaxIDs=combined$TaxIDs, means=rowMeans(combined[,colnames(combined)!= "TaxIDs"])))
  out_df$means <- as.numeric(out_df$means)
  out_df <- out_df[order(out_df$means, decreasing = T),]
  
  out_df <- out_df[!is.na(as.numeric(out_df$TaxIDs)),]
  out_df <- out_df[out_df$TaxIDs != "9606",]
  
  out_df <- out_df[1:1000,1]
  
  write.table(out_df, paste0("profiles/soil/known_species_list.tsv"), sep="\t", row.names = F, col.names = F)
}

make_ocean_profile <- function() {
  data <- preprocess("TaraPolar", 7)
  
  # Normalize and remove columns with 0 abundance
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    tmp <- tmp[!is.na(as.numeric(tmp$TaxIDs)),]
    tmp[nrow(tmp) + 1,] <- c(0, 100 - colSums(tmp[,-which(colnames(tmp)=="TaxIDs")]))
    tmp[nrow(tmp),which(colnames(tmp)=="TaxIDs")] <- "UNCLASSIFIED"
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    tmp[,-1] <- tmp[,-1][,colSums(tmp[,-1], na.rm = T)==100]
    data[[i]] <- tmp
  }
  
  dfs <- list()
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    tmp$Taxa <- NULL
    dfs[[i]] <- data.frame(cbind(means=rowMeans(tmp[,colnames(tmp)!= "TaxIDs"]),
                                 TaxIDs=tmp$TaxIDs))
    colnames(dfs[[i]]) <- c(paste0("means.", i), "TaxIDs")
  }
  
  combined <- dfs[[1]]
  for (i in 2:length(data)) {
    combined <- merge(combined, dfs[[i]], by=("TaxIDs"), all=T)
  }
  
  combined[is.na(combined)] <- 0
  combined[,colnames(combined)!= "TaxIDs"] <- lapply(combined[,colnames(combined)!= "TaxIDs"], as.numeric)
  
  out_df <- data.frame(cbind(TaxIDs=combined$TaxIDs, means=rowMeans(combined[,colnames(combined)!= "TaxIDs"])))
  out_df$means <- as.numeric(out_df$means)
  out_df <- out_df[order(out_df$means, decreasing = T),]
  
  out_df <- out_df[!is.na(as.numeric(out_df$TaxIDs)),]
  out_df <- out_df[out_df$TaxIDs != "9606",]
  
  out_df <- out_df[1:1000,1]
  
  write.table(out_df, paste0("profiles/ocean/known_species_list.tsv"), sep="\t", row.names = F, col.names = F)
}

make_animal_gut_profile <- function() {
  data1 <- preprocess("WildGut", 7)
  data2 <- preprocess("Cats", 7)
  data3 <- preprocess("Dogs", 7)
  
  data = list()
  for (i in 1:length(data1)) {
    data[[i]] = merge(data1[[i]], data2[[i]], by="TaxIDs") %>% merge(data3[[i]], by="TaxIDs")
  }
  
  # Normalize and remove columns with 0 abundance
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    tmp <- tmp[!is.na(as.numeric(tmp$TaxIDs)),]
    tmp[nrow(tmp) + 1,] <- c(0, 100 - colSums(tmp[,-which(colnames(tmp)=="TaxIDs")]))
    tmp[nrow(tmp),which(colnames(tmp)=="TaxIDs")] <- "UNCLASSIFIED"
    classified <- (100 - as.numeric(tmp[tmp$TaxIDs=="UNCLASSIFIED",-1])) / 100
    tmp <- tmp[tmp$TaxIDs!="UNCLASSIFIED",]
    tmp[,-1] <- mapply(renormalize, tmp[,-1], 1/classified)
    tmp[,-1] <- tmp[,-1][,colSums(tmp[,-1], na.rm = T)==100]
    data[[i]] <- tmp
  }
  
  dfs <- list()
  for (i in 1:length(data)) {
    tmp <- data[[i]]
    tmp$Taxa <- NULL
    dfs[[i]] <- data.frame(cbind(means=rowMeans(tmp[,colnames(tmp)!= "TaxIDs"]),
                                 TaxIDs=tmp$TaxIDs))
    colnames(dfs[[i]]) <- c(paste0("means.", i), "TaxIDs")
  }
  
  combined <- dfs[[1]]
  for (i in 2:length(data)) {
    combined <- merge(combined, dfs[[i]], by=("TaxIDs"), all=T)
  }
  
  combined[is.na(combined)] <- 0
  combined[,colnames(combined)!= "TaxIDs"] <- lapply(combined[,colnames(combined)!= "TaxIDs"], as.numeric)
  
  out_df <- data.frame(cbind(TaxIDs=combined$TaxIDs, means=rowMeans(combined[,colnames(combined)!= "TaxIDs"])))
  out_df$means <- as.numeric(out_df$means)
  out_df <- out_df[order(out_df$means, decreasing = T),]
  
  out_df <- out_df[!is.na(as.numeric(out_df$TaxIDs)),]
  out_df <- out_df[out_df$TaxIDs != "9606",]
  
  out_df <- out_df[1:1000,1]
  
  write.table(out_df, paste0("profiles/animal_gut/known_species_list.tsv"), sep="\t", row.names = F, col.names = F)
}

