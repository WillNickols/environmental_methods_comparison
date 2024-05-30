#! usr/bin/env Rscript

rm(list = ls())
library(docopt)
library(dplyr)
library(data.table)
library(stringr)

'Usage:
   parse_and_merge.R [-i <input> --abundance_dir <abundance_by_> --mash <mash> -o <output> --profiler <profiler>]

Options:
   -i simulations directory
   --abundance_dir abundance_by_
   --mash mash/ 
   -o output.tsv
   --profiler <profiler>

' -> doc 

opts <- docopt(doc)

# Read folder of abundances
in_dir <- gsub("/$", "", opts$abundance_dir)
files <- list.files(path=in_dir, full.names=TRUE)

read_nums <- files[grep("read_num",files)]

abundance_files <- files[grepl("\\.abundance\\.",files) & !grepl("done", files)]

# Combine abundance files
abundances <- fread(abundance_files[1])[,c(1,6)]
colnames(abundances) <- c("ID", (str_split(abundance_files[1], "/")[[1]] %>% tail(1) %>%
                             str_split(pattern = "\\.abundance"))[[1]][1])

for (abundance_file in abundance_files[-1]) {
  abundances$ID = as.character(abundances$ID)
  data <- fread(abundance_file)[,c(1,6)]
  colnames(data) <- c("ID", (str_split(abundance_file, "/")[[1]] %>% tail(1) %>%
                        str_split(pattern = "\\.abundance"))[[1]][1])
  if (nrow(data) > 0) {
    abundances <- full_join(abundances, data, by="ID")
  } else {
    abundances[,(str_split(abundance_file, "/")[[1]] %>% tail(1) %>%
                  str_split(pattern = "\\.abundance"))[[1]][1]] <- NA
  }
}

abundances <- data.frame("ID" = abundances$ID, "abundance" = coalesce(!!!abundances[,-1]))

# Merge in bin taxonomies
if (opts$profiler == "phylophlan") {
  map <- fread(gsub("/$", "", paste0(opts$i, "phylophlan/phylophlan_relab.tsv")))
  map <- map[,c(2,7)]
  colnames(map) <- c("ID", "Taxonomy")
} else {
  map <- fread(gsub("/$", "", paste0(gsub("simulations\\/assembly.*", "", opts$i), "simulations/gtdbtk/relab/gtdbtk_relab.tsv")))
  colnames(map) <- c("ID", "Taxonomy")
}

abundance <- left_join(abundances, map, by = "ID")

# Keep only medium and high quality
qa <- fread(gsub("/$", "", paste0(opts$i, "checkm/qa/checkm_qa_and_n50.tsv")))
qa <- qa[,c(1,5)]

abundance <- left_join(abundance, qa, by=c("ID"="bin_id"))
abundance <- abundance[abundance$keep == "keep" & !is.na(abundance$keep),]

out_df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(out_df) <- c("sample", "SGB", "bin_name", "real_tax", "reconstructed_tax", "abundance_real", "abundance_reconstructed")

sim_dir_tmp <- str_split(opts$i, "/")[[1]]
sim_dir <- paste0(paste0(sim_dir_tmp[1:(length(sim_dir_tmp) - 4)], collapse = "/"), "/")
profile <- read.csv(paste0(sim_dir, "profile.tsv"), sep="\t")
profile$filepath <- gsub(".*\\/", "", profile$filepath)

for (file in list.files(path = opts$mash)[!grepl("\\.txt|\\.msh", list.files(path = opts$mash))]) {
  mash_dist = read.csv(paste0(opts$mash, file, "/mash_dist_out.tsv"), sep="\t")
  rec_names <- gsub(".*\\.bins\\.", "", colnames(mash_dist)[-1]) %>% gsub(pattern="\\.fa", replacement = "")
  names(rec_names) <- gsub(".*\\/", "", mash_dist[,1][apply(as.matrix(mash_dist[,-1]), 2, which.min)]) %>% gsub(pattern="\\.fa", replacement = "")  %>% gsub(pattern="_genomic$", replacement = "")
  
  tmp_df <- data.frame(matrix(nrow = length(mash_dist[,1]), ncol = 7))
  colnames(tmp_df) <- c("sample", "SGB", "bin_name", "real_tax", "reconstructed_tax", "abundance_real", "abundance_reconstructed")
  tmp_df$sample <- file
  tmp_df$SGB <- gsub(".*\\/", "", mash_dist[,1]) %>% gsub(pattern="\\.fa", replacement = "") %>% gsub(pattern="_genomic$", replacement = "")
  tmp_df$bin_name <- rec_names[tmp_df$SGB]
  
  cur_profile <- read.csv(paste0(sim_dir, "simulations/reorganized/tax/", file, ".txt"), sep="\t", skip = 3)
  cur_profile <- cur_profile[cur_profile$RANK == "strain",]
  cur_profile <- left_join(cur_profile, profile, by=c("X_CAMI_genomeID"="genome_ID"))
  cut_first <- function(x) {
    tmp <- str_split(x, "\\.")[[1]]
    return(paste0(tmp[-1], collapse = "."))
  }
  tmp_cami_id <- sapply(cur_profile$X_CAMI_genomeID, cut_first)
  cur_profile$SGB <- ifelse(grepl("BIN_UNKNOWN_", cur_profile$X_CAMI_genomeID), gsub("\\.fa", "", cur_profile$filepath), tmp_cami_id)
  tmp_df$real_tax <- cur_profile$TAXPATH[match(tmp_df$SGB, cur_profile$SGB)]
  tmp_df$real_tax <- gsub("\\|\\|.*", "", tmp_df$real_tax)
  tmp_df$abundance_real <- cur_profile$PERCENTAGE[match(tmp_df$SGB, cur_profile$SGB)]
  
  tmp_df$reconstructed_tax <- abundance$Taxonomy[match(rec_names[tmp_df$SGB], abundance$ID)]
  tmp_df$abundance_reconstructed <- abundance$abundance[match(rec_names[tmp_df$SGB], abundance$ID)]
  tmp_df <- tmp_df[!is.na(tmp_df$real_tax),]
  
  out_df <- rbind(out_df,tmp_df)
}

write.table(out_df, opts$o, quote=FALSE, sep="\t", row.names=FALSE)
