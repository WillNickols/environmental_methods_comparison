#!/usr/bin/env Rscript

packages <- c("tidyverse", "data.table", "readxl", "docopt", "rjson")

install.packages(setdiff(packages, rownames(installed.packages())), repos="http://cran.us.r-project.org") 

# arguments

library(docopt)

'Usage:
   get_ena_file_reports.R [-a <accessions>]

Options:
   -a input table containing project accessions

' -> doc 

opts <- docopt(doc)

# load libraries 

library(tidyverse)
library(data.table)
library(readxl)
library(rjson)

#########################
# get accession numbers #
#########################

data <- opts$a

ena <- fread(data) %>%
  select(accession) %>%
  mutate(accession = strsplit(as.character(accession), ",")) %>% 
  unnest(accession) %>%
  mutate(accession = gsub(" ", "", gsub(".*=", "", accession))) %>%
  mutate(link = paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", accession, "&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_platform,library_name,library_strategy,fastq_ftp,submitted_ftp,sra_ftp,sample_alias,sample_title&format=tsv&download=true"))

#############################
# download the file reports #
#############################

file_report_dir <- "./output/file_reports/"
dir.create(file_report_dir, recursive=TRUE)

accessions <- levels(as.factor(ena$accession))

for (i in accessions) {
  
  accession_url <- as.data.frame(ena %>% filter(accession == i))[1,2]
  
  download.file(accession_url, paste0(file_report_dir, "/", i, "_file_report.tsv"))
  
}


############################
# combine the file reports #
############################

reports <- list.files(path=file_report_dir, pattern="_file_report.tsv", full.names=TRUE)

combined_reports <- list()

for (i in reports) {
  
  report <- fread(i)
  
  combined_reports <- rbind(combined_reports, report)
  
}

###############################################
# extract the ftp links from the file reports #
###############################################

out_dir <- "./output/reads/"
dir.create(out_dir, recursive=TRUE)

ftp <- combined_reports %>%
  mutate(fastq_ftp = strsplit(as.character(fastq_ftp), ";")) %>%
  unnest(fastq_ftp) %>%
  select(fastq_ftp) %>%
  mutate(reads = case_when(
    grepl("_.\\.fastq", fastq_ftp) ~ "paired",
    !grepl("_.\\.fastq", fastq_ftp) ~ "unpaired"
  )) %>%
  mutate(fastq_ftp = paste0("https://", fastq_ftp))

ftp_paired <- ftp %>%
  filter(reads == "paired") %>%
  select(fastq_ftp)

if (nrow(ftp_paired) > 0) {
  out_dir_p <- paste0(out_dir, "/paired")
  dir.create(out_dir_p, recursive=TRUE)
  write.table(ftp_paired, paste0(out_dir_p, "/ena_fastq_ftp_paired.txt"), sep="\t", quote=F, col.names=F, row.names=F)
}

ftp_unpaired <- ftp %>%
  filter(reads == "unpaired") %>%
  select(fastq_ftp)

if (nrow(ftp_unpaired) > 0) {
  out_dir_u <- paste0(out_dir, "/unpaired")
  dir.create(out_dir_u, recursive=TRUE)
  write.table(ftp_unpaired, paste0(out_dir_u, "/ena_fastq_ftp_unpaired.txt"), sep="\t", quote=F, col.names=F, row.names=F) 
}

############################
# download sample metadata #
############################

metadata_dir <- "./output/metadata/"
dir.create(metadata_dir, recursive=TRUE)

sample_ids <- levels(as.factor(combined_reports$sample_accession))

sample_to_study <- combined_reports %>%
  select(study_accession, sample_accession) %>%
  distinct()

# these variables will be removed from the metadata

variables_to_remove <- c(
  "domain",
  "geographic location (latitude)",
  "geographic location (longitude)",
  "lat lon",
  "human gut environmental package",
  "relationships.source",
  "relationships.type",
  "relationships.target",
  "externalReferences.url",
  "releaseDate",
  "updateDate",
  "submittedVia",
  "create",
  "_links.self.href",
  "_links.curationDomain.href",
  "_links.curationDomain.templated",
  "_links.curationLinks.href",
  "_links.curationLink.href",
  "_links.curationLink.templated"
)

combined_metadata <- c()

for (i in sample_ids) {
  
  json_file <- fromJSON(file = paste0("https://www.ebi.ac.uk/biosamples/samples/", i, ".json")) %>%
    unlist()
  
  sample_metadata <- json_file %>% 
    as.data.frame() %>%
    mutate(metadata = names(json_file)) %>%
    rename(value=1, metadata=2) %>%
    select(2, 1) %>%
    filter(!(metadata %in% variables_to_remove)) %>%
    mutate(metadata = gsub("characteristics\\.|host |\\.ontologyTerms$|\\.text$|\\.unit$", "", metadata)) %>%
    filter(!grepl("\\.tag", metadata)) %>%
    distinct() %>%
    rename(value = 2) %>%
    group_by(metadata) %>%
    summarise(value = paste0(value, collapse=" ")) %>%
    ungroup() %>%
    mutate(sample_accession = i)
  
  combined_metadata <- rbind(combined_metadata, sample_metadata)
  
}

combo_wide <- combined_metadata %>%
  merge(., sample_to_study, by="sample_accession") %>%
  distinct() %>%
  pivot_wider(names_from=metadata, values_from=value)

write.table(combo_wide, "biosamples_metadata.tsv", quote=FALSE, sep="\t", row.names=FALSE)

# end