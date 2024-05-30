remove(list=ls())
library(reshape2)
library(dplyr)

soil_conditions <- read.csv("analysis/metadata/forest_soil/41396_2018_279_MOESM2_ESM.csv", header = F, skip = 1)
soil_conditions <- soil_conditions[-1,c(3,4,13,14)]
colnames(soil_conditions) <- c("Location", "Sample", "pH_organic", "pH_mineral")
soil_conditions$pH_organic <- as.numeric(gsub(" .*", "", iconv(soil_conditions$pH_organic, from = "ISO-8859-1", to = "UTF-8")))
soil_conditions$pH_mineral <- as.numeric(gsub(" .*", "", iconv(soil_conditions$pH_mineral, from = "ISO-8859-1", to = "UTF-8")))
soil_conditions <- melt(soil_conditions,id.vars=c("Location", "Sample"))
soil_conditions$variable <- ifelse(soil_conditions$variable == "pH_mineral", "Mineral", "Organic")

read_info <- read.csv("analysis/metadata/forest_soil/41396_2018_279_MOESM4_ESM.csv", header = T, skip = 1)[,c(1,4)]

read_identifiers <- read.csv("analysis/metadata/forest_soil/filereport_read_run_PRJEB12502.tsv", sep = "\t")
read_identifiers <- read_identifiers[grepl(";", read_identifiers$fastq_ftp),]
read_identifiers$id = gsub(".*/", "", read_identifiers$submitted_ftp) %>% gsub(pattern="\\..*", replacement="")

read_identifiers <- read_identifiers[,c("run_accession", "id")]

read_identifiers <- left_join(read_identifiers, read_info, by=c("id"="Library.Name"))

read_identifiers$Sample = case_when(grepl("^A|^J|^O|^B|^L|^D", read_identifiers$id) ~ substr(read_identifiers$id, 1,2),
                                    grepl("^TX", read_identifiers$id) ~ substr(read_identifiers$id, 1,3))

merged_metadata = full_join(read_identifiers, soil_conditions, by=c("Sample", "Soil.Layer" = "variable"))
merged_metadata = merged_metadata[!is.na(merged_metadata$Soil.Layer) & !is.na(merged_metadata$run_accession) & !(is.na(merged_metadata$Location) & is.na(merged_metadata$value)),]
merged_metadata = merged_metadata %>% distinct()
merged_metadata <- merged_metadata[,c(1, 3, 5, 6)]

colnames(merged_metadata) <- c("sample_id", "soil_layer", "location", "pH")
write.table(merged_metadata, "analysis/metadata/forest_soil/merged_metadata.tsv", sep="\t", row.names = F)





