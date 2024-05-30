remove(list=ls())
library(reshape2)
library(dplyr)

hmp2 <- read.csv("analysis/metadata/human/hmp2_metadata.csv", check.names = F)
hmp2 <- hmp2[,c("External ID", "Participant ID", "data_type", "consent_age", "BMI", "sex")]
for (subject in unique(hmp2$`Participant ID`)) {
  hmp2$BMI[hmp2$`Participant ID` == subject] = mean(hmp2$BMI[hmp2$`Participant ID` == subject], na.rm = T)
}
hmp2 <- hmp2[hmp2$data_type == "metagenomics",]
accessions <- read.csv("analysis/metadata/sample_info/read_info.tsv", sep = "\t")
hmp2 <- hmp2[hmp2$`External ID` %in% accessions$ID,]
hmp2$data_type <- NULL
colnames(hmp2) <- c("sample_id", "participant", "age", "bmi", "sex")

write.table(hmp2, "analysis/metadata/human/merged_metadata.tsv", sep="\t", row.names = F)
