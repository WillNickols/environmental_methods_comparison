remove(list = ls())
library(dplyr)

# Read cleaned read numbers
df <- data.frame(matrix(nrow = 0, ncol = 2))
for (file in paste0("analysis/metadata/sample_info/read_nums/", list.files("analysis/metadata/sample_info/read_nums"))) {
  df_append <- read.csv(file, sep = "\t")
  df_append$study <- gsub(".*read_nums\\/", "", file) %>% gsub(pattern="_read_num.tsv", replacement="")
  df <- rbind(df, df_append)
}
df$ID <- gsub("^re", "", df$ID)
df$read_count <- rep(NA, nrow(df))

for (file in paste0("analysis/metadata/sample_info/ena_reports/", list.files("analysis/metadata/sample_info/ena_reports"))) {
  tmp <- read.csv(file, sep="\t")
  if (ncol(tmp) == 4) {
    tmp$library_name <- gsub("_MGX", "", tmp$library_name)
    df <- left_join(df, tmp, by=c("ID" = "library_name")) %>% 
      mutate(read_count = coalesce(read_count.x, read_count.y)) %>% 
      dplyr::select(-read_count.x, -read_count.y)
    df$ID[df$ID == df$library_name] <- df$run_accession
  } else {
    df <- left_join(df, tmp, by=c("ID" = "run_accession")) %>% 
      mutate(read_count = coalesce(read_count.x, read_count.y)) %>% 
      dplyr::select(-read_count.x, -read_count.y)
  }
}

df <- df[,c("ID", "num", "study", "read_count")]
# Paired end
df$read_count[df$study != "gator_soil"] <- 2*df$read_count[df$study != "gator_soil"]

if (sum(df$read_count < df$num) > 0) {
  stop("Cleaned reads should be less than original reads")
}

df$PMID = case_when(df$study == "acid_mine" ~ "32785287",
                     df$study == "animal_gut" ~ "36443458",
                     df$study == "cat_gut" ~ "35467389",
                     df$study == "coastal_sediment" ~ "28793929",
                     df$study == "dog_gut" ~ "32303546",
                     df$study == "forest_soil" ~ "30258172",
                     df$study == "gator_soil" ~ "32424717",
                     df$study == "human" ~ "31142855",
                     df$study == "saltmarsh" ~ "36443458",
                     df$study == "tara_polar" ~ "26029378",)

df$Project_ID = case_when(df$study == "acid_mine" ~ "PRJNA540505",
                    df$study == "animal_gut" ~ "PRJEB42019",
                    df$study == "cat_gut" ~ "PRJNA758898",
                    df$study == "coastal_sediment" ~ "PRJNA322450",
                    df$study == "dog_gut" ~ "PRJEB34360",
                    df$study == "forest_soil" ~ "PRJEB12502",
                    df$study == "gator_soil" ~ "PRJNA554694",
                    df$study == "human" ~ "PRJNA398089",
                    df$study == "saltmarsh" ~ "PRJEB42019",
                    df$study == "tara_polar" ~ "PRJEB9740",)

df$Pairing = case_when(df$study == "gator_soil" ~ "Unpaired",
                       TRUE ~ "Paired")

df$Study = case_when(df$study == "acid_mine" ~ "Acid Mine Drainage",
                     df$study == "animal_gut" ~ "Wild Animal Gut",
                     df$study == "cat_gut" ~ "Cat Gut",
                     df$study == "coastal_sediment" ~ "Coastal Sediment",
                     df$study == "dog_gut" ~ "Dog Gut",
                     df$study == "forest_soil" ~ "Forest Soil",
                     df$study == "gator_soil" ~ "Gator Nest",
                     df$study == "human" ~ "Human Gut",
                     df$study == "saltmarsh" ~ "Salt Marsh",
                     df$study == "tara_polar" ~ "Tara Polar",)

df <- df[,c("ID", "Study", "PMID", "Project_ID", "Pairing", "read_count", "num")]
colnames(df) <- c("ID", "Study", "PMID", "Project ID", "Pairing", "Raw read number", "Cleaned read number")
df <- df[order(df$Study),]

write.table(df, "analysis/metadata/sample_info/read_info.tsv", sep="\t", row.names = F)
