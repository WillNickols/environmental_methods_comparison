#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   phylophlan_tax_wrangling.R [-i <input> --tax <tax> --qa <combined> -o <output>]

Options:
   -i checkm abundance directory
   --tax phylophlan tax
   --qa checkm qa combined
   -o output

' -> doc 

opts <- docopt(doc)

# Read folder of abundances
in_dir <- gsub("/$", "", opts$i)
files <- list.files(path=in_dir, full.names=TRUE)

read_nums <- files[grep("read_num",files)]
taxonomies <- files[grep("taxonomy",files)]

# Combine abundance files
checkm <- c()
for (taxonomy in taxonomies) {
  data <- fread(taxonomy)[,c(1,6)]
  colnames(data) <- c("ID", "Abundance")
  checkm <- rbind(checkm, data)
}

# Merge in bin taxonomies
map <- fread(gsub("/$", "", opts$tax), skip = 3)
map <- map[,-5]
colnames(map) <- c("ID", "Species", "Genus", "Family")
map$Taxonomy <- case_when(sapply(str_split(map$Species, ":"), `[`, 4) %>% as.numeric <= 0.05 ~ sapply(str_split(map$Species, ":"), `[`, 3),
                          sapply(str_split(map$Genus, ":"), `[`, 4) %>% as.numeric <= 0.15 ~ sapply(str_split(map$Genus, ":"), `[`, 3),
                          sapply(str_split(map$Family, ":"), `[`, 4) %>% as.numeric <= 0.3 ~ sapply(str_split(map$Family, ":"), `[`, 3),
                          TRUE ~ NA_character_)
map$Taxonomy <- gsub("(.*):(.*):(.*):(.*)", "\\3", map$Taxonomy)
map <- data.frame(ID=map$ID, Taxonomy=map$Taxonomy)
abundance <- left_join(checkm, map, by = "ID")

# Keep only medium and high quality
qa <- fread(gsub("/$", "", opts$qa))
qa <- qa[,c(1,2)]
abundance <- left_join(abundance, qa, by=c("ID"="bin_id"))

# Reformat
abundance$Taxonomy <- gsub(pattern=" \\(root\\)", replacement="", abundance$Taxonomy)
abundance <- mutate(abundance, Taxonomy = case_when(
  grepl("unbinned", ID) ~ "UNKNOWN",
  is.na(Taxonomy) ~ "UNKNOWN",
  grepl("low_quality", quality) ~ "UNKNOWN",
  TRUE ~ Taxonomy
))

abundance$ID <- gsub("\\..*", "", abundance$ID)

samples <- unique(abundance$ID)
total_profile <- data.frame(as.character(matrix(nrow = 0, ncol = 1)))
colnames(total_profile) <- c("Taxonomy")
for (sample in samples) {
  # Keep only the current abundances
  profile <- abundance[abundance$ID==sample,c(2,3)]
  
  # Combine same taxonomic groups
  profile <- aggregate(.~Taxonomy,data=profile,FUN=sum)
  
  # Upcycle taxonomic groups 
  for (name in profile$Taxonomy) {
    abun <- profile$Abundance[profile$Taxonomy==name]
    newTax <- name
    while (grepl("\\|", newTax)) {
      newTax <- gsub("(\\|[^\\|]*)$", "", newTax)
      profile[nrow(profile) + 1,] <- c(newTax, abun)
    }
  }
  profile$Abundance <- as.numeric(profile$Abundance)
  
  # Combine same taxonomic groups again
  profile <- aggregate(.~Taxonomy,data=profile,FUN=sum)
  
  # Join to the full table
  colnames(profile) <- c("Taxonomy", sample)
  total_profile <- full_join(total_profile, profile)
}

# Alphabetize
total_profile <- total_profile[order(total_profile$Taxonomy),]

# Replace NA with 0
total_profile[is.na(total_profile)] <- 0

# Read in unmapped reads
mapped_props <- c()
for (read_num in read_nums) {
  text <- readChar(read_num, file.info(read_num)$size)
  temp <- strsplit(text, "\n")[[1]] %>% as.numeric()
  mapped <- temp[1]/temp[2]
  sample <- gsub(".*checkm\\/", "", read_num) %>% gsub(pattern="\\..*", replacement="")
  mapped_props <- rbind(mapped_props, c(sample, mapped))
}

colnames(mapped_props) <- c("ID", "prop")
mapped_props <- as.data.frame(mapped_props)
mapped_props$prop <- as.numeric(levels(mapped_props$prop))[mapped_props$prop]

for (ID in mapped_props$ID) {
  if (mapped_props$prop[mapped_props$ID==ID] == 0 ) {
    total_profile[,'new'] <- rep(0, nrow(total_profile))
    colnames(total_profile)[colnames(total_profile)=="new"] <- ID
    }
  total_profile[,ID] <- total_profile[,ID] * mapped_props$prop[mapped_props$ID==ID]
  total_profile[total_profile$Taxonomy=="UNKNOWN", ID] <- total_profile[total_profile$Taxonomy=="UNKNOWN", ID] + 
    100 * (1-mapped_props$prop[mapped_props$ID==ID])
}

write.table(total_profile, opts$o, quote=FALSE, sep="\t", row.names=FALSE)

#
