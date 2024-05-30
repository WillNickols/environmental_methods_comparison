remove(list=ls())
library(reshape2)
library(dplyr)

capitalize <- function(x) {
  paste(toupper(substring(x, 1,1)), substring(x, 2),
        sep="", collapse=" ")
}

biosamples_metadata <- read.csv("analysis/metadata/biosamples_metadata.tsv", sep = "\t", check.names = F)
combined_file_report <- read.csv("analysis/metadata/output_file_reports/combined_reports.tsv", sep = "\t")[, c("sample_accession", "run_accession")]
biosamples_metadata <- left_join(biosamples_metadata, combined_file_report)

# Acid mine
acid_mine <- biosamples_metadata[biosamples_metadata$study_accession == "PRJNA540505",]
acid_mine$season <- gsub(" .*", "", acid_mine$season)
acid_mine$source <- ifelse(grepl("sediment", acid_mine$`isolation source`), "sediment", "water")
# From https://doi.org/10.1371/journal.pone.0237599.s002:
acid_mine$pH <- case_when(acid_mine$season == "Summer" & acid_mine$source == "sediment" ~ 3.78,
                          acid_mine$season == "Winter" & acid_mine$source == "sediment" ~ 3.56,
                          acid_mine$season == "Summer" & acid_mine$source == "water" ~ 3.36,
                          acid_mine$season == "Winter" & acid_mine$source == "water" ~ 3.86)
acid_mine <- acid_mine[,c("run_accession", "season", "source", "pH")]
colnames(acid_mine) <- c("sample_id", "season", "source", "pH")
acid_mine$source <- ifelse(acid_mine$source == "sediment", "Sediment", "Water")

dir.create(file.path("analysis/metadata/acid_mine/"), showWarnings = FALSE)
write.table(acid_mine, "analysis/metadata/acid_mine/merged_metadata.tsv", sep="\t", row.names = F)

# Gator soil
gators <- biosamples_metadata[biosamples_metadata$study_accession == "PRJNA554694",]
gators$nest_location <- case_when(gators$`lat lon` == "31.96 N 85.10 W" ~ "Eufaula National Wildlife Refuge",
                                  gators$`lat lon` == "30.69 N 87.97 W" | gators$`lat lon` == "30.67 N 87.95 W" ~ "Mobile-Tensaw Delta",
                                  gators$`lat lon` == "29.89 N 92.96 W" | gators$`lat lon` == "29.88 N 92.96 W" | gators$`lat lon` == "29.85 N 92.98 W" | gators$`lat lon` == "29.85 N 92.97 W" ~ "Cameron Parish",
                                  TRUE ~ "J.D. Murphree Wildlife Management Area")
gators$nest_sublocation <- gsub(".*[0-9]", "", gators$name) %>% gsub(pattern="_.*", replacement="")
gators$nest_sublocation <- ifelse(gators$nest_sublocation == "A", "Nest exterior surface", "Egg-chamber after 50% of eggs removed")
gators <- gators[,c("run_accession", "nest_location", "nest_sublocation")]
colnames(gators) <- c("sample_id", "nest_location", "nest_sublocation")

dir.create(file.path("analysis/metadata/gator_soil/"), showWarnings = FALSE)
write.table(gators, "analysis/metadata/gator_soil/merged_metadata.tsv", sep="\t", row.names = F)

# Salt marsh
saltmarsh <- biosamples_metadata[biosamples_metadata$study_accession == "PRJEB42019" & biosamples_metadata$`env biome` == "marine salt marsh biome",]
saltmarsh$location <- saltmarsh$emp500_title
saltmarsh$location <- ifelse(saltmarsh$location == "Massachusetts salt marsh sediments", "Massachusetts", "Spain")
saltmarsh <- saltmarsh[,c("run_accession", "location")]
colnames(saltmarsh) <- c("sample_id", "location")

dir.create(file.path("analysis/metadata/saltmarsh/"), showWarnings = FALSE)
write.table(saltmarsh, "analysis/metadata/saltmarsh/merged_metadata.tsv", sep="\t", row.names = F)

# Tara oceans
tara <- biosamples_metadata[biosamples_metadata$study_accession == "PRJEB9740",]
tara$chlorophyll <- as.numeric(gsub(" .*", "", tara$`Chlorophyll Sensor`))
tara$chlorophyll <- ifelse(tara$chlorophyll == 99999, NA, tara$chlorophyll)
tara$nitrate <- as.numeric(gsub(" .*", "", tara$`nitrate sensor`))
tara$nitrate <- ifelse(tara$nitrate == 99999, NA, tara$nitrate)
tara$oxygen <- as.numeric(gsub(" .*", "", tara$`oxygen sensor`))
tara$oxygen <- ifelse(tara$oxygen == 99999, NA, tara$oxygen)
tara$depth <- ifelse(tara$depth == "10-100 m", "55 m", tara$depth)
tara$depth <- as.numeric(gsub(" .*", "", tara$depth))

tara <- tara[,c("run_accession", "chlorophyll", "nitrate", "oxygen", "depth")]
colnames(tara) <- c("sample_id", "chlorophyll", "nitrate", "oxygen", "depth")

dir.create(file.path("analysis/metadata/tara_polar/"), showWarnings = FALSE)
write.table(tara, "analysis/metadata/tara_polar/merged_metadata.tsv", sep="\t", row.names = F)

# Animal guts

animal <- biosamples_metadata[biosamples_metadata$study_accession == "PRJEB42019" & biosamples_metadata$`env biome` != "marine salt marsh biome",]
animal$diet <- case_when(animal$diet == "frugivore" ~ "herbivore",
                         animal$diet == "herbivore/insectivore" ~ "herbivore",
                         TRUE ~ animal$diet)

fly_conv <- animal$host_flight
names(fly_conv) <- animal$`common name`
fly_conv <- fly_conv[!is.na(fly_conv)]

animal$flight <- fly_conv[animal$`common name`]
animal$flight <- case_when(animal$`common name` == "common quail" ~ "terrestrial, flight",
                         animal$`common name` == "Raccoon dog" ~ "terrestrial, flightless",
                         animal$`common name` == "Sloth bear" ~ "terrestrial, flightless",
                         animal$`common name` == "Dwarf mongoose" ~ "terrestrial, flightless",
                         animal$`common name` == "Sand cat" ~ "terrestrial, flightless",
                         animal$`common name` == "Southern three-banded armadillo" ~ "terrestrial, flightless",
                         animal$`common name` == "Banded mongoose" ~ "terrestrial, flightless",
                         animal$`common name` == "Giant anteater" ~ "terrestrial, flightless",
                         animal$`common name` == "Aldabra giant tortoise" ~ "terrestrial, flightless",
                         animal$`common name` == "Glass lizard" ~ "terrestrial, flightless",
                         animal$`common name` == "Spiny-tailed monitor" ~ "terrestrial, flightless",
                         animal$`common name` == "Spider tortoise" ~ "terrestrial, flightless",
                         animal$`common name` == "Green anaconda" ~ "terrestrial, flightless",
                         animal$`common name` == "Paraguay horned frog" ~ "terrestrial, flightless",
                         animal$`common name` == "Misson golden-eyed tree frog" ~ "terrestrial, flightless",
                         animal$`common name` == "Emu" ~ "terrestrial, flightless",
                         animal$`common name` == "Rhea" ~ "terrestrial, flightless",
                         animal$`common name` == "Lappet faced vulture" ~ "terrestrial, flight",
                         animal$`common name` == "Vulturine guineafowl" ~ "terrestrial, flight",
                         animal$`common name` == "Brown kiwi" ~ "terrestrial, flightless",
                         TRUE ~ animal$flight
                         )

animal$flight <- gsub(",", "", animal$flight)
animal$flight <- case_when(animal$flight == "terrestrial flightless" ~ "Non-flying land animal",
                           animal$flight == "terrestrial flight" ~ "Flying land animal",
                           animal$flight == "marine" ~ "Marine animal")

animal$status <- gsub(" feces", "", animal$emp500_title)
animal$class <- gsub(" .*", "", animal$class)

animal <- animal[,c("run_accession", "diet", "flight", "status", "class")]
colnames(animal) <- c("sample_id", "diet", "flight", "habitat", "class")

animal$habitat <- ifelse(animal$habitat == "Whale", "Ocean", animal$habitat)
animal$habitat <- ifelse(animal$habitat == "Livestock", "Farm", animal$habitat)

animal$diet <- sapply(animal$diet, capitalize)

dir.create(file.path("analysis/metadata/animal_gut/"), showWarnings = FALSE)
write.table(animal, "analysis/metadata/animal_gut/merged_metadata.tsv", sep="\t", row.names = F)

