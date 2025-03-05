library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

setwd("2_data/14_LIPIDMAPS")
lipid_data <- read.csv("lipid_result.csv", header = TRUE, stringsAsFactors = FALSE)

lipid_data <- 
  lipid_data %>%
  mutate(across(all_of(c("superkingdom", "phylum", "class", "order", "family", "genus", "species")), ~ na_if(.x, "")))

lipid_dataset <- 
  lipid_data %>% 
  mutate(
    from_human = ifelse(coalesce(species, "") == "Homo sapiens", "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = ifelse(coalesce(superkingdom, "") == "Bacteria", "Yes", "Unknown"),
    from_which_bacteria = ifelse(from_bacteria == "Yes", sub(" \\(#\\d+\\)", "", Curated.from), "Unknown"),
    from_fungi = ifelse(
      coalesce(Curated.from, "") == "Fungi (#4751)" | coalesce(phylum, "") %in% c("Ascomycota", "Mucoromycota", "Basidiomycota", "Zoopagomycota"), 
      "Yes", 
      "Unknown"
    ),
    from_which_fungi = ifelse(
      from_fungi == "Yes", 
      sub(" \\(#\\d+\\)", "", coalesce(Curated.from, "Unknown")), 
      "Unknown"
    ),
    from_archaea = ifelse(coalesce(superkingdom, "") == "Archaea", "Yes", "Unknown"),
    from_which_archaea = ifelse(from_archaea == "Yes", sub(" \\(#\\d+\\)", "", Curated.from), "Unknown"),
    from_plant = ifelse(
      coalesce(phylum, "") %in% c("Chlorophyta", "Streptophyta", "Rhodophyta"), 
      "Yes", 
      "Unknown"
    ),
    from_which_plant = ifelse(from_plant == "Yes", sub(" \\(#\\d+\\)", "", Curated.from), "Unknown"),
    from_animal = ifelse(
      coalesce(phylum, "") %in% c("Arthropoda", "Chordata", "Porifera", "Mollusca", "Nematoda", "Echinodermata", "Cnidaria", "Bryozoa") & from_human != "Yes", 
      "Yes", 
      "Unknown"
    ),
    from_which_animal = ifelse(from_animal == "Yes", sub(" \\(#\\d+\\)", "", Curated.from), "Unknown"),
    from_environment = "Unknown",
    from_which_environment = "Unkonwn",
    from_virus = ifelse(coalesce(superkingdom, "") == "Viruses", "Yes", "Unknown"),
    from_which_virus = ifelse(from_virus == "Yes", sub(" \\(#\\d+\\)", "", Curated.from), "Unknown"),
    from_protist = ifelse(
      coalesce(phylum, "") %in% c("Ciliophora", "Evosea", "Foraminifera", "Bacillariophyta", "Haptophyta", "Endomyxa", "Apicomplexa", "Euglenozoa"), 
      "Yes", 
      "Unknown"
    ),
    from_which_protist = ifelse(from_protist == "Yes", sub(" \\(#\\d+\\)", "", Curated.from), "Unknown"),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",  
  )

lipid_dataset <- 
  lipid_dataset %>% 
  mutate(
    from_synthesis = ifelse(coalesce(Curated.from, "") == "synthetic construct (#32630)", "Yes", "No"),
  )


dir.create("2_data/14_LIPIDMAPS/lipidmaps_final_dataset", showWarnings = FALSE)
setwd("2_data/14_LIPIDMAPS/lipidmaps_final_dataset")

save(lipid_dataset, file = "lipid_final_dataset2.rda")
