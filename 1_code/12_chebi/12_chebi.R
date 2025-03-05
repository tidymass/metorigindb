library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(tidyverse)

chebi_data <- 
  read_csv("2_data/12_CHEBI/chebi_origin.csv")

filter_eukaryota <- 
  chebi_data %>% 
  dplyr::filter(superkingdom == "Eukaryota")

unique(filter_eukaryota$phylum)

chebi_database <- 
  chebi_data %>% 
  mutate(
    from_human = ifelse(coalesce(species, "") == "Homo sapiens", "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = ifelse(coalesce(superkingdom, "") == "Bacteria", "Yes", "Unknown"),
    from_which_bacteria = ifelse(from_bacteria == "Yes", species, "Unknown"),
    from_fungi = ifelse(
      coalesce(phylum, "") %in% c("Ascomycota", "Basidiomycota", "Mucoromycota"), 
      "Yes", 
      "Unknown"
    ),
    from_which_fungi = ifelse(from_fungi == "Yes", species, "Unknown"),
    from_archaea = ifelse(coalesce(superkingdom, "") == "Archaea", "Yes", "Unknown"),
    from_which_archaea = ifelse(from_archaea == "Yes", species, "Unknown"),
    from_plant = ifelse(
      coalesce(phylum, "") %in% c("Chlorophyta", "Streptophyta"), 
      "Yes", 
      "Unknown"
    ),
    from_which_plant = ifelse(from_plant == "Yes", species, "Unknown"),
    from_animal = ifelse(
      coalesce(phylum, "") %in% c("Chordata", "Arthropoda", "Cnidaria", "Platyhelminthes", "Nematoda", "Porifera", "Mollusca", "Echinodermata", "Bryozoa", "Annelida") & coalesce(species, "") != "Homo sapiens", 
      "Yes", 
      "Unknown"
    ),
    from_which_animal = ifelse(from_animal == "Yes", species, "Unknown"),
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(
      coalesce(phylum, "") %in% c("Euglenozoa", "Bacillariophyta", "Haptophyta", "Evosea", "Rhodophyta"), 
      "Yes", 
      "Unknown"
    ),
    from_which_protist = ifelse(from_protist == "Yes", species, "Unknown"),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown"
  )

# 假设 chebi_database 是已经处理过的数据框
chebi_merged <- 
  chebi_database %>%
  group_by(COMPOUND_ID) %>%
  summarize(
    #across(everything(), ~ .[sample(length(.), 1)]),
    across(
      !starts_with("from_"), 
      ~ .[sample.int(n(), 1)]
    ),
    from_human = ifelse(any(from_human == "Yes"), "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = ifelse(any(from_bacteria == "Yes"), "Yes", "Unknown"),
    from_which_bacteria = ifelse(all(from_which_bacteria == "Unknown"), 
                                 "Unknown", 
                                 paste(unique(na.omit(from_which_bacteria[from_which_bacteria != "Unknown"])), collapse = "{}")),
    from_fungi = ifelse(any(from_fungi == "Yes"), "Yes", "Unknown"),
    from_which_fungi = ifelse(all(from_which_fungi == "Unknown"), 
                              "Unknown", 
                              paste(unique(na.omit(from_which_fungi[from_which_fungi != "Unknown"])), collapse = "{}")),
    from_archaea = ifelse(any(from_archaea == "Yes"), "Yes", "Unknown"),
    from_which_archaea = ifelse(all(from_which_archaea == "Unknown"), 
                                "Unknown", 
                                paste(unique(na.omit(from_which_archaea[from_which_archaea != "Unknown"])), collapse = "{}")),
    from_plant = ifelse(any(from_plant == "Yes"), "Yes", "Unknown"),
    from_which_plant = ifelse(all(from_which_plant == "Unknown"), 
                              "Unknown", 
                              paste(unique(na.omit(from_which_plant[from_which_plant != "Unknown"])), collapse = "{}")),
    from_animal = ifelse(any(from_animal == "Yes"), "Yes", "Unknown"),
    from_which_animal = ifelse(all(from_which_animal == "Unknown"), 
                               "Unknown", 
                               paste(unique(na.omit(from_which_animal[from_which_animal != "Unknown"])), collapse = "{}")),
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(any(from_protist == "Yes"), "Yes", "Unknown"),
    from_which_protist = ifelse(all(from_which_protist == "Unknown"), 
                                "Unknown", 
                                paste(unique(na.omit(from_which_protist[from_which_protist != "Unknown"])), collapse = "{}")),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown"
  ) %>%
  ungroup()



dir.create("2_data/12_CHEBI/chebi_database", showWarnings = FALSE)
setwd("2_data/12_CHEBI/chebi_database")

save(chebi_merged, file = "chebi_database3.rda")



