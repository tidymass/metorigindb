library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

### to get KEGG database
library(KEGGgraph)
library(KEGGREST)
library(tidyverse)

load("3_data_analysis/1_KEGG/metabolites/metabolite_info_final.rda")
load("3_data_analysis/2_Organism_origin/organism_origin.rda")

get_from_which <- function(data, which_col) {
  values <- 
    data[[which_col]][data[[which_col]] != "No"] %>% 
    unique() %>% 
    sort()
  
  if (length(values) == 0) {
    "Unknown"  
  } else {
    paste(values, collapse = "{}")
  }
}


Metabolite_origin <- 
  metabolite_database_final %>%
  mutate(Name = str_split(organism, "\\{\\}")) %>%
  unnest(Name) %>% 
  left_join(organism_origin, by = "Name") %>% 
  group_by(KEGG.ID, Compound.name) %>%
  summarize(
    from_human = if_else(any(human == "Yes"), "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = if_else(any(bacteria == "Yes"), "Yes", "Unknown"),
    from_which_bacteria = get_from_which(cur_data(), "which_bacteria"),
    from_fungi = if_else(any(fungi == "Yes"), "Yes", "Unknown"),
    from_which_fungi = get_from_which(cur_data(), "which_fungi"),
    from_archaea = if_else(any(archaea == "Yes"), "Yes", "Unknown"),
    from_which_archaea = get_from_which(cur_data(), "which_archaea"),
    from_plant = if_else(any(plant == "Yes"), "Yes", "Unknown"),
    from_which_plant = get_from_which(cur_data(), "which_plant"),
    from_animal = if_else(any(animal == "Yes"), "Yes", "Unknown"),
    from_which_animal = get_from_which(cur_data(), "which_animal"),
    from_environment = "Unknown",
    from_which_environment = "Unkonwn",
    from_virus = if_else(any(virus == "Yes"), "Yes", "Unknown"),
    from_which_virus = get_from_which(cur_data(), "which_virus"),
    from_protist = if_else(any(protist == "Yes"), "Yes", "Unknown"),
    from_which_protist = get_from_which(cur_data(), "which_protist"),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",
    .groups = "drop"
  )

dir.create("2_data/3_Metabolite_origin", showWarnings = FALSE)
setwd("2_data/3_Metabolite_origin")

save(Metabolite_origin, file = "Metabolite_origin.rda")

  

  
