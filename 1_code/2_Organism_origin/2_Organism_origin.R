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

# get all organism info
orgs <- keggList("organism")
orgs <- as.data.frame(orgs)
colnames(orgs) <- c("KEGGid", "Name", "Species", "Phylogeny")

# get organism name list from metabolite table
load("2_data/1_KEGG/metabolite_info/metabolite_info_final.rda")
organism_list <- 
  metabolite_database_final %>% 
  pull(organism) %>% 
  str_split(pattern = "\\{\\}") %>%
  unlist() %>%
  unique() %>% 
  sort()

# get origin table for organisms in the list
organism_origin <-
  orgs %>% 
  filter(Name %in% organism_list) %>% 
  mutate(
    human = if_else(str_detect(Species, "Homo sapiens"), "Yes", "No"),
    bacteria = if_else(str_detect(Phylogeny, "Bacteria"), "Yes", "No"),
    which_bacteria = if_else(bacteria == "Yes", Species, "No"),
    fungi = if_else(str_detect(Phylogeny, "Fungi"), "Yes", "No"),
    which_fungi = if_else(fungi == "Yes", Species, "No"),
    archaea = if_else(str_detect(Phylogeny, "Archaea"), "Yes", "No"),
    which_archaea = if_else(archaea == "Yes", Species, "No"),
    plant = if_else(str_detect(Phylogeny, "Plants"), "Yes", "No"),
    which_plant = if_else(plant == "Yes", Species, "No"),
    animal = if_else(str_detect(Phylogeny, "Animals") & !str_detect(Species, "Homo sapiens"), "Yes", "No"),
    which_animal = if_else(animal == "Yes", Species, "No"),
    virus = if_else(str_detect(Phylogeny, "Virus"), "Yes", "No"),
    which_virus = if_else(virus == "Yes", Species, "No"),
    protist = if_else(str_detect(Phylogeny, "Protists"), "Yes", "No"),
    which_protist = if_else(protist == "Yes", Species, "No")
  ) %>%
  select(-Species, -Phylogeny) %>% 
  arrange(Name)

dir.create("3_data_analysis/2_Organism_origin", showWarnings = FALSE)
setwd("3_data_analysis/2_Organism_origin")

save(organism_origin, file = "organism_origin.rda")



  
