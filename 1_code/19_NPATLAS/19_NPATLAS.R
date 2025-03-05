library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

setwd("2_data/19_NPATLAS")
NPA_data <- readxl::read_excel("NPAtlas_download_2024_09.xlsx")

NPA_database <- 
  NPA_data %>% 
  mutate(
    from_human = "Unknown",
    from_which_part = "Unknown",
    from_bacteria = ifelse(origin_type == "Bacterium", "Yes", "Unknown"),
    from_which_bacteria = ifelse(from_bacteria == "Yes", origin_species, "Unknown"),
    from_fungi = ifelse(origin_type == "Fungus", "Yes", "Unknown"),
    from_which_fungi = ifelse(from_fungi == "Yes", origin_species, "Unknown"),
    from_archaea = "Unknown",
    from_which_archaea = "Unknown",
    from_plant = "Unknown",
    from_which_plant = "Unknown",
    from_animal = "Unknown",
    from_which_animal = "Unknown",
    from_environment = "Unknown",
    from_which_environment = "Unkonwn",
    from_virus = "Unkonwn",
    from_which_virus = "Unkonwn",
    from_protist = "Unkonwn",
    from_which_protist = "Unkonwn",
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",  
  )

dir.create("2_data/19_NPATLAS/NPATLAS_dataset", showWarnings = FALSE)
setwd("2_data/19_NPATLAS/NPATLAS_dataset")

save(NPA_database, file = "NPA_database.rda")
