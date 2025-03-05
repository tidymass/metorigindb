library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(jsonlite)

setwd("2_data/21_YMDB")
YMDB_data <- 
  read_csv("yeast-compounds-YMDB-2025-02-19.csv")

YMDB_database <- 
  YMDB_data %>% 
  mutate(
    from_human = "Unknown",
    from_which_part = "Unknown",
    from_bacteria = "Unknown",
    from_which_bacteria = "Unknown",
    from_fungi = "Yes",
    from_which_fungi = "Saccharomyces cerevisiae",
    from_archaea = "Unknown",
    from_which_archaea = "Unknown",
    from_plant = "Unknown",
    from_which_plant = "Unknown",
    from_animal = "Unknown",
    from_which_animal = "Unknown",
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = "Unknown",
    from_which_protist = "Unknown",
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",
  )

dir.create("2_data/21_YMDB/YMDB_database", showWarnings = FALSE)
setwd("2_data/21_YMDB/YMDB_database")

save(YMDB_database, file = "YMDB_database.rda")

