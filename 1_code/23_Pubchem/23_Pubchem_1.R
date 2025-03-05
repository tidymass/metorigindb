library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

setwd("2_data/23_Pubchem/Pubchem_compounds")

agrochemical_data <- read.csv("Agrochemical Information.csv", header = TRUE)
drug_data <- read.csv("Drug.csv", header = TRUE)
food_data <- read.csv("Food.csv", header = TRUE)

agrochemical_database <- 
  agrochemical_data %>% 
  mutate(
    from_environment = "Yes",
    from_drug = "Unknown",
    from_food = "Unknown"
  )

drug_database <- 
  drug_data %>% 
  mutate(
    from_environment = "Unknown",
    from_drug = "Yes",
    from_food = "Unknown"
  )

food_database <- 
  food_data %>% 
  mutate(
    from_environment = "Unknown",
    from_drug = "Unknown",
    from_food = "Yes"
  )

pubchem_combinded <- 
  rbind(agrochemical_database, drug_database, food_database)

pubchem_combinded_cleaned <- 
  pubchem_combinded %>% 
  group_by(cid) %>% 
  summarise(
    across(
      !starts_with("from_"), 
      ~ .[sample.int(n(), 1)]
    ),
    from_environment = ifelse(any(from_environment == "Yes"), "Yes", "Unknown"),
    from_drug = ifelse(any(from_drug == "Yes"), "Yes", "Unknown"),
    from_food = ifelse(any(from_food == "Yes"), "Yes", "Unknown")
  ) %>% 
  ungroup() %>% 
  arrange(cid)

dir.create("2_data/23_Pubchem/Pubchem_database", showWarnings = FALSE)
setwd("2_data/23_Pubchem/Pubchem_database")

save(pubchem_combinded_cleaned, file = "pubchem_database.rda")










