library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

setwd("2_data/18_T3DB")

t3db_data <- read.csv("toxins.csv", header = TRUE)

t3db_database <-
  t3db_data %>%
  mutate(
    from_human = "Unknown",
    from_which_human = "Unknown",
    from_bacteria = "Unknown",
    from_which_bacteria = "Unknown",
    from_fungi = "Unknown",
    from_which_fungi = "Unknown",
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
    from_which_food = "Unknown"
  )

t3db_database <- 
  t3db_database %>%
  mutate(
    from_drug = ifelse(str_detect(Categories, "\\bDrug\\b"), "Yes", from_drug),
    from_food = ifelse(str_detect(Categories, "\\bFood Toxin\\b"), "Yes", from_food),
    from_plant = ifelse(str_detect(Categories, "\\bPlant Toxin\\b"), "Yes", from_plant),
    from_environment = ifelse(
      str_detect(Categories, "\\bPollutant\\b|\\bAirborne Pollutant\\b|\\bHousehold Toxin\\b|\\bIndustrial/Workplace Toxin\\b|\\bCigarette Toxin\\b|\\bSynthetic Toxin\\b|\\bPesticide\\b"),
      "Yes", 
      from_environment
    ),
    from_animal = ifelse(str_detect(Categories, "\\bAnimal Toxin\\b"), "Yes", from_animal)
  )

dir.create("2_data/18_T3DB/T3DB_database", showWarnings = FALSE)
setwd("2_data/18_T3DB/T3DB_database")

save(t3db_database, file = "t3db_database.rda")

test <- t3db_ms1@spectra.info
