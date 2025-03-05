library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

load("2_data/17_RECON3/recon3_compound_data")

recon3d_database <- 
  recon3_compound_data %>% 
  mutate(
    from_human = "Yes",
    from_which_part = "Unknown",
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
    from_which_food = "Unknown",
    
  )

dir.create("2_data/17_RECON3/recon3d_database", showWarnings = FALSE)
setwd("2_data/17_RECON3/recon3d_database")

save(recon3d_database, file = "recon3d_database.rda")

