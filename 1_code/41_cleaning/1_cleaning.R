library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

load("2_data/41_MetOriginDB/MetOriginDB.rda")

colnames(MetOriginDB)

MetOriginDB_clean <- MetOriginDB %>%
  select(-from_fungi, -from_which_fungi, -from_archaea, -from_which_archaea, -from_virus, 
         -from_which_virus, -from_protist, -from_which_protist)

colnames(MetOriginDB_clean)
save(MetOriginDB_clean, file = "2_data/41_MetOriginDB/MetOriginDB_clean.rda")

