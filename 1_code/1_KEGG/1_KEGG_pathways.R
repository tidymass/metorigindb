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


####specific organism pathways
###get the organism list
organism_list <-
  keggList(database = "organism") %>%
  as.data.frame()

existing_list <-
  dir("2_data/1_KEGG/pathway/organism_pathways/")


setdiff(organism_list$organism, existing_list)
