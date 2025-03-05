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

load("2_data/23_Pubchem/Pubchem_database/pubchem_database.rda")
load("2_data/23_Pubchem/Pubchem_database/pubchem_tax_database.rda")

# 先确保所有列名一致
common_cols <- intersect(colnames(pubchem_combinded_cleaned), colnames(pubtax_database_mergy_filter))
extra_cols <- setdiff(colnames(pubtax_database_mergy_filter), colnames(pubchem_combinded_cleaned))

# 为 pubchem_combinded_cleaned 添加缺失列，并赋值为 "Unknown"
for (col in extra_cols) {
  pubchem_combinded_cleaned[[col]] <- "Unknown"
}

# 进行合并
pubchem_data <- bind_rows(pubtax_database_mergy_filter, pubchem_combinded_cleaned)

pubchem_data_merged <- 
  pubchem_data %>%
  group_by(cid) %>%
  summarize(
    across(!starts_with("from_"), ~ .[sample(length(.), 1)]),
  
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
    from_environment = ifelse(any(from_environment == "Yes"), "Yes", "Unknown"),
    from_which_environment = ifelse(all(from_which_environment == "Unknown"), 
                                    "Unknown", 
                                    paste(unique(na.omit(from_which_environment[from_which_environment != "Unknown"])), collapse = "{}")),
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(any(from_protist == "Yes"), "Yes", "Unknown"),
    from_which_protist = ifelse(all(from_which_protist == "Unknown"), 
                                "Unknown", 
                                paste(unique(na.omit(from_which_protist[from_which_protist != "Unknown"])), collapse = "{}")),
    from_drug = ifelse(any(from_drug == "Yes"), "Yes", "Unknown"),
    from_which_drug = ifelse(all(from_which_drug == "Unknown"), 
                             "Unknown", 
                             paste(unique(na.omit(from_which_drug[from_which_drug != "Unknown"])), collapse = "{}")),
    from_food = ifelse(any(from_food == "Yes"), "Yes", "Unknown"),
    from_which_food = ifelse(all(from_which_food == "Unknown"), 
                             "Unknown", 
                             paste(unique(na.omit(from_which_food[from_which_food != "Unknown"])), collapse = "{}"))
  ) %>%
  ungroup() %>% 
  arrange(cid)

setwd("2_data/23_Pubchem/Pubchem_database")

save(pubchem_data_merged, file = "pubchem_final.rda")


