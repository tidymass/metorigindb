library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

setwd("2_data/22_Knapsack")

load("knapsack_ms1.rda")

knapsack_data_expanded <- 
  knapsack_data %>% 
  mutate(across(c(Kingdom, Family, Species), ~ gsub("\\{\\}", "|", .))) %>%  
  separate_rows(Kingdom, Family, Species, sep = "\\|")

knapsack_dataset <- 
  knapsack_data_expanded %>% 
  mutate(
    from_human = ifelse(
      grepl("Homo sapiens", Species, ignore.case = TRUE),
      "Yes",
      "Unknown"
    ),
    from_which_part = ifelse(from_human == "Yes", Species, "Unknown"),
    from_bacteria = ifelse(Kingdom %in% c("Bacteria", "Firmicutes"), "Yes", "Unknown"),
    from_which_bacteria = ifelse(from_bacteria == "Yes", Species, "Unknown"),
    from_fungi = ifelse(Kingdom == "Fungi", "Yes", "Unknown"),
    from_which_fungi = ifelse(from_fungi == "Yes", Species, "Unknown"),
    from_archaea = ifelse(Kingdom %in% c("Euryarchaeota", "Crenarchaeota", "Archaea"), "Yes", "Unknown"),
    from_which_archaea = ifelse(from_archaea == "Yes", Species, "Unknown"),
    from_plant = ifelse(Kingdom %in% c("Plantae", "Plant", "Viridiplantae"), "Yes", "Unknown"),
    from_which_plant = ifelse(from_plant == "Yes", Species, "Unknown"),
    from_animal = ifelse(Kingdom %in% c("Animalia", "Sea fans") & from_human != "Yes", "Yes", "Unknown"),
    from_which_animal = ifelse(from_animal == "Yes", Species, "Unknown"),
    from_environment = "Unknown",
    from_which_environment = "Unkonwn",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(Kingdom %in% c("Chromalveolata", "Excavata", "Chromista", "Protozoa", "Amoebozoa"), "Yes", "Unknown"),
    from_which_protist = ifelse(from_protist == "Yes", Species, "Unknown"),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",  
  )

# 假设 chebi_database 是已经处理过的数据框
knapsack_dataset_cleaned <- 
  knapsack_dataset %>%
  group_by(C_ID) %>%
  summarize(
    across(
      !starts_with("from_"), 
      ~ .[sample.int(n(), 1)]
    ),
    from_human = ifelse(any(from_human == "Yes"), "Yes", "Unknown"),
    from_which_part = ifelse(all(from_which_part == "Unknown"), 
                             "Unknown", 
                             paste(unique(na.omit(from_which_part[from_which_part != "Unknown"])), collapse = "{}")),
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
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(any(from_protist == "Yes"), "Yes", "Unknown"),
    from_which_protist = ifelse(all(from_which_protist == "Unknown"), 
                                "Unknown", 
                                paste(unique(na.omit(from_which_protist[from_which_protist != "Unknown"])), collapse = "{}")),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown"
  ) %>%
  ungroup() %>% 
  arrange(C_ID)

# 1. 记录 df2 的原始列名（顺序）
orig_cols <- names(knapsack_dataset_cleaned)

# 2. 移除 df2 中待替换的列，然后通过 C_ID 从 df1 拼回新的列
knapsack_dataset_cleaned <- 
  knapsack_dataset_cleaned %>%
  select(-Kingdom, -Family, -Species) %>%
  left_join(
    knapsack_data %>% select(C_ID, Kingdom, Family, Species),
    by = "C_ID"
  )

# 3. 按照 df2 原有列的顺序重新排列
knapsack_dataset_cleaned <- 
  knapsack_dataset_cleaned %>%
  select(all_of(orig_cols))

dir.create("2_data/22_Knapsack/knapsack_database", showWarnings = FALSE)
setwd("2_data/22_Knapsack/knapsack_database")

save(knapsack_dataset_cleaned, file = "knapsack_database.rda")




