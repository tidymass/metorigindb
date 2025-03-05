library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

setwd("2_data/24_Reactome")

# 读取txt文件，假设文件路径为"data.txt"
ChEBI2Reactome_data <- read.table("ChEBI2Reactome.txt", sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE, comment.char = "")

# 生成列名 A, B, C, ..., 依赖数据框的列数
colnames(ChEBI2Reactome_data ) <- c("ChEBI_identifier", "Reactome_Physical_Entity_Stable_Identifier", "Reactome_Physical_Entity_Name", "Reactome_Pathway_Stable_identifier", "URL", "Event_Name", "Evidence_Code", "Species")

unique(ChEBI2Reactome_data$Species)

# 定义各分类下的物种名称向量
human_list <- c(
  "Homo sapiens" 
)

bacteria_list <- c(
  "Mycobacterium tuberculosis"
)


fungi_list <- c(
  "Saccharomyces cerevisiae",  
  "Schizosaccharomyces pombe"
)

archaea_list <- c()

plant_list <- c()

animal_list <- c(
  "Bos taurus",          
  "Canis familiaris",    
  "Drosophila melanogaster", 
  "Mus musculus",       
  "Rattus norvegicus",  
  "Sus scrofa",          
  "Xenopus tropicalis",  
  "Caenorhabditis elegans", 
  "Danio rerio",         
  "Gallus gallus"
)

environment_list <- c(
  # 题目中并没有明确环境来源示例，因此这里先留空
)

virus_list <- c(
  # 题目中没有出现病毒示例，留空
)

protist_list <- c(
  "Dictyostelium discoideum",   
  "Plasmodium falciparum"
)

drug_list <- c(
  # 如果有药物来源可以在此列出
)

food_list <- c(
  # 如果有食物来源可以在此列出
)

ChEBI2Reactome_database <- 
  ChEBI2Reactome_data %>%
  mutate(
    from_human = ifelse(Species %in% human_list, "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = ifelse(Species %in% bacteria_list, "Yes", "Unknown"),
    from_which_bacteria = ifelse(from_bacteria == "Yes", Species, "Unknown"),
    from_fungi = ifelse(Species %in% fungi_list, "Yes", "Unknown"),
    from_which_fungi = ifelse(from_fungi == "Yes", Species, "Unknown"),
    from_archaea = "Unknown",
    from_which_archaea = "Unknown",
    from_plant = "Unknown",
    from_which_plant = "Unknown",
    from_animal = ifelse(Species %in% animal_list, "Yes", "Unknown"),
    from_which_animal = ifelse(from_animal == "Yes", Species, "Unknown"),
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(Species %in% protist_list, "Yes", "Unknown"),
    from_which_protist = ifelse(from_protist == "Yes", Species, "Unknown"),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown"
  )

ChEBI2Reactome_database_cleaned <- 
  ChEBI2Reactome_database %>%
  group_by(ChEBI_identifier) %>%
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
  arrange(ChEBI_identifier)

ChEBI2Reactome_database_cleaned <- 
  ChEBI2Reactome_database %>%
  group_by(ChEBI_identifier) %>%
  summarize(
    # 1) For these columns, collapse all unique values separated by {}
    Reactome_Physical_Entity_Stable_Identifier = paste(unique(Reactome_Physical_Entity_Stable_Identifier), collapse = "{}"),
    Reactome_Physical_Entity_Name = paste(unique(Reactome_Physical_Entity_Name), collapse = "{}"),
    Evidence_Code = paste(unique(Evidence_Code), collapse = "{}"),
    Species = paste(unique(Species), collapse = "{}"),
    
    # 2) For these columns, collapse all values separated by {}
    Reactome_Pathway_Stable_identifier = paste(Reactome_Pathway_Stable_identifier, collapse = "{}"),
    URL = paste(URL, collapse = "{}"),
    Event_Name = paste(Event_Name, collapse = "{}"),
    
    from_human = ifelse(any(from_human == "Yes"), "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = ifelse(any(from_bacteria == "Yes"), "Yes", "Unknown"),
    from_which_bacteria = ifelse(
      all(from_which_bacteria == "Unknown"), 
      "Unknown", 
      paste(unique(na.omit(from_which_bacteria[from_which_bacteria != "Unknown"])), collapse = "{}")
    ),
    from_fungi = ifelse(any(from_fungi == "Yes"), "Yes", "Unknown"),
    from_which_fungi = ifelse(
      all(from_which_fungi == "Unknown"), 
      "Unknown", 
      paste(unique(na.omit(from_which_fungi[from_which_fungi != "Unknown"])), collapse = "{}")
    ),
    from_archaea = "Unknown",
    from_which_archaea = "Unknown",
    from_plant = "Unknown",
    from_which_plant = "Unknown",
    from_animal = ifelse(any(from_animal == "Yes"), "Yes", "Unknown"),
    from_which_animal = ifelse(
      all(from_which_animal == "Unknown"), 
      "Unknown", 
      paste(unique(na.omit(from_which_animal[from_which_animal != "Unknown"])), collapse = "{}")
    ),
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(any(from_protist == "Yes"), "Yes", "Unknown"),
    from_which_protist = ifelse(
      all(from_which_protist == "Unknown"), 
      "Unknown", 
      paste(unique(na.omit(from_which_protist[from_which_protist != "Unknown"])), collapse = "{}")
    ),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",
    
    .groups = "drop"  # ensure ungrouping after summarize
  ) %>%
  arrange(ChEBI_identifier)

dir.create("2_data/22_Knapsack/knapsack_database", showWarnings = FALSE)
setwd("2_data/24_Reactome/Reactome_database")

save(ChEBI2Reactome_database_cleaned, file = "ChEBI2Reactome_database.rda")




