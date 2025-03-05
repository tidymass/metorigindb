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

MIMEDB <- 
  read_csv("/Users/ejo/Projects/tidymass/metabolite_database_/2_data/20_MIMEDB/mimedb_metabolites_v1.csv")

# Define keywords for each category
keywords <- list(
  human = c("human", "humans", "endogenous", "body", "urine", "blood", "produced by the body"),
  bacteria = c("bacteria", "microbial", "microflora", "gut flora", "Escherichia", "Clostridium", "Lactobacillus", "Bifidobacterium"),
  fungi = c("fungi", "fungus", "yeast", "Candida", "Aspergillus", "Saccharomyces"),
  archaea = c("archaea", "methanogen", "haloarchaea"),
  plant = c("plant", "vegetable", "fruit", "herb", "leaf", "root", "seed"),
  animal = c("animal", "cow", "milk", "meat", "livestock", "fish", "egg"),
  environment = c("environment", "soil", "water", "air", "marine", "sediment"),
  virus = c("virus", "viral", "HIV", "influenza"),
  protist = c("protist", "protozoa", "algae", "Plasmodium"),
  drug = c("drug", "pharmaceutical", "medication", "antibiotic", "antiviral"),
  food = c("food", "diet", "consumption", "beverage", "milk", "cheese", "vegetables")
)

# 定义分类函数
classify_source <- function(description) {
  result <- list(
    from_human = "Unknown",
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
    from_which_food = "Unknown"
  )
  
  # 遍历关键词列表，检查分类
  for (category in names(keywords)) {
    if (any(str_detect(tolower(description), regex(paste(keywords[[category]], collapse = "|"), ignore_case = TRUE)))) {
      result[[paste0("from_", category)]] <- "Yes"
    }
  }
  
  return(result)
}

# 应用分类到 MIMEDB 数据库
MIMEDB_database <- 
  MIMEDB %>%
  rowwise() %>%
  mutate(
    classification = list(classify_source(description))
  ) %>%
  unnest_wider(classification)


dir.create("2_data/20_MIMEDB/MIMEDB_database", showWarnings = FALSE)
setwd("2_data/20_MIMEDB/MIMEDB_database")

save(MIMEDB_database, file = "MIMEDB_database.rda")







