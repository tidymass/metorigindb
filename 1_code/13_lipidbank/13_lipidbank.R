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

load("2_data/13_LIPIDBANK/lipidbank_ms1.rda")
lipibank_data <- 
  lipidbank_ms1@spectra.info 

lipidbank_database <- 
  lipibank_data %>% 
  dplyr::select(-`source[, -remove_idx]`, -From_animal, -From_microbiota, -From_archaea, -From_bacteria, -From_fungi, -From_food, -From_plant,
                -From_eukaryota, -From_other, -From_environment, -From_human, -From_drug)

lipidbank_database


filter1 <- 
  lipidbank_database %>%
  dplyr::filter(!is.na(Source))

unique_source <- 
  unique(filter1$Source)



library(stringr)

# 假设你的数据存储在lipidbank_data中
source_texts <- filter1$Source

# 使用正则表达式去除<< Ref. ... >>
cleaned_texts <- str_replace_all(source_texts, "<< Ref\\. .*?>>", "")

# 按句号拆分，并去掉前后空格
split_texts <- unlist(str_split(cleaned_texts, "\\.", simplify = FALSE))
split_texts <- trimws(split_texts)

# 去除空字符串并去重
unique_results <- unique(split_texts[split_texts != ""])

unique_results

# 假设lipidbank_data为你的数据框
source_texts <- lipidbank_database$Source

# 初始化所有列
cols <- c(
  "from_human", "from_which_part", "from_bacteria", "from_which_bacteria",
  "from_fungi", "from_which_fungi", "from_archaea", "from_which_archaea",
  "from_plant", "from_which_plant", "from_animal", "from_which_animal",
  "from_environment", "from_which_environment", "from_virus", "from_which_virus",
  "from_protist", "from_which_protist", "from_drug", "from_which_drug",
  "from_food", "from_which_food"
)

lipidbank_database[cols] <- "Unknown"

# 分类关键词定义
category_keywords <- list(
  from_human = "human|patient|infant|child|adult|newborn",
  from_which_part = "urine|serum|bile|feces|liver|plasma|blood|meconium|amniotic fluid|duodenal bile",
  from_bacteria = "bacteria|bacterial|microbial|microbiological|Pseudomonas|Mycobacterium|Streptomyces|Corynebacterium",
  from_which_bacteria = "Pseudomonas.*?|Mycobacterium.*?|Streptomyces.*?|Corynebacterium.*?|Achromobacter.*?|Serratia.*?|Rhodopseudomonas.*?",
  from_fungi = "fungus|fungi|Penicillium|Cunninghamella",
  from_which_fungi = "Penicillium.*?|Cunninghamella.*?",
  from_animal = "animal|mammal|seal|sealions|walrus|snake|viper|bird|owl|heron|pelican|python|monkey|rat|dog|manatee|shark|ray|fish|frog|toad|bullfrog|coelacanth|rabbit|alligator|hagfish|sturgeon|lamprey|turtle",
  from_which_animal = paste(
    c("seal", "snake", "viper", "manatee", "shark", "ray", "rabbit", "coelacanth", "bullfrog", "frog", "toad", "alligator", "hagfish", "sturgeon", "lamprey", "turtle", "Rana.*?", "Bufo.*?", "Latimeria.*?", "Vipera.*?", "Bitis.*?", "Hydrurga.*?", "Zalophus.*?", "Trichechus.*?", "Heptatretus.*?"), collapse = "|")
)


# 提取信息
for (i in seq_along(source_texts)) {
  text <- source_texts[i]
  
  for (cat in names(category_keywords)) {
    pattern <- category_keywords[[cat]]
    matches <- unique(str_extract_all(text, regex(pattern, ignore_case = TRUE))[[1]])
    
    if (length(matches) > 0) {
      if (str_detect(cat, "from_which")) {
        lipidbank_database[i, cat] <- paste(unique(matches), collapse = "{}")
      } else {
        lipidbank_database[i, cat] <- "Yes"
      }
    } else {
      lipidbank_database[i, cat] <- "Unknown"
    }
  }
}

# 规则：当from_*列为Unknown时，相应的from_which_*列也必须为Unknown
mapping <- list(
  "from_human" = "from_which_part",
  "from_bacteria" = "from_which_bacteria",
  "from_fungi" = "from_which_fungi",
  "from_archaea" = "from_which_archaea",
  "from_plant" = "from_which_plant",
  "from_animal" = "from_which_animal",
  "from_environment" = "from_which_environment",
  "from_virus" = "from_which_virus",
  "from_protist" = "from_which_protist",
  "from_drug" = "from_which_drug",
  "from_food" = "from_which_food"
)

for (main_col in names(mapping)) {
  which_col <- mapping[[main_col]]
  lipidbank_database[[which_col]][lipidbank_database[[main_col]] == "Unknown"] <- "Unknown"
}

lipidbank_database[lipidbank_database == "NA"] <- "Unknown"

dir.create("2_data/13_LIPIDBANK/lipibank_database", showWarnings = FALSE)
setwd("2_data/13_LIPIDBANK/lipibank_database")

save(lipidbank_database, file = "lipibank_database2.rda")
