library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

Lotus <- read_csv("2_data/16_Lotus/lotus.csv")

Lotus_data <- 
  Lotus %>% 
  dplyr::select('_id', lotus_id, wikidata_id, smiles, inchi, inchikey, inchi2D, inchikey2D, smiles2D, traditional_name, iupac_name, molecular_formula, molecular_weight, allTaxa)

setwd("2_data/16_Lotus")

save(Lotus_data, file = "Lotus_data.rda")

Lotus_data$allTaxa.head(10)
head(Lotus_data$allTaxa, 300)

# 创建新的分类列，初始化为"Unknown"
new_cols <- c("from_human", "from_which_part", "from_bacteria", "from_which_bacteria",
              "from_fungi", "from_which_fungi", "from_archaea", "from_which_archaea",
              "from_plant", "from_which_plant", "from_animal", "from_which_animal", 
              "from_environment", "from_which_environment",
              "from_virus", "from_which_virus", "from_protist", "from_which_protist",
              "from_drug", "from_which_drug", "from_food", "from_which_food")

Lotus_dataset <- Lotus_data
Lotus_dataset[, new_cols] <- "Unknown"

# 定义分类的优先级和关键词
categories <- list(
  human = list(
    keywords = c("homo", "human", "homo sapiens"),
    from_col = "from_human",
    which_col = "from_which_part",
    exclude = c("bacteria","cyanobacteria","proteobacteria","actinobacteria","gammaproteobacteria",
                "fungi","ascomycota","basidiomycota","zygomycota","mucoromycota",
                "archaea",
                "plantae","viridiplantae","streptophyta","magnoliopsida","liliopsida","equisetopsida",
                "tracheophyta","rhodophyta","chlorophyta","archaeplastida",
                "animalia","metazoa","chordata","porifera","cnidaria","arthropoda","mollusca",
                "virus","protist","chromista","ochrophyta")
  ),
  bacteria = list(
    keywords = c("bacteria","cyanobacteria","proteobacteria","actinobacteria","gammaproteobacteria"),
    from_col = "from_bacteria",
    which_col = "from_which_bacteria",
    exclude = c("fungi","eukaryota","ascomycota","basidiomycota","zygomycota","mucoromycota",
                "archaea","sorangium cellulosum",
                "plantae","viridiplantae","streptophyta","magnoliopsida","liliopsida","equisetopsida",
                "tracheophyta","rhodophyta","chlorophyta","archaeplastida",
                "animalia","metazoa","chordata","porifera","cnidaria","arthropoda","mollusca",
                "virus","protist","chromista","ochrophyta",
                "homo","human","homo sapiens")
  ),
  fungi = list(
    keywords = c("fungi","ascomycota","basidiomycota","zygomycota","mucoromycota",
                 "sordariomycetes","eurotiomycetes","agaricomycetes","dothideomycetes"),
    from_col = "from_fungi",
    which_col = "from_which_fungi",
    exclude = c("bacteria","cyanobacteria","proteobacteria","actinobacteria","gammaproteobacteria",
                "archaea",
                "plantae","viridiplantae","streptophyta","magnoliopsida","liliopsida","equisetopsida",
                "tracheophyta","rhodophyta","chlorophyta","archaeplastida",
                "animalia","metazoa","chordata","porifera","cnidaria","arthropoda","mollusca",
                "virus","protist","chromista","ochrophyta",
                "homo","human","homo sapiens")
  ),
  archaea = list(
    keywords = c("archaea"),
    from_col = "from_archaea",
    which_col = "from_which_archaea",
    exclude = c("bacteria","fungi","plantae","viridiplantae","streptophyta","magnoliopsida",
                "liliopsida","equisetopsida","tracheophyta","rhodophyta","chlorophyta","archaeplastida",
                "animalia","metazoa","virus","protist","chromista","ochrophyta",
                "homo","human","homo sapiens")
  ),
  plant = list(
    keywords = c("plantae","viridiplantae","streptophyta","magnoliopsida","liliopsida",
                 "equisetopsida","tracheophyta","rhodophyta","chlorophyta","archaeplastida"),
    from_col = "from_plant",
    which_col = "from_which_plant",
    exclude = c("bacteria","fungi","archaea",
                "animalia","metazoa","chordata","porifera","cnidaria","arthropoda","mollusca",
                "virus","protist","chromista","ochrophyta",
                "homo","human","homo sapiens")
  ),
  animal = list(
    keywords = c("animalia","metazoa","chordata","porifera","cnidaria","arthropoda","mollusca",
                 "bryozoa","echinodermata","tunicate","gymnolaemata","ascidiacea",
                 "gastropoda","demospongiae","anthozoa"),
    from_col = "from_animal",
    which_col = "from_which_animal",
    exclude = c("bacteria","fungi","archaea",
                "plantae","viridiplantae","streptophyta","magnoliopsida","liliopsida","equisetopsida",
                "tracheophyta","rhodophyta","chlorophyta","archaeplastida",
                "virus","protist","chromista","ochrophyta",
                "homo","human","homo sapiens")
  ),
  virus = list(
    keywords = c("virus"),
    from_col = "from_virus",
    which_col = "from_which_virus",
    exclude = c("bacteria","fungi","archaea",
                "plantae","viridiplantae","streptophyta","magnoliopsida","liliopsida","equisetopsida",
                "tracheophyta","rhodophyta","chlorophyta","archaeplastida",
                "animalia","metazoa","protist","chromista","ochrophyta",
                "homo","human","homo sapiens")
  ),
  protist = list(
    keywords = c("protist","chromista","ochrophyta","amoebozoa","dinoflagellata","stramenopiles"),
    from_col = "from_protist",
    which_col = "from_which_protist",
    exclude = c("bacteria","fungi","archaea",
                "plantae","viridiplantae","streptophyta","magnoliopsida","liliopsida","equisetopsida",
                "tracheophyta","rhodophyta","chlorophyta","archaeplastida",
                "animalia","metazoa","virus",
                "homo","human","homo sapiens")
  )
)
# 处理每一行的函数
process_row <- function(taxa) {
  taxa_clean <- gsub("\\[|'|\\]", "", taxa)
  taxa_list <- strsplit(taxa_clean, ", ")[[1]]
  taxa_list <- tolower(taxa_list)
  
  result <- list()
  for (col in new_cols) {
    result[[col]] <- "Unknown"
  }
  
  for (cat in categories) {
    if (any(taxa_list %in% cat$keywords)) {
      result[[cat$from_col]] <- "Yes"
      
      candidate_names <- taxa_list[!taxa_list %in% cat$exclude]
      candidate_names <- candidate_names[!candidate_names %in% cat$keywords]
      
      species_names <- grep(" ", candidate_names, value = TRUE, ignore.case = TRUE)
      if (length(species_names) > 0) {
        result[[cat$which_col]] <- species_names[1]
      } else if (length(candidate_names) > 0) {
        result[[cat$which_col]] <- candidate_names[1]
      } else {
        result[[cat$which_col]] <- "Unknown"
      }
      break
    }
  }
  return(result)
}

# 应用处理函数到每一行
processed_data <- lapply(Lotus_dataset$allTaxa, process_row)

# 将结果合并到原数据框
for (col in new_cols) {
  Lotus_dataset[[col]] <- sapply(processed_data, function(x) x[[col]])
}

# 处理空白格（即 ""）为 "Unknown"
Lotus_dataset[, grep("from_which_", names(Lotus_dataset))] <- lapply(
  Lotus_dataset[, grep("from_which_", names(Lotus_dataset))], 
  function(x) sapply(x, function(y) ifelse(y == "" | y == " ", "Unknown", y))
)

# 将结果转换为首字母大写
capitalize <- function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

Lotus_dataset[, grep("from_which_", names(Lotus_dataset))] <- lapply(
  Lotus_dataset[, grep("from_which_", names(Lotus_dataset))], 
  function(x) sapply(x, function(y) ifelse(y == "Unknown", "Unknown", capitalize(y)))
)

Lotus_dataset$from_which_part <- "Unknown"

dir.create("2_data/16_Lotus/Lotus_dataset", showWarnings = FALSE)
setwd("2_data/16_Lotus/Lotus_dataset")

save(Lotus_dataset, file = "Lotus_dataset.rda")

library(stringr)
library(jsonlite)

# 去除单引号，转换为 R 列表
extract_unique_names <- function(source_col) {
  unique(unlist(lapply(source_col, function(x) fromJSON(gsub("'", "\"", x)))))
}

# 提取唯一的名称
unique_list <- extract_unique_names(Lotus_data$allTaxa)
head(unique_list, 3)
unique_list






