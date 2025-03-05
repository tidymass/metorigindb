library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)
library(dplyr)


input_dir <- "2_data/14_LIPIDMAPS/lipidmaps_clawer/lipidmaps_all"
lipidmap_files <- list.files(input_dir, pattern = "\\.rda$", 
                             full.names = TRUE)

lipidmap_list <- list()

for (i in 1:length(lipidmap_files)) {

  lipidmapx = load(lipidmap_files[i])
  
  lipidmap_list[[i]] <- get(lipidmapx)
  
  rm(list = lipidmapx)
}

all_colnames <- unique(unlist(lapply(lipidmap_list, colnames)))

lipidmap_dataset <- do.call(rbind, lapply(lipidmap_list, function(df) {
  # 找到当前dataframe缺少的列
  missing_cols <- setdiff(all_colnames, colnames(df))
  
  # 为缺少的列添加NA值
  if (length(missing_cols) > 0) {
    df[missing_cols] <- NA
  }
  
  # 确保列的顺序一致
  df <- df[all_colnames]
  
  return(df)
}))

rownames(lipidmap_dataset) <- seq_len(nrow(lipidmap_dataset))

dir.create("2_data/14_LIPIDMAPS/lipidmaps_dataset", showWarnings = FALSE)
setwd("2_data/14_LIPIDMAPS/lipidmaps_dataset")

save(lipidmap_dataset, file = "lipidmap_data_all.rda")



