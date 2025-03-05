library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(dplyr)
library(tidyverse)

setwd("2_data/22_Knapsack/K_combine")

rda_files <- list.files(pattern = "\\.rda$")

data_list <- lapply(rda_files, function(file) {
  load(file)  # 载入 .rda 文件
  x           # 直接使用 x 变量
})

# 合并所有 dataframe
merged_data <- bind_rows(data_list)
# 移除所有列均为 NA 的行
knapsack_data <- 
  merged_data %>%
  dplyr::filter(!if_all(everything(), is.na))

setwd("2_data/22_Knapsack")

save(knapsack_data, file = "knapsack_ms1.rda")
