library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(ggvenn)
library(tidyverse)

df <- ms1_database

# 构建用于 Venn 图的数据列表
venn_data <- list(
  Human = which(df$from_human == "Yes"),
  Environment = which(df$from_environment == "Yes")
)

# 颜色设置
metabolite_source_color <- c(
  "Human" = "#2c61a1",
  "Plant" = "#78938a",
  "Food" = "#f5eddc",
  "Bacteria" = "#0f1c5c",
  "Animal" = "#d2b48c",
  "Environment" = "#8f354b", 
  "Drug" = "#000000"
)

# 提取颜色配置
fill_colors <- c(
  metabolite_source_color["Human"] %>% unname(),
  metabolite_source_color["Environment"] %>% unname()
)

# 画图
venn_plot <- 
ggvenn(
    venn_data,
    show_percentage = TRUE,
    fill_color = fill_colors,
    stroke_size = 0.5,
    set_name_size = 5,
    text_size = 5
  ) +
  labs(title = "Metabolite Sources") +
  theme(plot.title = element_text(hjust = 0.5))


