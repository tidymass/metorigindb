library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(ComplexUpset)
library(ggplot2)
library(UpSetR)
library(tidyverse)

load("2_data/35_Combine_LOTUS/ms1_database4.rda")

df <- ms1_database

df_bin <- data.frame(
  from_human = as.numeric(df$from_human    == "Yes"),
  from_bacteria = as.numeric(df$from_bacteria == "Yes"),
  from_plant = as.numeric(df$from_plant    == "Yes"),
  from_animal = as.numeric(df$from_animal    == "Yes"),
  from_environment = as.numeric(df$from_environment    == "Yes"),
  from_drug = as.numeric(df$from_drug    == "Yes"),
  from_food = as.numeric(df$from_food     == "Yes")
)

# 使用 UpSetR 的 upset() 函数绘制 UpSet plot
plot <-
UpSetR::upset(
  df_bin,
  sets = c("from_human", "from_bacteria", 
           "from_plant", "from_animal", "from_environment",
           "from_drug", "from_food")
)


UpSetR::upset(df_bin, 
      sets = c("from_human", "from_bacteria", "from_plant", 
               "from_animal", "from_environment", "from_drug", "from_food"),
      sets.bar.color = "#56B4E9")

UpSetR::upset(df_bin, 
      sets = c("from_human", "from_bacteria", "from_plant", 
               "from_animal", "from_environment", "from_drug", "from_food"),
      sets.bar.color = "#56B4E9",
      #order.by = "freq",
      keep.order = FALSE)

# 设置列顺序与颜色一一对应
set_order <- c("from_plant", "from_human","from_food", 
               "from_bacteria", "from_animal", 
                "from_drug", "from_environment")

set_colors <- c("#78938a", "#2c61a1", "#f5eddc", 
                "#0f1c5c", "#d2b48c", "#000000", "#8f354b")

P <- 
UpSetR::upset(df_bin, 
      sets = set_order,
      sets.bar.color = set_colors,
      main.bar.color = "gray30",
      #order.by = "freq",
      keep.order = FALSE)


pdf(file = "/Users/yijiang/Desktop/test.pdf", width = 14, height = 5)
P
dev.off()





