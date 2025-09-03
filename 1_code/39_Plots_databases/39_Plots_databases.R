library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(UpSetR)
library(ggvenn)
library(tidyverse)
library(ComplexUpset)
library(ggupset)
library(ggplot2)


# color for sources
metabolite_source_color <- c(
  "Human" = "#2c61a1",
  "Plant" = "#78938a",
  "Food" = "#f5eddc",
  "Bacteria" = "#0f1c5c",
  "Animal" = "#d2b48c",
  "Environment" = "#8f354b", 
  "Drug" = "#000000"
)

load("2_data/10_BIGG/BIGG_dataset/BIGG_database2.rda")
load("2_data/12_CHEBI/chebi_database/chebi_database3.rda")
load("2_data/7_HMDB/HMDB_Origin/HMDB_Origin_final.rda")
load("2_data/6_Combine/Real_final_combined/kegg_database.rda")
load("2_data/16_Lotus/Lotus_dataset/Lotus_dataset3.rda")
load("2_data/20_MIMEDB/MIMEDB_database/MIMEDB_database3.rda")
load("2_data/24_Reactome/Reactome_database/ChEBI2Reactome_database.rda")
load("2_data/18_T3DB/T3DB_database/t3db_database.rda")


# 1. Venn plot for BIGG

setwd("2_data/39_Plots_databases")

df <- BIGG_database

venn_data <- list(
  Human = which(df$from_human == "Yes"),
  Bacteria = which(df$from_bacteria == "Yes"),
  Animal = which(df$from_animal == "Yes")
)

# color settings
fill_colors <- c(
  metabolite_source_color["Human"] %>% unname(),
  metabolite_source_color["Bacteria"] %>% unname(),
  metabolite_source_color["Animal"] %>% unname()
)


P1 <- 
ggvenn(
    venn_data,
    show_percentage = TRUE,
    fill_color = fill_colors,
    stroke_size = 0.5,
    set_name_size = 5,
    text_size = 5
  ) +
  labs(title = "VennPlot for BIGG") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("BIGG.pdf", plot = P1, device = "pdf", width = 13, height = 8)

# 2. UpsetPlot for CHEBI
df <- chebi_merged

df_bin <- data.frame(
  from_human = as.numeric(df$from_human    == "Yes"),
  from_bacteria = as.numeric(df$from_bacteria == "Yes"),
  from_plant = as.numeric(df$from_plant    == "Yes"),
  from_animal = as.numeric(df$from_animal    == "Yes")
  #from_environment = as.numeric(df$from_environment    == "Yes"),
  #from_drug = as.numeric(df$from_drug    == "Yes"),
  #from_food = as.numeric(df$from_food     == "Yes")
)

set_colors <- c(metabolite_source_color["Human"], metabolite_source_color["Plant"], metabolite_source_color["Animal"], 
                metabolite_source_color["Bacteria"])

P2 <- 
UpSetR::upset(
  df_bin,
  sets = c("from_human", "from_bacteria", 
           "from_plant", "from_animal"),
  sets.bar.color = set_colors,
  main.bar.color = "gray30"
  # number.angles = 15,   # 交集中 metabolite 数字的角度
  # point.size = 3.5,     # 点的大小
  # line.size = 3,        # 交集线的粗细
  # order.by = "freq",    # 以频次进行排序
  # decreasing = TRUE     # 降序排列
  # # text.scale = c(3, 3, 3, 1.4, 3, 1.2)
)

pdf(file = "CHEBI.pdf", width = 13, height = 8)
P2
dev.off()

# 3. Upset Plot for HMDB
df <- hmdb_metabolite

df_bin <- data.frame(
  from_human = as.numeric(df$from_human    == "Yes"),
  from_bacteria = as.numeric(df$from_bacteria == "Yes"),
  from_plant = as.numeric(df$from_plant    == "Yes"),
  from_animal = as.numeric(df$from_animal    == "Yes"),
  from_environment = as.numeric(df$from_environment    == "Yes"),
  from_drug = as.numeric(df$from_drug    == "Yes"),
  from_food = as.numeric(df$from_food     == "Yes")
)

set_colors <- c(metabolite_source_color["Human"], 
                metabolite_source_color["Food"], 
                metabolite_source_color["Drug"], 
                #metabolite_source_color["Bacteria"],
                metabolite_source_color["Bacteria"],
                metabolite_source_color["Environment"],
                metabolite_source_color["Plant"]
                )

P3 <- 
  UpSetR::upset(
    df_bin,
    sets = c("from_human", "from_bacteria", 
             "from_plant", "from_environment",
             "from_drug", "from_food"),
    sets.bar.color = set_colors,
    main.bar.color = "gray30"
    # number.angles = 15,   # 交集中 metabolite 数字的角度
    # point.size = 3.5,     # 点的大小
    # line.size = 3,        # 交集线的粗细
    # order.by = "freq",    # 以频次进行排序
    # decreasing = TRUE     # 降序排列
    # # text.scale = c(3, 3, 3, 1.4, 3, 1.2)
  )

pdf(file = "HMDB.pdf", width = 13, height = 8)
P3
dev.off()
# 4. Upset Plot for KEGG

df <- kegg_database

df_bin <- data.frame(
  from_human = as.numeric(df$from_human    == "Yes"),
  from_bacteria = as.numeric(df$from_bacteria == "Yes"),
  from_plant = as.numeric(df$from_plant    == "Yes"),
  from_animal = as.numeric(df$from_animal    == "Yes"),
  from_environment = as.numeric(df$from_environment    == "Yes"),
  from_drug = as.numeric(df$from_drug    == "Yes"),
  from_food = as.numeric(df$from_food     == "Yes")
)

set_colors <- c(metabolite_source_color["Drug"], 
                metabolite_source_color["Bacteria"], 
                metabolite_source_color["Plant"], 
                metabolite_source_color["Animal"],
                metabolite_source_color["Human"]
)

P4 <- 
  UpSetR::upset(
    df_bin,
    sets = c("from_human", "from_bacteria", 
             "from_plant", "from_animal",
             "from_drug"),
    sets.bar.color = set_colors,
    main.bar.color = "gray30"
    # number.angles = 15,   # 交集中 metabolite 数字的角度
    # point.size = 3.5,     # 点的大小
    # line.size = 3,        # 交集线的粗细
    # order.by = "freq",    # 以频次进行排序
    # decreasing = TRUE     # 降序排列
    # # text.scale = c(3, 3, 3, 1.4, 3, 1.2)
  )

pdf(file = "KEGG.pdf", width = 13, height = 8)
P4
dev.off()

# 5. Upset Plot for Lotus
df <- Lotus_dataset

df_bin <- data.frame(
  from_human = as.numeric(df$from_human    == "Yes"),
  from_bacteria = as.numeric(df$from_bacteria == "Yes"),
  from_plant = as.numeric(df$from_plant    == "Yes"),
  from_animal = as.numeric(df$from_animal    == "Yes"),
  from_environment = as.numeric(df$from_environment    == "Yes"),
  from_drug = as.numeric(df$from_drug    == "Yes"),
  from_food = as.numeric(df$from_food     == "Yes")
)

set_colors <- c(metabolite_source_color["Plant"], 
                metabolite_source_color["Animal"], 
                metabolite_source_color["Bacteria"], 
                metabolite_source_color["Human"]
)

P5 <- 
  UpSetR::upset(
    df_bin,
    sets = c("from_human", "from_bacteria", 
             "from_plant", "from_animal"
             ),
    sets.bar.color = set_colors,
    main.bar.color = "gray30"
    # number.angles = 15,   # 交集中 metabolite 数字的角度
    # point.size = 3.5,     # 点的大小
    # line.size = 3,        # 交集线的粗细
    # order.by = "freq",    # 以频次进行排序
    # decreasing = TRUE     # 降序排列
    # # text.scale = c(3, 3, 3, 1.4, 3, 1.2)
  )

pdf(file = "LOTUS.pdf", width = 13, height = 8)
P5
dev.off()


# 6. Upset Plot for MIMEDB
df <- MIMEDB_database

df_bin <- data.frame(
  from_human = as.numeric(df$from_human    == "Yes"),
  from_bacteria = as.numeric(df$from_bacteria == "Yes"),
  from_plant = as.numeric(df$from_plant    == "Yes"),
  from_animal = as.numeric(df$from_animal    == "Yes"),
  from_environment = as.numeric(df$from_environment    == "Yes"),
  from_drug = as.numeric(df$from_drug    == "Yes"),
  from_food = as.numeric(df$from_food     == "Yes")
)

set_colors <- c(metabolite_source_color["Bacteria"], 
                metabolite_source_color["Animal"], 
                metabolite_source_color["Plant"], 
                metabolite_source_color["Human"],
                metabolite_source_color["Food"],
                metabolite_source_color["Environment"],
                metabolite_source_color["Drug"]
)

P6 <- 
  UpSetR::upset(
    df_bin,
    sets = c("from_human", "from_bacteria", "from_environment", "from_drug", "from_food",
             "from_plant", "from_animal"
    ),
    sets.bar.color = set_colors,
    main.bar.color = "gray30"
    # number.angles = 15,   # 交集中 metabolite 数字的角度
    # point.size = 3.5,     # 点的大小
    # line.size = 3,        # 交集线的粗细
    # order.by = "freq",    # 以频次进行排序
    # decreasing = TRUE     # 降序排列
    # # text.scale = c(3, 3, 3, 1.4, 3, 1.2)
  )

pdf(file = "MIMEDB.pdf", width = 13, height = 8)
P6
dev.off()


# 7. Venn Plot for REACTOME

df <- ChEBI2Reactome_database_cleaned

venn_data <- list(
  Human = which(df$from_human == "Yes"),
  Bacteria = which(df$from_bacteria == "Yes"),
  Animal = which(df$from_animal == "Yes")
)

# color settings
fill_colors <- c(
  metabolite_source_color["Human"] %>% unname(),
  metabolite_source_color["Bacteria"] %>% unname(),
  metabolite_source_color["Animal"] %>% unname()
)


P7 <- 
  ggvenn(
    venn_data,
    show_percentage = TRUE,
    fill_color = fill_colors,
    stroke_size = 0.5,
    set_name_size = 5,
    text_size = 5
  ) +
  labs(title = "VennPlot for REACTOME") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("REACTOME.pdf", plot = P7, device = "pdf", width = 13, height = 8)


# 8. Upset Plot for T3DB
df <- t3db_database

df_bin <- data.frame(
  from_human = as.numeric(df$from_human    == "Yes"),
  from_bacteria = as.numeric(df$from_bacteria == "Yes"),
  from_plant = as.numeric(df$from_plant    == "Yes"),
  from_animal = as.numeric(df$from_animal    == "Yes"),
  from_environment = as.numeric(df$from_environment    == "Yes"),
  from_drug = as.numeric(df$from_drug    == "Yes"),
  from_food = as.numeric(df$from_food     == "Yes")
)

set_colors <- c(metabolite_source_color["Environment"], 
                metabolite_source_color["Food"], 
                metabolite_source_color["Drug"], 
                metabolite_source_color["Animal"],
                metabolite_source_color["Plant"]
)

P8 <- 
  UpSetR::upset(
    df_bin,
    sets = c("from_environment", "from_drug", "from_food",
             "from_plant", "from_animal"
    ),
    sets.bar.color = set_colors,
    main.bar.color = "gray30"
    # number.angles = 15,   # 交集中 metabolite 数字的角度
    # point.size = 3.5,     # 点的大小
    # line.size = 3,        # 交集线的粗细
    # order.by = "freq",    # 以频次进行排序
    # decreasing = TRUE     # 降序排列
    # # text.scale = c(3, 3, 3, 1.4, 3, 1.2)
  )

pdf(file = "T3DB.pdf", width = 13, height = 8)
P8
dev.off()


