library(r4projects)
library(patchwork)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

load("2_data/6_Combine/Real_final_combined/kegg_database.rda")      #KEGG
load("2_data/7_HMDB/HMDB_Origin/HMDB_Origin_final.rda")             #HMDB
load("2_data/8_FOODB/FOODB_database/FOODB_database.rda")            #FOODB
load("2_data/9_DRUGBANK/DRUGBANK_database/drugbank_database.rda")   #DRUGBANK
load("2_data/12_CHEBI/chebi_database/chebi_database3.rda")          #CHEBI
load("2_data/10_BIGG/BIGG_dataset/BIGG_database2.rda")              #BIGG
load("2_data/20_MIMEDB/MIMEDB_database/MIMEDB_database3.rda")       #MIMEDB
load("2_data/24_Reactome/Reactome_database/ChEBI2Reactome_database.rda") #REACTOME
load("2_data/18_T3DB/T3DB_database/t3db_database.rda")              #T3DB
load("2_data/16_Lotus/Lotus_dataset/Lotus_dataset3.rda")            #LOTUS
load("2_data/15_MODELSEED/modelseed_dataset/modelseed_database_final.rda") #MODELSEED


# Define the color palette
metabolite_source_color <- c(
  "Human" = "#2c61a1",
  "Plant" = "#78938a", 
  "Food" = "#f5eddc",
  "Bacteria" = "#0f1c5c",
  "Animal" = "#d2b48c", 
  "Environment" = "#8f354b",
  "Drug" = "#000000"
)


# Function to process each database and count metabolites by source
process_database <- function(data, db_name) {
  # Select only the relevant columns (exclude other from_which* columns)
  relevant_cols <- c("from_human", "from_bacteria", "from_plant", 
                     "from_animal", "from_drug", "from_environment", "from_food")
  
  # Make sure we only select columns that actually exist in the data
  cols_to_use <- relevant_cols[relevant_cols %in% names(data)]
  
  # Filter to retain only the needed columns
  df <- data %>% select(all_of(cols_to_use))
  
  # Count number of 'Yes' values for each source
  counts <- sapply(df, function(x) sum(x == "Yes", na.rm = TRUE))
  
  # Create a data frame with the counts
  result <- data.frame(
    database = db_name,
    source = c("Human", "Bacteria", "Plant", "Animal", "Drug", "Environment", "Food"),
    count = c(
      counts["from_human"], 
      counts["from_bacteria"], 
      counts["from_plant"],
      counts["from_animal"], 
      counts["from_drug"], 
      counts["from_environment"], 
      counts["from_food"]
    )
  )
  return(result)
}

results_foodb <- process_database(foodb_database, "FooDB")
results_hmdb <- process_database(hmdb_metabolite, "HMDB")
results_kegg <- process_database(kegg_database, "KEGG")
results_drugbank <- process_database(drugbank_database, "DRUGBANK")
results_chebi <- process_database(chebi_merged, "CHEBI")
results_bigg <- process_database(BIGG_database, "BIGG")
results_mimedb <- process_database(MIMEDB_database, "MIMEDB")
results_reactome <- process_database(ChEBI2Reactome_database_cleaned, "REACTOME")
results_t3db <- process_database(t3db_database, "T3DB")
results_lotus <- process_database(Lotus_dataset, "LOTUS")
results_modelseed <- process_database(modelseed_database_final, "MODELSEED")
  
  

# Combine all results
all_results <- bind_rows(
  results_foodb,
  results_hmdb,
  results_kegg,
  results_drugbank,
  results_chebi,
  results_bigg,
  results_mimedb,
  results_reactome,
  results_t3db,
  results_lotus,
  results_modelseed
) %>%
  mutate(
    source = factor(
      source,
      levels = c("Environment", "Drug", "Animal", "Bacteria", "Food", "Human", "Plant")
    )
  )

# 仅仅总数
P1 <- 
  ggplot(all_results_total, aes(x = database, y = total_count)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Metabolite Total Counts by Database", 
    x = "Database", 
    y = "Number of Metabolites"
  ) +
  theme_minimal()

# 百分比图
P2 <- ggplot(all_results, aes(x = database, y = count, fill = source)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = metabolite_source_color, name = "Metabolite Source") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Metabolite Sources (Proportions)", x = "Database", y = "Proportion") +
  theme_bw()


combined_plot <- 
P1 / P2

setwd("2_data/37_BarPlots")
# 保存为 PDF
ggsave("metabolite_source_summary3.pdf", plot = combined_plot, width = 12, height = 10)






# 总数图
# P1 <- ggplot(all_results, aes(x = database, y = count, fill = source)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = metabolite_source_color, name = "Metabolite Source") +
#   labs(title = "Metabolite Sources (Total Counts)", x = "Database", y = "Number of Metabolites") +
#   theme_minimal()

# 计算每个数据库的总数（即原始数据框的行数）
all_results_total <- data.frame(
  database = c("FOODB", "HMDB", "KEGG", "DRUGBANK", "CHEBI", 
               "BIGG", "MIMEDB", "REACTOME", "T3DB", "LOTUS", "MODELSEED"),
  total_count = c(
    nrow(foodb_database),
    nrow(hmdb_metabolite),
    nrow(kegg_database),
    nrow(drugbank_database),
    nrow(chebi_merged),
    nrow(BIGG_database),
    nrow(MIMEDB_database),
    nrow(ChEBI2Reactome_database_cleaned),
    nrow(t3db_database),
    nrow(Lotus_dataset),
    nrow(modelseed_database_final)
  )
)

# 按字母顺序设置数据库顺序
all_results_total$database <- factor(all_results_total$database, levels = sort(unique(all_results_total$database)))

# 绘图
P1 <- 
ggplot(all_results_total, aes(x = database, y = total_count)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = total_count), vjust = -0.5, size = 3.5) +
  labs(
    title = "Total Number of Metabolites per Database", 
    x = "Database", 
    y = "Number of Metabolites"
  ) +
  scale_y_continuous(labels = scales::label_comma()) +  # 使用完整数字格式
  theme_bw()



P1 /P2

