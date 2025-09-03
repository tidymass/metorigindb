library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

load("2_data/25_Combine_KEGG_HMDB/kegg_database_new3.rda")
load("2_data/25_Combine_KEGG_HMDB/hmdb_database_new4.rda")

source("1_code/merge_functions.R")

hmdb_database_new <- 
  hmdb_database_new %>% 
  mutate(
    from_which_environment = ifelse(from_which_environment == "Unkonwn", "Unknown", from_which_environment),
    from_virus = ifelse(from_virus == "Unkonwn", "Unknown", from_virus),
    from_which_virus = ifelse(from_which_virus == "Unkonwn", "Unknown", from_which_virus),
    from_protist = ifelse(from_protist == "Unkonwn", "Unknown", from_protist),
    from_which_protist = ifelse(from_which_protist == "Unkonwn", "Unknown", from_which_protist)
  )





# 1. 基于kegg id匹配
kegg_match_kegg <- 
  kegg_database_new %>%
  filter(!is.na(KEGG_ID) & KEGG_ID %in% hmdb_database_new$KEGG_ID)

kegg_remaining <- 
  kegg_database_new %>%
  filter(!(KEGG_ID %in% kegg_match_kegg$KEGG_ID))

#save(kegg_remaining, file = "2_data/26_Combine_KEGG_HMDB2/remaining_data/kegg_remaining.rda")

# 将所有数值列转换为字符串
hmdb_database_new <- 
  hmdb_database_new %>%
  mutate(
    Monoisotopic_weight = as.character(Monoisotopic_weight),
    Average_weight = as.character(Average_weight)
  )


df1 <- hmdb_database_new
df2 <- kegg_match_kegg
match_col <- "KEGG_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  #for(i in 1:10) {
  cat(i, " ")
  # 若a列为NA或找不到匹配值则该行保持不变
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    # 在 df2 中找到对应的行
    df2_match <- 
      df2 %>%
      filter(!!sym(match_col) == df1[[match_col]][i])
    
    # 将 df1 的该行与 df2 的匹配行组合成一个临时 dataframe
    df_temp <- bind_rows(df1[i, ], df2_match)
    # 此时 df_temp 有两行，它们的列名都相同
    
    # 使用 summarise + across 对指定列做合并，并只保留一行
    merged_row <- 
      df_temp %>%
      group_by(!!sym(match_col)) %>% 
      summarise(
        
        # 
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, HMDB_ID, CAS_ID, INCHI_ID, INCHIKEY_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # 
        across(c(Synonyms, Formula_all, Compound_description, Database_source, CAS_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_kegg_id <- bind_rows(res_list)
#combined_kegg_id <- do.call(rbind, res_list)



save(combined_kegg_id, file = "2_data/26_Combine_KEGG_HMDB2/combined_data/combined_kegg_id2.rda")




# kegg_id <- 
#   kegg_remaining %>%
#   pull(KEGG_ID) %>%
#   na.omit() %>% 
#   as.list()
# 
# save(kegg_id, file = "2_data/25_Combine_KEGG_HMDB/kegg_id.rda")
# 
# # 假设 kegg_id 已经是一个 list
# kegg_id <- list("C00001", "C00002", "C00030")  # 示例数据
# 
# # 初始化结果列表
# results_list <- vector("list", length(kegg_id))
# 
# 
# # 遍历每个 KEGG ID 并获取 HMDB ID
# for (i in seq_along(kegg_id)) {
#   cat(i, " ")
#   result <- masstools::convert_metabolite_id(
#     query = kegg_id[[i]],
#     from = "KEGG",
#     to = "Human Metabolome Database",
#     top = 1,
#     server = c("cts.fiehnlab", "chemspider", "openai"),
#     chemspider_apikey = "",  # 如有 API Key，请填写
#     openai_apikey = ""
#   )
#   
#   # 存入列表
#   results_list[[i]] <- result
#   
#   Sys.sleep(1)
# }
# 
# # 合并所有 DataFrame 结果
# final_results <- do.call(rbind, results_list)
# 
# # 初始化结果列表
# results_list <- vector("list", length(kegg_id))
# 
# 
# # 设置结果保存的文件夹
# results_dir <- "2_data/25_Combine_KEGG_HMDB/kegg_to_hmdb"  # 你可以修改为你的目标文件夹
# 
# 
# # 遍历每个 KEGG ID 并获取 HMDB ID
# for (i in seq_along(kegg_id)) {
#   cat(i, " ")
#   
#   # 定义保存的文件路径
#   result_file <- file.path(results_dir, paste0(kegg_id[[i]], ".rda"))
#   
#   # 检查文件是否已存在，若存在则跳过
#   if (file.exists(result_file)) {
#     cat("跳过", kegg_id[[i]], "（文件已存在）\n")
#     next
#   }
#   
#   # 执行转换
#   result <- masstools::convert_metabolite_id(
#     query = kegg_id[[i]],
#     from = "KEGG",
#     to = "Human Metabolome Database",
#     top = 1,
#     server = c("cts.fiehnlab", "chemspider", "openai"),
#     chemspider_apikey = "",  # 如有 API Key，请填写
#     openai_apikey = ""
#   )
#   
#   # 保存结果到 .rda 文件
#   save(result, file = result_file)
#   
#   Sys.sleep(1)
# }

results_dir <- "2_data/25_Combine_KEGG_HMDB/kegg_to_hmdb"
# 获取指定目录下所有 .rda 文件
files <- list.files(path = results_dir, pattern = "\\.rda$", full.names = TRUE)

# 初始化空列表
all_results <- list()

# 遍历文件并加载数据
for (file in files) {
  load(file)  # 变量 result 被加载
  all_results[[file]] <- result  # 存入列表
}


# 如果所有 result 均是相同列结构的数据框，可直接使用 bind_rows
final_df <- bind_rows(all_results)

final_df <- 
  final_df %>% 
  rename(KEGG_ID = KEGG)

kegg_database_new2 <- 
  kegg_remaining %>% 
  left_join(final_df, by = "KEGG_ID")

kegg_database_new2 <- 
  kegg_database_new2 %>% 
  select(-HMDB_ID) %>% 
  rename(HMDB_ID = `Human Metabolome Database`)


# 2. 基于hmdb id匹配
kegg_match_hmdb <- 
  kegg_database_new2 %>%
  filter(!is.na(HMDB_ID) & HMDB_ID %in% combined_kegg_id$HMDB_ID)

kegg_remaining2 <- 
  kegg_database_new2 %>%
  filter(!(HMDB_ID %in% kegg_match_hmdb$HMDB_ID))

save(kegg_remaining2, file = "2_data/26_Combine_KEGG_HMDB2/remaining_data/kegg_remaining2.rda")


df1 <- combined_kegg_id
df2 <- kegg_match_hmdb
match_col <- "HMDB_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  #for(i in 1:10) {
  cat(i, " ")
  # 若a列为NA或找不到匹配值则该行保持不变
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    # 在 df2 中找到对应的行
    df2_match <- 
      df2 %>%
      filter(!!sym(match_col) == df1[[match_col]][i])
    
    # 将 df1 的该行与 df2 的匹配行组合成一个临时 dataframe
    df_temp <- bind_rows(df1[i, ], df2_match)
    # 此时 df_temp 有两行，它们的列名都相同
    
    # 使用 summarise + across 对指定列做合并，并只保留一行
    merged_row <- 
      df_temp %>%
      group_by(!!sym(match_col)) %>% 
      summarise(
        
        # 
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, CAS_ID, INCHI_ID, INCHIKEY_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # 
        across(c(Synonyms, Formula_all, Compound_description, Database_source, KEGG_ID, CAS_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_hmdb_id <- bind_rows(res_list)

#save(combined_hmdb_id, file = "2_data/26_Combine_KEGG_HMDB2/combined_data/combined_hmdb_id.rda")


# kegg_id <- 
#   kegg_remaining2 %>%
#   pull(KEGG_ID) %>%
#   na.omit() %>% 
#   as.list()
# 
# save(kegg_id, file = "2_data/25_Combine_KEGG_HMDB/kegg_id2.rda")
# 
# # 假设 kegg_id 已经是一个 list
# kegg_id <- list("C00001", "C00002", "C00030")  # 示例数据
# 
# # 设置结果保存的文件夹
# results_dir <- "2_data/25_Combine_KEGG_HMDB/kegg_to_cas"  # 你可以修改为你的目标文件夹
# 
# # 遍历每个 KEGG ID 并获取 HMDB ID
# for (i in seq_along(kegg_id)) {
#   cat(i, " ")
#   
#   # 定义保存的文件路径
#   result_file <- file.path(results_dir, paste0(kegg_id[[i]], ".rda"))
#   
#   # 检查文件是否已存在，若存在则跳过
#   if (file.exists(result_file)) {
#     cat("跳过", kegg_id[[i]], "（文件已存在）\n")
#     next
#   }
#   
#   # 执行转换
#   result <- masstools::convert_metabolite_id(
#     query = kegg_id[[i]],
#     from = "KEGG",
#     to = "CAS",
#     top = 1,
#     server = c("cts.fiehnlab", "chemspider", "openai"),
#     chemspider_apikey = "",  # 如有 API Key，请填写
#     openai_apikey = ""
#   )
#   
#   # 保存结果到 .rda 文件
#   save(result, file = result_file)
#   
#   Sys.sleep(1)
# }

results_dir <- "2_data/25_Combine_KEGG_HMDB/kegg_to_cas" 

# 获取指定目录下所有 .rda 文件
files <- list.files(path = results_dir, pattern = "\\.rda$", full.names = TRUE)

# 初始化空列表
all_results <- list()

# 遍历文件并加载数据
for (file in files) {
  load(file)  # 变量 result 被加载
  all_results[[file]] <- result  # 存入列表
}


# 如果所有 result 均是相同列结构的数据框，可直接使用 bind_rows
final_df <- bind_rows(all_results)

final_df <- 
  final_df %>% 
  rename(KEGG_ID = KEGG)

kegg_database_new3 <- 
  kegg_remaining2 %>% 
  left_join(final_df, by = "KEGG_ID")

# 
kegg_database_new3 <- 
  kegg_database_new3 %>%
  mutate(
    CAS_ID = coalesce(CAS_ID, CAS)
    ) %>%
  
  select(-CAS)  # 删除原来的 CAS 列


# 2. 基于cas id匹配
kegg_match_cas <- 
  kegg_database_new3 %>%
  filter(!is.na(CAS_ID) & CAS_ID %in% combined_hmdb_id$CAS_ID)

kegg_remaining3 <- 
  kegg_database_new3 %>%
  filter(!(CAS_ID %in% kegg_match_cas$CAS_ID))

#save(kegg_remaining3, file = "2_data/25_Combine_KEGG_HMDB/kegg_remaining3.rda")

save(kegg_remaining3, file = "2_data/26_Combine_KEGG_HMDB2/remaining_data/kegg_remaining3.rda")




df1 <- combined_hmdb_id
df2 <- kegg_match_cas
match_col <- "CAS_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  #for(i in 1:10) {
  cat(i, " ")
  # 若a列为NA或找不到匹配值则该行保持不变
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    # 在 df2 中找到对应的行
    df2_match <- 
      df2 %>%
      filter(!!sym(match_col) == df1[[match_col]][i])
    
    # 将 df1 的该行与 df2 的匹配行组合成一个临时 dataframe
    df_temp <- bind_rows(df1[i, ], df2_match)
    # 此时 df_temp 有两行，它们的列名都相同
    
    # 使用 summarise + across 对指定列做合并，并只保留一行
    merged_row <- 
      df_temp %>%
      group_by(!!sym(match_col)) %>% 
      summarise(
        
        # change here
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, HMDB_ID, INCHI_ID, INCHIKEY_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, Formula_all, Compound_description, Database_source, KEGG_ID, CAS_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_cas_id <- bind_rows(res_list)

#save(combined_hmdb_id, file = "2_data/25_Combine_KEGG_HMDB/HMDB_KEGG_Combined/combined_hmdb_id.rda")
save(combined_cas_id, file = "2_data/26_Combine_KEGG_HMDB2/combined_data/combined_cas_id.rda")


# kegg_id <- 
#   kegg_remaining3 %>%
#   pull(KEGG_ID) %>%
#   na.omit() %>% 
#   as.list()
# 
# save(kegg_id, file = "2_data/26_Combine_KEGG_HMDB2/Kegg_id/kegg_id.rda")
# length(unique(kegg_match_cas$KEGG_ID))
# 
# 
# kegg_id <- kegg_id[1:3000]
# # 假设 kegg_id 已经是一个 list
# kegg_id <- list("C00001", "C00002", "C00030")  # 示例数据
# 
# # 设置结果保存的文件夹
# results_dir <- "2_data/26_Combine_KEGG_HMDB2/kegg_to_inchi"  # 你可以修改为你的目标文件夹
# 
# # 遍历每个 KEGG ID 并获取 HMDB ID
# for (i in seq_along(kegg_id)) {
#   cat(i, " ")
#   
#   # 定义保存的文件路径
#   result_file <- file.path(results_dir, paste0(kegg_id[[i]], ".rda"))
#   
#   # 检查文件是否已存在，若存在则跳过
#   if (file.exists(result_file)) {
#     cat("跳过", kegg_id[[i]], "（文件已存在）\n")
#     next
#   }
#   
#   # 执行转换
#   result <- masstools::convert_metabolite_id(
#     query = kegg_id[[i]],
#     from = "KEGG",
#     to = "InChIKey",
#     top = 1,
#     server = c("cts.fiehnlab", "chemspider", "openai"),
#     chemspider_apikey = "",  # 如有 API Key，请填写
#     openai_apikey = ""
#   )
#   
#   # 保存结果到 .rda 文件
#   save(result, file = result_file)
#   
#   Sys.sleep(1)
# }


results_dir <- "2_data/26_Combine_KEGG_HMDB2/kegg_to_inchi" 

# 获取指定目录下所有 .rda 文件
files <- list.files(path = results_dir, pattern = "\\.rda$", full.names = TRUE)

# 初始化空列表
all_results <- list()

# 遍历文件并加载数据
for (file in files) {
  load(file)  # 变量 result 被加载
  all_results[[file]] <- result  # 存入列表
}


# 如果所有 result 均是相同列结构的数据框，可直接使用 bind_rows
final_df <- bind_rows(all_results)


final_df <- 
  final_df %>% 
  rename(KEGG_ID = KEGG)

kegg_database_new4 <- 
  kegg_remaining3 %>% 
  left_join(final_df, by = "KEGG_ID")

kegg_database_new4 <- 
  kegg_database_new4 %>% 
  select(-INCHIKEY_ID) %>% 
  rename(INCHIKEY_ID = InChIKey)


# 2. 基于inchikey_id匹配
kegg_match_inchikey <- 
  kegg_database_new4 %>%
  filter(!is.na(INCHIKEY_ID) & INCHIKEY_ID %in% combined_cas_id$INCHIKEY_ID)

kegg_remaining4 <- 
  kegg_database_new4 %>%
  filter(!(INCHIKEY_ID %in% kegg_match_inchikey$INCHIKEY_ID))

save(kegg_remaining4, file = "2_data/26_Combine_KEGG_HMDB2/remaining_data/kegg_remaining4.rda")


df1 <- combined_cas_id
df2 <- kegg_match_inchikey
match_col <- "INCHIKEY_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  #for(i in 1:10) {
  cat(i, " ")
  # 若a列为NA或找不到匹配值则该行保持不变
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    # 在 df2 中找到对应的行
    df2_match <- 
      df2 %>%
      filter(!!sym(match_col) == df1[[match_col]][i])
    
    # 将 df1 的该行与 df2 的匹配行组合成一个临时 dataframe
    df_temp <- bind_rows(df1[i, ], df2_match)
    # 此时 df_temp 有两行，它们的列名都相同
    
    # 使用 summarise + across 对指定列做合并，并只保留一行
    merged_row <- 
      df_temp %>%
      group_by(!!sym(match_col)) %>% 
      summarise(
        
        # change here
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, CAS_ID, HMDB_ID, INCHI_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, Formula_all, Compound_description, Database_source, KEGG_ID, CAS_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_inchikey_id <- bind_rows(res_list)
save(combined_inchikey_id, file = "2_data/26_Combine_KEGG_HMDB2/combined_data/combined_inchikey_id.rda")

KEGG_HMDB_Combined_Final <- rbind(combined_inchikey_id, kegg_remaining4)
save(KEGG_HMDB_Combined_Final, file = "2_data/26_Combine_KEGG_HMDB2/combined_data/KEGG_HMDB_Combined_Final2.rda")



