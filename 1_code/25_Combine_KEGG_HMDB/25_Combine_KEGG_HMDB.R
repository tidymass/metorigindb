library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

load("2_data/6_Combine/Real_final_combined/kegg_database.rda")
load("2_data/7_HMDB/HMDB_Origin/HMDB_Origin_final.rda")

hmdb_database <- hmdb_metabolite


count_braces <- str_count(kegg_database$PubChem_ID, "\\{")
rows_with_multiple_braces <- which(count_braces > 1)

any(is.na(kegg_database$PubChem_ID))

# # 按 {} 分割，最多保留前两个值
# split_values <- str_split_fixed(kegg_database$PubChem_ID, "\\{\\}", 2)
# 
# # 添加新列，并将空字符串 ("") 替换为 NA
# kegg_database <- 
#   kegg_database %>%
#   mutate(
#     PUBCHEM_COMPOUND_ID = na_if(split_values[, 1], ""),
#     PUBCHEM_SUBSTANCE_ID = na_if(split_values[, 2], "")
#   )

# 使用strsplit按"{}"分割，然后取第1个、2个元素
split_list <- strsplit(kegg_database$PubChem_ID, "\\{\\}")

# 提取第一段
kegg_database$PUBCHEM_COMPOUND_ID <- sapply(split_list, function(x) {
  if (length(x) >= 1) x[1] else NA
})

# 提取第二段
kegg_database$PUBCHEM_SUBSTANCE_ID <- sapply(split_list, function(x) {
  if (length(x) >= 2) x[2] else NA
})

kegg_database <- 
  kegg_database %>% 
  dplyr::rename(Synonyms = NAME)

# 按{}分割后，仅保留第一部分
kegg_database$NAME <- sapply(strsplit(kegg_database$Synonyms, "\\{\\}"), function(x) x[1])

any(str_detect(kegg_database$Formula, "\\{\\}"))

# split formula
kegg_database <- 
  kegg_database %>% 
  dplyr::rename(Formula_all = FORMULA)

kegg_database <- 
  kegg_database %>% 
  dplyr::rename(CAS_ID_all = CAS_ID)

kegg_database$Formula <- sapply(strsplit(kegg_database$Formula_all, "\\{\\}"), function(x) x[1])

save(kegg_database, file = "2_data/25_Combine_KEGG_HMDB/kegg_database.rda")


kegg_database_new <- 
  kegg_database %>%
  mutate(
    Compound_name         = NAME,
    Monoisotopic_weight   = EXACT_MASS,
    Average_weight        = MOL_WEIGHT,
    Compound_description  = COMMENT,
    SMILES_ID             = NA,
    SMILES_ID_all         = NA,
    INCHI_ID              = NA,
    INCHIKEY_ID           = NA,
    CAS_ID                = word(CAS_ID_all, 1, sep = "\\{\\}"),
    KEGG_ID               = COMPOUND_ID,
    KEGG_DRUG_ID          = DRUG_ID,
    HMDB_ID               = NA,
    CHEMSPIDER_ID         = NA,
    DRUGBANK_ID           = NA,  
    FOODB_ID              = NA,
    CHEBI_ID              = ChEBI_ID,
    CHEMBL_ID             = ChEMBL_ID,
    PDB_CCD_ID            = `PDB-CCD_ID`,
    KNAPSACK_ID           = KNApSAcK_ID,
    LIPIDMAPS_ID          = LIPIDMAPS_ID,
    LIPIDBANK_ID          = LipidBank_ID,
    BIOCYC_ID             = NA,
    BIGG_ID               = NA,
    WIKIPEDIA_ID          = NA,
    METLIN_ID             = NA,
    Database_source       = "KEGG"
  ) %>%
  dplyr::select(
    Compound_name,
    Synonyms,
    Formula,
    Formula_all,
    Monoisotopic_weight,
    Average_weight,
    Compound_description,
    HMDB_ID,
    KEGG_ID,
    CAS_ID,
    CAS_ID_all,
    INCHI_ID,
    INCHIKEY_ID,
    SMILES_ID,
    SMILES_ID_all,
    KEGG_DRUG_ID,
    CHEMSPIDER_ID,
    DRUGBANK_ID,
    FOODB_ID,
    PUBCHEM_COMPOUND_ID,
    PUBCHEM_SUBSTANCE_ID,
    CHEBI_ID,
    CHEMBL_ID,
    PDB_CCD_ID,
    `3DMET_ID`,
    `NIKKAJI_ID`,
    `KNAPSACK_ID`,
    `LIPIDMAPS_ID`,
    `LIPIDBANK_ID`,
    `BIOCYC_ID`,
    `BIGG_ID`,
    `WIKIPEDIA_ID`,
    `METLIN_ID`,
    from_human,
    from_which_part,
    from_bacteria,
    from_which_bacteria,
    from_fungi,
    from_which_fungi,
    from_archaea,
    from_which_archaea,
    from_plant,
    from_which_plant,
    from_animal,
    from_which_animal,
    from_environment,
    from_which_environment,
    from_virus,
    from_which_virus,
    from_protist,
    from_which_protist,
    from_drug,
    from_which_drug,
    from_food,
    from_which_food,
    Database_source
  )

hmdb_database <- 
  hmdb_database %>%
  mutate(kegg_id = str_trim(kegg_id))

save(hmdb_database, file = "2_data/25_Combine_KEGG_HMDB/hmdb_database.rda")



hmdb_database_new <-
  hmdb_database %>%
  mutate(
    Compound_name        = name,
    Synonyms             = synonyms,
    Formula              = chemical_formula,
    Formula_all          = chemical_formula,
    Monoisotopic_weight  = monisotopic_molecular_weight,
    Average_weight       = average_molecular_weight,
    Compound_description = description,
    SMILES_ID            = smiles,
    SMILES_ID_all        = smiles,
    INCHI_ID             = inchi,
    INCHIKEY_ID          = inchikey,
    KEGG_ID = if_else(str_detect(kegg_id, "C"), kegg_id, NA),
    KEGG_DRUG_ID = if_else(str_detect(kegg_id, "D"), kegg_id, NA),
    HMDB_ID              = accession,
    CAS_ID               = cas_registry_number,
    CAS_ID_all           = cas_registry_number,
    CHEMSPIDER_ID        = chemspider_id,
    DRUGBANK_ID          = drugbank_id,
    FOODB_ID             = foodb_id,
    PUBCHEM_COMPOUND_ID  = pubchem_compound_id,
    PUBCHEM_SUBSTANCE_ID = NA,
    CHEBI_ID             = chebi_id,
    CHEMBL_ID            = NA,
    PDB_CCD_ID           = NA,
    `3DMET_ID`           = NA,
    NIKKAJI_ID           = NA,
    KNAPSACK_ID          = NA,
    LIPIDMAPS_ID         = NA,
    LIPIDBANK_ID         = NA,
    BIOCYC_ID            = biocyc_id,
    BIGG_ID              = bigg_id,
    WIKIPEDIA_ID         = wikipedia_id,
    METLIN_ID            = metlin_id,
    Database_source      = "HMDB"
  ) %>%
  dplyr::select(
    Compound_name,
    Synonyms,
    Formula,
    Formula_all,
    Monoisotopic_weight,
    Average_weight,
    Compound_description,
    HMDB_ID,
    KEGG_ID,
    CAS_ID,
    CAS_ID_all,
    INCHI_ID,
    INCHIKEY_ID,
    SMILES_ID,
    SMILES_ID_all,
    KEGG_DRUG_ID,
    CHEMSPIDER_ID,
    DRUGBANK_ID,
    FOODB_ID,
    PUBCHEM_COMPOUND_ID,
    PUBCHEM_SUBSTANCE_ID,
    CHEBI_ID,
    CHEMBL_ID,
    PDB_CCD_ID,
    `3DMET_ID`,
    `NIKKAJI_ID`,
    `KNAPSACK_ID`,
    `LIPIDMAPS_ID`,
    `LIPIDBANK_ID`,
    `BIOCYC_ID`,
    `BIGG_ID`,
    `WIKIPEDIA_ID`,
    `METLIN_ID`,
    from_human,
    from_which_part,
    from_bacteria,
    from_which_bacteria,
    from_fungi,
    from_which_fungi,
    from_archaea,
    from_which_archaea,
    from_plant,
    from_which_plant,
    from_animal,
    from_which_animal,
    from_environment,
    from_which_environment,
    from_virus,
    from_which_virus,
    from_protist,
    from_which_protist,
    from_drug,
    from_which_drug,
    from_food,
    from_which_food,
    Database_source
  )

  
# 将所有数值列转换为字符串
hmdb_database_new <- 
  hmdb_database_new %>%
  mutate(
    Monoisotopic_weight = as.character(Monoisotopic_weight),
    Average_weight = as.character(Average_weight)
  )

save(kegg_database_new, file = "2_data/25_Combine_KEGG_HMDB/kegg_database_new3.rda")
save(hmdb_database_new, file = "2_data/25_Combine_KEGG_HMDB/hmdb_database_new4.rda")


# # 1. 基于kegg id匹配
# kegg_match_kegg <- 
#   kegg_database_new %>%
#   filter(!is.na(KEGG_ID) & KEGG_ID %in% hmdb_database_new$KEGG_ID)
# 
# kegg_remaining <- 
#   kegg_database_new %>%
#   filter(!(KEGG_ID %in% kegg_match_kegg$KEGG_ID))









# # (a) 花括号去重合并
# merge_braces <- function(x) {
#   x <- x[!is.na(x)]
# 
#   if (length(x) == 0) return(NA)
# 
#   # 2. 以 {} 进行分割 -> 打平为字符向量 -> 去重 -> 重新拼接
#   x %>%
#     str_split("\\{\\}") %>%    # 拆分成列表
#     unlist() %>%          # 转换为字符向量
#     unique() %>%               # 去重
#     paste(collapse = "{}")     # 用 {} 重新拼接成字符串
# }
# # 
# # (b) 根据数据来源优先选择
# #    - 默认优先 "HMDB"
# take_first_source <- function(x, preferred_source) {
#   # 在 summarise 的分组环境中，获取当前分组下所有列
#   all_data <- cur_data()
# 
#   # 取出数据库来源列（确保你的数据中有名为 "Database_source" 的列）
#   src <- all_data$Database_source
# 
#   # 只保留 x 非 NA 的行
#   valid_idx <- !is.na(x)
#   x_valid   <- x[valid_idx]
#   src_valid <- src[valid_idx]
# 
#   # 如果全是 NA，直接返回 NA（真正的 NA，而不是 "NA" 字符串）
#   if (length(x_valid) == 0) {
#     return(NA)
#   }
# 
#   # 找到首选数据库所在行的索引
#   preferred_idx <- which(src_valid == preferred_source)
#   if (length(preferred_idx) > 0) {
#     # 如果有来自首选数据库的值，返回该组中【第一个】非 NA 值
#     return(x_valid[preferred_idx[1]])
#   } else {
#     # 否则，返回该组的第一个非 NA 值
#     return(x_valid[1])
#   }
# }
# 
# 
# 
# merge_from_column <- function(x, col_name) {
# 
#   # (a) 如果是 from_which_ 开头
#   if (str_starts(col_name, "from_which_")) {
#     x_no_na <- na.omit(x)
#     x_no_unknown <- x_no_na[x_no_na != "Unknown"]
#     if (length(x_no_unknown) == 0) {
#       return("Unknown")
#     } else {
#       return(
#         x_no_unknown %>%
#           str_split("\\{\\}") %>%    # 拆分成列表
#           unlist() %>%          # 转换为字符向量
#           unique() %>%               # 去重
#           paste(collapse = "{}")
#       )
#     }
#   }
# 
#   # (b) 否则当作 from_xxx 列
#   #     若任何一条记录为 "Yes" 则返回 "Yes"，否则返回 "Unknown"
#   if (any(x == "Yes", na.rm = TRUE)) {
#     return("Yes")
#   } else {
#     return("Unknown")
#   }
}
# 
# any(str_detect(combined_hmdb_kegg$Compound_name, "\\{\\}"))
# 
# 
# 
# # 将 KEGG_ID 为 NA 和非 NA 的行拆分
# combined_hmdb_kegg_na <- combined_hmdb_kegg %>%
#   filter(is.na(KEGG_ID))
# 
# combined_hmdb_kegg_non_na <- combined_hmdb_kegg %>%
#   filter(!is.na(KEGG_ID))
# 
# 
# # combined_hmdb_kegg_non_na_grouped <-
# #   combined_hmdb_kegg_non_na %>%
# #   group_by(KEGG_ID) %>%
# #   summarize(
# #     # 用 take_first_source
# #     across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight), ~ take_first_source(.x, "HMDB")),
# # 
# #     # 用 merge_braces
# #     across(c(Synonyms, Formula_all, Compound_description, Database_source, SMILES_ID_all), merge_braces),
# # 
# #     # 对所有 _ID 结尾的列用 merge_braces
# #     across(ends_with("_ID"), merge_braces),
# # 
# #     # 对所有 from_ 开头的列用 merge_from_column
# #     across(starts_with("from_"), ~ merge_from_column(.x, cur_column())),
# # 
# #     .groups = "drop"
# #   ) #%>%
# # dplyr::select(all_of(names(combined_hmdb_kegg))) %>%
# # arrange(KEGG_ID)
# 
# 
# combined_hmdb_kegg_non_na_grouped <- 
#   combined_hmdb_kegg_non_na %>%
#   group_by(KEGG_ID) %>%
#   summarize(
#     # 1) Columns that should use take_first_source
#     across(
#       c(Compound_name, Formula, Monoisotopic_weight, Average_weight,
#         SMILES_ID, INCHI_ID, INCHIKEY_ID),
#       ~ take_first_source(.x, "HMDB")
#     ),
#     
#     # 2) Columns that should use merge_braces
#     across(
#       c(Synonyms, Formula_all, Compound_description, Database_source, SMILES_ID_all),
#       merge_braces
#     ),
#     
#     # 3) All other _ID columns should use merge_braces
#     across(
#       ends_with("_ID") & !colnames(.) %in% c("SMILES_ID", "INCHI_ID", "INCHIKEY_ID"),
#       merge_braces
#     ),
#     
#     # 4) Columns starting with "from_" should use merge_from_column
#     across(
#       starts_with("from_"),
#       ~ merge_from_column(.x, cur_column())
#     ),
#     
#     .groups = "drop"
#   )
# 
# combined_hmdb_kegg_final <- bind_rows(
#   combined_hmdb_kegg_non_na_grouped,
#   combined_hmdb_kegg_na
# )
# 
# which(str_detect(combined_hmdb_kegg_final$Monoisotopic_weight, "\\{\\}"))
# 
# colnames(combined_hmdb_kegg_final)
# 
# # choose only one as representation
# combined_hmdb_kegg_final$Monoisotopic_weight <- word(combined_hmdb_kegg_final$Monoisotopic_weight, 1, sep = "\\{\\}")
# combined_hmdb_kegg_final$Average_weight <- word(combined_hmdb_kegg_final$Average_weight, 1, sep = "\\{\\}")
# 
# 
# which(str_detect(kegg_database_new$SMILES_ID, "\\{\\}"))
# 
# 
# 
# combined_hmdb_kegg_test <- 
#   combined_hmdb_kegg %>% 
#   dplyr::select(Compound_name, Formula, Average_weight, KEGG_ID, from_human, from_which_part, from_which_bacteria, Database_source) %>% 
#   head(10)
# 
# #######
# 
# combined_hmdb_kegg$KEGG_ID <- trimws(combined_hmdb_kegg$KEGG_ID)
# 
# 
# combined_hmdb_kegg_na <- combined_hmdb_kegg %>%
#   filter(is.na(KEGG_ID))
# 
# combined_hmdb_kegg_non_na <- combined_hmdb_kegg %>%
#   filter(!is.na(KEGG_ID))
# 
# # 对 KEGG_ID 非 NA 的行进行分组汇总
# combined_hmdb_kegg_non_na_grouped <- combined_hmdb_kegg_non_na %>%
#   group_by(KEGG_ID) %>%
#   summarize(
#     across(!starts_with("from_"), merge_braces),
#     across(starts_with("from_"), ~ merge_from_column(.x, cur_column())),
#     .groups = "drop"
#   )
# 
# # 将处理完的非 NA 数据与 NA 数据再合并
# combined_hmdb_kegg_final <- bind_rows(
#   combined_hmdb_kegg_non_na_grouped,
#   combined_hmdb_kegg_na
# )
# 
# combined_hmdb_kegg_final_ <- 
#   combined_hmdb_kegg_final %>% 
#   group_by(KEGG_ID)
# 
# 
# 
# 
# 




