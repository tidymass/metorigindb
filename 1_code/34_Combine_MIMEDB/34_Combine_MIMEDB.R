library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

source("1_code/merge_functions.R")

load("2_data/33_Combine_MODELSEED/MODELSEED_Combined_Final3")
load("2_data/20_MIMEDB/MIMEDB_database/MIMEDB_database3.rda")

MIMEDB_database[] <- lapply(MIMEDB_database, function(x) {
  if (is.character(x)) {
    x[x == "NULL"] <- NA  # 替换 "NULL" 为 NA
  }
  return(x)  # 保持原数据类型
})

#
MIMEDB_database_new <- 
  MIMEDB_database %>%
  mutate(
    Compound_name = name,
    Synonyms = name,
    Formula = moldb_formula,
    Formula_all = moldb_formula,
    
    Monoisotopic_weight = moldb_mono_mass,
    Average_weight = moldb_average_mass,
    
    Compound_description = description,
    HMDB_ID = hmdb_id,
    HMDB_ID_all = hmdb_id,
    KEGG_ID = NA,
    KEGG_ID_all = NA,
    
    CAS_ID = cas,
    CAS_ID_all = cas,
    
    INCHI_ID = moldb_inchi,
    INCHI_ID_all = moldb_inchi,
    INCHIKEY_ID = moldb_inchikey,
    INCHIKEY_ID_all = moldb_inchikey,
    
    SMILES_ID = moldb_smiles,
    SMILES_ID_all = moldb_smiles,
    
    KEGG_DRUG_ID = NA,
    CHEMSPIDER_ID = NA,
    
    DRUGBANK_ID = NA,
    FOODB_ID = NA,
    PUBCHEM_COMPOUND_ID = NA,
    PUBCHEM_SUBSTANCE_ID = NA,
    
    CHEBI_ID = NA,
    CHEBI_ID_all = NA,
    CHEMBL_ID = NA,
    PDB_CCD_ID = NA,
    `3DMET_ID` = NA,
    `NIKKAJI_ID` = NA,
    `KNAPSACK_ID` = NA,
    `LIPIDMAPS_ID` = NA,
    `LIPIDBANK_ID` = NA,
    `BIOCYC_ID` = NA,
    `BIGG_ID` = NA,
    BIGG_IDENTIFIER_ID = NA,
    `WIKIPEDIA_ID` = NA,
    `METLIN_ID` = NA,
    T3DB_ID = NA, 
    REACTOME_ID = NA,
    MODELSEED_ID = NA,
    MIMEDB_ID = mime_id,
    Database_source = "MIMEDB"
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
    HMDB_ID_all,
    KEGG_ID,
    KEGG_ID_all,
    CAS_ID,
    CAS_ID_all,
    INCHI_ID,
    INCHI_ID_all,
    INCHIKEY_ID,
    INCHIKEY_ID_all,
    SMILES_ID,
    SMILES_ID_all,
    KEGG_DRUG_ID,
    CHEMSPIDER_ID,
    DRUGBANK_ID,
    FOODB_ID,
    PUBCHEM_COMPOUND_ID,
    PUBCHEM_SUBSTANCE_ID,
    CHEBI_ID,
    CHEBI_ID_all,
    CHEMBL_ID,
    PDB_CCD_ID,
    `3DMET_ID`,
    `NIKKAJI_ID`,
    `KNAPSACK_ID`,
    `LIPIDMAPS_ID`,
    `LIPIDBANK_ID`,
    `BIOCYC_ID`,
    `BIGG_ID`,
    BIGG_IDENTIFIER_ID,
    `WIKIPEDIA_ID`,
    `METLIN_ID`,
    T3DB_ID,
    REACTOME_ID,
    MODELSEED_ID,
    MIMEDB_ID,
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

MODELSEED_Combined_Final <- 
  MODELSEED_Combined_Final %>% 
  mutate(
    MIMEDB_ID = NA
  )

MIMEDB_database_new <- 
  MIMEDB_database_new  %>%
  mutate(
    Monoisotopic_weight = as.character(Monoisotopic_weight),
    Average_weight = as.character(Average_weight),
    PUBCHEM_COMPOUND_ID = as.character(PUBCHEM_COMPOUND_ID),
    CHEBI_ID = as.character(CHEBI_ID),
    CHEBI_ID_all = as.character(CHEBI_ID_all)
  )

database_new <- MIMEDB_database_new

database_new <- 
  database_new %>% 
  mutate(
    from_bacteria = "Yes"
  )

# 1. 基于hmdb匹配
# 注意是否是有多的id
hmdb_match <- 
  database_new %>%
  filter(!is.na(HMDB_ID) & HMDB_ID %in% MODELSEED_Combined_Final$HMDB_ID)

remaining1 <- 
  database_new %>%
  filter(!(HMDB_ID %in% hmdb_match$HMDB_ID))

# change here
df1 <- MODELSEED_Combined_Final
df2 <- hmdb_match

# change here
match_col <- "HMDB_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  
  cat(i, " ")
  
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    
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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, CHEBI_ID, KEGG_ID, CAS_ID, INCHI_ID, INCHIKEY_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, MIMEDB_ID, MODELSEED_ID, REACTOME_ID, T3DB_ID, Formula_all, Compound_description, INCHI_ID_all, INCHIKEY_ID_all, Database_source, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID_all, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combine_hmdb_id <- bind_rows(res_list)

# 2. 基于inchikey匹配
# 注意是否是有多的id
inchikey_match <- 
  remaining1 %>%
  filter(!is.na(INCHIKEY_ID) & INCHIKEY_ID %in% combine_hmdb_id$INCHIKEY_ID)

remaining2 <- 
  remaining1 %>%
  filter(!(INCHIKEY_ID %in% inchikey_match$INCHIKEY_ID))

# change here
df1 <- combine_hmdb_id
df2 <- inchikey_match

# change here
match_col <- "INCHIKEY_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  
  cat(i, " ")
  
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    
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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, CHEBI_ID, HMDB_ID, CAS_ID, INCHI_ID, KEGG_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, MIMEDB_ID, MODELSEED_ID, REACTOME_ID, T3DB_ID, Formula_all, Compound_description, INCHI_ID_all, INCHIKEY_ID_all, Database_source, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID_all, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combine_inchikey_id <- bind_rows(res_list)


# 3. 基于inchi匹配
# 注意是否是有多的id
inchi_match <- 
  remaining2 %>%
  filter(!is.na(INCHI_ID) & INCHI_ID %in% combine_inchikey_id$INCHI_ID)

remaining3 <- 
  remaining2 %>%
  filter(!(INCHI_ID %in% inchi_match$INCHI_ID))

# change here
df1 <- combine_inchikey_id
df2 <- inchi_match

# change here
match_col <- "INCHI_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  
  cat(i, " ")
  
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    
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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, CHEBI_ID, HMDB_ID, CAS_ID, INCHIKEY_ID, KEGG_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, MIMEDB_ID, MODELSEED_ID, REACTOME_ID, T3DB_ID, Formula_all, Compound_description, INCHI_ID_all, INCHIKEY_ID_all, Database_source, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID_all, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combine_inchi_id <- bind_rows(res_list)



# 4. 基于cas匹配
# 注意是否是有多的id
cas_match <- 
  remaining3 %>%
  filter(!is.na(CAS_ID) & CAS_ID %in% combine_inchi_id$CAS_ID)

remaining4 <- 
  remaining3 %>%
  filter(!(CAS_ID %in% cas_match$CAS_ID))

# change here
df1 <- combine_inchi_id
df2 <- cas_match

# change here
match_col <- "CAS_ID"

res_list <- vector("list", length = nrow(df1))

## 循环
for(i in seq_len(nrow(df1))) {
  
  cat(i, " ")
  
  if (is.na(df1[[match_col]][i]) || !(df1[[match_col]][i] %in% df2[[match_col]])) {
    res_list[[i]] <- df1[i, ]
  } else {
    
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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, CHEBI_ID, HMDB_ID, INCHI_ID, INCHIKEY_ID, KEGG_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, MIMEDB_ID, MODELSEED_ID, REACTOME_ID, T3DB_ID, Formula_all, Compound_description, INCHI_ID_all, INCHIKEY_ID_all, Database_source, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID_all, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combine_cas_id <- bind_rows(res_list)

MIMEDB_Combined_Final <- rbind(combine_cas_id, remaining4)
save(MIMEDB_Combined_Final, file = "2_data/34_Combine_MIMEDB/MIMEDB_Combined_Final_test4.rda")
