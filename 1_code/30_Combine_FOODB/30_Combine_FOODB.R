library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

source("1_code/merge_functions.R")

load("2_data/29_Combine_CHEBI/CHEBI_Combined_Final2.rda")
load("2_data/8_FOODB/FOODB_database/FOODB_database.rda")

test <- 
  CHEBI_Combined_Final %>% 
  dplyr::filter(!is.na(HMDB_ID) & HMDB_ID %in% foodb_database$HMDB.ID)

sum(test$from_human == "Yes")


foodb_database <- 
  foodb_database %>%
  mutate(HMDB.ID = ifelse(
    str_detect(HMDB.ID, "^HMDB\\d+$"),  # 确保是 HMDB+数字格式
    str_replace(HMDB.ID, "HMDB(\\d+)", function(x) {
      paste0("HMDB", str_pad(str_extract(x, "\\d+"), width = 7, side = "left", pad = "0"))
    }),
    HMDB.ID  # 对于NA或者其他非符合条件的值，保持不变
  ))


# FOODB UNKNOWN
foodb_database <- 
  foodb_database %>%
  mutate(across(c(from_which_environment, from_virus, from_which_virus, from_protist, from_which_protist), ~ ifelse(. == "Unkonwn", "Unknown", .)))


FOODB_database_new <- 
  foodb_database %>%
  mutate(
    Compound_name = Compound.name,
    Synonyms = Synonyms,
    Formula = Formula,
    Formula_all = Formula,
    Monoisotopic_weight = mz,
    Average_weight = Average.mass,
    Compound_description = Description,
    HMDB_ID = HMDB.ID,
    HMDB_ID_all = HMDB.ID,
    KEGG_ID = KEGG.ID,
    KEGG_ID_all = KEGG.ID,
    CAS_ID = CAS.ID,
    CAS_ID_all = CAS.ID,
    INCHI_ID = INCHI.ID,
    INCHI_ID_all = INCHI.ID,
    INCHIKEY_ID = INCHIKEY.ID,
    INCHIKEY_ID_all = INCHIKEY.ID,
    SMILES_ID = SMILES.ID,
    SMILES_ID_all = SMILES.ID,
    KEGG_DRUG_ID = NA,
    CHEMSPIDER_ID = CHEMSPIDER.ID,
    DRUGBANK_ID = NA,
    FOODB_ID = FOODB.ID,
    PUBCHEM_COMPOUND_ID = PUBCHEM.ID,
    PUBCHEM_SUBSTANCE_ID = NA,
    CHEBI_ID = CHEBI.ID,
    CHEMBL_ID = NA,
    PDB_CCD_ID = NA,
    `3DMET_ID` = NA,
    `NIKKAJI_ID` = NA,
    `KNAPSACK_ID` = NA,
    `LIPIDMAPS_ID` = NA,
    `LIPIDBANK_ID` = NA,
    `BIOCYC_ID` = BIOCYC.ID,
    `BIGG_ID` = NA,
    BIGG_IDENTIFIER_ID = NA,
    `WIKIPEDIA_ID` = WIKIPEDIA.ID,
    `METLIN_ID` = NA,
    Database_source = "FOODB"
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


database_new <- 
  FOODB_database_new %>%
  mutate(
    Monoisotopic_weight = as.character(Monoisotopic_weight),
    Average_weight = as.character(Average_weight),
    PUBCHEM_COMPOUND_ID = as.character(PUBCHEM_COMPOUND_ID),
    CHEBI_ID = as.character(CHEBI_ID),
  )



# 1. 基于kegg id匹配 change here
kegg_match <- 
  database_new %>%
  filter(!is.na(KEGG_ID) & KEGG_ID %in% CHEBI_Combined_Final$KEGG_ID)

remaining1 <- 
  database_new %>%
  filter(!(KEGG_ID %in% kegg_match$KEGG_ID))

# change here
df1 <- CHEBI_Combined_Final
df2 <- kegg_match

# change here
match_col <- "KEGG_ID"

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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, HMDB_ID, CAS_ID, INCHI_ID, INCHIKEY_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, Formula_all, Compound_description, INCHI_ID_all, INCHIKEY_ID_all, Database_source, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column())),
        
        from_human = "Unknown"
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_kegg_id <- bind_rows(res_list)


# 2. 基于HMDB id匹配
hmdb_match <- 
  remaining1 %>%
  filter(!is.na(HMDB_ID) & HMDB_ID %in% combined_kegg_id$HMDB_ID)

remaining2 <- 
  remaining1 %>%
  filter(!(HMDB_ID %in% hmdb_match$HMDB_ID))

df1 <- combined_kegg_id
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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, KEGG_ID, CAS_ID, INCHI_ID, INCHIKEY_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, Formula_all, Compound_description, INCHI_ID_all, INCHIKEY_ID_all, Database_source, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column())),
        
        from_human = "Unknown"
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_hmdb_id <- bind_rows(res_list)

# 3. 基于cas id匹配
cas_match <- 
  remaining2 %>%
  filter(!is.na(CAS_ID) & CAS_ID %in% combined_hmdb_id$CAS_ID)

remaining3 <- 
  remaining2 %>%
  filter(!(CAS_ID %in% cas_match$CAS_ID))

df1 <- combined_hmdb_id
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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, KEGG_ID, INCHIKEY_ID, INCHI_ID, HMDB_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, Formula_all, Compound_description, Database_source, INCHI_ID_all, INCHIKEY_ID_all, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column())),
        
        from_human = "Unknown"
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_cas_id <- bind_rows(res_list)


# 4. 基于inchi id匹配
inchi_match <- 
  remaining3 %>%
  filter(!is.na(INCHI_ID) & INCHI_ID %in% combined_cas_id$INCHI_ID)

remaining4 <- 
  remaining3 %>%
  filter(!(INCHI_ID %in% inchi_match$INCHI_ID))

df1 <- combined_cas_id
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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, KEGG_ID, INCHIKEY_ID, CAS_ID, HMDB_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, Formula_all, Compound_description, Database_source, INCHI_ID_all, INCHIKEY_ID_all, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column())),
        
        from_human = "Unknown"
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_inchi_id <- bind_rows(res_list)

FOODB_Combined_Final <- rbind(combined_inchi_id, remaining4)
save(FOODB_Combined_Final, file = "2_data/30_Combine_FOODB/FOODB_Combined_Final4.rda")

inchikey_match <- 
  remaining4 %>%
  filter(!is.na(INCHIKEY_ID) & INCHIKEY_ID %in% combined_inchi_id$INCHIKEY_ID)

