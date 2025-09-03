library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

source("1_code/merge_functions.R")

load("2_data/26_Combine_KEGG_HMDB2/combined_data/KEGG_HMDB_Combined_Final2.rda")
load("2_data/10_BIGG/BIGG_dataset/BIGG_database2.rda")


BIGG_database_new <- 
  BIGG_database %>%
  mutate(
    Compound_name         = Compound.name,
    Synonyms              = Synonyms,
    Formula               = Formula,
    Formula_all           = Formula,
    Monoisotopic_weight   = mz,
    Average_weight        = NA,
    Compound_description  = NA,
    SMILES_ID             = NA,
    SMILES_ID_all         = NA,
    INCHI_ID              = NA,
    INCHI_ID_all          = NA,
    INCHIKEY_ID           = INCHIKEY.ID,
    INCHIKEY_ID_all       = INCHIKEY.ID,
    CAS_ID                = NA,
    CAS_ID_all            = NA,
    KEGG_ID               = word(KEGG.ID, 1, sep = "\\{\\}"),
    KEGG_ID_all           = KEGG.ID,
    KEGG_DRUG_ID          = KEGG_DRUG.ID,
    HMDB_ID               = word(HMDB.ID, 1, sep = "\\{\\}"),
    HMDB_ID_all           = HMDB.ID,
    CHEMSPIDER_ID         = NA,
    DRUGBANK_ID           = NA,  
    FOODB_ID              = NA,
    PUBCHEM_COMPOUND_ID   = NA,
    PUBCHEM_SUBSTANCE_ID  = NA,
    CHEBI_ID              = CHEBI.ID,
    CHEMBL_ID             = NA,
    PDB_CCD_ID            = NA,
    `3DMET_ID`            = NA,
    `NIKKAJI_ID`          = NA,
    KNAPSACK_ID           = NA,
    LIPIDMAPS_ID          = LIPIDMAPS.ID,
    LIPIDBANK_ID          = NA,
    BIOCYC_ID             = BIOCYC.ID,
    BIGG_ID               = NA,
    BIGG_IDENTIFIER_ID    = BIGG.ID,
    WIKIPEDIA_ID          = NA,
    METLIN_ID             = NA,
    Database_source       = "BIGG"
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




KEGG_HMDB_Combined <- 
  KEGG_HMDB_Combined_Final %>%
  mutate(
    HMDB_ID_all = HMDB_ID,
    KEGG_ID_all = KEGG_ID,
    KEGG_ID = word(KEGG_ID, 1, sep = "\\{\\}"),
    INCHI_ID_all = INCHI_ID,
    INCHIKEY_ID_all = INCHIKEY_ID,
    BIGG_IDENTIFIER_ID = NA
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

BIGG_database_new <- 
  BIGG_database_new %>%
  mutate(
    Monoisotopic_weight = as.character(Monoisotopic_weight),
    Average_weight = as.character(Average_weight)
  )

# 1. 基于kegg id匹配
bigg_match_kegg <- 
  BIGG_database_new %>%
  filter(!is.na(KEGG_ID) & KEGG_ID %in% KEGG_HMDB_Combined$KEGG_ID)

bigg_remaining <- 
  BIGG_database_new %>%
  filter(!(KEGG_ID %in% bigg_match_kegg$KEGG_ID))

df1 <- KEGG_HMDB_Combined
df2 <- bigg_match_kegg

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
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_kegg_id <- bind_rows(res_list)



# 2. 基于HMDB id匹配
bigg_match_hmdb <- 
  bigg_remaining %>%
  filter(!is.na(HMDB_ID) & HMDB_ID %in% combined_kegg_id$HMDB_ID)

bigg_remaining2 <- 
  bigg_remaining %>%
  filter(!(HMDB_ID %in% bigg_match_hmdb$HMDB_ID))

df1 <- combined_kegg_id
df2 <- bigg_match_hmdb

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
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_hmdb_id <- bind_rows(res_list)



# 3. 基于INCHIKEY id匹配
bigg_match_inchikey <- 
  bigg_remaining2 %>%
  filter(!is.na(INCHIKEY_ID) & INCHIKEY_ID %in% combined_hmdb_id$INCHIKEY_ID)

bigg_remaining3 <- 
  bigg_remaining2 %>%
  filter(!(INCHIKEY_ID %in% bigg_match_inchikey$INCHIKEY_ID))

df1 <- combined_hmdb_id
df2 <- bigg_match_inchikey

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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, KEGG_ID, CAS_ID, INCHI_ID, HMDB_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, Formula_all, Compound_description, Database_source, INCHI_ID_all, INCHIKEY_ID_all, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combined_inchikey_id <- bind_rows(res_list)

BIGG_Combined_Final <- rbind(combined_inchikey_id, bigg_remaining3)
save(BIGG_Combined_Final, file = "2_data/27_Combine_BIGG/combined_data/BIGG_Combined_Final3.rda")

BIGG_Combined_Final <- 
  BIGG_Combined_Final %>% 
  mutate(
    from_which_environment = "Unknown"
  )













