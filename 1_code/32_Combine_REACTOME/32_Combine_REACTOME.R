library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

source("1_code/merge_functions.R")

load("2_data/31_Combine_T3DB/T3DB_Combined_Final3.rda")
load("2_data/24_Reactome/Reactome_database/ChEBI2Reactome_database.rda")



#
REACTOME_database_new <- 
  ChEBI2Reactome_database_cleaned %>%
  mutate(
    Compound_name = word(Reactome_Physical_Entity_Name, 1, sep = "\\{\\}"),
    Synonyms = Reactome_Physical_Entity_Name,
    Formula = NA,
    Formula_all = NA,
    
    Monoisotopic_weight = NA,
    Average_weight = NA,
    
    Compound_description = NA,
    HMDB_ID = NA,
    HMDB_ID_all = NA,
    KEGG_ID = NA,
    KEGG_ID_all = NA,
    
    CAS_ID = NA,
    CAS_ID_all = NA,
    
    INCHI_ID = NA,
    INCHI_ID_all = NA,
    INCHIKEY_ID = NA,
    INCHIKEY_ID_all = NA,
    
    SMILES_ID = NA,
    SMILES_ID_all = NA,
    
    KEGG_DRUG_ID = NA,
    CHEMSPIDER_ID = NA,
    
    DRUGBANK_ID = NA,
    FOODB_ID = NA,
    PUBCHEM_COMPOUND_ID = NA,
    PUBCHEM_SUBSTANCE_ID = NA,
    
    CHEBI_ID = word(ChEBI_identifier, 1, sep = "\\{\\}"),
    CHEBI_ID_all = ChEBI_identifier,
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
    REACTOME_ID = Reactome_Physical_Entity_Stable_Identifier,
    Database_source = "REACTOME"
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

T3DB_Combined_Final <- 
  T3DB_Combined_Final %>% 
  mutate(
    CHEBI_ID_all = CHEBI_ID,
    CHEBI_ID = word(CHEBI_ID, 1 ,sep = "\\{\\}"),
    REACTOME_ID = NA
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

REACTOME_database_new <- 
  REACTOME_database_new  %>%
  mutate(
    Monoisotopic_weight = as.character(Monoisotopic_weight),
    Average_weight = as.character(Average_weight),
    PUBCHEM_COMPOUND_ID = as.character(PUBCHEM_COMPOUND_ID),
    CHEBI_ID = as.character(CHEBI_ID),
    CHEBI_ID_all = as.character(CHEBI_ID_all)
  )

# 1. 基于chebi匹配
# 注意是否是有多的id
chebi_match <- 
  REACTOME_database_new %>%
  filter(!is.na(CHEBI_ID) & CHEBI_ID %in% T3DB_Combined_Final$CHEBI_ID)

remaining1 <- 
  REACTOME_database_new %>%
  filter(!(CHEBI_ID %in% chebi_match$CHEBI_ID))

# change here
df1 <- T3DB_Combined_Final
df2 <- chebi_match

# change here
match_col <- "CHEBI_ID"

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
        across(c(Compound_name, Formula, Monoisotopic_weight, Average_weight, KEGG_ID, HMDB_ID, CAS_ID, INCHI_ID, INCHIKEY_ID, SMILES_ID), ~ take_first_source(.x, "HMDB")),
        
        # change here
        across(c(Synonyms, REACTOME_ID, T3DB_ID, Formula_all, Compound_description, INCHI_ID_all, INCHIKEY_ID_all, Database_source, KEGG_ID_all, BIGG_IDENTIFIER_ID, CAS_ID_all, HMDB_ID_all, SMILES_ID_all, KEGG_DRUG_ID, CHEMSPIDER_ID, DRUGBANK_ID, FOODB_ID, PUBCHEM_COMPOUND_ID, PUBCHEM_SUBSTANCE_ID, CHEBI_ID_all, CHEMBL_ID, PDB_CCD_ID, '3DMET_ID', NIKKAJI_ID, KNAPSACK_ID, LIPIDMAPS_ID, LIPIDBANK_ID, BIOCYC_ID, BIGG_ID, WIKIPEDIA_ID, METLIN_ID), merge_braces),
        
        #
        across(starts_with("from_"), ~ merge_from_column(.x, cur_column()))
      )
    
    # merged_row 现在是一行，作为合并后的结果
    res_list[[i]] <- merged_row
  }
}

combine_chebi_id <- bind_rows(res_list)

REACTOME_Combined_Final <- rbind(combine_chebi_id, remaining1)
save(REACTOME_Combined_Final, file = "2_data/32_Combine_REACTOME/REACTOME_Combined_Final3.rda")




