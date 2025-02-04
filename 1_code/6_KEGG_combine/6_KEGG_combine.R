library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

### to get KEGG database
library(KEGGgraph)
library(KEGGREST)
library(tidyverse)

load("2_data/3_Metabolite_origin/Metabolite_origin.rda")
load("2_data/4_KEGG_all/kegg_compound_dataframe/kegg_compound_dataframe.rda")
load("2_data/5_DRUG/kegg_drug_dataframe/kegg_drug_dataframe.rda")


Combined_dataframe <- 
  bind_rows(kegg_compound_dataframe, kegg_drug_dataframe) %>% 
  dplyr::select("COMPOUND_ID", "DRUG_ID", "NAME", "FORMULA", "EXACT_MASS", "MOL_WEIGHT", 
         "REACTION", "INTERACTION", "CLASS", "METABOLISM", "PATHWAY_ID", "PATHWAY_NAME", 
         "MODULE_ID", "MODULE_NAME", "TARGET_TARGET", "TARGET_PATHWAY", "ENZYME", "BRITE", 
         "SOURCE", "CAS_ID", "PubChem_ID", "ChEBI_ID", "ChEMBL_ID", "PDB-CCD_ID", "3DMET_ID", 
         "NIKKAJI_ID", "KNApSAcK_ID", "LIPIDMAPS_ID", "LipidBank_ID", "LigandBox_ID", "DrugBank_ID", 
         "ATOM", "BOND", "COMMENT")

Combined_dataframe$compound_id_first = ifelse(Combined_dataframe$COMPOUND_ID == "NA", NA, strsplit(Combined_dataframe$COMPOUND_ID, " ") %>% 
                                                sapply(`[[`, 1))

Compound_combined_result <- 
  Combined_dataframe %>%
  filter(!is.na(compound_id_first)) %>%  
  group_by(compound_id_first) %>%
  summarise(
    
    COMPOUND_ID = ifelse(
      all(is.na(COMPOUND_ID)), NA, 
      paste(unique(unlist(strsplit(na.omit(COMPOUND_ID), " "))), collapse = " ")
    ),
    
    DRUG_ID = ifelse(
      all(is.na(DRUG_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(DRUG_ID), " "))), collapse = " ")
    ),
    
    NAME = ifelse(
      all(is.na(NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    FORMULA = ifelse(
      all(is.na(FORMULA)), NA,
      paste(unique(na.omit(FORMULA)), collapse = "{}")
    ),

    EXACT_MASS = ifelse(
      all(is.na(EXACT_MASS)), NA,
      paste(unique(na.omit(EXACT_MASS)), collapse = "{}")
    ),    
    
    MOL_WEIGHT = ifelse(
      all(is.na(MOL_WEIGHT)), NA,
      paste(unique(na.omit(MOL_WEIGHT)), collapse = "{}")
    ),   
    
    REACTION = ifelse(
      all(is.na(REACTION)), NA,  
      paste(unique(unlist(strsplit(na.omit(REACTION), " "))), collapse = " ")
    ),
    
    INTERACTION = ifelse(
      all(is.na(INTERACTION)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(INTERACTION), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    CLASS = ifelse(
      all(is.na(CLASS)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(CLASS), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    METABOLISM = ifelse(
      all(is.na(METABOLISM)), NA,  
      paste(na.omit(METABOLISM), collapse = "{}")
    ),
    
    PATHWAY_ID = ifelse(
      all(is.na(PATHWAY_ID)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(PATHWAY_ID), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    PATHWAY_NAME = ifelse(
      all(is.na(PATHWAY_NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(PATHWAY_NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    MODULE_ID = ifelse(
      all(is.na(MODULE_ID)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(MODULE_ID), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    MODULE_NAME = ifelse(
      all(is.na(MODULE_NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(MODULE_NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    TARGET_TARGET = ifelse(
      all(is.na(TARGET_TARGET)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(TARGET_TARGET), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    TARGET_PATHWAY = ifelse(
      all(is.na(TARGET_PATHWAY)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(TARGET_PATHWAY), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    ENZYME = ifelse(
      all(is.na(ENZYME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(ENZYME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    NAME = ifelse(
      all(is.na(NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    BRITE = ifelse(
      all(is.na(BRITE)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(BRITE), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    SOURCE = ifelse(
      all(is.na(SOURCE)), NA,  
      paste(na.omit(SOURCE), collapse = "{}")
    ),
    
    CAS_ID = ifelse(
      all(is.na(CAS_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(CAS_ID), " "))), collapse = " ")
    ),
    
    PubChem_ID = ifelse(
      all(is.na(PubChem_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(PubChem_ID), " "))), collapse = " ")
    ),
    
    ChEBI_ID = ifelse(
      all(is.na(ChEBI_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(ChEBI_ID), " "))), collapse = " ")
    ),
    
    ChEMBL_ID = ifelse(
      all(is.na(ChEMBL_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(ChEMBL_ID), " "))), collapse = " ")
    ),
    
    `PDB-CCD_ID` = ifelse(
      all(is.na(`PDB-CCD_ID`)), NA,  
      paste(unique(unlist(strsplit(na.omit(`PDB-CCD_ID`), " "))), collapse = " ")
    ),
    
    `3DMET_ID` = ifelse(
      all(is.na(`3DMET_ID`)), NA,  
      paste(unique(unlist(strsplit(na.omit(`3DMET_ID`), " "))), collapse = " ")
    ),
    
    NIKKAJI_ID = ifelse(
      all(is.na(NIKKAJI_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(NIKKAJI_ID), " "))), collapse = " ")
    ),
    
    KNApSAcK_ID = ifelse(
      all(is.na(KNApSAcK_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(KNApSAcK_ID), " "))), collapse = " ")
    ),
    
    LIPIDMAPS_ID = ifelse(
      all(is.na(LIPIDMAPS_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(LIPIDMAPS_ID), " "))), collapse = " ")
    ),
    
    LipidBank_ID = ifelse(
      all(is.na(LipidBank_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(LipidBank_ID), " "))), collapse = " ")
    ),
    
    LigandBox_ID = ifelse(
      all(is.na(LigandBox_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(LigandBox_ID), " "))), collapse = " ")
    ),
    
    DrugBank_ID = ifelse(
      all(is.na(DrugBank_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(DrugBank_ID), " "))), collapse = " ")
    ),
    
    ATOM = ifelse(
      all(is.na(ATOM)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(ATOM), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    BOND = ifelse(
      all(is.na(BOND)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(BOND), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    COMMENT = ifelse(
      all(is.na(COMMENT)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(COMMENT), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    )
    
  ) %>%
  ungroup() %>%
  bind_rows(Combined_dataframe %>% filter(is.na(compound_id_first)))

Compound_combined_result <- 
  Compound_combined_result %>% 
  dplyr::select(-compound_id_first)

Compound_combined_result$drug_id_first = ifelse(Compound_combined_result$DRUG_ID == "NA", NA, strsplit(Compound_combined_result$DRUG_ID, " ") %>% 
                                                sapply(`[[`, 1))

Final_combined_result <- 
  Compound_combined_result %>%
  filter(!is.na(drug_id_first)) %>%  
  group_by(drug_id_first) %>%
  summarise(
    
    COMPOUND_ID = ifelse(
      all(is.na(COMPOUND_ID)), NA, 
      paste(unique(unlist(strsplit(na.omit(COMPOUND_ID), " "))), collapse = " ")
    ),
    
    DRUG_ID = ifelse(
      all(is.na(DRUG_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(DRUG_ID), " "))), collapse = " ")
    ),
    
    NAME = ifelse(
      all(is.na(NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    FORMULA = ifelse(
      all(is.na(FORMULA)), NA,
      paste(unique(na.omit(FORMULA)), collapse = "{}")
    ),
    
    EXACT_MASS = ifelse(
      all(is.na(EXACT_MASS)), NA,
      paste(unique(na.omit(EXACT_MASS)), collapse = "{}")
    ),    
    
    MOL_WEIGHT = ifelse(
      all(is.na(MOL_WEIGHT)), NA,
      paste(unique(na.omit(MOL_WEIGHT)), collapse = "{}")
    ),   
    
    REACTION = ifelse(
      all(is.na(REACTION)), NA,  
      paste(unique(unlist(strsplit(na.omit(REACTION), " "))), collapse = " ")
    ),
    
    INTERACTION = ifelse(
      all(is.na(INTERACTION)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(INTERACTION), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    CLASS = ifelse(
      all(is.na(CLASS)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(CLASS), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    METABOLISM = ifelse(
      all(is.na(METABOLISM)), NA,  
      paste(na.omit(METABOLISM), collapse = "{}")
    ),
    
    PATHWAY_ID = ifelse(
      all(is.na(PATHWAY_ID)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(PATHWAY_ID), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    PATHWAY_NAME = ifelse(
      all(is.na(PATHWAY_NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(PATHWAY_NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    MODULE_ID = ifelse(
      all(is.na(MODULE_ID)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(MODULE_ID), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    MODULE_NAME = ifelse(
      all(is.na(MODULE_NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(MODULE_NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    TARGET_TARGET = ifelse(
      all(is.na(TARGET_TARGET)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(TARGET_TARGET), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    TARGET_PATHWAY = ifelse(
      all(is.na(TARGET_PATHWAY)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(TARGET_PATHWAY), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    ENZYME = ifelse(
      all(is.na(ENZYME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(ENZYME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    NAME = ifelse(
      all(is.na(NAME)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(NAME), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    BRITE = ifelse(
      all(is.na(BRITE)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(BRITE), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    SOURCE = ifelse(
      all(is.na(SOURCE)), NA,  
      paste(na.omit(SOURCE), collapse = "{}")
    ),
    
    CAS_ID = ifelse(
      all(is.na(CAS_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(CAS_ID), " "))), collapse = " ")
    ),
    
    PubChem_ID = ifelse(
      all(is.na(PubChem_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(PubChem_ID), " "))), collapse = " ")
    ),
    
    ChEBI_ID = ifelse(
      all(is.na(ChEBI_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(ChEBI_ID), " "))), collapse = " ")
    ),
    
    ChEMBL_ID = ifelse(
      all(is.na(ChEMBL_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(ChEMBL_ID), " "))), collapse = " ")
    ),
    
    `PDB-CCD_ID` = ifelse(
      all(is.na(`PDB-CCD_ID`)), NA,  
      paste(unique(unlist(strsplit(na.omit(`PDB-CCD_ID`), " "))), collapse = " ")
    ),
    
    `3DMET_ID` = ifelse(
      all(is.na(`3DMET_ID`)), NA,  
      paste(unique(unlist(strsplit(na.omit(`3DMET_ID`), " "))), collapse = " ")
    ),
    
    NIKKAJI_ID = ifelse(
      all(is.na(NIKKAJI_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(NIKKAJI_ID), " "))), collapse = " ")
    ),
    
    KNApSAcK_ID = ifelse(
      all(is.na(KNApSAcK_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(KNApSAcK_ID), " "))), collapse = " ")
    ),
    
    LIPIDMAPS_ID = ifelse(
      all(is.na(LIPIDMAPS_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(LIPIDMAPS_ID), " "))), collapse = " ")
    ),
    
    LipidBank_ID = ifelse(
      all(is.na(LipidBank_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(LipidBank_ID), " "))), collapse = " ")
    ),
    
    LigandBox_ID = ifelse(
      all(is.na(LigandBox_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(LigandBox_ID), " "))), collapse = " ")
    ),
    
    DrugBank_ID = ifelse(
      all(is.na(DrugBank_ID)), NA,  
      paste(unique(unlist(strsplit(na.omit(DrugBank_ID), " "))), collapse = " ")
    ),
    
    ATOM = ifelse(
      all(is.na(ATOM)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(ATOM), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    BOND = ifelse(
      all(is.na(BOND)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(BOND), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    ),
    
    COMMENT = ifelse(
      all(is.na(COMMENT)), NA, 
      paste(unique(unlist(strsplit(paste(na.omit(COMMENT), collapse = "{}"), "\\{\\}"))), collapse = "{}")
    )
    
  ) %>%
  ungroup() %>%
  bind_rows(Compound_combined_result %>% filter(is.na(drug_id_first)))

Final_combined_result <- 
  Final_combined_result %>% 
  dplyr::select(-drug_id_first)

dir.create("2_data/6_Combine/final_combine_result", showWarnings = FALSE)
setwd("2_data/6_Combine/final_combine_result")

save(Final_combined_result, file = "Final_combined_result.rda")

load("2_data/6_Combine/final_combine_result/Final_combined_result.rda")


Metabolite_origin <- 
  Metabolite_origin %>% 
  dplyr::select(-Compound.name)

Metabolite_origin_compound <- 
  Metabolite_origin %>% 
  filter(grepl("C", KEGG.ID)) 

Metabolite_origin_compound <- 
  Metabolite_origin_compound %>%
  dplyr::rename(COMPOUND_ID = KEGG.ID)

Metabolite_origin_drug <- 
  Metabolite_origin %>% 
  filter(!grepl("C", KEGG.ID)) 

Metabolite_origin_drug <- 
  Metabolite_origin_drug %>%
  dplyr::rename(DRUG_ID = KEGG.ID)

join_with_compound <- 
  Final_combined_result %>%
  left_join(Metabolite_origin_compound, by = "COMPOUND_ID")

join_with_compound_missing <- 
  join_with_compound %>% filter(is.na(from_human))

join_with_drug <- 
  join_with_compound_missing %>%
  dplyr::select("COMPOUND_ID", "DRUG_ID", "NAME", "FORMULA", "EXACT_MASS", "MOL_WEIGHT", 
                "REACTION", "INTERACTION", "CLASS", "METABOLISM", "PATHWAY_ID", "PATHWAY_NAME", 
                "MODULE_ID", "MODULE_NAME", "TARGET_TARGET", "TARGET_PATHWAY", "ENZYME", "BRITE", 
                "SOURCE", "CAS_ID", "PubChem_ID", "ChEBI_ID", "ChEMBL_ID", "PDB-CCD_ID", "3DMET_ID", 
                "NIKKAJI_ID", "KNApSAcK_ID", "LIPIDMAPS_ID", "LipidBank_ID", "LigandBox_ID", "DrugBank_ID", 
                "ATOM", "BOND", "COMMENT") %>%  
  left_join(Metabolite_origin_drug, by = "DRUG_ID")

Real_final_combined <- 
  bind_rows(
    join_with_compound %>% filter(!is.na(from_human)),
    join_with_drug 
  )

Real_final_combined <- 
  Real_final_combined %>% 
  arrange(COMPOUND_ID, DRUG_ID)

Real_final_combined <- 
  Real_final_combined %>%
  mutate(across(c(from_human, from_which_part, from_bacteria, from_which_bacteria, from_fungi, from_which_fungi, from_archaea, 
                  from_which_archaea, from_plant, from_which_plant, from_animal, from_which_animal, from_environment, 
                  from_which_environment, from_virus, from_which_virus, from_protist, from_which_protist, from_drug, from_which_drug, 
                  from_food, from_which_food), ~ replace_na(.x, "Unknown")))

dir.create("2_data/6_Combine/Real_final_combined", showWarnings = FALSE)
setwd("2_data/6_Combine/Real_final_combined")

save(Real_final_combined, file = "Real_final_combined.rda")
