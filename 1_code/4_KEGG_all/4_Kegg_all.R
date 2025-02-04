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

load("2_data/4_KEGG_all/kegg_compound_database")

kegg_compound_database2 <- vector("list", length = length(kegg_compound_database))

id_systems <-
  kegg_compound_database %>%
  purrr::map(function(x){
    x$DBLINKS
  }) %>%
  unlist() %>%
  stringr::str_split(": ") %>%
  lapply(function(x){
    x[1]
  }) %>%
  unlist() %>%
  unique()

names_systems <-
  kegg_compound_database %>%
  purrr::map(function(x){
    names(x)
  }) %>%
  unlist() %>%
  unique()

template <-
  matrix(NA, nrow = 1, ncol = 26) %>%
  as.data.frame()

colnames(template) <- c("ENTRY", "NAME", "FORMULA", "EXACT_MASS", "MOL_WEIGHT", "REMARK", "REACTION", "PATHWAY_ID", "PATHWAY_NAME", "MODULE_ID", "MODULE_NAME", "ENZYME",
                        "BRITE", paste(id_systems, "ID", sep = "_"), "ATOM", "BOND", "COMMENT")



for(i in 1:length(kegg_compound_database)){
#for(i in 1:12){
  
  cat(i, " ")
  
  x = kegg_compound_database[[i]]
  
  y = template
  
  y$ENTRY <- unname(x$ENTRY)
  
  y$NAME <-
    
    paste(stringr::str_replace(x$NAME, "\\;$", ""), collapse = "{}")
  
  y$FORMULA <- ifelse(is.null(x$FORMULA), NA, x$FORMULA)
  
  y$EXACT_MASS <- ifelse(is.null(x$EXACT_MASS), NA, x$EXACT_MASS)
  
  y$MOL_WEIGHT <- ifelse(is.null(x$MOL_WEIGHT), NA, x$MOL_WEIGHT)

  y$REMARK <- ifelse(is.null(x$REMARK), NA, x$REMARK)
  
  y$REACTION <- ifelse(is.null(x$REACTION), NA, paste(x$REACTION, collapse = "{}"))
  
  y$PATHWAY_ID <- ifelse(is.null(x$PATHWAY), NA, paste(names(x$PATHWAY), collapse = "{}"))
  
  y$PATHWAY_NAME <- ifelse(is.null(x$PATHWAY), NA, paste(x$PATHWAY, collapse = "{}"))
  
  y$MODULE_ID <- ifelse(is.null(x$MODULE), NA, paste(names(x$MODULE), collapse = "{}"))
  
  y$MODULE_NAME <- ifelse(is.null(x$MODULE), NA, paste(x$MODULE, collapse = "{}"))
  
  y$ENZYME <- ifelse(is.null(x$ENZYME), NA, paste(stringr::str_trim(x$ENZYME, side = "both"), collapse = "{}"))
  
  y$BRITE <- ifelse(is.null(x$BRITE), NA, paste(stringr::str_trim(x$BRITE, side = "both"), collapse = "{}"))
  
  
  if(any(stringr::str_detect(x$DBLINKS, "CAS"))){
    
    y$CAS_ID <- stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "PubChem"))){
    
    y$PubChem_ID <- stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "ChEBI"))){
    
    y$ChEBI_ID <-
      
      stringr::str_replace(grep("ChEBI", x$DBLINKS, value = TRUE), "ChEBI: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "ChEMBL"))){
    
    y$ChEMBL_ID <-
      
      stringr::str_replace(grep("ChEMBL", x$DBLINKS, value = TRUE), "ChEMBL: ", "") %>%
      
      stringr::str_trim(side = "both")
    
    
    
  }
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "PDB-CCD"))){
    
    y$`PDB-CCD_ID` <-
      
      stringr::str_replace(grep("PDB-CCD", x$DBLINKS, value = TRUE), "PDB-CCD: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "3DMET"))){
    
    y$`3DMET_ID` <-
      
      stringr::str_replace(grep("3DMET", x$DBLINKS, value = TRUE), "3DMET: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "NIKKAJI"))){
    
    y$NIKKAJI_ID <-
      
      stringr::str_replace(grep("NIKKAJI", x$DBLINKS, value = TRUE), "NIKKAJI: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "KNApSAcK"))){
    
    y$KNApSAcK_ID <-
      
      stringr::str_replace(grep("KNApSAcK", x$DBLINKS, value = TRUE), "KNApSAcK: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "LIPIDMAPS" ))){
    
    y$LIPIDMAPS_ID <-
      
      stringr::str_replace(grep("LIPIDMAPS", x$DBLINKS, value = TRUE), "LIPIDMAPS: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "LipidBank"))){
    
    y$LipidBank_ID <-
      
      stringr::str_replace(grep("LipidBank", x$DBLINKS, value = TRUE), "LipidBank: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }

    
  y$ATOM <- ifelse(is.null(x$ATOM), NA, paste(stringr::str_trim(x$ATOM, side = "both"), collapse = "{}"))  
  
  y$BOND <- ifelse(is.null(x$BOND), NA, paste(stringr::str_trim(x$BOND, side = "both"), collapse = "{}"))

  y$COMMENT <- ifelse(is.null(x$COMMENT), NA, paste(stringr::str_trim(x$COMMENT, side = "both"), collapse = "{}"))
    
  kegg_compound_database2[[i]] <- y
}

kegg_compound_dataframe <- 
  kegg_compound_database2%>%
  do.call(rbind, .) %>% 
  as.data.frame()

kegg_compound_dataframe <- 
  kegg_compound_dataframe %>% 
  dplyr::rename(
    DRUG_ID = REMARK,
    COMPOUND_ID = ENTRY
    )
kegg_compound_dataframe$DRUG_ID <- ifelse(!is.na(kegg_compound_dataframe$DRUG_ID), stringr::str_replace(kegg_compound_dataframe$DRUG_ID, ".*: ", ""), NA)
kegg_compound_dataframe$DRUG_ID <- sapply(kegg_compound_dataframe$DRUG_ID , function(x) {
  if (is.na(x)) {
    return(NA) 
  }
  if (grepl("^Same as:", x)) { 
    sub("^Same as: ?", "", x) 
  } else {
    NA 
  }
})

dir.create("2_data/4_KEGG_all/kegg_compound_dataframe", showWarnings = FALSE)
setwd("2_data/4_KEGG_all/kegg_compound_dataframe")

save(kegg_compound_dataframe, file = "kegg_compound_dataframe.rda")


