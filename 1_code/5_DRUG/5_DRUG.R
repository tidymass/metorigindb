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

load("2_data/5_DRUG/kegg_drug_database")

kegg_drug_database2 <- vector("list", length = length(kegg_drug_database))

id_systems <-
  kegg_drug_database %>%
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
  kegg_drug_database %>%
  purrr::map(function(x){
    names(x)
  }) %>%
  unlist() %>%
  unique()

template <-
  matrix(NA, nrow = 1, ncol = 24) %>%
  as.data.frame()

get_list_depth <- function(x) {
  if (!is.list(x)) {
    return(0)  
  } else {
    return(1 + max(sapply(x, get_list_depth, simplify = TRUE), na.rm = TRUE))
  }
}

colnames(template) <- c("ENTRY", "NAME", "FORMULA", "EXACT_MASS", "MOL_WEIGHT", "CLASS", "METABOLISM", "REMARK", "INTERACTION", "TARGET_TARGET", "TARGET_PATHWAY", "SOURCE",
                        "BRITE", paste(id_systems, "ID", sep = "_"), "ATOM", "BOND", "COMMENT")


for(i in 1:length(kegg_drug_database)){
#for(i in 1:100){
  
  cat(i, " ")
  
  x = kegg_drug_database[[i]]
  
  y = template
  
  y$ENTRY <- unname(x$ENTRY)
  
  y$NAME <-
    
    paste(stringr::str_replace(x$NAME, "\\;$", ""), collapse = "{}")
  
  y$FORMULA <- ifelse(is.null(x$FORMULA), NA, x$FORMULA)
  
  y$EXACT_MASS <- ifelse(is.null(x$EXACT_MASS), NA, x$EXACT_MASS)
  
  y$MOL_WEIGHT <- ifelse(is.null(x$MOL_WEIGHT), NA, x$MOL_WEIGHT)
  
  y$CLASS <- ifelse(is.null(x$CLASS), NA, paste(stringr::str_trim(x$CLASS, side = "both"), collapse = "{}"))

  y$METABOLISM <- ifelse(is.null(x$METABOLISM), NA, x$METABOLISM)

  y$REMARK <- ifelse(is.null(x$REMARK), NA, paste(x$REMARK, collapse = "{}"))  
    
  y$INTERACTION <- 
    if (length(x$INTERACTION) != 0) {
      if (length(x$INTERACTION) == 1){
        ifelse(x$INTERACTION == "", NA, paste(x$INTERACTION, collapse = "{}"))
      } else{
        paste(x$INTERACTION, collapse = "{}")
      }
    } else {
      NA
    }
    
    
      
  y$TARGET_TARGET <- 
    if (get_list_depth(x$TARGET) == 1) {
      ifelse(is.null(x$TARGET$TARGET), NA, paste(x$TARGET$TARGET, collapse = "{}"))
    } else {
      ifelse(is.null(x$TARGET), NA, paste(x$TARGET, collapse = "{}"))
    }
  
  y$TARGET_PATHWAY <- 
    if (get_list_depth(x$TARGET) == 1) {
      ifelse(is.null(x$TARGET$PATHWAY), NA, paste(x$TARGET$PATHWAY, collapse = "{}"))
    } else {
      NA
    }
  
  y$SOURCE <- ifelse(is.null(x$SOURCE), NA, paste(x$SOURCE, collapse = "{}"))  
  
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
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "LigandBox"))){
    
    y$LigandBox_ID <-
      
      stringr::str_replace(grep("LigandBox", x$DBLINKS, value = TRUE), "LigandBox: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "NIKKAJI"))){
    
    y$NIKKAJI_ID <-
      
      stringr::str_replace(grep("NIKKAJI", x$DBLINKS, value = TRUE), "NIKKAJI: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  
  
  
  
  
  
  if(any(stringr::str_detect(x$DBLINKS, "DrugBank"))){
    
    y$DrugBank_ID <-
      
      stringr::str_replace(grep("DrugBank", x$DBLINKS, value = TRUE), "DrugBank: ", "") %>%
      
      stringr::str_trim(side = "both")
    
  }
  
  y$ATOM <- ifelse(is.null(x$ATOM), NA, paste(stringr::str_trim(x$ATOM, side = "both"), collapse = "{}"))  
  
  y$BOND <- ifelse(is.null(x$BOND), NA, paste(stringr::str_trim(x$BOND, side = "both"), collapse = "{}"))
  
  y$COMMENT <- ifelse(is.null(x$COMMENT), NA, paste(stringr::str_trim(x$COMMENT, side = "both"), collapse = "{}")) 
  
  kegg_drug_database2[[i]] <- y
}

kegg_drug_dataframe <- 
  kegg_drug_database2%>%
  do.call(rbind, .) %>% 
  as.data.frame()

kegg_drug_dataframe <- 
  kegg_drug_dataframe %>% 
  dplyr::rename(
    DRUG_ID = ENTRY,
    COMPOUND_ID = REMARK
  )

kegg_drug_dataframe$COMPOUND_ID <- sapply(kegg_drug_dataframe$COMPOUND_ID, function(x) {
  if (is.na(x)) {
    return(NA) 
  }
  first_part <- strsplit(x, "\\{\\}")[[1]][1]
  if (grepl("^Same as:", first_part)) { 
    sub("^Same as: ?", "", first_part)
  } else {
    NA 
  }
})

dir.create("2_data/5_DRUG/kegg_drug_dataframe", showWarnings = FALSE)
setwd("2_data/5_DRUG/kegg_drug_dataframe")

save(kegg_drug_dataframe, file = "kegg_drug_dataframe.rda")

