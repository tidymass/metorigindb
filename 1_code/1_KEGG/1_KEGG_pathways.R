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

#####extract metabolites from KEGG pathways

files <- dir("2_data/1_KEGG/pathway/organism_pathways/")

metabolite_database <- vector("list", length(files))

for(i in 1:length(files)){
  cat(i, " ")
  load(file.path(paste0("2_data/1_KEGG/pathway/organism_pathways/",files[i]), "pathway.rda"))
  metabolite_info <-
    do.call(rbind, pathway@compound_list) %>% 
    as.data.frame() %>% 
    dplyr::distinct(KEGG.ID, .keep_all = TRUE)
  metabolite_info$organism <- files[i]
  metabolite_info
  metabolite_database[[i]] <- metabolite_info
}


metabolite_database <-
metabolite_database %>% 
  do.call(rbind, .) %>% 
  as.data.frame()


metabolite_database %>% 
  dplyr::arrange(KEGG.ID)


unique_KEGG_ID <- 
  unique(metabolite_database$KEGG.ID)


metabolite_database_final <-
purrr::map(unique_KEGG_ID, function(x){
  cat(x, " ")
  temp_organism <-
  metabolite_database %>% 
    dplyr::filter(KEGG.ID == x) %>% 
    pull(organism) %>% 
    paste0(collapse = "{}")
  
  temp_metablite_info <-
  metabolite_database %>% 
    dplyr::filter(KEGG.ID == x) %>% 
    dplyr::distinct(KEGG.ID, .keep_all = TRUE)
  
  temp_metablite_info$organism <- temp_organism
  temp_metablite_info
})


metabolite_database_final <-
  metabolite_database_final %>% 
  do.call(rbind, .) %>%
  as.data.frame()


dir.create("3_data_analysis/1_KEGG/metabolites", showWarnings = FALSE)
setwd("3_data_analysis/1_KEGG/metabolites")

save(metabolite_info, file = "metabolite_info.rda")
