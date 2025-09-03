####HMDB database: https://hmdb.ca/

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



###download the metabolite data here: https://hmdb.ca/downloads

###unzip the file and put it in the folder "2_data/HMDB/MS1/"



hmdb <- read_xml("hmdb_metabolites.xml")



hmdb <- as_list(hmdb)



x = hmdb$hmdb[[1]]



hmdb_metabolite <- vector(mode = "list", length = length(hmdb$hmdb))



for(i in 1:length(hmdb$hmdb)){
  
  cat(i, " ")
  
  x <- hmdb$hmdb[[i]]
  
  biospecimen_locations <-
    
    paste(stringr::str_trim(unname(unlist(x$biological_properties$biospecimen_locations))),
          
          collapse = "{}")
  
  cellular_locations <-
    
    paste(stringr::str_trim(unlist(x$biological_properties$cellular_locations)),
          
          collapse = "{}")
  
  tissue_locations <-
    
    paste(stringr::str_trim(unlist(x$biological_properties$tissue_locations)),
          
          collapse = "{}")
  
  # pathway <-
  
  # lapply(x$biological_properties$pathways, function(y){
  
  #   pathway_name <- unlist(unlist(y$name))
  
  #   pathway_smpdb_id <- unlist(unlist(y$smpdb_id))
  
  #   pathway_kegg_id <- unlist(unlist(y$kegg_map_id))
  
  #   pathway_name <-
  
  #     ifelse(is.null(pathway_name), NA, pathway_name)
  
  #   pathway_smpdb_id <-
  
  #     ifelse(is.null(pathway_smpdb_id), NA, pathway_smpdb_id)
  
  #   pathway_kegg_id <-
  
  #     ifelse(is.null(pathway_kegg_id), NA, pathway_kegg_id)
  
  #   c(pathway_name = pathway_name,
  
  #     pathway_smpdb_id = pathway_smpdb_id,
  
  #     pathway_kegg_id = pathway_kegg_id)
  
  # }) %>%
  
  #   do.call(rbind, .) %>%
  
  #   as.data.frame()
  
  #
  
  # pathway$pathway_name <-
  
  #   paste(pathway$pathway_name, collapse = "{}")
  
  # pathway$pathway_smpdb_id <-
  
  #   paste(pathway$pathway_smpdb_id[!is.na(pathway$pathway_smpdb_id)], collapse = "{}")
  
  # pathway$pathway_kegg_id <-
  
  #   paste(pathway$pathway_kegg_id[!is.na(pathway$pathway_kegg_id)], collapse = "{}")
  
  #
  
  # pathway <-
  
  # pathway %>%
  
  #   dplyr::distinct(pathway_name, .keep_all = TRUE)
  
  #
  
  source <-
    
    tryCatch(
      
      lapply(x$ontology[[2]]$descendants[[1]]$descendants, function(z){
        
        unname(unlist(z$term[[1]]))
        
      }) %>%
        
        unlist() %>%
        
        unname() %>%
        
        paste0(collapse = "{}"),
      
      error = function(e) NA
      
    )
  
  
  
  
  
  average_molecular_weight = ifelse(is.null(unlist(x$average_molecular_weight)),
                                    
                                    NA,
                                    
                                    as.numeric(unlist(x$average_molecular_weight)))
  
  monisotopic_molecular_weight = ifelse(is.null(unlist(x$monisotopic_molecular_weight)),
                                        
                                        NA,
                                        
                                        as.numeric(unlist(x$monisotopic_molecular_weight)))
  
  
  
  secondary_accessions = paste(stringr::str_trim(unlist(x$secondary_accessions)),
                               
                               collapse = "{}")
  
  
  
  synonyms = paste(stringr::str_trim(unname(unlist(x$synonyms))), collapse = "{}")
  
  
  
  hmdb_metabolite[[i]] <-
    
    data.frame(version = ifelse(is.null(unlist(x$version)), NA, unlist(x$version)),
               
               creation_date = ifelse(is.null(unlist(x$creation_date)), NA, unlist(x$creation_date)),
               
               update_date = ifelse(is.null(unlist(x$update_date)), NA, unlist(x$update_date)),
               
               accession = ifelse(is.null(unlist(x$accession)), NA, unlist(x$accession)),
               
               status = ifelse(is.null(unlist(x$status)), NA, unlist(x$status)),
               
               secondary_accessions = secondary_accessions,
               
               name = unlist(x$name),
               
               description = ifelse(is.null(unlist(x$description)), NA, unlist(x$description)),
               
               synonyms = synonyms,
               
               chemical_formula = ifelse(is.null(unlist(x$chemical_formula)), NA, unlist(x$chemical_formula)),
               
               average_molecular_weight = average_molecular_weight,
               
               monisotopic_molecular_weight = monisotopic_molecular_weight,
               
               iupac_name = ifelse(is.null(unlist(x$iupac_name)), NA, unlist(x$iupac_name)),
               
               traditional_iupac = ifelse(is.null(unlist(x$traditional_iupac)), NA, unlist(x$traditional_iupac)),
               
               cas_registry_number = ifelse(is.null(unlist(x$cas_registry_number)), NA, unlist(x$cas_registry_number)),
               
               smiles = ifelse(is.null(unlist(x$smiles)), NA, unlist(x$smiles)),
               
               inchi = ifelse(is.null(unlist(x$inchi)), NA, unlist(x$inchi)),
               
               inchikey = ifelse(is.null(unlist(x$inchikey)), NA, unlist(x$inchikey)),
               
               kingdom = ifelse(is.null(unlist(x$taxonomy$kingdom)),
                                
                                NA, unlist(x$taxonomy$kingdom)),
               
               super_class = ifelse(is.null(unlist(x$taxonomy$super_class)),
                                    
                                    NA, unlist(x$taxonomy$super_class)),
               
               class = ifelse(is.null(unlist(x$taxonomy$class)),
                              
                              NA, unlist(x$taxonomy$class)),
               
               sub_class = ifelse(is.null(unlist(x$taxonomy$sub_class)),
                                  
                                  NA, unlist(x$taxonomy$sub_class)),
               
               state = ifelse(is.null(unlist(x$state)), NA, unlist(x$state)),
               
               source = source,
               
               biospecimen_locations = biospecimen_locations,
               
               cellular_locations = cellular_locations,
               
               tissue_locations = tissue_locations,
               
               chemspider_id = ifelse(is.null(unlist(x$chemspider_id)), NA, unlist(x$chemspider_id)),
               
               drugbank_id = ifelse(is.null(unlist(x$drugbank_id)), NA, unlist(x$drugbank_id)),
               
               foodb_id = ifelse(is.null(unlist(x$foodb_id)), NA, unlist(x$foodb_id)),
               
               pubchem_compound_id = ifelse(is.null(unlist(x$pubchem_compound_id)), NA, unlist(x$pubchem_compound_id)),
               
               chebi_id = ifelse(is.null(unlist(x$chebi_id)), NA, unlist(x$chebi_id)),
               
               kegg_id = ifelse(is.null(unlist(x$kegg_id)), NA, unlist(x$kegg_id)),
               
               biocyc_id = ifelse(is.null(unlist(x$biocyc_id)), NA, unlist(x$biocyc_id)),
               
               bigg_id = ifelse(is.null(unlist(x$bigg_id)), NA, unlist(x$bigg_id)),
               
               wikipedia_id = ifelse(is.null(unlist(x$wikipedia_id)), NA, unlist(x$wikipedia_id)),
               
               metlin_id = ifelse(is.null(unlist(x$metlin_id)), NA, unlist(x$metlin_id))
               
    )
  
}



hmdb_metabolite <-
  
  hmdb_metabolite %>%
  
  dplyr::bind_rows() %>%
  
  as.data.frame()


setwd("2_data/7_HMDB")
save(hmdb_metabolite, file = "hmdb_metabolite.rda", compress = "xz")

load("hmdb_metabolite.rda")

hmdb_metabolite$biospecimen_locations[hmdb_metabolite$biospecimen_locations == ""] <- NA
hmdb_metabolite$tissue_locations[hmdb_metabolite$tissue_locations == ""] <- NA


hmdb_metabolite <- 
  hmdb_metabolite %>% 
  rowwise() %>% 
  mutate(
    from_human = "Yes",
    from_which_part = ifelse(
      is.na(biospecimen_locations) & is.na(tissue_locations), "Unknown", 
      paste(unique(na.omit(unlist(strsplit(paste(na.omit(c(biospecimen_locations, tissue_locations))), "\\{\\}")))), 
            collapse = "{}")
    ),
    from_bacteria = "Unknown",
    from_which_bacteria = "Unknown",
    from_fungi = "Unknown",
    from_which_fungi = "Unknown",
    from_archaea = "Unknown",
    from_which_archaea = "Unknown",
    from_plant = "Unknown",
    from_which_plant = "Unknown",
    from_animal = "Unknown",
    from_which_animal = "Unknown",
    from_environment = "Unknown",
    from_which_environment = "Unkonwn",
    from_virus = "Unkonwn",
    from_which_virus = "Unkonwn",
    from_protist = "Unkonwn",
    from_which_protist = "Unkonwn",
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",
    
  )

#
setwd("2_data/7_HMDB/HMDB_files")

#
drug_orgin <- read.csv("Drug", sep = ",", header = TRUE, stringsAsFactors = FALSE)
food_orgin <- read.csv("Food", sep = ",", header = TRUE, stringsAsFactors = FALSE)
environment_orgin <- read.csv("Environment", sep = ",", header = TRUE, stringsAsFactors = FALSE)
Mic_orgin <- read.csv("Microbial", sep = ",", header = TRUE, stringsAsFactors = FALSE)
plant_orgin <- read.csv("Plant", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# 1. Drug
idx <- which(hmdb_metabolite$accession %in% drug_orgin$HMDB_ID)

# 修改对应行的列值
hmdb_metabolite$from_drug[idx]        <- "Yes"
hmdb_metabolite$from_human[idx]       <- "Unknown"
hmdb_metabolite$from_which_part[idx]  <- "Unknown"


# 2. Food
idx <- which(hmdb_metabolite$accession %in% food_orgin$HMDB_ID)

# 修改对应行的列值
hmdb_metabolite$from_food[idx]        <- "Yes"
hmdb_metabolite$from_human[idx]       <- "Unknown"
hmdb_metabolite$from_which_part[idx]  <- "Unknown"


# 3. Environment
idx <- which(hmdb_metabolite$accession %in% environment_orgin$HMDB_ID)

# 修改对应行的列值
hmdb_metabolite$from_environment[idx]        <- "Yes"
hmdb_metabolite$from_human[idx]       <- "Unknown"
hmdb_metabolite$from_which_part[idx]  <- "Unknown"



# 4. Mic
idx <- which(hmdb_metabolite$accession %in% Mic_orgin$HMDB_ID)

# 修改对应行的列值
hmdb_metabolite$from_bacteria[idx]        <- "Yes"
hmdb_metabolite$from_human[idx]       <- "Unknown"
hmdb_metabolite$from_which_part[idx]  <- "Unknown"



# 5. Plant
idx <- which(hmdb_metabolite$accession %in% plant_orgin$HMDB_ID)

# 修改对应行的列值
hmdb_metabolite$from_plant[idx]        <- "Yes"
hmdb_metabolite$from_human[idx]       <- "Unknown"
hmdb_metabolite$from_which_part[idx]  <- "Unknown"

dir.create("2_data/7_HMDB/HMDB_Origin", showWarnings = FALSE)
setwd("2_data/7_HMDB/HMDB_Origin")

save(hmdb_metabolite, file = "~/tidymass/metabolite_database/2_data/7_HMDB/HMDB_Origin/HMDB_Origin_final.rda")

