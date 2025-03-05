###

no_source()

setwd(masstools::get_project_wd())

source("1_code/3_utils.R")

rm(list = ls())

setwd('2_data/DRUGBANK/')



library(dbparser)

library(dplyr)

library(ggplot2)

library(XML)



library(MetaDBparse)



# # Description <- NULL

# # base.loc <- file.path("database", "drugbank_source")

# # if (!dir.exists(base.loc)) {

# #   dir.create(base.loc)

# # }

# #

# # # zip.file <- file.path(base.loc, "drugbank.zip")

# #

# # # utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

# # input <- file.path("full database.xml")

# # header <- readLines(input, n = 10)

# # hasInfo <-

# #   grep(

# #     x = header,

# #     pattern = "version",

# #     value = TRUE,

# #     perl = TRUE

# #   )[2]

# # version <-

# #   stringr::str_match(string = hasInfo, pattern = "version=\"(.*)\" exported")[, 2]

# # theurl <-

# #   RCurl::getURL("https://go.drugbank.com/stats",

# #                 .opts = list(ssl.verifypeer = FALSE))

# # tables <- XML::readHTMLTable(theurl, header = FALSE)

# # stats <- data.table::as.data.table(tables[[1]])

# # colnames(stats) <- c("Description", "Count")

# # n <-

# #   as.numeric(as.character(gsub(

# #     x = stats[Description == "Total Number of Drugs"]$Count,

# #     pattern = ",",

# #     replacement = ""

# #   )))

# # envir <- environment()

# # envir$db.formatted <-

# #   data.frame(

# #     compoundname = rep(NA, n),

# #     baseformula = rep(NA, n),

# #     identifier = rep(NA, n),

# #     structure = rep(NA, n),

# #     charge = rep(NA, n),

# #     description = rep(NA, n)

# #   )

# # envir$pb <- pbapply::startpb(min = 0, max = n)

# # envir$idx <- 0

# # metabolite <- function(currNode, currEnvir = envir) {

# #   if (currEnvir$idx %% 10 == 0) {

# #     pbapply::setpb(currEnvir$pb, currEnvir$idx)

# #   }

# #   currEnvir$idx <- currEnvir$idx + 1

# #   properties <- currNode[["calculated-properties"]]

# #   if (is.null(properties)) {

# #     properties <- currNode[["experimental-properties"]]

# #   }

# #   proplist <- XML::xmlToList(properties)

# #   if (length(proplist) == 0) {

# #     return(NULL)

# #   }

# #   which.form <- which(sapply(proplist, function(x) {

# #     if ("kind" %in% names(x)) {

# #       res <- x[["kind"]] == "Molecular Formula"

# #     }

# #     else {

# #       res <- FALSE

# #     }

# #     res

# #   }))

# #   which.struc <- which(sapply(proplist, function(x) {

# #     if ("kind" %in% names(x)) {

# #       res <- x[["kind"]] == "SMILES"

# #     }

# #     else {

# #       res <- FALSE

# #     }

# #     res

# #   }))

# #   which.charge <- which(sapply(proplist, function(x) {

# #     if ("kind" %in% names(x)) {

# #       res <- x[["kind"]] == "Physiological Charge"

# #     }

# #     else {

# #       res <- FALSE

# #     }

# #     res

# #   }))

# #   if (length(which.form) == 0 & length(which.struc) == 0) {

# #     return(NULL)

# #   }

# #   currEnvir$db.formatted[currEnvir$idx, "compoundname"] <-

# #     XML::xmlValue(currNode[["name"]])

# #   currEnvir$db.formatted[currEnvir$idx, "identifier"] <-

# #     XML::xmlValue(currNode[["drugbank-id"]])

# #   currEnvir$db.formatted[currEnvir$idx, "baseformula"] <-

# #     proplist[[which.form]][["value"]]

# #   currEnvir$db.formatted[currEnvir$idx, "structure"] <-

# #     if (length(which.struc) > 0) {

# #       proplist[[which.struc]][["value"]]

# #     }

# #   else {

# #     ""

# #   }

# #   currEnvir$db.formatted[currEnvir$idx, "description"] <-

# #     XML::xmlValue(currNode[["description"]])

# #   currEnvir$db.formatted[currEnvir$idx, "charge"] <-

# #     if (length(which.charge) > 0) {

# #       proplist[[which.charge]][["value"]]

# #     }

# #   else {

# #     0

# #   }

# # }

# #

# # res <-

# #   XML::xmlEventParse(

# #     file = input,

# #     branches = list(drug = metabolite, `drugbank-metabolite-id-value` = print)

# #   )

# #

# # envir$db.formatted <- envir$db.formatted[-1,]

# # drugbank2 = envir$db.formatted

# # save(drugbank2, file = "drugbank2")

# load("drugbank2")

#

# ## parse data from XML and save it to memory

# # read_drugbank_xml_db("full database.xml")

# #

# # ## load drugs data

# # drugs <- drugs()

# #

# # ## load drug groups data

# # drug_groups <- drug_groups()

# #

# # ## load drug targets actions data

# # drug_targets_actions <- targets_actions()

# # save(drugs, file = "drugs")

# #

# # load("drugs")

# #

# #

# # ###only remain small molecules

# # general_information =

# #   drugs$general_information %>%

# #   dplyr::filter(type == "small molecule") %>%

# #   dplyr::rename(Lab.ID = primary_key) %>%

# #   dplyr::select(Lab.ID, everything()) %>%

# #   dplyr::mutate(drugbank.ID = Lab.ID) %>%

# #   dplyr::rename(Compound.name = name,

# #                 CAS.ID = cas_number,

# #                 mz = monoisotopic_mass,

# #                 UNII.ID = unii,

# #                 Description = description,

# #                 Create_date = created,

# #                 Updated_date = updated,

# #                 Average.mass = average_mass,

# #                 State = state,

# #   ) %>%

# #   dplyr::select(-c(synthesis_reference, fda_label, msds))

# #

# # drug_classification <-

# #   drugs$drug_classification %>%

# #   dplyr::filter(drugbank_id %in% general_information$Lab.ID) %>%

# #   dplyr::rename(Lab.ID = drugbank_id) %>%

# #   dplyr::select(Lab.ID, everything()) %>%

# #   dplyr::select(-c(description, direct_parent, alternative_parents, substituents)) %>%

# #   dplyr::rename(Kingdom = kingdom,

# #                 Super_class = superclass,

# #                 Class = class,

# #                 Sub_class = class)

# #

# # library(plyr)

# #

# # synonyms =

# #   drugs$synonyms %>%

# #   dplyr::filter(`drugbank-id` %in% general_information$Lab.ID) %>%

# #   plyr::dlply(.variables = .(`drugbank-id`)) %>%

# #   purrr::map(function(x) {

# #     if (nrow(x) == 1) {

# #       return(x)

# #     } else{

# #       x$synonym = paste(x$synonym, collapse = "{}")

# #       x$language = paste(x$language, collapse = "{}")

# #       x$coder = paste(x$coder, collapse = "{}")

# #       x[1, , drop = FALSE]

# #     }

# #   }) %>%

# #   do.call(rbind, .) %>%

# #   as.data.frame() %>%

# #   dplyr::rename(Lab.ID = `drugbank-id`) %>%

# #   dplyr::select(Lab.ID, everything()) %>%

# #   dplyr::select(-c(coder, language))

# #

# # drugbank <-

# #   general_information %>%

# #   dplyr::left_join(drug_classification, by = c("Lab.ID")) %>%

# #   dplyr::left_join(synonyms, by = c("Lab.ID"))

# #

# # drugbank$RT = NA

# # drugbank$HMDB.ID = NA

# # drugbank$KEGG.ID = NA

# # drugbank$mz.pos = NA

# # drugbank$mz.neg = NA

# # drugbank$Submitter = "DRUGBANK"

# #

# # drugbank$From_drug = TRUE

# #

# # save(drugbank, file = "drugbank")

# load("drugbank")

#

# drugbank2 <-

#   drugbank2 %>%

#   dplyr::select(baseformula, identifier) %>%

#   dplyr::filter(!is.na(identifier)) %>%

#   dplyr::rename(Lab.ID = identifier,

#                 Formula = baseformula)

#

# drugbank <-

#   drugbank %>%

#   dplyr::left_join(drugbank2, by = "Lab.ID")

#

# drugbank =

#   drugbank %>%

#   dplyr::select(

#     Lab.ID,

#     mz,

#     RT,

#     CAS.ID,

#     HMDB.ID,

#     KEGG.ID,

#     Formula,

#     mz.pos,

#     mz.neg,

#     Submitter,

#     everything()

#   ) %>%

#   dplyr::filter(!is.na(Formula) & !is.na(mz))

#

# drugbank$mz = as.numeric(drugbank$mz)

#

# drugbank <- as.data.frame(drugbank)

#

# drugbank <-

# drugbank %>%

#   dplyr::rename(Synonym = synonym) %>%

#   dplyr::mutate(From_drug = "Yes")

#

#

# load(here::here("2_data/HMDB/MS1/hmdb_ms1.rda"))

#

# idx <-

# match(drugbank$CAS.ID, hmdb_ms1@spectra.info$CAS.ID)

#

# colnames(drugbank)

# colnames(hmdb_ms1@spectra.info)

#

# drugbank$HMDB.ID <- hmdb_ms1@spectra.info$HMDB.ID[idx]

# drugbank$KEGG.ID <- hmdb_ms1@spectra.info$KEGG.ID[idx]

#

# drugbank$From_human <- NA

# drugbank$From_human <- hmdb_ms1@spectra.info$From_human[idx]

#

# drugbank$Average.mass <- as.numeric(drugbank$Average.mass)

#

# drugbank$Super_class %>% head()

# hmdb_ms1@spectra.info$Super_class %>% head()

#

# drugbank$SMILES.ID <- NA

# drugbank$SMILES.ID <- hmdb_ms1@spectra.info$SMILES.ID[idx]

#

# drugbank$INCHI.ID <- NA

# drugbank$INCHI.ID <- hmdb_ms1@spectra.info$INCHI.ID[idx]

#

# drugbank$INCHIKEY.ID <- NA

# drugbank$INCHIKEY.ID <- hmdb_ms1@spectra.info$INCHIKEY.ID[idx]

#

# drugbank[which(drugbank == "", arr.ind = TRUE)] <- NA

#

# openxlsx::write.xlsx(drugbank, file = "drugbank.xlsx",

#                      asTable = TRUE, overwrite = TRUE)



drugbank_ms1 <-
  
  construct_database(
    
    path = ".",
    
    version = "2022-04-21",
    
    metabolite.info.name = "drugbank.xlsx",
    
    source = "drugbank",
    
    link = "https://go.drugbank.com/",
    
    creater = "Xiaotao Shen",
    
    email = "shenxt@stanford.edu",
    
    rt = FALSE,
    
    threads = 3
    
  )



load("../HMDB/MS1/hmdb_ms1.rda")

load("../KEGG/kegg_ms1.rda")



source(here::here("1_code/3_utils.R"))



intersect(colnames(drugbank_ms1@spectra.info),
          
          colnames(hmdb_ms1@spectra.info))



setdiff(colnames(hmdb_ms1@spectra.info),
        
        colnames(drugbank_ms1@spectra.info))



drugbank_ms1 <-
  
  update_metid_database_info(
    
    database = drugbank_ms1,
    
    ref_database = hmdb_ms1,
    
    by = c(
      
      "CAS.ID",
      
      "HMDB.ID",
      
      "KEGG.ID",
      
      "Formula",
      
      "mz.pos",
      
      "mz.neg",
      
      "Compound.name",
      
      "SMILES.ID",
      
      "INCHI.ID",
      
      "INCHIKEY.ID"
      
    ),
    
    combine_columns = c(
      
      "CAS.ID",
      
      "HMDB.ID",
      
      "KEGG.ID",
      
      "Formula",
      
      "mz.pos",
      
      "mz.neg",
      
      "Compound.name",
      
      "SMILES.ID",
      
      "INCHI.ID",
      
      "INCHIKEY.ID",
      
      "Kingdom",
      
      "Super_class",
      
      "Sub_class"
      
    ),
    
    new_columns = c(
      
      "status",
      
      "Synonyms",
      
      "IUPAC_name",
      
      "Traditional_IUPAC_name",
      
      "Class",
      
      "Biospecimen_locations",
      
      "Cellular_locations",
      
      "Tissue_locations",
      
      "CHEMSPIDER.ID",
      
      "DRUGBANK.ID",
      
      "FOODB.ID",
      
      "PUBCHEM.ID",
      
      "CHEBI.ID",
      
      "BIOCYC.ID",
      
      "BIGG.ID",
      
      "WIKIPEDIA.ID",
      
      "METLIN.ID"
      
    )
    
  )



intersect(colnames(drugbank_ms1@spectra.info),
          
          colnames(kegg_ms1@spectra.info))



setdiff(colnames(kegg_ms1@spectra.info),
        
        colnames(drugbank_ms1@spectra.info))



source(here::here("1_code/3_utils.R"))



drugbank_ms1 <-
  
  update_metid_database_info(
    
    database = drugbank_ms1,
    
    ref_database = kegg_ms1,
    
    by = c(
      
      "Compound.name",
      
      "CAS.ID",
      
      "HMDB.ID",
      
      "KEGG.ID",
      
      "INCHI.ID",
      
      "FOODB.ID",
      
      "WIKIPEDIA.ID",
      
      "DRUGBANK.ID",
      
      "PUBCHEM.ID",
      
      "CHEBI.ID",
      
      "INCHI.ID",
      
      "INCHIKEY.ID"
      
    ),
    
    combine_columns = c(
      
      "Compound.name",
      
      "CAS.ID",
      
      "HMDB.ID",
      
      "KEGG.ID",
      
      "INCHI.ID",
      
      "FOODB.ID",
      
      "WIKIPEDIA.ID",
      
      "DRUGBANK.ID",
      
      "PUBCHEM.ID",
      
      "CHEBI.ID",
      
      "INCHI.ID",
      
      "INCHIKEY.ID"
      
    ),
    
    new_columns = c(
      
      "CHEMBL.ID",
      
      "LIPIDMAPS.ID",
      
      "LIPIDBANK.ID",
      
      "KEGG_DRUG.ID"
      
    )
    
  )



load(here::here("2_data/source_system/source_system.rda"))



library(tidyverse)

library(tidyselect)

library(metid)



drugbank_ms1 <-
  
  update_metid_database_source_system(
    
    database = drugbank_ms1,
    
    source_system = source_system,
    
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    
    prefer = "database"
    
  )
install.packages("openxlsx")
library(openxlsx)


save(drugbank_ms1, file = "drugbank_ms1.rda")

setwd("2_data/9_DRUGBANK")
drugbank_database <- 
  read.xlsx("drugbank.xlsx")

drugbank_database <- 
  drugbank_database %>% 
  mutate(
    from_human = ifelse(is.na(From_human), "Unknown", From_human),
    from_which_part = "Unknown",
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
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = "Unknown",
    from_which_protist = "Unknown",
    from_drug = ifelse(is.na(From_drug), "Unknown", From_drug),
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown",
    
  )

dir.create("2_data/9_DRUGBANK/DRUGBANK_database", showWarnings = FALSE)
setwd("DRUGBANK_database")

save(drugbank_database, file = "drugbank_database.rda")



