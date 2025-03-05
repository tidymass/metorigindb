library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(tidyverse)

pubtax_data <- read.csv("2_data/23_Pubchem/Pubchem_compounds/pubtax_data_origin.csv", stringsAsFactors = FALSE)


filter_eukaryota <- 
  pubtax_data %>% 
  dplyr::filter(superkingdom == "Eukaryota")

unique(filter_eukaryota$phylum)

pubtax_database <- 
  pubtax_data %>% 
  mutate(
    from_human = ifelse(species == "Homo sapiens", "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = ifelse(superkingdom == "Bacteria", "Yes", "Unknown"),
    from_which_bacteria = ifelse(from_bacteria == "Yes", taxname, "Unknown"),
    from_fungi = ifelse(
      phylum %in% c("Ascomycota", "Basidiomycota", "Mucoromycota", "Microsporidia", 
                    "Blastocladiomycota", "Zoopagomycota"), 
      "Yes", 
      "Unknown"
    ),
    from_which_fungi = ifelse(from_fungi == "Yes", taxname, "Unknown"),
    from_archaea = ifelse(superkingdom == "Archaea", "Yes", "Unknown"),
    from_which_archaea = ifelse(from_archaea == "Yes", taxname, "Unknown"),
    from_plant = ifelse(
      phylum %in% c("Streptophyta", "Chlorophyta", "Rhodophyta"), 
      "Yes", 
      "Unknown"
    ),
    from_which_plant = ifelse(from_plant == "Yes", taxname, "Unknown"),
    from_animal = ifelse(
      phylum %in% c("Arthropoda", "Chordata", "Nematoda", "Mollusca", "Platyhelminthes", 
                    "Porifera", "Cnidaria", "Echinodermata", "Annelida", "Bryozoa", 
                    "Hemichordata", "Nemertea", "Brachiopoda") & from_human != "Yes", 
      "Yes", 
      "Unknown"
    ),
    from_which_animal = ifelse(from_animal == "Yes", taxname, "Unknown"),
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(
      phylum %in% c("Euglenozoa", "Bacillariophyta", "Oomycota", "Evosea", "Haptophyta", 
                    "Ciliophora", "Discosea", "Parabasalia", "Foraminifera", "Apicomplexa", 
                    "Endomyxa"), 
      "Yes", 
      "Unknown"
    ),
    from_which_protist = ifelse(from_protist == "Yes", taxname, "Unknown"),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown"
  )

taxonomy <- read.csv("2_data/23_Pubchem/Pubchem_compounds/taxonomy.csv", stringsAsFactors = FALSE)
taxonomy <- 
  taxonomy %>% 
  dplyr::rename(cid = X.cid)

# 假设 df1 和 df2 都包含 "cid" 列，df2 包含 from_ 相关列
pubtax_database_mergy <- 
  pubtax_database %>%
  left_join(taxonomy, by = "cid")

pubtax_database_mergy_filter <- 
  pubtax_database_mergy %>%
  select("cid", "cmpdname.y", "cmpdsynonym", "mw", 
         "mf", "polararea", "complexity", "xlogp", 
         "heavycnt", "hbonddonor", "hbondacc", "rotbonds", 
         "inchi", "smiles", "inchikey", "iupacname", 
         "exactmass", "monoisotopicmass", "charge", "covalentunitcnt", 
         "isotopeatomcnt", "totalatomstereocnt", "definedatomstereocnt", "undefinedatomstereocnt",
         "totalbondstereocnt", "definedbondstereocnt", "undefinedbondstereocnt", "pclidcnt", 
         "gpidcnt", "gpfamilycnt", "meshheadings", "annothits", 
         "annothitcnt", "aids", "cidcdate", "sidsrcname", 
         "depcatg", "annotation", "from_human", 
         "from_which_part", "from_bacteria", "from_which_bacteria", "from_fungi", 
         "from_which_fungi", "from_archaea", "from_which_archaea", "from_plant", 
         "from_which_plant", "from_animal", "from_which_animal", "from_environment", 
         "from_which_environment", "from_virus", "from_which_virus", "from_protist", 
         "from_which_protist", "from_drug", "from_which_drug", "from_food", 
         "from_which_food")

pubtax_database_mergy_filter <- 
  pubtax_database_mergy_filter %>% 
  dplyr::rename(cmpdname = cmpdname.y)

setwd("2_data/23_Pubchem/Pubchem_database")

save(pubtax_database_mergy_filter, file = "pubchem_tax_database.rda")














