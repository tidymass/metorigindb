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

load("2_data/10_BIGG/metabolite/bigg_ms1.RData")
BIGG_data <- 
  bigg_ms1@spectra.info

BIGG_data <- 
  BIGG_data %>% 
  select(-c(from_human, from_which_part, from_bacteria, from_which_bacteria,
            from_fungi, from_which_fungi, from_eukaryote, from_which_eukaryote,
            from_archaea, from_which_archaea, from_plant, from_which_plant,
            from_animal, from_which_animal, from_environment, from_which_environment,
            from_virus, from_which_virus, from_drug, from_which_drug,
            from_food, from_which_food))

unique_organisms <- 
  BIGG_data %>% 
  separate_rows(organism, sep = "\\{\\}") %>% 
  distinct(organism) %>% 
  pull(organism)

# 定义各分类下的物种名称向量
human_list <- c(
  "Homo sapiens" 
)

bacteria_list <- c(
  "Escherichia coli str. K-12 substr. MG1655",
  "Staphylococcus aureus subsp. aureus N315",
  "Escherichia coli APEC O1",
  "Yersinia pestis CO92",
  "Shigella flexneri 2a str. 301",
  "Pseudomonas putida KT2440",
  "Helicobacter pylori 26695",
  "Mycobacterium tuberculosis H37Rv",
  "Escherichia coli BW2952",
  "Escherichia coli BL21(DE3)",
  "Escherichia coli CFT073",
  "Escherichia coli O127:H6 str. E2348/69",
  "Escherichia coli 042",
  "Escherichia coli 55989",
  "Escherichia coli ABU 83972",
  "Escherichia coli B str. REL606",
  "Escherichia coli 'BL21-Gold(DE3)pLysS AG'",
  "Escherichia coli DH1",
  "Escherichia coli str. K-12 substr. DH10B",
  "Escherichia coli O139:H28 str. E24377A",
  "Escherichia coli ED1a",
  "Escherichia coli O157:H7 str. EC4115",
  "Escherichia coli HS",
  "Escherichia coli IAI1",
  "Escherichia coli IAI39",
  "Escherichia coli NA114",
  "Escherichia coli O103:H2 str. 12009",
  "Escherichia coli O111:H- str. 11128",
  "Escherichia coli O26:H11 str. 11368",
  "Escherichia coli IHE3034",
  "Escherichia coli ATCC 8739",
  "Escherichia coli 536",
  "Escherichia coli S88",
  "Escherichia coli O157:H7 str. Sakai",
  "Escherichia coli SE11",
  "Escherichia coli SE15",
  "Escherichia coli SMS-3-5",
  "Escherichia coli O157:H7 str. TW14359",
  "Escherichia coli UMN026",
  "Escherichia coli W",
  "Escherichia coli KO11FL",
  "Escherichia coli ETEC H10407",
  "Escherichia coli O55:H7 str. CB9615",
  "Escherichia coli LF82",
  "Escherichia coli O83:H1 str. NRG 857C",
  "Shigella flexneri 2a str. 2457T",
  "Shigella boydii CDC 3083-94",
  "Shigella boydii Sb227",
  "Shigella dysenteriae Sd197",
  "Shigella flexneri 5 str. 8401",
  "Shigella flexneri 2002017",
  "Shigella sonnei Ss046",
  "Escherichia coli UM146",
  "Escherichia coli UMNK88",
  "Escherichia coli UTI89",
  "Escherichia coli O157:H7 str. EDL933",
  "Klebsiella pneumoniae subsp. pneumoniae MGH 78578",
  "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2",
  "Geobacter metallireducens GS-15",
  "Synechocystis sp. PCC 6803",
  "Thermotoga maritima MSB8",
  "Clostridium ljungdahlii DSM 13528",
  "Escherichia coli str. K-12 substr. W3110",
  "Bacillus subtilis subsp. subtilis str. 168",
  "Lactococcus lactis subsp. cremoris MG1363",
  "Synechococcus elongatus PCC 7942",
  "Staphylococcus aureus subsp. aureus USA300_TCH1516",
  "Acinetobacter baumannii AYE",
  "Escherichia coli DH5[alpha",
  "Escherichia coli C",
  "Clostridioides difficile 630",
  "Salmonella pan-reactome",
  "Escherichia coli DH5[alpha]"
)


fungi_list <- c(
  "Saccharomyces cerevisiae S288C"
  # 如果有其他真菌也可以补充
)

archaea_list <- c(
  "Methanosarcina barkeri str. Fusaro"
  # 题目示例中只有这个古菌
)

plant_list <- c(
  # 如果将某些藻类视为植物，可以添加。但题目中的绿藻更倾向于放在 protist
  # 这里暂时留空
)

animal_list <- c(
  "Mus musculus",
  "Cricetulus griseus"
  # 题目中的小鼠、仓鼠等
)

environment_list <- c(
  # 题目中并没有明确环境来源示例，因此这里先留空
)

virus_list <- c(
  # 题目中没有出现病毒示例，留空
)

protist_list <- c(
  "Chlamydomonas reinhardtii",
  "Phaeodactylum tricornutum CCAP 1055/1",
  "Plasmodium falciparum 3D7",
  "Plasmodium vivax Sal-1",
  "Plasmodium berghei",
  "Plasmodium cynomolgi strain B",
  "Plasmodium knowlesi strain H",
  "Trypanosoma cruzi Dm28c"
)

drug_list <- c(
  # 如果有药物来源可以在此列出
)

food_list <- c(
  # 如果有食物来源可以在此列出
)

# 处理 DataFrame
# 处理 DataFrame
BIGG_database <- 
  BIGG_data %>%
  mutate(
    from_human = ifelse(str_detect(organism, str_c(human_list, collapse = "|")), "Yes", "Unknown"),
    from_which_part = "Unknown",
    from_bacteria = ifelse(str_detect(organism, str_c(str_replace_all(bacteria_list, "(\\[|\\])", "\\\\\\1"), collapse = "|")), "Yes", "Unknown"),
    from_which_bacteria = sapply(
      str_extract_all(organism, str_c(str_replace_all(bacteria_list, "(\\[|\\])", "\\\\\\1"), collapse = "|")), 
      function(x) ifelse(length(x) > 0, paste(unique(x), collapse = "{}"), "Unknown")
    ),
    from_fungi = ifelse(str_detect(organism, str_c(fungi_list, collapse = "|")), "Yes", "Unknown"),
    from_which_fungi = sapply(str_extract_all(organism, str_c(fungi_list, collapse = "|")), function(x) ifelse(length(x) > 0, paste(unique(x), collapse = "{}"), "Unknown")),
    from_archaea = ifelse(str_detect(organism, str_c(archaea_list, collapse = "|")), "Yes", "Unknown"),
    from_which_archaea = sapply(str_extract_all(organism, str_c(archaea_list, collapse = "|")), function(x) ifelse(length(x) > 0, paste(unique(x), collapse = "{}"), "Unknown")),
    from_plant = "Unknown",
    from_which_plant = "Unknown",
    from_animal = ifelse(str_detect(organism, str_c(animal_list, collapse = "|")), "Yes", "Unknown"),
    from_which_animal = sapply(str_extract_all(organism, str_c(animal_list, collapse = "|")), function(x) ifelse(length(x) > 0, paste(unique(x), collapse = "{}"), "Unknown")),
    from_environment = "Unknown",
    from_which_environment = "Unknown",
    from_virus = "Unknown",
    from_which_virus = "Unknown",
    from_protist = ifelse(str_detect(organism, str_c(protist_list, collapse = "|")), "Yes", "Unknown"),
    from_which_protist = sapply(str_extract_all(organism, str_c(protist_list, collapse = "|")), function(x) ifelse(length(x) > 0, paste(unique(x), collapse = "{}"), "Unknown")),
    from_drug = "Unknown",
    from_which_drug = "Unknown",
    from_food = "Unknown",
    from_which_food = "Unknown"
  )

dir.create("2_data/10_BIGG/BIGG_dataset", showWarnings = FALSE)
setwd("2_data/10_BIGG/BIGG_dataset")

save(BIGG_database, file = "BIGG_database2.rda")



test <- bigg_ms1@spectra.info
