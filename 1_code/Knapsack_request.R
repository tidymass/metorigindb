library(rvest)
library(httr)
library(xml2)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# 函数：request_knapsack_compound
request_knapsack_compound <- function(C_ID = "C00000001") {
  
  # 构建目标 URL
  url_template <- paste0("https://www.knapsackfamily.com/knapsack_core/information.php?word=", C_ID)
  
  # 尝试读取网页
  page <- tryCatch(
    read_html(url_template),
    error = function(e) {
      message("读取网页出错: ", e)
      return(NULL)
    }
  )
  
  # 如果无法读取页面，则返回 NULL
  if (is.null(page)) {
    return(NULL)
  }
  
  # 通用解析字段的函数
  extract_field <- function(page, field_name) {
    node <- page %>%
      html_node(xpath = paste0("//th[@class='inf' and contains(., '", field_name, "')]/following-sibling::td[1]"))
    
    if (is.na(node)) {
      return("")
    } else {
      raw_text <- as.character(node)
      clean_text <- raw_text %>%
        # 把 <br> 标签替换为 {} ，再去除所有HTML标签并压缩空格
        str_replace_all("(?i)<br\\s*/?>", "{}") %>%
        str_replace_all("<[^>]+>", "") %>%
        str_squish()
      return(clean_text)
    }
  }
  
  # 提取字段
  clean_name  <- extract_field(page, "Name")
  formula     <- extract_field(page, "Formula")
  mw          <- extract_field(page, "Mw")
  cas_rn      <- extract_field(page, "CAS RN")
  
  # 提取 C_ID 字段
  c_id_text <- extract_field(page, "C_ID")
  
  # 仅保留第一个形如 CXXXXXXXX 的片段（C 开头 + 8 位数字）
  c_id <- c_id_text %>%
    # 清理一些可能的多余文本
    str_replace_all("//-->", "") %>%
    str_replace_all("\\{\\}", "") %>%
    str_extract("C\\d{8}") %>%
    coalesce("")
  
  inchi_key   <- extract_field(page, "InChIKey")
  inchi_code  <- extract_field(page, "InChICode")
  smiles      <- extract_field(page, "SMILES")
  start_substs <- extract_field(page, "Start Substs in Alk. Biosynthesis Prediction")
  
  # 解析 Organism 表格 (table.org)
  organism_table <- page %>% html_node("table.org")
  
  if (is.na(organism_table)) {
    organism_data <- tibble(Kingdom = character(),
                            Family  = character(),
                            Species = character())
  } else {
    org_trs <- organism_table %>% html_nodes("tr")
    if (length(org_trs) > 0) {
      organism_data <- org_trs %>%
        map_dfr(~ {
          tds <- .x %>% html_nodes("td")
          # 对于标题行或内容行长度不够的情况，做个保护
          if (length(tds) < 3) {
            tibble(
              Kingdom = NA_character_,
              Family  = NA_character_,
              Species = NA_character_
            )
          } else {
            tibble(
              Kingdom = tds[1] %>% html_text(trim = TRUE),
              Family  = tds[2] %>% html_text(trim = TRUE),
              Species = tds[3] %>% html_text(trim = TRUE)
            )
          }
        })
    } else {
      organism_data <- tibble(Kingdom = character(),
                              Family  = character(),
                              Species = character())
    }
  }
  
  # 合并物种信息
  if (nrow(organism_data) == 0) {
    combined_kingdom <- ""
    combined_family  <- ""
    combined_species <- ""
  } else {
    combined_kingdom <- paste(organism_data$Kingdom,  collapse = "{}")
    combined_family  <- paste(organism_data$Family,   collapse = "{}")
    combined_species <- paste(organism_data$Species,  collapse = "{}")
  }
  
  # 组合所有字段
  x <- tibble(
    Name    = clean_name,
    Formula = formula,
    Mw      = mw,
    CAS_RN  = cas_rn,
    C_ID    = c_id,
    InChIKey = inchi_key,
    InChICode = inchi_code,
    SMILES = smiles,
    Start_Substs_in_Alk._Biosynthesis_Prediction_ = start_substs,
    Kingdom = combined_kingdom,
    Family  = combined_family,
    Species = combined_species
  )
  
  # 将所有空字符串 ("") 统一替换为 NA
  x <- x %>%
    mutate(across(everything(), ~ na_if(.x, "")))
  
  return(x)
}

# 函数调用示例：
# x <- request_knapsack_compound(C_ID = "C00000100")
# print(x)

