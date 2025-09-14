library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(taxonomizr)
# sql_db <- "/Users/yijiang/Tidymass/taxadb/nameNode.sqlite"
# prepareDatabase(sqlFile = sql_db)

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

split_bacteria_cell_raw <- function(x) {
  if (is.na(x) || !nzchar(x) || trimws(x) == "Unknown") {
    return("Unknown")
  }
  parts <- unlist(strsplit(x, "\\{\\}"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return("Unknown")
  parts
}

## 你已有：
split_bacteria_cell_raw <- function(x) {
  if (is.na(x) || !nzchar(x) || trimws(x) == "Unknown") {
    return("Unknown")
  }
  parts <- unlist(strsplit(x, "\\{\\}"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return("Unknown")
  parts
}


SQL_DB_DEFAULT <- "~/tidymass/metorigindb/2_data/42_SQL/nameNode.sqlite"
RANKS_OUT  <- c("phylum","class","order","family","genus","species")
RANKS_TAKE <- RANKS_OUT
## - augment_with_taxonomy_offline_rawmatch()（我们会在里面改“解析 taxid”的那一行）

## 1) 名称规范化：去首尾空格 & 合并多空格
.norm_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- stringr::str_squish(x)
  # 空串统一设为 NA，避免进映射表
  x[!nzchar(x)] <- NA_character_
  x
}
## 2) 用手工表构建 override map：名称(规范化) -> ID（字符）
make_override_map <- function(unknown_bacteria_id) {
  unknown_bacteria_id %>%
    dplyr::transmute(
      key = .norm_name(Bacteria),
      val = as.character(NCBI_ID)
    ) %>%
    dplyr::filter(!is.na(key), !is.na(val), nzchar(val), key != "Unknown") %>%
    dplyr::distinct(key, .keep_all = TRUE) %>%
    tibble::deframe()   # 得到一个命名向量 / list： override_map[["规范化名称"]] = "ID"
}

## 3) 解析函数（先查 override，再走原逻辑）
resolve_taxid_with_override <- function(name, sql_db, override_map = NULL) {
  if (is.na(name) || !nzchar(name) || identical(name, "Unknown")) return("Unknown")
  
  nm <- .norm_name(name)[1]
  if (is.na(nm) || !nzchar(nm)) return("Unknown")
  
  # override 优先 —— 用 %in% 判断
  if (!is.null(override_map) && nm %in% names(override_map)) {
    return(as.character(override_map[[nm]]))
  }
  
  # 原逐级截短匹配
  tokens <- strsplit(nm, "\\s+")[[1]]
  if (length(tokens) == 0) return("Unknown")
  for (k in seq.int(length(tokens), 1)) {
    cand <- paste(tokens[1:k], collapse = " ")
    id <- suppressWarnings(taxonomizr::getId(cand, sqlFile = sql_db))
    if (!is.na(id)) return(as.character(id))
  }
  "Unknown"
}

## 4) 在你的 augment_with_taxonomy_offline_rawmatch() 里，增加一个 override_map 参数，
##    并仅改“解析 taxid”的那一行。其余保持原样即可。
augment_with_taxonomy_offline_rawmatch <- function(df,
                                                   col = "from_which_bacteria",
                                                   sql_db = SQL_DB_DEFAULT,
                                                   override_map = NULL) {
  if (!col %in% names(df)) {
    stop(sprintf("列 '%s' 不存在于提供的数据框中。", col))
  }
  # 1) 行内拆分
  name_lists <- df[[col]] %>% purrr::map(split_bacteria_cell_raw)
  # 2) 唯一名称全集
  all_names_raw <- unique(unlist(name_lists, use.names = FALSE))
  # 3) 解析 taxid —— 只改这一段：调用 resolve_taxid_with_override()
  resolved_ids <- vapply(
    all_names_raw,
    function(x) resolve_taxid_with_override(x, sql_db, override_map),
    FUN.VALUE = character(1)
  )
  
  # 4) 批量 getTaxonomy（只对匹配到的 id）
  hit_idx <- which(resolved_ids != "Unknown")
  tax_df <- NULL
  if (length(hit_idx) > 0) {
    tax_mat <- taxonomizr::getTaxonomy(
      as.integer(resolved_ids[hit_idx]),
      sqlFile = sql_db,
      desiredTaxa = RANKS_TAKE,
      getNames  = TRUE
    )
    tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
    tax_df[is.na(tax_df)] <- ""
    row_for_i <- rep(NA_integer_, length(all_names_raw))
    row_for_i[hit_idx] <- seq_along(hit_idx)
  } else {
    row_for_i <- rep(NA_integer_, length(all_names_raw))
    tax_df <- data.frame()
  }
  
  .pick_safe <- function(x) {
    if (is.null(x) || length(x) == 0) return("Unknown")
    x <- x[1]; if (is.na(x) || !nzchar(x)) return("Unknown"); x
  }
  
  # 5) 名称 -> (taxid, 各级)
  lookup <- vector("list", length(all_names_raw))
  names(lookup) <- all_names_raw
  for (i in seq_along(all_names_raw)) {
    if (identical(resolved_ids[i], "Unknown") || is.na(row_for_i[i])) {
      lookup[[i]] <- c(
        list(taxid = "Unknown"),
        as.list(setNames(rep("Unknown", length(RANKS_OUT)), RANKS_OUT))
      ); next
    }
    this <- tax_df[row_for_i[i], , drop = FALSE]
    lookup[[i]] <- list(
      taxid   = as.character(resolved_ids[i]),
      phylum  = .pick_safe(this$phylum),
      class   = .pick_safe(this$class),
      order   = .pick_safe(this$order),
      family  = .pick_safe(this$family),
      genus   = .pick_safe(this$genus),
      species = .pick_safe(this$species)
    )
  }
  
  # 6) 折叠回去
  pick_and_collapse <- function(vec_names, field) {
    if (length(vec_names) == 1 && identical(vec_names, "Unknown")) return("Unknown")
    vals <- vapply(vec_names, function(nn) lookup[[nn]][[field]], FUN.VALUE = character(1))
    vals[is.na(vals) | !nzchar(vals)] <- "Unknown"
    paste(vals, collapse = "{}")
  }
  new_cols <- tibble::tibble(
    bacteria_ncbi_id = purrr::map_chr(name_lists, ~ pick_and_collapse(.x, "taxid")),
    bacteria_phylum  = purrr::map_chr(name_lists, ~ pick_and_collapse(.x, "phylum")),
    bacteria_class   = purrr::map_chr(name_lists, ~ pick_and_collapse(.x, "class")),
    bacteria_order   = purrr::map_chr(name_lists, ~ pick_and_collapse(.x, "order")),
    bacteria_family  = purrr::map_chr(name_lists, ~ pick_and_collapse(.x, "family")),
    bacteria_genus   = purrr::map_chr(name_lists, ~ pick_and_collapse(.x, "genus")),
    bacteria_species = purrr::map_chr(name_lists, ~ pick_and_collapse(.x, "species"))
  )
  dplyr::bind_cols(df, new_cols)
}

## 5) 使用方式（从头跑一遍，但遇到手工表名称就用其ID）
# 读取初始数据
load("2_data/41_MetOriginDB/MetOriginDB_clean.rda")
# 读取手工表（两列：Bacteria, NCBI_ID）
unknown_bacteria_id <- readr::read_csv("2_data/41_MetOriginDB/unknown_bacteria_id.csv", show_col_types = FALSE)

# 构建 override_map
override_map <- make_override_map(unknown_bacteria_id)

# 跑增强（与原来一致，但多传了 override_map）
sql_db <- "~/tidymass/metorigindb/2_data/42_SQL/nameNode.sqlite"
metorigindb_database <- augment_with_taxonomy_offline_rawmatch(
  df      = MetOriginDB_clean,
  col     = "from_which_bacteria",
  sql_db  = sql_db,
  override_map = override_map
)

# 可选：保存结果
save(metorigindb_database, file = "~/tidymass/metorigindb/2_data/41_MetOriginDB/metorigindb_database_filling.rda")
# write.csv(metorigindb_database, "2_data/41_MetOriginDB/MetOriginDB_clean_aug.csv", row.names = FALSE)



metorigindb_database_count <- 
  metorigindb_database %>% 
  mutate(
    bacteria_count = ifelse(is.na(from_which_bacteria) | from_which_bacteria == "Unknown", 
                            0, 
                            str_count(from_which_bacteria, "\\{\\}") + 1),
    Unknown_count = str_count(bacteria_ncbi_id, "Unknown"))

real_test <- metorigindb_database_count %>% 
  dplyr::select(bacteria_count, Unknown_count)

unknown_bacteria <- 
  metorigindb_database %>%
  mutate(
    from_list = strsplit(from_which_bacteria, "\\{\\}"),
    id_list   = strsplit(bacteria_ncbi_id, "\\{\\}")
  ) %>%
  mutate(
    unknown_names = map2(from_list, id_list, ~ .x[.y == "Unknown"])
  ) %>%
  pull(unknown_names) %>%
  unlist() %>%
  trimws() %>%
  setdiff("Unknown") %>%
  unique()



