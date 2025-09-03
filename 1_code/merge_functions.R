# (a) 花括号去重合并
merge_braces <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(NA)
  
  # 2. 以 {} 进行分割 -> 打平为字符向量 -> 去重 -> 重新拼接
  x %>%
    str_split("\\{\\}") %>%    # 拆分成列表
    unlist() %>%          # 转换为字符向量
    unique() %>%               # 去重
    paste(collapse = "{}")     # 用 {} 重新拼接成字符串
}
# 
# (b) 根据数据来源优先选择
#    - 默认优先 "HMDB"
take_first_source <- function(x, preferred_source) {
  # 在 summarise 的分组环境中，获取当前分组下所有列
  all_data <- cur_data()
  
  # 取出数据库来源列（确保你的数据中有名为 "Database_source" 的列）
  src <- all_data$Database_source
  
  # 只保留 x 非 NA 的行
  valid_idx <- !is.na(x)
  x_valid   <- x[valid_idx]
  src_valid <- src[valid_idx]
  
  # 如果全是 NA，直接返回 NA（真正的 NA，而不是 "NA" 字符串）
  if (length(x_valid) == 0) {
    return(NA)
  }
  
  # 找到首选数据库所在行的索引
  preferred_idx <- which(src_valid == preferred_source)
  if (length(preferred_idx) > 0) {
    # 如果有来自首选数据库的值，返回该组中【第一个】非 NA 值
    return(x_valid[preferred_idx[1]])
  } else {
    # 否则，返回该组的第一个非 NA 值
    return(x_valid[1])
  }
}



merge_from_column <- function(x, col_name) {
  
  # (a) 如果是 from_which_ 开头
  if (str_starts(col_name, "from_which_")) {
    x_no_na <- na.omit(x)
    x_no_unknown <- x_no_na[x_no_na != "Unknown"]
    if (length(x_no_unknown) == 0) {
      return("Unknown")
    } else {
      return(
        x_no_unknown %>%
          str_split("\\{\\}") %>%    # 拆分成列表
          unlist() %>%          # 转换为字符向量
          unique() %>%               # 去重
          paste(collapse = "{}")
      )
    }
  }
  
  # (b) 否则当作 from_xxx 列
  #     若任何一条记录为 "Yes" 则返回 "Yes"，否则返回 "Unknown"
  if (any(x == "Yes", na.rm = TRUE)) {
    return("Yes")
  } else {
    return("Unknown")
  }
}