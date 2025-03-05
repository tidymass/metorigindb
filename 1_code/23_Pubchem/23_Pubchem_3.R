library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)


library(jsonlite)
library(dplyr)

# 设定存放 .json 文件的文件夹路径
json_folder <- "C:\\Users\\Ejo\\Desktop\\test"

# 获取文件夹下所有 .json 文件的完整路径
json_files <- list.files(json_folder, pattern = "\\.json$", full.names = TRUE)

# 读取并合并
df_list <- lapply(json_files, function(file) {
  # 从当前文件读取 JSON
  dat <- fromJSON(file)
  
  # 转为数据框后，只取需要的列
  # 如果文件确实都包含 cid, taxid, taxname 这三列，可直接这样取
  # 若字段名称可能不一致或有缺失，需要做相应的错误处理
  dat_df <- data.frame(
    cid = dat$cid,
    taxid = dat$taxid,
    taxname = dat$taxname,
    stringsAsFactors = FALSE
  )
  
  return(dat_df)
})

# 将所有文件的结果合并为一个 data.frame
final_df <- bind_rows(df_list)


