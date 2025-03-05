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

pubchem_tax_data <- read.csv("2_data/23_Pubchem/Pubchem_compounds/taxonomy.csv", header = TRUE)

cid_vector <- 
  pubchem_tax_data %>% 
  dplyr::pull(X.cid)

# 设置工作目录（可选）
setwd("/Users/ejo/Desktop")

# 指定保存文件夹路径
save_folder <- "/Users/ejo/Desktop/test_files"


# 获取可用的核心数，并行计算使用 (留 1 个核心空闲)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 使用 foreach 并行下载
foreach(cid = cid_vector, .packages = "utils") %dopar% {
  
  # 目标文件路径
  file_name <- paste0(cid, ".json")
  file_path <- file.path(save_folder, file_name)
  
  # 如果文件已存在，跳过下载
  if (!file.exists(file_path)) {
    
    # 构造下载链接
    download_url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?",
      "infmt=json&outfmt=json&query={",
      "\"download\":\"*\",",
      "\"collection\":\"consolidatedcompoundtaxonomy\",",
      "\"order\":[\"cid,asc\"],",
      "\"start\":1,",
      "\"limit\":10000000,",
      "\"downloadfilename\":\"pubchem_cid_", cid, "_consolidatedcompoundtaxonomy\",",
      "\"where\":{\"ands\":[{\"cid\":\"", cid, "\"}]}",
      "}"
    )
    
    # 打印下载信息（仅在主进程中显示）
    message("正在下载 cid = ", cid, " 的文件...")
    
    # 下载文件
    tryCatch({
      download.file(url = download_url,
                    destfile = file_path,
                    mode = "wb")
      message("下载完成: ", file_name)
    }, error = function(e) {
      warning("下载失败: ", file_name, " - ", e$message)
    })
    
  } else {
    return(NULL)
  }
}

# 关闭并行计算
stopCluster(cl)

