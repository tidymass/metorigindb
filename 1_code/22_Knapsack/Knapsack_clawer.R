library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source('1_code/Knapsack_request.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

C_IDs <- sprintf("C%08d", 1:64652)

output_dir <- "2_data/22_Knapsack/K_combine"

# download the files from 1 to 10000 in the list
for (i in 1:length(C_IDs)) {
  cat(i, " ")
  C_id <- C_IDs[i]
  
  # construct the file path and check whether it exists
  output_file <- file.path(output_dir, paste0(C_id, ".rda"))
  if (file.exists(output_file)) {
    next()
  }
  
  # request the data
  x <- tryCatch({
    request_knapsack_compound(C_ID = C_id)
  }, error = function(e) {
    cat("Error for C_ID:", C_id, "-", e$message, "\n")
    NULL
  })
  
  # save as .rda file
  if (!is.null(x)) {
    save(x, file = output_file)
  } else {
    cat("Failed to retrieve data for C_ID:", C_id, "\n")
  }
  
  # add rate limits
  Sys.sleep(1)
}













