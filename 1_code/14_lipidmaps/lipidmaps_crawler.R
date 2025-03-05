library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

library(cinf)

library(ChemmineR)

library(massdatabase)

setwd('2_data/14_LIPIDMAPS')

lipidmaps <-
  read.SDFset(sdfstr = "structures.sdf", skipErrors = TRUE)

lipidmaps_result <-
  1:length(lipidmaps) %>%
  purrr::map(function(i){
    cat(i, " ")
    x <- lipidmaps[[i]]
    result <-
      tryCatch(matrix(x[[4]], nrow = 1) %>%
                 as.data.frame(), error = NULL)
    if(is.null(result)){
      return(NULL)
    }
    colnames(result) <- names(x[[4]])
    result
  })

all_column_name <-
  lipidmaps_result %>%
  lapply(colnames) %>%
  unlist() %>%
  unique()


lipidmaps_result <-
  lipidmaps_result <-
  1:length(lipidmaps_result) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <- lipidmaps_result[[i]]
    
    diff_names <- setdiff(all_column_name, colnames(x))
    if (length(diff_names) == 0) {
      return(x[, all_column_name])
    }
    add_info <-
      matrix(NA, nrow = 1, ncol = length(diff_names)) %>%
      as.data.frame()
    colnames(add_info) <- diff_names
    cbind(x, add_info) %>%
      dplyr::select(all_column_name)
  })


lipidmaps_result <-
  lipidmaps_result %>%
  dplyr::bind_rows() %>%
  as.data.frame()

library(dplyr)

# load the lipid bank ids
load("lm_ids.rda")

# set the path to store outputs
output_dir <- ""

# download the files from 1 to 300 in the list
for (i in 1:300) {
  cat(i, " ")
  lm_id <- lm_ids[i]
  
  # construct the file path and check whether it exists
  output_file <- file.path(output_dir, paste0(lm_id, ".rda"))
  if (file.exists(output_file)) {
    next()
  }
  
  # request the data
  x <- tryCatch({
    request_lipidmaps_lipid(lipid_id = lm_id)
  }, error = function(e) {
    cat("Error for LM_ID:", lm_id, "-", e$message, "\n")
    NULL
  })
  
  # save as .rda file
  if (!is.null(x)) {
    save(x, file = output_file)
  } else {
    cat("Failed to retrieve data for LM_ID:", lm_id, "\n")
  }
  
  # add rate limits
  Sys.sleep(1)
}
