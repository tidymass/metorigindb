library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

load("2_data/0_MS1/ms1_database2.rda")

# select all the columns start with "from_"
from_cols <- grep("^from_", names(ms1_database), value = TRUE)

# keep rows with at least one "Yes"
ms1_database_filter <- ms1_database[apply(ms1_database[from_cols], 1, function(row) any(row == "Yes")), ]

MetOriginDB <- ms1_database_filter

save(MetOriginDB, file = "2_data/41_MetOriginDB/MetOriginDB.rda")

