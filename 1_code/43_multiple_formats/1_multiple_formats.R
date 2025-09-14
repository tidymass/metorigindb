library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

# load the metorigin database
load("2_data/41_MetOriginDB/metorigindb_database_xz.rda")

colnames(metorigindb_database)

# put bacteria together
# original columns
cols <- colnames(metorigindb_database)

# columns to move
cols_to_move <- c("bacteria_ncbi_id", "bacteria_phylum", "bacteria_class", 
                  "bacteria_order", "bacteria_family", "bacteria_genus", "bacteria_species")

# find target position
pos <- match("from_which_bacteria", cols)

# new order of columns
new_order <- append(cols[-match(cols_to_move, cols)], cols_to_move, after = pos)

metorigindb_database_final <- metorigindb_database[, new_order]

# save as rda
save(metorigindb_database_final, file = "2_data/41_MetOriginDB/metorigindb_database_final.rda", compress = "xz")
# save as csv
write.csv(metorigindb_database_final, file = "2_data/41_MetOriginDB/metorigindb_database.csv", row.names = FALSE)
# save as tsv

library(data.table)
fwrite(
  metorigindb_database_final,
  file = "2_data/41_MetOriginDB/metorigindb_database.tsv",
  sep = "\t",
  quote = TRUE,
  na = ""
)

# ways to read
data_csv <- read.csv("~/Tidymass/MetOriginDB/2_data/41_MetOriginDB/metorigindb_database.csv", header = TRUE, sep = ",")
# 
data_tsv <- read.delim(
  "2_data/41_MetOriginDB/metorigindb_database.tsv",
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,           
  blank.lines.skip = TRUE 
)
