# create_tsvs.R
#
# Ch.1 (Aggregation & QC) — dump per-object metadata to TSV for a directory
# of ingested Breastscape RDS objects (used as an input to rbind_tsvs_breastscape.R).
#
# Inputs:
#   BREASTSCAPE_OBJECTS_DIR - directory of ingested Breastscape *.RDS objects
# Outputs:
#   <BREASTSCAPE_OBJECTS_DIR>/<project_name>.tsv per object
#
# Usage:
#   BREASTSCAPE_OBJECTS_DIR=/path/to/temp_objects/breastscape Rscript create_tsvs.R

source("common/functions.R")
library(dplyr)

input_dir <- Sys.getenv("BREASTSCAPE_OBJECTS_DIR", unset = "data/temp_objects/breastscape")

files_in <- list.files(
  path = input_dir,
  pattern = "\\.RDS$",
  full.names = TRUE
)

for (f in files_in) {
  print("loading")
  obj <- readRDS(f)
  project_name <- obj@project.name
  out_path <- paste0(input_dir, "/", project_name, ".tsv")
  print("saving")
  save_metadata_tsv_2(obj, out_path)
  print(paste0("finished ", project_name))
}
