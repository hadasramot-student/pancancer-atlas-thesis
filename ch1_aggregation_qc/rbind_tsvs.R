# rbind_tsvs.R
#
# Ch.1 (Aggregation & QC) — concatenate all per-object metadata TSVs
# (written during ingestion/annotation) into one atlas-wide metadata table.
# Used as the source table for the atlas composition figures/tables.
#
# Inputs:
#   ATLAS_DATA_DIR - directory containing per-object *.tsv metadata dumps
# Outputs:
#   <ATLAS_DATA_DIR>/all_metadata_merged.tsv
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript rbind_tsvs.R

library(readr)
library(dplyr)
library(purrr)
library(tools)

input_dir  <- Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects")
output_tsv <- file.path(input_dir, "all_metadata_merged.tsv")

# 1) list all .tsv files
tsv_files <- list.files(
  path       = input_dir,
  pattern    = "\\.tsv$",
  full.names = TRUE
)

# 2) read and add a column with source file (optional but useful)
metadata_list <- map(tsv_files, function(path) {
  dat <- read_tsv(path, show_col_types = FALSE)
  dat$source_object <- file_path_sans_ext(basename(path))

  # harmonize types only if columns exist
  cols_char <- c(
    "Meta.SingleR",
    "Meta.Sctype",
    "Meta.scimilarity",
    "patient",
    "Consensus_Cell_Type"
  )

  for (cn in cols_char) {
    if (cn %in% names(dat)) {
      dat[[cn]] <- as.character(dat[[cn]])
    }
  }

  dat
})

all_meta <- bind_rows(metadata_list)

write_tsv(all_meta, output_tsv)
