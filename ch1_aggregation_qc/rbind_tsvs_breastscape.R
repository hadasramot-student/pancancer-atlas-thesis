# rbind_tsvs_breastscape.R
#
# Ch.1 (Aggregation & QC) — same as rbind_tsvs.R, scoped to the Breastscape
# source subdirectory (kept separate since Breastscape objects carry a
# predicted_doublet column the other sources don't).
#
# Inputs:
#   BREASTSCAPE_OBJECTS_DIR - directory of per-object *.tsv metadata dumps
# Outputs:
#   <BREASTSCAPE_OBJECTS_DIR>/all_metadata_merged.tsv
#
# Usage:
#   BREASTSCAPE_OBJECTS_DIR=/path/to/temp_objects/breastscape Rscript rbind_tsvs_breastscape.R

library(readr)
library(dplyr)
library(purrr)
library(tools)

input_dir  <- Sys.getenv("BREASTSCAPE_OBJECTS_DIR", unset = "data/temp_objects/breastscape")
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
    "predicted_doublet",
    "Meta.SingleR",
    "Meta.Sctype",
    "Meta.scimilarity",
    "patient"
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
