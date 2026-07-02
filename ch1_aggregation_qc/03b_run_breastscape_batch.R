# 03b_run_breastscape_batch.R
#
# Ch.1 (Aggregation & QC) — batch driver for the Breastscape source.
#
# Calls 00_main_dispatch.R (--ptype breastscape) once per .h5ad file found
# in the raw Breastscape directory.
#
# Inputs:
#   RAW_BREASTSCAPE_DIR - directory of raw per-project Breastscape .h5ad files
# Outputs:
#   Invokes 00_main_dispatch.R --ptype breastscape --h5ad_path <file> once
#   per project (each call writes its own ingested/QC'd RDS via
#   03_ingest_breastscape.R).
#
# Usage:
#   RAW_BREASTSCAPE_DIR=/path/to/breastscape_h5ads Rscript 03b_run_breastscape_batch.R

input_dir <- Sys.getenv("RAW_BREASTSCAPE_DIR", unset = "data/raw/breastscape")

# all .h5ad files in that folder, with full paths
files_in <- list.files(
  path = input_dir,
  pattern = "\\.h5ad$",
  full.names = TRUE
)

for (f in files_in) {
  cmd <- sprintf(
    "Rscript 00_main_dispatch.R --ptype breastscape --h5ad_path %s",
    shQuote(f)
  )
  message("Running: ", cmd)
  system(cmd)
}
