# 01_ingest_scbasecount.R
#
# Ch.1 (Aggregation & QC) — batch driver for the ScBaseCount source.
#
# Scans a directory of per-project ScBaseCount RDS files and calls
# 00_main_dispatch.R (--ptype scbasecount) on each one that hasn't already
# been ingested (i.e. isn't present, post-QC, in the output directory).
#
# Inputs:
#   RAW_SCBASECOUNT_DIR - directory of raw per-project ScBaseCount .RDS files
#   ATLAS_DATA_DIR       - output directory checked for already-ingested
#                           "<project>_step4.RDS" files (skip list)
# Outputs:
#   Invokes 00_main_dispatch.R --ptype scbasecount --rds <file> once per
#   unprocessed project (each call writes its own "<project>_step4.RDS").
#
# Usage:
#   RAW_SCBASECOUNT_DIR=/path/to/scbasecount_rds \
#   ATLAS_DATA_DIR=/path/to/temp_objects \
#     Rscript 01_ingest_scbasecount.R

input_dir   <- Sys.getenv("RAW_SCBASECOUNT_DIR", unset = "data/raw/scbasecount")
exclude_dir <- Sys.getenv("ATLAS_DATA_DIR",       unset = "data/temp_objects")

# all RDS in main
main_files <- list.files(input_dir, pattern = "\\.RDS$", full.names = TRUE)

# all RDS already ingested (take only the names)
exclude_names <- basename(list.files(exclude_dir, pattern = "\\.RDS$", full.names = TRUE))
exclude_names_normalized <- sub("(?i)_step4(?=\\.rds$)", "", exclude_names, perl = TRUE)

# keep only those not already ingested
files_to_process <- main_files[!basename(main_files) %in% exclude_names_normalized]

for (f in files_to_process) {
  cmd <- sprintf(
    "Rscript 00_main_dispatch.R --ptype scbasecount --rds %s",
    shQuote(f)
  )
  message("Running: ", cmd)
  system(cmd)
}
