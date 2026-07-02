#!/usr/bin/env Rscript
# dump_session_info.R
#
# Run this inside the R environment used for the pipeline (e.g. the Rstudio
# apptainer container) to capture exact package versions for the thesis
# appendix / reproducibility record. Writes sessionInfo.txt next to this
# script. If `renv` is available and this repo has been renv::init()'d, also
# refreshes renv.lock.
#
# Usage:
#   Rscript dump_session_info.R

for (pkg in readLines("r_packages.txt")) {
  pkg <- trimws(sub("#.*$", "", pkg))
  if (!nzchar(pkg) || grepl("/", pkg)) next   # skip comments/blank + GitHub-only refs
  suppressWarnings(suppressMessages(try(library(pkg, character.only = TRUE), silent = TRUE)))
}

sink("sessionInfo.txt")
print(sessionInfo())
sink()
message("Wrote sessionInfo.txt")

if (requireNamespace("renv", quietly = TRUE) && file.exists("renv.lock")) {
  renv::snapshot()
  message("Updated renv.lock")
}
