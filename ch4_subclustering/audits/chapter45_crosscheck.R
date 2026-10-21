#!/usr/bin/env Rscript
# chapter45_crosscheck.R - sum cell counts across annotated_lineages/*.RDS
# (Ch.4 lineage x cancer-type post-annotation objects) as a sanity check
# against the atlas total. Loads one file at a time to bound memory.

suppressPackageStartupMessages(library(Seurat))

TO <- "/ofir_grp/users/hadasra/temp_objects/annotated_lineages"
files <- list.files(TO, pattern = "_annotated\\.RDS$", full.names = TRUE)
cat("Found", length(files), "annotated_lineages files\n\n")

total <- 0
results <- data.frame(file = character(0), n_cells = integer(0))
for (f in files) {
  n <- tryCatch({
    obj <- readRDS(f)
    nc <- ncol(obj)
    rm(obj); gc(verbose = FALSE)
    nc
  }, error = function(e) {
    cat("  ERROR reading", basename(f), ":", conditionMessage(e), "\n")
    NA
  })
  cat(sprintf("%-55s %10s\n", basename(f), ifelse(is.na(n), "ERROR", format(n, big.mark=","))))
  if (!is.na(n)) {
    total <- total + n
    results <- rbind(results, data.frame(file = basename(f), n_cells = n))
  }
}

cat("\n=== TOTAL across", nrow(results), "successfully-read annotated_lineages files:", format(total, big.mark=","), "cells ===\n")
write.csv(results, "/ofir_grp/users/hadasra/temp_objects/chapter4_crosscheck_counts.csv", row.names = FALSE)
