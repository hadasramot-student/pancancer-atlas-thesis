#!/usr/bin/env Rscript
# check_batch_counts.R
#
# Ch.3 (Batch correction) — reads every {CellType}_{CancerType}_merged.RDS
# and reports how many unique batches (project_name) each object contains.
# Used to decide which lineage x cancer-type objects need batch correction
# at all (single-batch objects are skipped).
#
# Inputs:  ATLAS_DATA_DIR/merged_by_query/<CancerType>/<CellType>_<CancerType>_merged.RDS
# Outputs: ATLAS_DATA_DIR/merged_by_query/batch_counts.csv
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript check_batch_counts.R

suppressPackageStartupMessages(library(Seurat))

BASE       <- file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "merged_by_query")
OUTPUT_CSV <- file.path(BASE, "batch_counts.csv")

CANCER_TYPES <- c(
  "Brain", "Colorectal", "Head_and_Neck", "Kidney", "Liver_Biliary",
  "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"
)

CELL_TYPES <- c(
  "B-cell", "Endothelial", "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells",
  "Myeloid", "NK-cells", "T-cells", "Epithelial", "other"
)

results <- vector("list", length(CANCER_TYPES) * length(CELL_TYPES))
idx <- 0L

for (ct in CANCER_TYPES) {
  for (cel in CELL_TYPES) {
    idx <- idx + 1L
    rds_path <- file.path(BASE, ct, paste0(cel, "_", ct, "_merged.RDS"))

    if (!file.exists(rds_path)) {
      cat(sprintf("[MISSING] %s\n", rds_path))
      results[[idx]] <- data.frame(
        cell_type   = cel,
        cancer_type = ct,
        rds_path    = rds_path,
        n_batches   = NA_integer_,
        status      = "missing",
        stringsAsFactors = FALSE
      )
      next
    }

    cat(sprintf("[READ] %s\n", rds_path))
    tryCatch({
      obj       <- readRDS(rds_path)
      n_batches <- length(unique(obj$project_name))
      status    <- if (n_batches > 1L) "multi_batch" else "single_batch"
      cat(sprintf("  -> %d batch(es) [%s]\n", n_batches, status))
      results[[idx]] <- data.frame(
        cell_type   = cel,
        cancer_type = ct,
        rds_path    = rds_path,
        n_batches   = n_batches,
        status      = status,
        stringsAsFactors = FALSE
      )
      rm(obj); gc(verbose = FALSE)
    }, error = function(e) {
      cat(sprintf("  -> ERROR: %s\n", conditionMessage(e)))
      results[[idx]] <<- data.frame(
        cell_type   = cel,
        cancer_type = ct,
        rds_path    = rds_path,
        n_batches   = NA_integer_,
        status      = paste0("error: ", conditionMessage(e)),
        stringsAsFactors = FALSE
      )
    })
  }
}

out <- do.call(rbind, results)
write.csv(out, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("\n=== Batch count summary ===\n"))
cat(sprintf("Output written to: %s\n", OUTPUT_CSV))
cat(sprintf("  Multi-batch  (will run scANVI): %d\n",
            sum(out$status == "multi_batch",  na.rm = TRUE)))
cat(sprintf("  Single-batch (will skip):       %d\n",
            sum(out$status == "single_batch", na.rm = TRUE)))
cat(sprintf("  Missing RDS  (will skip):       %d\n",
            sum(out$status == "missing",      na.rm = TRUE)))
cat(sprintf("  Errors:                         %d\n",
            sum(startsWith(out$status, "error"), na.rm = TRUE)))
