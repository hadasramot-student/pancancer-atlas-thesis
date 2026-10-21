#!/usr/bin/env Rscript
# extract_breast_chunk_metadata.R
#
# Loads all 13 Breast Epithelial chunk RDS files, extracts per-cell metadata,
# and saves as a CSV to be joined onto the h5ad.

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
suppressPackageStartupMessages(library(Seurat))

CHUNK_DIR <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/tmp_merge/Epithelial"
OUT_CSV   <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/breast_epithelial_metadata.csv"

chunk_files <- sort(list.files(CHUNK_DIR, pattern = "chunk_.*\\.RDS$", full.names = TRUE))
log_msg(sprintf("Found %d chunk files", length(chunk_files)))

all_meta <- list()
all_cols <- NULL  # union of all columns seen across chunks

for (i in seq_along(chunk_files)) {
  f <- chunk_files[i]
  log_msg(sprintf("Loading chunk %d/%d: %s", i, length(chunk_files), basename(f)))
  obj <- readRDS(f)
  log_msg(sprintf("  %d cells, %d metadata columns", ncol(obj), ncol(obj@meta.data)))

  meta <- obj@meta.data
  meta$cell_barcode <- rownames(meta)

  # Track union of all columns
  all_cols <- union(all_cols, colnames(meta))

  all_meta[[i]] <- meta
  rm(obj); gc()
}

# Align all chunks to the full column union (fill missing cols with NA)
log_msg(sprintf("Total metadata columns across all chunks: %d", length(all_cols)))
all_meta <- lapply(all_meta, function(df) {
  missing <- setdiff(all_cols, colnames(df))
  if (length(missing) > 0) {
    df[missing] <- NA
  }
  df[, all_cols, drop = FALSE]
})

log_msg("Combining all chunks...")
combined <- do.call(rbind, all_meta)
log_msg(sprintf("Total cells: %d", nrow(combined)))

if ("scimilarity_celltype" %in% colnames(combined)) {
  log_msg("scimilarity_celltype distribution (top 15):")
  top <- sort(table(combined$scimilarity_celltype), decreasing = TRUE)[1:15]
  print(top)
}

log_msg(sprintf("Writing to: %s", OUT_CSV))
write.csv(combined, OUT_CSV, row.names = FALSE)
log_msg("Done.")
