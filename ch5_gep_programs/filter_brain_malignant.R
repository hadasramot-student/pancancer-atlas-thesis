#!/usr/bin/env Rscript
# filter_brain_malignant.R
#
# Ch.5 (GEP programs) — tumor/malignant cell selection, Brain special case.
# Unlike the other 9 cancer types (where the Ch.4 "Epithelial_<Cancer>"
# lineage-split object already *is* the tumor compartment), Brain tumors
# are not epithelial, so they land in the "other_Brain" lineage object
# under a "Malignant" cell_type label. This script filters that object down
# to just the Malignant cells and saves it under an "Epithelial_Brain_..."
# name so it loads identically to the other cancers' cNMF/GeneNMF inputs.
#
# Inputs:  ATLAS_DATA_DIR/merged_by_query/Brain/other_Brain_merged_full_scanvi_final.RDS
# Outputs: ATLAS_DATA_DIR/merged_by_query/Brain/Epithelial_Brain_merged_full_scanvi_final.RDS
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript filter_brain_malignant.R

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
.libPaths(c(userlib[nzchar(userlib)], "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

atlas_data_dir <- Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects")
INPUT  <- file.path(atlas_data_dir, "merged_by_query/Brain/other_Brain_merged_full_scanvi_final.RDS")
OUTPUT <- file.path(atlas_data_dir, "merged_by_query/Brain/Epithelial_Brain_merged_full_scanvi_final.RDS")

log_msg("Loading other_Brain_merged_full_scanvi_final.RDS...")
obj <- readRDS(INPUT)
log_msg(sprintf("Loaded: %d cells x %d genes", ncol(obj), nrow(obj)))

log_msg("Cell type distribution:")
print(sort(table(obj$cell_type), decreasing = TRUE))

log_msg("Filtering to Malignant cells...")
malignant <- subset(obj, subset = cell_type == "Malignant")
log_msg(sprintf("Malignant cells: %d", ncol(malignant)))

if (ncol(malignant) == 0) {
  stop("No Malignant cells found — check cell_type column values above")
}

# Tag cancer type for consistency with other objects
malignant$cancer_type <- "Brain"

log_msg(sprintf("Saving to: %s", OUTPUT))
saveRDS(malignant, OUTPUT)
log_msg("Done!")
