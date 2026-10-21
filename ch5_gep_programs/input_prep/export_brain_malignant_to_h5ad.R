#!/usr/bin/env Rscript
log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

RDS  <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Brain/Epithelial_Brain_merged_full_scanvi_final.RDS"
OUT  <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Brain/Epithelial_Brain_merged_full.h5ad"

log_msg("Loading Brain Malignant RDS...")
obj <- readRDS(RDS)
log_msg(sprintf("Cells: %d  Genes: %d", ncol(obj), nrow(obj)))

# Export via SeuratDisk
suppressPackageStartupMessages(library(SeuratDisk))

tmp_h5seurat <- sub("\\.h5ad$", ".h5seurat", OUT)
log_msg("Converting to h5seurat...")
SaveH5Seurat(obj, filename = tmp_h5seurat, overwrite = TRUE)

log_msg("Converting to h5ad...")
Convert(tmp_h5seurat, dest = "h5ad", overwrite = TRUE)

file.remove(tmp_h5seurat)
log_msg(sprintf("Done. Saved: %s", OUT))
