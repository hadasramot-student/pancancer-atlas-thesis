#!/usr/bin/env Rscript
# batch4_harmony.R
#
# Ch.3 (Batch correction) — Harmony integration.
# Inputs:  --input_rds  (merged Seurat object, raw counts)
# Outputs: --output_rds (adds harmony reduction + umap.harmony + harmony_clusters)
#
# Usage:
#   Rscript batch4_harmony.R --input_rds obj_merged.RDS \
#     --output_rds obj_merged_harmony.RDS --batch_col project_name

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
})

option_list <- list(
  make_option(c("--input_rds"),  type="character", default=NULL),
  make_option(c("--output_rds"), type="character", default=NULL),
  make_option(c("--batch_col"),  type="character", default="project_name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_rds)) stop("Please provide --input_rds")
if (is.null(opt$output_rds)) {
  opt$output_rds <- sub("\\.RDS$", "_harmony.RDS", opt$input_rds)
}

batch_col <- opt$batch_col

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

log_msg(paste("Loading:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Loaded: %d cells, %d genes", ncol(obj), nrow(obj)))

# --------------------------------------------------
# Harmony requires PCA to already exist
# --------------------------------------------------
DefaultAssay(obj) <- "RNA"

log_msg("Normalizing and finding variable features")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)

max_pcs <- min(50L, floor(ncol(obj) / 2L))
log_msg(sprintf("Running PCA with %d PCs", max_pcs))
obj <- RunPCA(obj, npcs = max_pcs, verbose = FALSE)

# --------------------------------------------------
# Harmony integration
# --------------------------------------------------
log_msg("Running Harmony integration")
obj <- RunHarmony(
  obj,
  group.by.vars  = batch_col,
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = TRUE
)

use_dims <- min(30L, max_pcs)

log_msg("Running FindNeighbors on Harmony embedding")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:use_dims)

log_msg("Running FindClusters")
obj <- FindClusters(obj, resolution = 0.6, cluster.name = "harmony_clusters")

log_msg("Running UMAP on Harmony embedding")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:use_dims,
               reduction.name = "umap.harmony")

log_msg("Finished Harmony integration")

log_msg(paste("Saving:", opt$output_rds))
saveRDS(obj, opt$output_rds)
log_msg("Done")