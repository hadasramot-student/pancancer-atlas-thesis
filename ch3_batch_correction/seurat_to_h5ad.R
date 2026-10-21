#!/usr/bin/env Rscript
# seurat_to_h5ad.R
#
# Ch.3 (Batch correction) — converts a merged Seurat object to h5ad so it
# can be handed to scANVI (scanvi_batch_correction.py / export_scanvi_to_csv.py).
#
# Inputs:  --input_rds  (merged Seurat object, raw counts)
# Outputs: --output_h5ad (AnnData with counts layer + obs metadata)
#
# Usage:
#   Rscript seurat_to_h5ad.R --input_rds obj_merged.RDS --output_h5ad obj_merged.h5ad

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(anndata)
})

# --------------------------------------------------
# CLI arguments
# --------------------------------------------------
option_list <- list(
  make_option(c("--input_rds"), type = "character", default = NULL,
              help = "Path to input Seurat RDS file (required)", metavar = "character"),
  
  make_option(c("--output_h5ad"), type = "character", default = NULL,
              help = "Path to output h5ad file (required)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_rds))  stop("Please provide --input_rds")
if (is.null(opt$output_h5ad)) stop("Please provide --output_h5ad")

# --------------------------------------------------
# Logger
# --------------------------------------------------
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

# --------------------------------------------------
# Load Seurat object
# --------------------------------------------------
log_msg(paste("Loading Seurat object from:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Loaded object: %d cells x %d genes", ncol(obj), nrow(obj)))

# --------------------------------------------------
# Ensure RNA assay exists and join layers (Seurat v5)
# --------------------------------------------------
if (!("RNA" %in% Assays(obj))) {
  stop("RNA assay not found in Seurat object")
}
DefaultAssay(obj) <- "RNA"

# Always attempt JoinLayers — safe to call even on already-joined objects
log_msg("Joining RNA layers (safe to call even if already joined)")
obj <- tryCatch({
  JoinLayers(obj)
}, error = function(e) {
  log_msg(sprintf("JoinLayers failed (may already be joined): %s", e$message))
  obj  # return unchanged if it fails
})

log_msg("RNA layers ready for conversion")

# --------------------------------------------------
# Extract counts and log-normalised data
# --------------------------------------------------
log_msg("Extracting counts (RNA, layer = 'counts')")
counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")

log_msg("Extracting log-normalised data (RNA, layer = 'data')")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

# Defensive: for single-project objects, the RNA layer can still be named
# "counts.<project>" after JoinLayers (nothing to join with only one layer),
# and Seurat's internal layer resolver can return a matrix whose rownames
# don't survive that path -- always re-stamp gene names from the object's
# own authoritative Features(), not whatever the fetched matrix happened to
# carry, so var_names can never silently end up as a numeric placeholder.
stopifnot(nrow(counts_mat) == length(Features(obj)))
rownames(counts_mat) <- Features(obj)
rownames(data_mat)   <- Features(obj)

# AnnData expects cells x genes (Seurat is genes x cells)
counts_mat <- t(counts_mat)
data_mat   <- t(data_mat)

# --------------------------------------------------
# Build and sanitize cell metadata (obs)
# --------------------------------------------------
log_msg("Preparing cell metadata (obs)")
meta <- obj@meta.data

# Drop clearly problematic or unused columns (example: empty `patient`)
if ("patient" %in% colnames(meta)) {
  log_msg("Dropping 'patient' column from metadata (empty / problematic)")
  meta$patient <- NULL
}

if ("paper" %in% colnames(meta)) {
  log_msg("Dropping 'paper' column from metadata (problematic for h5ad)")
  meta$paper <- NULL
}

# Drop all-NA columns outright (carry zero information, and casting
# them to character has been observed to intermittently fail during
# the reticulate/anndata write step -- e.g. TypeError: "Can't
# implicitly convert non-string objects to strings" on an all-None
# object column -- so remove the ambiguity instead of working around it)
all_na_cols <- vapply(meta, function(x) all(is.na(x)), logical(1))
if (any(all_na_cols)) {
  log_msg(paste(
    "Dropping all-NA columns (empty / would risk a flaky anndata write):",
    paste(names(which(all_na_cols)), collapse = ", ")
  ))
  meta <- meta[, !all_na_cols, drop = FALSE]
}

# Force each column to be numeric, logical, or character
for (col in colnames(meta)) {
  x <- meta[[col]]
  
  if (is.factor(x)) {
    meta[[col]] <- as.character(x)
  } else if (is.list(x) || is.data.frame(x)) {
    meta[[col]] <- sapply(x, function(v) paste(unlist(v), collapse = ";"))
  } else if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    meta[[col]] <- as.character(x)
  } else if (!is.numeric(x) && !is.logical(x) && !is.character(x)) {
    meta[[col]] <- as.character(x)
  }
}

# Final enforcement & attribute cleanup
for (col in colnames(meta)) {
  x <- meta[[col]]
  if (!is.numeric(x) && !is.logical(x) && !is.character(x)) {
    meta[[col]] <- as.character(x)
  }
  attr(meta[[col]], "class") <- class(meta[[col]])[1]
}

log_msg("Metadata column classes after sanitization:")
print(sapply(meta, class))

# --------------------------------------------------
# Build AnnData object
# --------------------------------------------------
log_msg("Building AnnData object")

cell_names <- rownames(meta)
gene_names <- colnames(counts_mat)

var_df <- data.frame(
  gene_name = gene_names,
  row.names = gene_names
)

# Create AnnData with log-normalized data in X and raw counts as a layer
# counts layer is required for seurat_v3 HVG selection in scVI/scANVI
ad <- AnnData(
  X      = data_mat,                  # log-normalised expression
  obs    = meta,                      # cell metadata
  var    = var_df,                    # gene metadata (minimal)
  layers = list(counts = counts_mat)  # raw counts for HVG selection
)

ad$obs_names <- cell_names
ad$var_names <- gene_names

log_msg("AnnData object ready (X = log-normalized data, layers['counts'] = raw counts)")

# --------------------------------------------------
# Write h5ad
# --------------------------------------------------
log_msg(paste("Saving h5ad to:", opt$output_h5ad))
ad$write_h5ad(opt$output_h5ad)

log_msg(paste("Saved:", opt$output_h5ad))
log_msg("Done")