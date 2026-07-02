#!/usr/bin/env Rscript
# scanvi_into_seurat.R
#
# Ch.3 (Batch correction) — re-imports the scANVI latent embedding + labels
# (from export_scanvi_to_csv.py) back into the Seurat object as a
# DimReduc, alongside the CCA/MNN/ComBat/Harmony reductions.
#
# Inputs:  --seurat (RDS), --latent_csv, --labels_csv (from export_scanvi_to_csv.py)
# Outputs: --output (RDS with a "scanvi" DimReduc + scanvi_pred metadata column)
#
# Usage:
#   Rscript scanvi_into_seurat.R --seurat obj.RDS \
#     --latent_csv obj_scanvi_latent.csv --labels_csv obj_scanvi_labels.csv \
#     --output obj_with_scanvi.RDS

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option(c("-s", "--seurat"), type = "character",
              help = "Input Seurat .RDS file"),
  make_option(c("-l", "--latent_csv"), type = "character",
              help = "CSV with scANVI latent, cols: cell, SCANVI_1..N"),
  make_option(c("-a", "--labels_csv"), type = "character",
              help = "CSV with scANVI labels, cols: cell, label_col"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output Seurat .RDS file"),
  make_option(c("--label_col"), type = "character", default = "scanvi_pred",
              help = "Column in labels_csv with labels [default %default]"),
  make_option(c("--meta_col"), type = "character", default = "scanvi_pred",
              help = "Metadata column name to create in Seurat [default %default]"),
  make_option(c("--reduction"), type = "character", default = "scanvi",
              help = "Name of DimReduc to create [default %default]"),
  make_option(c("--overwrite"), action = "store_true", default = FALSE,
              help = "Overwrite output RDS if it exists")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$seurat) || is.null(opt$latent_csv) || is.null(opt$labels_csv) || is.null(opt$output)) {
  stop("Please provide --seurat, --latent_csv, --labels_csv and --output")
}

if (!file.exists(opt$seurat)) stop("Seurat file not found: ", opt$seurat)
if (!file.exists(opt$latent_csv)) stop("Latent CSV not found: ", opt$latent_csv)
if (!file.exists(opt$labels_csv)) stop("Labels CSV not found: ", opt$labels_csv)
if (file.exists(opt$output) && !opt$overwrite) {
  stop("Output file exists: ", opt$output, " (use --overwrite to replace)")
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

log_msg("Loading Seurat object: ", opt$seurat)
seu <- readRDS(opt$seurat)
cells_seu <- colnames(seu)

## --- Latent CSV ---
log_msg("Reading latent CSV: ", opt$latent_csv)
latent_df <- read.csv(opt$latent_csv, stringsAsFactors = FALSE)

if (!"cell" %in% colnames(latent_df)) {
  stop("Latent CSV must have a 'cell' column.")
}

# keep only overlapping cells
latent_df <- latent_df[latent_df$cell %in% cells_seu, , drop = FALSE]
if (nrow(latent_df) == 0) {
  stop("No overlapping cells between latent CSV and Seurat object.")
}

# numeric latent matrix (drop 'cell' column)
latent_mat <- as.matrix(latent_df[, setdiff(colnames(latent_df), "cell"), drop = FALSE])

# align rows to Seurat cell order
latent_mat <- latent_mat[match(cells_seu, latent_df$cell), , drop = FALSE]
rownames(latent_mat) <- cells_seu

# ensure reasonable column names
if (is.null(colnames(latent_mat)) || any(colnames(latent_mat) == "")) {
  colnames(latent_mat) <- paste0(toupper(opt$reduction), "_", seq_len(ncol(latent_mat)))
}

log_msg("Latent matrix dimensions: ", paste(dim(latent_mat), collapse = " x "))

## --- Create DimReduc ---
log_msg("Creating DimReduc '", opt$reduction, "'")
seu[[opt$reduction]] <- CreateDimReducObject(
  embeddings = latent_mat,
  key = paste0(toupper(opt$reduction), "_"),
  assay = DefaultAssay(seu)
)

## --- Labels CSV ---
log_msg("Reading labels CSV: ", opt$labels_csv)
labels_df <- read.csv(opt$labels_csv, stringsAsFactors = FALSE)

if (!"cell" %in% colnames(labels_df)) {
  stop("Labels CSV must have a 'cell' column.")
}
if (!(opt$label_col %in% colnames(labels_df))) {
  stop("Label column '", opt$label_col, "' not found in labels CSV.")
}

labels_vec <- labels_df[[opt$label_col]]
names(labels_vec) <- labels_df$cell

labels_aligned <- labels_vec[cells_seu]
if (any(is.na(labels_aligned))) {
  missing_n <- sum(is.na(labels_aligned))
  warning("Missing labels for ", missing_n, " cells; they will be NA.")
}

seu@meta.data[[opt$meta_col]] <- labels_aligned

log_msg("Added metadata column '", opt$meta_col, "' to Seurat object.")
log_msg("Saving updated Seurat object to: ", opt$output)
saveRDS(seu, file = opt$output)
log_msg("Done.")

# Rscript scanvi_into_seurat.R \
# --seurat $ATLAS_DATA_DIR/merged_by_query/Hematopoietic_Stem_and_Progenitor_cells_Breast_merged_batchco.RDS \
# --latent_csv $ATLAS_DATA_DIR/merged_by_query/Hematopoietic_Stem_and_Progenitor_cells_Breast_scanvi_latent.csv \
# --labels_csv $ATLAS_DATA_DIR/merged_by_query/Hematopoietic_Stem_and_Progenitor_cells_Breast_scanvi_labels.csv \
# --output $ATLAS_DATA_DIR/merged_by_query/Hematopoietic_Stem_and_Progenitor_cells_Breast_merged_batchco_with_scanvi.RDS \
# --label_col scanvi_pred \
# --meta_col scanvi_pred \
# --reduction scanvi \
# --overwrite