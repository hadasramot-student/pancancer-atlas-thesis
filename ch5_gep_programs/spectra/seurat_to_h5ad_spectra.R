#!/usr/bin/env Rscript
# seurat_to_h5ad_spectra.R
#
# Ch.5 (GEP programs) — Spectra pipeline step A. Convert ONE clustered immune
# Seurat (.RDS, from Ch.4 clustering) to .h5ad for Spectra.
#   X            = log-normalised expression (RNA 'data' layer)   [cells x genes]
#   layers$counts = raw counts (RNA 'counts' layer)
#   obs          = a curated subset of metadata (incl. cell_type_key)
#
# Adds (a) obs-column selection and (b) optional stratified per-file
# subsampling for the huge lineages (T-cells / Myeloid).
#
# Usage:
#   Rscript seurat_to_h5ad_spectra.R \
#     --input_rds  /path/to/clustered_lineages/T-cells_Lung_clustered.RDS \
#     --output_h5ad /path/to/clustered_lineages/h5ad/T-cells_Lung_clustered.h5ad \
#     --celltype_key scimilarity_celltype \
#     --max_cells 40000          # 0 = keep all cells in this file

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
if (nzchar(userlib)) .libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(anndata)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input_rds",   type = "character", default = NULL),
  make_option("--output_h5ad",  type = "character", default = NULL),
  make_option("--celltype_key", type = "character", default = "scimilarity_celltype"),
  make_option("--max_cells",    type = "integer",   default = 0L,
              help = "Stratified subsample target for THIS file (0 = keep all)"),
  make_option("--seed",         type = "integer",   default = 0L)
)))
if (is.null(opt$input_rds))   stop("Please provide --input_rds")
if (is.null(opt$output_h5ad)) stop("Please provide --output_h5ad")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

# obs columns to keep (only those present are exported)
keep_cols <- c(
  opt$celltype_key,
  "Consensus_Cell_Type", "cell_type", "scimilarity_celltype",
  "scimilarity_vs2nd", "project_name", "cancer_type", "sample",
  "orig.ident", "nCount_RNA", "nFeature_RNA", "pct_counts_mt",
  "seurat_clusters", "doublet_score"
)

log_msg(paste("Loading:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Loaded: %d cells x %d genes", ncol(obj), nrow(obj)))

if (!("RNA" %in% Assays(obj))) stop("RNA assay not found")
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  log_msg("JoinLayers on Assay5 RNA")
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

# --- optional stratified per-file subsample ---
if (!is.null(opt$max_cells) && opt$max_cells > 0 && ncol(obj) > opt$max_cells) {
  set.seed(opt$seed)
  ck <- opt$celltype_key
  if (!(ck %in% colnames(obj@meta.data))) {
    log_msg(paste0("celltype_key '", ck, "' absent; random subsample"))
    keep <- sample(colnames(obj), opt$max_cells)
  } else {
    grp  <- as.character(obj@meta.data[[ck]])
    frac <- opt$max_cells / ncol(obj)
    keep <- unlist(lapply(split(colnames(obj), grp), function(cells) {
      k <- max(1, round(length(cells) * frac)); k <- min(k, length(cells))
      sample(cells, k)
    }), use.names = FALSE)
  }
  log_msg(sprintf("Subsampled %d -> %d cells", ncol(obj), length(keep)))
  obj <- subset(obj, cells = keep)
}

log_msg("Extracting counts + log-normalised data")
counts_mat <- t(GetAssayData(obj, assay = "RNA", layer = "counts"))
data_mat   <- t(GetAssayData(obj, assay = "RNA", layer = "data"))

# --- obs ---
meta <- obj@meta.data
present <- intersect(keep_cols, colnames(meta))
meta <- meta[, present, drop = FALSE]
for (col in colnames(meta)) {
  x <- meta[[col]]
  if (is.factor(x) || inherits(x, "Date") || inherits(x, "POSIXt")) meta[[col]] <- as.character(x)
  else if (is.list(x)) meta[[col]] <- vapply(x, function(v) paste(unlist(v), collapse = ";"), character(1))
  else if (!is.numeric(x) && !is.logical(x) && !is.character(x)) meta[[col]] <- as.character(x)
  if (is.character(x)) meta[[col]][is.na(meta[[col]])] <- "Unknown"
}
log_msg("obs columns exported:")
print(sapply(meta, class))

gene_names <- colnames(counts_mat)
ad <- AnnData(
  X      = data_mat,
  obs    = meta,
  var    = data.frame(gene_name = gene_names, row.names = gene_names),
  layers = list(counts = counts_mat)
)
ad$obs_names <- rownames(meta)
ad$var_names <- gene_names

dir.create(dirname(opt$output_h5ad), recursive = TRUE, showWarnings = FALSE)
log_msg(paste("Writing:", opt$output_h5ad))
ad$write_h5ad(opt$output_h5ad)
log_msg("Done")
