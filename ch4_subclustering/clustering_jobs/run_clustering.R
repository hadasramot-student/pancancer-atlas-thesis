#!/usr/bin/env Rscript
# run_clustering.R
#
# Ch.4 (Subclustering) — per-lineage HVG selection -> PCA (30 PCs, default)
# -> FindNeighbors -> RunUMAP -> Louvain/Leiden clustering at 5 resolutions
# (0.25, 0.5, 0.75, 1.0, 2.0) for one lineage Seurat object.
#
# Inputs:  --rds (a merged, batch-corrected lineage object)
# Outputs: <outdir>/<lineage>_clustered.RDS, <lineage>_cluster_summary.txt
#
# Usage: Rscript run_clustering.R --rds <path> --lineage <name> [--outdir <dir>]

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(future)
})

option_list <- list(
  make_option("--rds",     type = "character", help = "Path to input .RDS Seurat object"),
  make_option("--lineage", type = "character", help = "Lineage name label (used in output filename)"),
  make_option("--outdir",  type = "character", default = NULL,
              help = "Output directory (default: same as input RDS)"),
  make_option("--n_pcs",   type = "integer",   default = 30L,
              help = "Number of PCs for PCA [default: 30]"),
  make_option("--min_dist", type = "double",   default = 0.3,
              help = "UMAP min_dist [default: 0.3]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$rds),     file.exists(opt$rds))
stopifnot(!is.null(opt$lineage))

outdir <- if (!is.null(opt$outdir)) opt$outdir else dirname(opt$rds)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(outdir, paste0(opt$lineage, "_clustered.RDS"))
out_summary <- file.path(outdir, paste0(opt$lineage, "_cluster_summary.txt"))

cat(sprintf("[%s] Loading: %s\n", Sys.time(), opt$rds))
obj <- readRDS(opt$rds)
cat(sprintf("[%s] Loaded: %d cells x %d genes\n", Sys.time(), ncol(obj), nrow(obj)))

# ── Join layers if Seurat v5 split layers ─────────────────────────────────
if ("RNA" %in% names(obj@assays)) {
  lyrs <- SeuratObject::Layers(obj[["RNA"]])
  if (length(lyrs) > 1 && any(grepl("^counts\\.", lyrs))) {
    cat(sprintf("[%s] Joining split RNA layers...\n", Sys.time()))
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  }
}
DefaultAssay(obj) <- "RNA"

# ── Normalise + HVGs (only if not already done) ───────────────────────────
if (is.null(VariableFeatures(obj)) || length(VariableFeatures(obj)) == 0) {
  cat(sprintf("[%s] NormalizeData...\n", Sys.time()))
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  cat(sprintf("[%s] FindVariableFeatures (10000)...\n", Sys.time()))
  obj <- FindVariableFeatures(obj, nfeatures = 10000, verbose = FALSE)
  cat(sprintf("[%s] ScaleData...\n", Sys.time()))
  obj <- ScaleData(obj, verbose = FALSE)
} else {
  cat(sprintf("[%s] HVGs already present (%d). Skipping Normalize/ScaleData.\n",
              Sys.time(), length(VariableFeatures(obj))))
  # Ensure data layer exists; if not, normalise
  if (!"data" %in% SeuratObject::Layers(obj[["RNA"]])) {
    cat(sprintf("[%s] No 'data' layer found — running NormalizeData + ScaleData.\n", Sys.time()))
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
  }
}

# ── PCA ───────────────────────────────────────────────────────────────────
cat(sprintf("[%s] RunPCA (%d PCs)...\n", Sys.time(), opt$n_pcs))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = opt$n_pcs, verbose = FALSE)

# ── Decide embedding for neighbors: scANVI latent if present, else PCA ───
use_reduction <- "pca"
use_dims      <- seq_len(opt$n_pcs)

scanvi_reductions <- c("scanvi", "scANVI", "SCANVI", "latent.scanvi")
for (red_name in scanvi_reductions) {
  if (red_name %in% Reductions(obj)) {
    n_latent <- ncol(Embeddings(obj, reduction = red_name))
    use_reduction <- red_name
    use_dims      <- seq_len(n_latent)
    cat(sprintf("[%s] Using scANVI embedding '%s' (%d dims) for neighbors.\n",
                Sys.time(), red_name, n_latent))
    break
  }
}
if (use_reduction == "pca") {
  cat(sprintf("[%s] No scANVI embedding found — using PCA (%d dims) for neighbors.\n",
              Sys.time(), opt$n_pcs))
}

# ── Neighbors ─────────────────────────────────────────────────────────────
options(future.globals.maxSize = 16 * 1024^3)
plan("sequential")
cat(sprintf("[%s] FindNeighbors (reduction=%s, dims=1:%d)...\n",
            Sys.time(), use_reduction, max(use_dims)))
obj <- FindNeighbors(obj, reduction = use_reduction, dims = use_dims, verbose = FALSE)

# ── UMAP ──────────────────────────────────────────────────────────────────
cat(sprintf("[%s] RunUMAP (n_components=2, min_dist=%.2f)...\n", Sys.time(), opt$min_dist))
obj <- RunUMAP(obj, reduction = use_reduction, dims = use_dims,
               n.components = 2L, min.dist = opt$min_dist, verbose = FALSE)

# ── Clustering at 5 resolutions ───────────────────────────────────────────
resolutions <- c(0.25, 0.5, 0.75, 1.0, 2.0)
cluster_counts <- integer(length(resolutions))
names(cluster_counts) <- as.character(resolutions)

for (i in seq_along(resolutions)) {
  res <- resolutions[i]
  col <- paste0("leiden_res", gsub("\\.", "_", as.character(res)))
  cat(sprintf("[%s] FindClusters (resolution=%.2f -> column '%s')...\n",
              Sys.time(), res, col))
  obj <- FindClusters(obj, resolution = res, cluster.name = col,
                      algorithm = 4,   # Leiden
                      verbose = FALSE)
  cluster_counts[i] <- length(unique(obj@meta.data[[col]]))
}

umap_ok <- "umap" %in% Reductions(obj)

# ── Save ──────────────────────────────────────────────────────────────────
cat(sprintf("[%s] Saving to %s...\n", Sys.time(), out_rds))
saveRDS(obj, out_rds)

# ── Summary ───────────────────────────────────────────────────────────────
summary_lines <- c(
  paste0("Lineage:         ", opt$lineage),
  paste0("Input RDS:       ", opt$rds),
  paste0("Cells:           ", ncol(obj)),
  paste0("Embedding used:  ", use_reduction),
  paste0("UMAP computed:   ", ifelse(umap_ok, "Yes", "No")),
  "Cluster counts by resolution:",
  paste(sprintf("  res=%.2f : %d clusters", resolutions, cluster_counts), collapse = "\n"),
  paste0("Output RDS:      ", out_rds)
)
writeLines(summary_lines, out_summary)
cat("\n===== SUMMARY =====\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("===================\n")
