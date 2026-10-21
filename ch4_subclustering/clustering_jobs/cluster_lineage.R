#!/usr/bin/env Rscript
# PCA + UMAP + multi-resolution clustering for one lineage object.
# Usage:
#   Rscript cluster_lineage.R --input <path.RDS> --output <path.RDS> \
#                             --lineage <name> [--n_pcs 30] [--min_dist 0.3]

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option("--input",   type = "character", help = "Input RDS path"),
  make_option("--output",  type = "character", help = "Output RDS path"),
  make_option("--lineage", type = "character", help = "Lineage name (for logging)"),
  make_option("--n_pcs",   type = "integer",   default = 30),
  make_option("--min_dist",type = "double",    default = 0.3),
  make_option("--default_res", type = "double", default = 0.50,
              help = "Default clustering resolution for Idents + UMAP [0.50]"),
  make_option("--summary_out", type = "character", default = NULL,
              help = "Optional path to write cluster-summary TSV")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat(sprintf("[%s] Loading %s\n", Sys.time(), opt$lineage))
obj <- readRDS(opt$input)
DefaultAssay(obj) <- "RNA"
cat(sprintf("[%s] Loaded: %d cells x %d genes\n", Sys.time(), ncol(obj), nrow(obj)))

# ── Join split layers (Seurat v5 merged objects) ──────────────────────────────
if (inherits(obj[["RNA"]], "Assay5")) {
  cat(sprintf("[%s] JoinLayers (Seurat v5 multi-layer object)\n", Sys.time()))
  obj <- JoinLayers(obj)
}

# ── Normalise + HVGs ──────────────────────────────────────────────────────────
cat(sprintf("[%s] NormalizeData + FindVariableFeatures (10k HVGs)\n", Sys.time()))
obj <- NormalizeData(obj, verbose = FALSE)
hvg <- tryCatch(
  FindVariableFeatures(obj, nfeatures = 10000, selection.method = "vst", verbose = FALSE),
  error = function(e) {
    cat(sprintf("[%s] VST failed (%s); falling back to mean.var.plot\n", Sys.time(), conditionMessage(e)))
    FindVariableFeatures(obj, nfeatures = 10000, selection.method = "mean.var.plot", verbose = FALSE)
  }
)
obj <- hvg

# ── PCA ───────────────────────────────────────────────────────────────────────
cat(sprintf("[%s] ScaleData + RunPCA (%d PCs)\n", Sys.time(), opt$n_pcs))
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, npcs = opt$n_pcs, verbose = FALSE)

# ── Choose embedding for neighbours: scANVI latent > Harmony > PCA ───────────
if ("scanvi" %in% names(obj@reductions)) {
  reduction_use <- "scanvi"
  dims_use      <- 1:ncol(Embeddings(obj, "scanvi"))
  cat(sprintf("[%s] Using scANVI latent embedding (%d dims)\n",
              Sys.time(), length(dims_use)))
} else if ("harmony" %in% names(obj@reductions)) {
  reduction_use <- "harmony"
  dims_use      <- 1:min(opt$n_pcs, ncol(Embeddings(obj, "harmony")))
  cat(sprintf("[%s] Using Harmony embedding (%d dims)\n",
              Sys.time(), length(dims_use)))
} else {
  reduction_use <- "pca"
  dims_use      <- 1:opt$n_pcs
  cat(sprintf("[%s] Using PCA embedding (%d dims)\n",
              Sys.time(), length(dims_use)))
}

# ── Neighbours + UMAP ─────────────────────────────────────────────────────────
cat(sprintf("[%s] FindNeighbors\n", Sys.time()))
obj <- FindNeighbors(obj, reduction = reduction_use, dims = dims_use, verbose = FALSE)

cat(sprintf("[%s] RunUMAP (n_components=2, min_dist=%.1f)\n", Sys.time(), opt$min_dist))
obj <- RunUMAP(obj, reduction = reduction_use, dims = dims_use,
               n.components = 2L, min.dist = opt$min_dist, verbose = FALSE)

# ── Leiden / Louvain clustering at 6 resolutions ──────────────────────────────
resolutions <- sort(unique(c(0.25, 0.30, 0.5, 0.75, 1.0, 2.0)))
cluster_cols <- character(length(resolutions))

for (i in seq_along(resolutions)) {
  res  <- resolutions[i]
  col  <- sprintf("louvain_res%.2f", res)
  cluster_cols[i] <- col
  cat(sprintf("[%s] FindClusters res=%.2f\n", Sys.time(), res))
  obj  <- FindClusters(obj, resolution = res, algorithm = 1,   # 1 = Louvain
                       cluster.name = col, verbose = FALSE)
}

# Set Idents to the requested default resolution
default_col <- sprintf("louvain_res%.2f", opt$default_res)
if (default_col %in% colnames(obj@meta.data)) {
  Idents(obj) <- default_col
  cat(sprintf("[%s] Idents set to %s\n", Sys.time(), default_col))
} else {
  warning("default_res column '", default_col, "' not found; Idents unchanged")
}

# ── Summary ───────────────────────────────────────────────────────────────────
n_cells   <- ncol(obj)
n_clusters <- sapply(cluster_cols, function(col) {
  if (col %in% colnames(obj@meta.data)) {
    nlevels(factor(obj@meta.data[[col]]))
  } else NA_integer_
})
umap_ok <- "umap" %in% names(obj@reductions)

cat("\n====== CLUSTER SUMMARY ======\n")
cat(sprintf("Lineage    : %s\n", opt$lineage))
cat(sprintf("Cells      : %s\n", format(n_cells, big.mark = ",")))
cat(sprintf("Embedding  : %s\n", reduction_use))
cat(sprintf("UMAP       : %s\n", if (umap_ok) "yes" else "NO"))
cat(sprintf("Default res: %.2f\n", opt$default_res))
for (i in seq_along(resolutions)) {
  marker <- if (resolutions[i] == opt$default_res) " <-- default" else ""
  cat(sprintf("  res=%.2f -> %d clusters (%s)%s\n",
              resolutions[i], n_clusters[i], cluster_cols[i], marker))
}
cat("=============================\n\n")

# Optionally write summary to TSV
if (!is.null(opt$summary_out)) {
  nc <- setNames(as.list(n_clusters), paste0("res_", sprintf("%.2f", resolutions)))
  df <- as.data.frame(c(
    list(lineage = opt$lineage, n_cells = n_cells,
         embedding = reduction_use, umap = umap_ok,
         default_res = opt$default_res),
    nc
  ), stringsAsFactors = FALSE)
  write.table(df, opt$summary_out, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("[%s] Summary written to %s\n", Sys.time(), opt$summary_out))
}

# ── UMAP plots ────────────────────────────────────────────────────────────────
plot_dir <- file.path(dirname(opt$output), "umap_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Detect a meta-type / lineage column to colour by
  meta_col <- intersect(
    c("Consensus_Cell_Type", "consensus_cell_type", "meta_type",
      "cell_type", "majority_cell_type", "louvain_res0.50"),
    colnames(obj@meta.data)
  )[1]

  top_cluster_col <- default_col   # use requested default resolution

  plots_to_save <- list(
    cluster   = list(col = top_cluster_col,
                     title = sprintf("%s — Louvain res %.2f", opt$lineage, opt$default_res)),
    meta_type = if (!is.na(meta_col))
                  list(col = meta_col,
                       title = sprintf("%s — %s", opt$lineage, meta_col))
                else NULL
  )

  for (pname in names(plots_to_save)) {
    pinfo <- plots_to_save[[pname]]
    if (is.null(pinfo)) next
    png_path <- file.path(plot_dir,
                          sprintf("%s_UMAP_%s.png", opt$lineage, pname))
    p <- DimPlot(obj, reduction = "umap", group.by = pinfo$col,
                 label = TRUE, repel = TRUE, pt.size = 0.1) +
         ggtitle(pinfo$title) +
         theme(legend.position = "right")
    ggsave(png_path, p, width = 10, height = 8, dpi = 150)
    cat(sprintf("[%s] Saved UMAP plot: %s\n", Sys.time(), png_path))
  }
} else {
  cat(sprintf("[%s] ggplot2 not available — skipping UMAP plots\n", Sys.time()))
}

# ── Save RDS ──────────────────────────────────────────────────────────────────
cat(sprintf("[%s] Saving to %s\n", Sys.time(), opt$output))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, opt$output)
cat(sprintf("[%s] Done.\n", Sys.time()))
