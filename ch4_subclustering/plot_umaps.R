#!/usr/bin/env Rscript
# plot_umaps.R
#
# Re-generate UMAP plots for a clustered lineage object.
# Overwrites existing PNGs in umap_plots/:
#
#   {lineage}_UMAP_cluster.png   — coloured by louvain_res0.50 cluster ID
#   {lineage}_UMAP_scimilarity.png — coloured by scimilarity_cluster_label
#                                    (majority-vote SCimilarity per cluster)
#
# Reads from: <input_dir>/<lineage>_clustered.RDS
# Writes to:  <input_dir>/umap_plots/

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option("--input_dir", type = "character", default = NULL,
    help = "Directory containing *_clustered.RDS files [required]"),
  make_option("--lineage", type = "character", default = NULL,
    help = "Lineage name (e.g. T-cells_Lung). Matches *_clustered.RDS [required]"),
  make_option("--cluster_col", type = "character", default = "louvain_res0.50",
    help = "Cluster column for _UMAP_cluster.png [louvain_res0.50]"),
  make_option("--scim_col", type = "character", default = "scimilarity_cluster_label",
    help = "SCimilarity majority-vote column for _UMAP_scimilarity.png [scimilarity_cluster_label]"),
  make_option("--overwrite", action = "store_true", default = TRUE,
    help = "Overwrite existing PNGs [default TRUE]"),
  make_option("--pt_size", type = "double", default = 0.1,
    help = "Point size for DimPlot [0.1]"),
  make_option("--width",  type = "integer", default = 12L,
    help = "Plot width in inches [12]"),
  make_option("--height", type = "integer", default = 9L,
    help = "Plot height in inches [9]"),
  make_option("--dpi", type = "integer", default = 150L,
    help = "Plot DPI [150]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_dir)) stop("--input_dir is required")
if (is.null(opt$lineage))   stop("--lineage is required")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

rds_path <- file.path(opt$input_dir, paste0(opt$lineage, "_clustered.RDS"))
if (!file.exists(rds_path)) stop("RDS not found: ", rds_path)

plot_dir <- file.path(opt$input_dir, "umap_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Loading: ", rds_path)
obj <- readRDS(rds_path)
log_msg("Loaded: ", ncol(obj), " cells")

# Verify UMAP exists
if (!"umap" %in% names(obj@reductions))
  stop("No 'umap' reduction found in object. Run cluster_lineage.R first.")

meta_cols <- colnames(obj@meta.data)

# ── Helper: save one DimPlot ──────────────────────────────────────────────────
save_umap <- function(col, png_name, title, n_colors = NULL) {
  if (!col %in% meta_cols) {
    log_msg("  WARNING: column '", col, "' not found — skipping ", png_name)
    return(invisible(NULL))
  }
  png_path <- file.path(plot_dir, png_name)
  log_msg("  Plotting ", png_name, " (group.by = '", col, "')")

  n_groups <- length(unique(obj@meta.data[[col]]))

  # Build a broad qualitative palette that scales with n_groups
  if (n_groups <= 8) {
    pal <- scales::hue_pal()(n_groups)
  } else if (n_groups <= 20) {
    pal <- c(RColorBrewer::brewer.pal(min(n_groups, 12), "Paired"),
             RColorBrewer::brewer.pal(min(n_groups - 12, 8), "Dark2"))[seq_len(n_groups)]
  } else {
    pal <- scales::hue_pal()(n_groups)
  }

  p <- DimPlot(
      obj,
      reduction  = "umap",
      group.by   = col,
      label      = TRUE,
      repel      = TRUE,
      pt.size    = opt$pt_size,
      label.size = 3
    ) +
    ggtitle(title) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "right",
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(size = 12, face = "bold")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1))

  ggsave(png_path, p,
         width  = opt$width,
         height = opt$height,
         dpi    = opt$dpi)
  log_msg("  Saved: ", png_path)
}

# ── Plot 1: cluster IDs ───────────────────────────────────────────────────────
cluster_col <- opt$cluster_col
if (!cluster_col %in% meta_cols) {
  # fallback to any louvain column
  cluster_col <- grep("^louvain_res", meta_cols, value = TRUE)[1]
  log_msg("'", opt$cluster_col, "' not found, using: ", cluster_col)
}

save_umap(
  col      = cluster_col,
  png_name = sprintf("%s_UMAP_cluster.png", opt$lineage),
  title    = sprintf("%s — Louvain clusters (%s)", opt$lineage, cluster_col)
)

# ── Plot 2: SCimilarity majority label ────────────────────────────────────────
scim_col <- opt$scim_col
if (!scim_col %in% meta_cols) {
  # fallback: try scimilarity_label
  for (candidate in c("scimilarity_label", "scimilarity_celltype",
                      "scanvi_pred", "Meta.scimilarity")) {
    if (candidate %in% meta_cols) { scim_col <- candidate; break }
  }
}

save_umap(
  col      = scim_col,
  png_name = sprintf("%s_UMAP_meta_type.png", opt$lineage),
  title    = sprintf("%s — SCimilarity cluster label (%s)", opt$lineage, scim_col)
)

log_msg("Done: ", opt$lineage)
