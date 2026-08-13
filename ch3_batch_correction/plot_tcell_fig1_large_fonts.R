#!/usr/bin/env Rscript
# plot_tcell_fig1_large_fonts.R
#
# One-off re-render of the T-cells "Batch correction overview" figure
# (thesis Fig. 15) with larger typography for print legibility at 6in
# display width. Derived from plot_publication_figures.R (Figure 1 logic
# only, T-cells only). Font sizes for panel titles, legend text/title,
# and figure title/subtitle are scaled ~1.5-2x relative to panel size;
# canvas is widened/heightened modestly so larger legends don't clip,
# keeping net relative increase in the 1.5-2x range. No data or plot-type
# changes vs. the original figure.
#
# Inputs:  ATLAS_DATA_DIR/merged_by_query/Breast/T-cells_Breast_merged_step4_harmony_dev015_nocca.RDS
#          (+ matching scANVI latent CSV)
# Outputs: ATLAS_DATA_DIR/merged_by_query/publication_figures/T-cells_Breast_fig1_overview_large_fonts.{pdf,png}
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript plot_tcell_fig1_large_fonts.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

BASE      <- file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "merged_by_query")
BASE_CELL <- file.path(BASE, "Breast")
OUT_DIR   <- file.path(BASE, "publication_figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cel <- "T-cells"
TOP_N_CELLTYPES <- 10L

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---------------------------------------------------------------
# Colour helpers (unchanged from original script)
# ---------------------------------------------------------------
discrete_colors <- function(n) {
  base_cols <- c(
    "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3",
    "#A65628","#F781BF","#66C2A5","#FC8D62","#8DA0CB",
    "#E78AC3","#A6D854","#FFD92F","#1B9E77","#D95F02",
    "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D",
    "#1F78B4","#B2DF8A","#33A02C","#FB9A99","#FDBF6F",
    "#CAB2D6","#6A3D9A","#B15928","#FF7F00","#FFFF99"
  )
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  colorRampPalette(base_cols)(n)
}

make_color_map <- function(vals) {
  lvls <- sort(unique(na.omit(as.character(vals))))
  setNames(discrete_colors(length(lvls)), lvls)
}

# ---------------------------------------------------------------
# Ensure UMAP; compute from reduction if not already present
# ---------------------------------------------------------------
ensure_umap <- function(obj, from_reduction, umap_name,
                        n_neighbors = 30L, min_dist = 0.3) {
  if (umap_name %in% names(obj@reductions)) return(obj)
  log_msg("  Computing UMAP '", umap_name, "' from '", from_reduction, "'...")
  ndim <- ncol(Embeddings(obj, from_reduction))
  obj  <- RunUMAP(obj,
                  reduction      = from_reduction,
                  dims           = seq_len(min(50L, ndim)),
                  reduction.name = umap_name,
                  reduction.key  = paste0(toupper(gsub("[^A-Za-z0-9]", "", umap_name)), "_"),
                  n.neighbors    = n_neighbors,
                  min.dist       = min_dist,
                  verbose        = FALSE)
  obj
}

# ---------------------------------------------------------------
# Single UMAP panel (raw ggplot2) — fonts scaled ~1.7-2x vs. original
# (base_size 9->18, panel title 10->18, legend.text 7->14,
#  legend.title 8->16, legend.key.size 0.3cm->0.55cm,
#  legend point size 3->5.5)
# ---------------------------------------------------------------
umap_panel <- function(coords, vals, color_map,
                       title        = "",
                       pt_size      = 0.3,
                       alpha        = 0.5,
                       guide_type   = "none",
                       legend_title = "",
                       legend_ncol  = 1L,
                       xlim_range   = NULL,
                       ylim_range   = NULL) {

  idx <- sample(nrow(coords))
  df  <- data.frame(
    x   = coords[idx, 1],
    y   = coords[idx, 2],
    col = factor(as.character(vals)[idx], levels = names(color_map)),
    stringsAsFactors = FALSE
  )

  leg <- if (guide_type == "right") {
    guide_legend(
      override.aes = list(size = 5, alpha = 1),
      ncol         = legend_ncol,
      keywidth     = unit(0.5, "cm"),
      keyheight    = unit(0.5, "cm"),
      title.position = "top"
    )
  } else "none"

  p <- ggplot(df, aes(x, y, colour = col)) +
    geom_point(size = pt_size, alpha = alpha, stroke = 0, shape = 16) +
    scale_colour_manual(values   = color_map,
                        name     = legend_title,
                        na.value = "grey80",
                        guide    = leg) +
    labs(title = title) +
    theme_void(base_size = 15) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 17,
                                     margin = margin(b = 5)),
      legend.position = if (guide_type == "right") "right" else "none",
      legend.text     = element_text(size = 12),
      legend.title    = element_text(size = 14, face = "bold"),
      legend.key.size = unit(0.5, "cm"),
      legend.spacing.y = unit(0.05, "cm"),
      plot.margin     = margin(4, 4, 4, 4)
    )

  if (!is.null(xlim_range)) p <- p + coord_cartesian(xlim = xlim_range, ylim = ylim_range)
  p
}

# ---------------------------------------------------------------
# Load T-cells step4 object + scANVI latent (current on-disk layout:
# per-cancer subdir Breast/, not flat merged_by_query/)
# ---------------------------------------------------------------
log_msg("=== ", cel, " ===")

step4_path <- file.path(BASE_CELL,
  paste0(cel, "_Breast_merged_step4_harmony_dev015_nocca.RDS"))
if (!file.exists(step4_path)) stop("step4 RDS not found: ", step4_path)

log_msg("  Loading step4 object...")
obj <- readRDS(step4_path)
log_msg("  Cells: ", ncol(obj), "  |  Reductions: ",
        paste(names(obj@reductions), collapse = ", "))

obj <- ensure_umap(obj, "pca",        "umap.pca")
obj <- ensure_umap(obj, "mnn",        "umap.mnn")
obj <- ensure_umap(obj, "pca.combat", "umap.combat")
obj <- ensure_umap(obj, "harmony",    "umap.harmony")

lat_csv <- file.path(BASE_CELL,
  paste0(cel, "_Breast_merged_step6_scanvi_dev015_nocca_scanvi_latent.csv"))

has_scanvi <- FALSE
if (!file.exists(lat_csv)) {
  log_msg("  scANVI latent CSV not found — 4-column figure")
} else {
  log_msg("  Loading scANVI latent...")
  lat_df  <- read.csv(lat_csv, row.names = 1, stringsAsFactors = FALSE)
  common  <- intersect(colnames(obj), rownames(lat_df))
  log_msg("  Cells in common: ", length(common), " / ", ncol(obj))

  if (length(common) < 0.85 * ncol(obj)) {
    log_msg("  Too few common cells — skipping scANVI")
  } else {
    obj     <- obj[, common]
    lat_mat <- as.matrix(lat_df[common, , drop = FALSE])
    colnames(lat_mat) <- paste0("SCANVI_", seq_len(ncol(lat_mat)))
    obj[["scanvi"]] <- CreateDimReducObject(
      embeddings = lat_mat,
      key        = "SCANVI_",
      assay      = DefaultAssay(obj)
    )
    obj <- ensure_umap(obj, "scanvi", "umap.scanvi")
    has_scanvi <- TRUE
  }
}

n_cells <- ncol(obj)
pt_size <- min(2.0, max(0.05, 120000 / n_cells))

methods <- list(
  list(umap = "umap.pca",     label = "Uncorrected"),
  list(umap = "umap.mnn",     label = "MNN"),
  list(umap = "umap.combat",  label = "ComBat"),
  list(umap = "umap.harmony", label = "Harmony")
)
if (has_scanvi) methods[[5]] <- list(umap = "umap.scanvi", label = "scANVI")
n_cols <- length(methods)

ct_broad  <- as.character(obj$Meta.scimilarity)
broad_map <- make_color_map(ct_broad)
n_broad   <- length(broad_map)

batch_vals <- as.character(obj$project_name)
batch_map  <- make_color_map(batch_vals)
n_studies  <- length(batch_map)
batch_ncol <- if (n_studies > 20) 2L else 1L

# figure canvas: height is sized dynamically off the *actual* number of
# Study-legend rows so the enlarged legend text/keys never clip, no matter
# how many studies are in the object. Constants below (per_row_in, header
# padding, safety margin) were back-calculated from a first attempt that
# clipped 3 of 22 legend rows at fig_h=8.5in with legend.text=14/key=0.55cm
# -- i.e. don't shrink this without re-validating against a render.
legend_rows      <- ceiling(n_studies / batch_ncol)
per_row_in       <- 0.20   # empirical, at legend.text=12pt / key=0.5cm
legend_block_in  <- legend_rows * per_row_in + 0.35   # + legend title/header
row_avail_in     <- legend_block_in * 1.25            # 25% safety margin
title_block_in   <- 1.8                               # title+subtitle+col headers

fig_w <- 3.2 * n_cols + 3.2
fig_h <- title_block_in + 2 * row_avail_in + 1.0       # +1in extra pad

log_msg("  Building Figure 1 (large fonts)...")

row1 <- vector("list", n_cols)
row2 <- vector("list", n_cols)

for (i in seq_len(n_cols)) {
  m      <- methods[[i]]
  coords <- Embeddings(obj, m$umap)
  xr     <- range(coords[, 1]) * c(1.02, 1.02)
  yr     <- range(coords[, 2]) * c(1.02, 1.02)

  row1[[i]] <- umap_panel(
    coords, ct_broad, broad_map,
    title        = m$label,
    pt_size      = pt_size,
    alpha        = 0.5,
    guide_type   = if (i == n_cols) "right" else "none",
    legend_title = "Compartment",
    legend_ncol  = 1L,
    xlim_range   = xr,
    ylim_range   = yr
  )

  row2[[i]] <- umap_panel(
    coords, batch_vals, batch_map,
    title        = "",
    pt_size      = pt_size,
    alpha        = 0.5,
    guide_type   = if (i == n_cols) "right" else "none",
    legend_title = "Study",
    legend_ncol  = batch_ncol,
    xlim_range   = xr,
    ylim_range   = yr
  )
}

fig1 <- wrap_plots(c(row1, row2), nrow = 2, ncol = n_cols,
                   widths = c(rep(1, n_cols - 1), 1)) +
  plot_annotation(
    title    = paste0(cel, " - Breast cancer | Batch correction overview"),
    subtitle = sprintf("%s cells  |  %d studies  |  %d compartment labels",
                       format(n_cells, big.mark = ","), n_studies, n_broad),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 22, hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5, colour = "grey45")
    )
  )

log_msg(sprintf("  n_studies=%d  legend_rows=%d  fig_w=%.2fin  fig_h=%.2fin",
                n_studies, legend_rows, fig_w, fig_h))

out_pdf1 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig1_overview_large_fonts.pdf"))
out_png1 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig1_overview_large_fonts.png"))
ggsave(out_pdf1, fig1, width = fig_w, height = fig_h,
       device = pdf, useDingbats = FALSE)
ggsave(out_png1, fig1, width = fig_w, height = fig_h,
       dpi = 300, device = png)
log_msg("  Saved: ", out_png1)
log_msg("Done.")
