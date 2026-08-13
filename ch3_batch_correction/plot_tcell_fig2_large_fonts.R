#!/usr/bin/env Rscript
# plot_tcell_fig2_large_fonts.R
#
# One-off re-render of the T-cells "Cell type subtypes" figure (thesis
# Fig. 16) with larger typography for print legibility at 6in display
# width. Derived from plot_publication_figures.R (Figure 2 logic only,
# T-cells only) and from plot_tcell_fig1_large_fonts.R (font-scaling
# convention used for Fig. 15). Font sizes for panel titles, legend
# text/title, and figure title/subtitle are scaled ~1.5-2x relative to
# panel size; canvas is widened/heightened modestly so the enlarged
# legend (10 subtypes + Other) doesn't clip, keeping the net relative
# increase in the 1.5-2x range. No data or plot-type changes vs. the
# original figure.
#
# Inputs:  ATLAS_DATA_DIR/merged_by_query/Breast/T-cells_Breast_merged_step4_harmony_dev015_nocca.RDS
#          (+ matching scANVI latent CSV)
# Outputs: ATLAS_DATA_DIR/merged_by_query/publication_figures/T-cells_Breast_fig2_subtypes_large_fonts.{pdf,png}
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript plot_tcell_fig2_large_fonts.R

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

recode_top_n <- function(vals, top_n = TOP_N_CELLTYPES) {
  counts <- sort(table(vals), decreasing = TRUE)
  keep   <- names(counts)[seq_len(min(top_n, length(counts)))]
  ifelse(vals %in% keep, vals, "Other")
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
# (base_size 9->15, panel title 10->17, legend.text 7->12,
#  legend.title 8->14, legend.key.size 0.3cm->0.5cm,
#  legend point size 3->5). Same scaling used for Fig. 15 (fig1).
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
      override.aes   = list(size = 5, alpha = 1),
      ncol           = legend_ncol,
      keywidth       = unit(0.5, "cm"),
      keyheight      = unit(0.5, "cm"),
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
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 17,
                                      margin = margin(b = 5)),
      legend.position  = if (guide_type == "right") "right" else "none",
      legend.text      = element_text(size = 12),
      legend.title     = element_text(size = 14, face = "bold"),
      legend.key.size  = unit(0.5, "cm"),
      legend.spacing.y = unit(0.05, "cm"),
      plot.margin      = margin(4, 4, 4, 4)
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

# ---------------------------------------------------------------
# Fine-grained cell-type subtype labels (top N + Other), matching
# the original Figure 2 logic in plot_publication_figures.R.
# T-cells has no DEMOTE_LABELS entry, so no demotion needed here.
# ---------------------------------------------------------------
ct_fine_raw  <- as.character(obj$scimilarity_celltype)
ct_fine      <- recode_top_n(ct_fine_raw, TOP_N_CELLTYPES)
n_shown      <- length(setdiff(unique(ct_fine), "Other"))
ct_fine_lvls <- c(sort(setdiff(unique(ct_fine), "Other")), "Other")
ct_fine_cols <- c(discrete_colors(n_shown), "#BBBBBB")
names(ct_fine_cols) <- ct_fine_lvls
fine_map <- ct_fine_cols

# ---------------------------------------------------------------
# Canvas sizing: same per-column deltas as the Fig.15 (fig1)
# large-fonts script (2.8->3.2in/col, 2.2->3.2in margin), since the
# panel aesthetics (umap_panel) are identical. Height is driven by
# the single subtype legend (up to TOP_N_CELLTYPES + "Other" rows,
# ncol=1), using the same empirical per-row constant (0.20in/row at
# legend.text=12pt / key=0.5cm) validated for Fig.15. A panel_min_in
# floor keeps the actual UMAP scatter row from being squeezed when
# the legend is short.
# ---------------------------------------------------------------
legend_rows     <- n_shown + 1L   # + "Other"
per_row_in      <- 0.20
legend_block_in <- legend_rows * per_row_in + 0.35
row_avail_in    <- legend_block_in * 1.25
title_block_in  <- 1.6
panel_min_in    <- 4.0

fig_w <- 3.2 * n_cols + 3.2
fig_h <- title_block_in + max(row_avail_in, panel_min_in) + 0.5

log_msg("  Building Figure 2 (large fonts)...")

row_f2 <- vector("list", n_cols)

for (i in seq_len(n_cols)) {
  m      <- methods[[i]]
  coords <- Embeddings(obj, m$umap)
  xr     <- range(coords[, 1]) * c(1.02, 1.02)
  yr     <- range(coords[, 2]) * c(1.02, 1.02)

  row_f2[[i]] <- umap_panel(
    coords, ct_fine, fine_map,
    title        = m$label,
    pt_size      = pt_size,
    alpha        = 0.5,
    guide_type   = if (i == n_cols) "right" else "none",
    legend_title = sprintf("Cell type\n(top %d + Other)", TOP_N_CELLTYPES),
    legend_ncol  = 1L,
    xlim_range   = xr,
    ylim_range   = yr
  )
}

fig2 <- wrap_plots(row_f2, nrow = 1, ncol = n_cols) +
  plot_annotation(
    title    = paste0(cel, " - Breast cancer | Cell type subtypes"),
    subtitle = sprintf("%s cells  |  %d subtype labels",
                       format(n_cells, big.mark = ","),
                       length(unique(ct_fine_raw))),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 22, hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5, colour = "grey45")
    )
  )

log_msg(sprintf("  n_shown=%d  legend_rows=%d  fig_w=%.2fin  fig_h=%.2fin",
                n_shown, legend_rows, fig_w, fig_h))

out_pdf2 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig2_subtypes_large_fonts.pdf"))
out_png2 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig2_subtypes_large_fonts.png"))
ggsave(out_pdf2, fig2, width = fig_w, height = fig_h,
       device = pdf, useDingbats = FALSE)
ggsave(out_png2, fig2, width = fig_w, height = fig_h,
       dpi = 300, device = png)
log_msg("  Saved: ", out_png2)
log_msg("Done.")
