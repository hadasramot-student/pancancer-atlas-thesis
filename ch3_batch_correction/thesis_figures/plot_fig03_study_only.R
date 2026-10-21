#!/usr/bin/env Rscript
# plot_fig03_study_only.R
#
# Standalone single-row version of the "study of origin" panel from
# plot_publication_figures.R Figure 1 (bottom row), for the
# T-cells - Breast cancer object only.
#
# Figure (1 x 5):
#   Cols   : Uncorrected | MNN | ComBat | Harmony | scANVI
#   Colour : project_name (study / batch)
#   Legend : right side
#
# Reuses the same discrete_colors()/make_color_map()/umap_panel()/ensure_umap()
# logic as plot_publication_figures.R so the visual style matches exactly;
# just drops the Meta.scimilarity compartment row and keeps only the
# batch/study row, rendered as its own complete figure with column headers.
#
# Run via R container:
#   sbatch --mem=128G ~/rstudio/rstudio.sbatch \
#     --script /ofir_grp/users/hadasra/create_sc_object_workflow/merge/plot_fig03_study_only.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

BASE       <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query"
INPUT_BASE <- file.path(BASE, "Breast")   # per-cell-type RDS/CSV live under Breast/ as of the 2026-07-06 reorg
OUT_DIR    <- file.path(BASE, "publication_figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CEL <- "T-cells"

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---------------------------------------------------------------
# Colour helpers (identical to plot_publication_figures.R)
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
# Single UMAP panel (raw ggplot2) - identical to plot_publication_figures.R
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
      override.aes = list(size = 3, alpha = 1),
      ncol         = legend_ncol,
      keywidth     = unit(0.35, "cm"),
      keyheight    = unit(0.35, "cm")
    )
  } else "none"

  p <- ggplot(df, aes(x, y, colour = col)) +
    geom_point(size = pt_size, alpha = alpha, stroke = 0, shape = 16) +
    scale_colour_manual(values   = color_map,
                        name     = legend_title,
                        na.value = "grey80",
                        guide    = leg) +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 10,
                                     margin = margin(b = 3)),
      legend.position = if (guide_type == "right") "right" else "none",
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.3, "cm"),
      plot.margin     = margin(4, 4, 4, 4)
    )

  if (!is.null(xlim_range)) p <- p + coord_cartesian(xlim = xlim_range, ylim = ylim_range)
  p
}

# ---------------------------------------------------------------
# Load object
# ---------------------------------------------------------------
log_msg("=== ", CEL, " ===")

step4_path <- file.path(INPUT_BASE, paste0(CEL, "_Breast_merged_step4_harmony_dev015_nocca.RDS"))
if (!file.exists(step4_path)) stop("step4 RDS not found: ", step4_path)

log_msg("Loading step4 object...")
obj <- readRDS(step4_path)
log_msg("Cells: ", ncol(obj), "  |  Reductions: ", paste(names(obj@reductions), collapse = ", "))

# ---- ensure method UMAPs (same as plot_publication_figures.R) ----
obj <- ensure_umap(obj, "pca",        "umap.pca")
obj <- ensure_umap(obj, "mnn",        "umap.mnn")
obj <- ensure_umap(obj, "pca.combat", "umap.combat")
obj <- ensure_umap(obj, "harmony",    "umap.harmony")

# ---- scANVI latent ----
lat_csv <- file.path(INPUT_BASE, paste0(CEL, "_Breast_merged_step6_scanvi_dev015_nocca_scanvi_latent.csv"))

has_scanvi <- FALSE
if (!file.exists(lat_csv)) {
  log_msg("scANVI latent CSV not found — 4-column figure")
} else {
  log_msg("Loading scANVI latent...")
  lat_df <- read.csv(lat_csv, row.names = 1, stringsAsFactors = FALSE)
  common <- intersect(colnames(obj), rownames(lat_df))
  log_msg("Cells in common: ", length(common), " / ", ncol(obj))

  if (length(common) < 0.85 * ncol(obj)) {
    log_msg("Too few common cells — skipping scANVI")
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

# ---- method columns ----
methods <- list(
  list(umap = "umap.pca",     label = "Uncorrected"),
  list(umap = "umap.mnn",     label = "MNN"),
  list(umap = "umap.combat",  label = "ComBat"),
  list(umap = "umap.harmony", label = "Harmony")
)
if (has_scanvi) methods[[5]] <- list(umap = "umap.scanvi", label = "scANVI")
n_cols <- length(methods)

# ---- colour map: study of origin ----
batch_vals <- as.character(obj$project_name)
batch_map  <- make_color_map(batch_vals)
n_studies  <- length(batch_map)
batch_ncol <- if (n_studies > 20) 2L else 1L

fig_w <- 2.8 * n_cols + 2.2

# =============================================================
# FIGURE 03 — 1 row: study of origin only, with column headers
# =============================================================
log_msg("Building single-row study-of-origin figure...")

row_study <- vector("list", n_cols)

for (i in seq_len(n_cols)) {
  m      <- methods[[i]]
  coords <- Embeddings(obj, m$umap)
  xr     <- range(coords[, 1]) * c(1.02, 1.02)
  yr     <- range(coords[, 2]) * c(1.02, 1.02)

  row_study[[i]] <- umap_panel(
    coords, batch_vals, batch_map,
    title        = m$label,                          # column header, now on the only row
    pt_size      = pt_size,
    alpha        = 0.5,
    guide_type   = if (i == n_cols) "right" else "none",
    legend_title = "Study",
    legend_ncol  = batch_ncol,
    xlim_range   = xr,
    ylim_range   = yr
  )
}

fig03 <- wrap_plots(row_study, nrow = 1, ncol = n_cols) +
  plot_annotation(
    title    = paste0(CEL, " - Breast cancer | Study of origin across batch correction methods"),
    subtitle = sprintf("%s cells  |  %d studies", format(n_cells, big.mark = ","), n_studies),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey45")
    )
  )

out_pdf <- file.path(OUT_DIR, paste0(CEL, "_Breast_fig03_study_only.pdf"))
out_png <- file.path(OUT_DIR, paste0(CEL, "_Breast_fig03_study_only.png"))

ggsave(out_pdf, fig03, width = fig_w, height = 4.5,
       device = pdf, useDingbats = FALSE)
ggsave(out_png, fig03, width = fig_w, height = 4.5,
       dpi = 300, device = png)

log_msg("Saved: ", out_pdf)
log_msg("Saved: ", out_png)
log_msg("Done.")
