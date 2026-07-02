#!/usr/bin/env Rscript
# plot_publication_figures.R
#
# Ch.3 (Batch correction) — method-comparison figures. For each Breast
# cell-type object, produces two publication-quality figures:
#
#   Figure 1 (2 x n_cols):
#     Row 1 : colored by Meta.scimilarity   (broad compartment labels)
#     Row 2 : colored by project_name       (study / batch)
#     Cols  : Uncorrected | MNN | ComBat | Harmony | scANVI
#     Legend: right side, separate per row
#
#   Figure 2 (1 x n_cols):
#     Colored by scimilarity_celltype (top 10 + Other)
#     Same 5 columns as Figure 1
#     Legend: right side
#
# Both figures saved as PDF and PNG (300 dpi).
#
# Inputs:  ATLAS_DATA_DIR/merged_by_query/<CellType>_Breast_merged_step4_harmony_dev015_nocca.RDS
#          (+ the scANVI latent CSV for the same object)
# Outputs: ATLAS_DATA_DIR/merged_by_query/publication_figures/*.{pdf,png}
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript plot_publication_figures.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

BASE    <- file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "merged_by_query")
OUT_DIR <- file.path(BASE, "publication_figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CELL_TYPES <- c(
  "B-cell", "Endothelial", "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells",
  "Myeloid", "NK-cells", "T-cells"
)

TOP_N_CELLTYPES <- 10L

# Generic / uninformative labels to force into "Other" before top-N selection.
# This ensures those colour slots go to actual subtypes instead.
DEMOTE_LABELS <- list(
  "Endothelial" = c("endothelial cell"),
  "Fibroblasts" = c("fibroblast")
)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---------------------------------------------------------------
# Colour helpers
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
# Single UMAP panel (raw ggplot2)
# guide_type : "none" | "right"
# legend_ncol: columns in legend (for large legends)
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
# Main loop
# ---------------------------------------------------------------
for (cel in CELL_TYPES) {

  log_msg("=== ", cel, " ===")

  step4_path <- file.path(BASE,
    paste0(cel, "_Breast_merged_step4_harmony_dev015_nocca.RDS"))

  if (!file.exists(step4_path)) {
    log_msg("  step4 RDS not found — skipping")
    next
  }

  log_msg("  Loading step4 object...")
  obj <- readRDS(step4_path)
  log_msg("  Cells: ", ncol(obj), "  |  Reductions: ",
          paste(names(obj@reductions), collapse = ", "))

  # ---- ensure method UMAPs ----
  obj <- ensure_umap(obj, "pca",        "umap.pca")
  obj <- ensure_umap(obj, "mnn",        "umap.mnn")
  obj <- ensure_umap(obj, "pca.combat", "umap.combat")
  obj <- ensure_umap(obj, "harmony",    "umap.harmony")

  # ---- scANVI latent ----
  lat_csv <- file.path(BASE,
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

  # ---- method columns ----
  methods <- list(
    list(umap = "umap.pca",     label = "Uncorrected"),
    list(umap = "umap.mnn",     label = "MNN"),
    list(umap = "umap.combat",  label = "ComBat"),
    list(umap = "umap.harmony", label = "Harmony")
  )
  if (has_scanvi) methods[[5]] <- list(umap = "umap.scanvi", label = "scANVI")
  n_cols <- length(methods)

  # ---- colour maps ----
  ct_broad  <- as.character(obj$Meta.scimilarity)
  broad_map <- make_color_map(ct_broad)
  n_broad   <- length(broad_map)

  batch_vals <- as.character(obj$project_name)
  batch_map  <- make_color_map(batch_vals)
  n_studies  <- length(batch_map)
  batch_ncol <- if (n_studies > 20) 2L else 1L

  ct_fine_raw  <- as.character(obj$scimilarity_celltype)
  # Demote generic labels to "Other" so top-N slots go to real subtypes
  demote <- DEMOTE_LABELS[[cel]]
  ct_fine_adj <- if (!is.null(demote)) {
    ifelse(tolower(ct_fine_raw) %in% tolower(demote), "Other", ct_fine_raw)
  } else {
    ct_fine_raw
  }
  ct_fine      <- recode_top_n(ct_fine_adj, TOP_N_CELLTYPES)
  n_shown      <- length(setdiff(unique(ct_fine), "Other"))
  ct_fine_lvls <- c(sort(setdiff(unique(ct_fine), "Other")), "Other")
  ct_fine_cols <- c(discrete_colors(n_shown), "#BBBBBB")
  names(ct_fine_cols) <- ct_fine_lvls
  fine_map <- ct_fine_cols

  # figure width: 2.8 per column + right margin for legend
  fig_w <- 2.8 * n_cols + 2.2

  # =============================================================
  # FIGURE 1 — 2 rows: broad compartment (top) | study/batch (bottom)
  # =============================================================
  log_msg("  Building Figure 1...")

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
        plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey45")
      )
    )

  out_pdf1 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig1_overview.pdf"))
  out_png1 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig1_overview.png"))
  ggsave(out_pdf1, fig1, width = fig_w, height = 7.5,
         device = pdf, useDingbats = FALSE)
  ggsave(out_png1, fig1, width = fig_w, height = 7.5,
         dpi = 300, device = png)
  log_msg("  Saved: ", out_pdf1)

  # =============================================================
  # FIGURE 2 — 1 row: fine-grained cell type subtypes
  # =============================================================
  log_msg("  Building Figure 2...")

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
        plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey45")
      )
    )

  out_pdf2 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig2_subtypes.pdf"))
  out_png2 <- file.path(OUT_DIR, paste0(cel, "_Breast_fig2_subtypes.png"))
  ggsave(out_pdf2, fig2, width = fig_w, height = 4.5,
         device = pdf, useDingbats = FALSE)
  ggsave(out_png2, fig2, width = fig_w, height = 4.5,
         dpi = 300, device = png)
  log_msg("  Saved: ", out_pdf2)

  rm(obj); gc(verbose = FALSE)
  cat("\n")
}

log_msg("Done. All figures in: ", OUT_DIR)
