#!/usr/bin/env Rscript
# plot_batch_correction_comparison.R
#
# Ch.3 (Batch correction) — method-comparison figures. For each Breast
# cell-type object, produces a 2x4 UMAP grid:
#   Columns  : MNN | ComBat | Harmony | scANVI
#   Top row  : colored by project_name  (batch mixing)
#   Bottom row: colored by scimilarity_celltype (biology preservation)
#
# All 4 methods use the dev015_nocca objects so comparisons are on the same cells.
# Data sources:
#   MNN / ComBat / Harmony  : step4_harmony_dev015_nocca.RDS (contains all 3 reductions+UMAPs)
#   scANVI                  : step6_scanvi_dev015_nocca_scanvi_latent.csv -> UMAP computed here
#
# Inputs:  ATLAS_DATA_DIR/merged_by_query/<CellType>_Breast_merged_step4_harmony_dev015_nocca.RDS
# Outputs: ATLAS_DATA_DIR/merged_by_query/batch_correction_comparison_plots/*.{pdf,png}
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript plot_batch_correction_comparison.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

BASE    <- file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "merged_by_query")
OUT_DIR <- file.path(BASE, "batch_correction_comparison_plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CELL_TYPES <- c(
  "B-cell",
  "Endothelial",
  "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells",
  "Myeloid",
  "NK-cells",
  "T-cells"
)

# Maximum labelled cell types in the bottom-row legend.
# Labels beyond this are collapsed to "Other".
TOP_N_CELLTYPES <- 10L

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---------------------------------------------------------------
# Discrete colour palette — scales to any number of levels.
# "Other" is always grey.
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

# Collapse rare cell-type labels into "Other".
# Returns a named character vector mapping original → display label.
recode_celltypes <- function(vals, top_n = TOP_N_CELLTYPES) {
  counts    <- sort(table(vals), decreasing = TRUE)
  keep      <- names(counts)[seq_len(min(top_n, length(counts)))]
  recoded   <- ifelse(vals %in% keep, vals, "Other")
  recoded
}

# ---------------------------------------------------------------
# Single UMAP panel — shuffled draw order so no batch sits on top.
# guide_type: "none" | "bottom"  (bottom = legend below the panel)
# ---------------------------------------------------------------
umap_panel <- function(coords, vals, color_map,
                       title      = "",
                       subtitle   = "",
                       pt_size    = 0.25,
                       guide_type = "none",
                       legend_title = "") {

  idx <- sample(length(vals))
  df  <- data.frame(
    x   = coords[idx, 1],
    y   = coords[idx, 2],
    col = factor(vals[idx], levels = names(color_map)),
    stringsAsFactors = FALSE
  )

  leg <- switch(guide_type,
    bottom = guide_legend(
      override.aes = list(size = 3, alpha = 1),
      ncol         = 4,
      keywidth     = unit(0.4, "cm"),
      keyheight    = unit(0.4, "cm")
    ),
    "none"
  )

  ggplot(df, aes(x, y, colour = col)) +
    geom_point(size = pt_size, alpha = 0.65, stroke = 0, shape = 16) +
    scale_colour_manual(values   = color_map,
                        name     = legend_title,
                        na.value = "grey80",
                        guide    = leg) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 9) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 10,
                                   margin = margin(b = 1)),
      plot.subtitle = element_text(hjust = 0.5, size = 7.5, colour = "grey50",
                                   margin = margin(b = 2)),
      legend.position  = if (guide_type == "bottom") "bottom" else "none",
      legend.text      = element_text(size = 7.5),
      legend.title     = element_text(size = 8, face = "bold"),
      legend.spacing.x = unit(0.15, "cm"),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# ---------------------------------------------------------------
# Ensure a UMAP is present; compute it if missing
# ---------------------------------------------------------------
ensure_umap <- function(obj, from_reduction, umap_name,
                        n_neighbors = 30, min_dist = 0.3) {
  if (umap_name %in% names(obj@reductions)) return(obj)
  log_msg("  Computing UMAP '", umap_name, "' from '", from_reduction, "'...")
  ndim <- ncol(Embeddings(obj, from_reduction))
  obj  <- RunUMAP(obj,
                  reduction      = from_reduction,
                  dims           = seq_len(min(50L, ndim)),
                  reduction.name = umap_name,
                  reduction.key  = paste0(toupper(gsub("\\.", "", umap_name)), "_"),
                  n.neighbors    = n_neighbors,
                  min.dist       = min_dist,
                  verbose        = FALSE)
  obj
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
  log_msg("  Cells: ", ncol(obj))

  # --- ensure all three method UMAPs exist ---
  obj <- ensure_umap(obj, "mnn",        "umap.mnn")
  obj <- ensure_umap(obj, "pca.combat", "umap.combat")
  obj <- ensure_umap(obj, "harmony",    "umap.harmony")

  # --- add scANVI latent as a reduction ---
  lat_csv <- file.path(BASE,
    paste0(cel, "_Breast_merged_step6_scanvi_dev015_nocca_scanvi_latent.csv"))

  if (!file.exists(lat_csv)) {
    log_msg("  scANVI latent CSV not found — skipping scANVI column")
    has_scanvi <- FALSE
  } else {
    log_msg("  Loading scANVI latent CSV...")
    lat_df <- read.csv(lat_csv, row.names = 1, stringsAsFactors = FALSE)
    common  <- intersect(colnames(obj), rownames(lat_df))
    log_msg("  Cells in common: ", length(common), " / ", ncol(obj))

    if (length(common) < 0.85 * ncol(obj)) {
      log_msg("  Too few common cells — skipping scANVI")
      has_scanvi <- FALSE
    } else {
      obj <- obj[, common]
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

  # --- point size scaled to cell count ---
  pt <- ifelse(ncol(obj) > 20000, 0.08, ifelse(ncol(obj) > 8000, 0.15, 0.35))

  # --- batch colour map (many studies; no legend shown) ---
  batch_map <- make_color_map(obj$project_name)
  n_studies <- length(unique(obj$project_name))

  # --- cell-type: keep top TOP_N_CELLTYPES, collapse rest to "Other" ---
  ct_recoded   <- recode_celltypes(obj$scimilarity_celltype)
  ct_levels    <- c(sort(setdiff(unique(ct_recoded), "Other")), "Other")
  n_shown      <- sum(ct_levels != "Other")
  ct_cols      <- c(discrete_colors(n_shown), "Other" = "#BBBBBB")
  names(ct_cols)[seq_len(n_shown)] <- ct_levels[seq_len(n_shown)]
  celltype_map <- ct_cols

  n_orig_labels <- length(unique(obj$scimilarity_celltype))

  # --- define columns ---
  methods <- list(
    list(umap = "umap.mnn",     label = "MNN"),
    list(umap = "umap.combat",  label = "ComBat"),
    list(umap = "umap.harmony", label = "Harmony")
  )
  if (has_scanvi) methods[[4]] <- list(umap = "umap.scanvi", label = "scANVI")
  n_cols <- length(methods)

  # --- build panels ---
  # Top row  : colored by batch, NO legend (just shows mixing visually)
  # Bottom row: colored by cell type, legend collected at bottom of figure
  top_row    <- vector("list", n_cols)
  bottom_row <- vector("list", n_cols)

  for (i in seq_len(n_cols)) {
    m            <- methods[[i]]
    coords       <- Embeddings(obj, m$umap)
    batch_vals   <- as.character(obj$project_name)
    ct_vals      <- ct_recoded

    # Top row: method name as title, study count as subtitle, NO legend
    top_row[[i]] <- umap_panel(
      coords, batch_vals, batch_map,
      title      = m$label,
      subtitle   = if (i == 1) sprintf("n = %d studies", n_studies) else "",
      pt_size    = pt,
      guide_type = "none"
    )

    # Bottom row: "Cell type" label only on first column, legend on ALL panels
    # (patchwork guide collection will deduplicate into one shared legend)
    bottom_row[[i]] <- umap_panel(
      coords, ct_vals, celltype_map,
      title        = if (i == 1) "Cell type" else "",
      pt_size      = pt,
      guide_type   = "bottom",
      legend_title = sprintf("Cell type  (top %d shown; %d total)",
                             n_shown, n_orig_labels)
    )
  }

  # --- assemble: top row then bottom row, collect legend at bottom ---
  fig <- wrap_plots(c(top_row, bottom_row), nrow = 2, ncol = n_cols) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  fig <- fig +
    plot_annotation(
      title    = paste0(cel, " - Breast cancer | Batch correction comparison"),
      subtitle = sprintf("%s cells  |  %d studies  |  %d cell-type labels",
                         format(ncol(obj), big.mark = ","),
                         n_studies, n_orig_labels),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey45")
      )
    )

  # --- save (extra height for bottom legend) ---
  out_pdf <- file.path(OUT_DIR, paste0(cel, "_Breast_batch_comparison_v2.pdf"))
  ggsave(out_pdf, fig,
         width  = 3.4 * n_cols,
         height = 9.0,
         device = pdf, useDingbats = FALSE)
  log_msg("  Saved: ", out_pdf)

  rm(obj); gc(verbose = FALSE)
  cat("\n")
}

log_msg("Done. All figures in: ", OUT_DIR)
