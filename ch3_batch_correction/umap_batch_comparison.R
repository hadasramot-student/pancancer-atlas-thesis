#!/usr/bin/env Rscript
# umap_batch_comparison.R
#
# Ch.3 (Batch correction) — 2x5 UMAP panel comparing batch correction
# methods for one lineage x cancer-type object (example: Myeloid, Breast).
# Row 1: coloured by study of origin; Row 2: coloured by cell-type subtype.
# Conditions: Uncorrected (PCA), ComBat, MNN, Harmony, scANVI.
#
# Inputs: --rds_path : a *_merged_full_scanvi_final.RDS with pca/combat/
#         mnn/harmony/scanvi reductions already computed (batch1-4 + scANVI)
# Outputs: --out_png : the comparison figure
#
# Usage:
#   Rscript umap_batch_comparison.R \
#     --rds_path data/temp_objects/merged_by_query/Breast/Myeloid_Breast_merged_full_scanvi_final.RDS \
#     --out_png outputs/figures/Myeloid_Breast_batch_correction_UMAP.png

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

if (nzchar(Sys.getenv("R_USER_LIBRARY"))) .libPaths(c(Sys.getenv("R_USER_LIBRARY"), .libPaths()))

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

option_list <- list(
  make_option("--rds_path", type = "character",
              default = "data/temp_objects/merged_by_query/Breast/Myeloid_Breast_merged_full_scanvi_final.RDS"),
  make_option("--out_png", type = "character",
              default = "outputs/figures/Myeloid_Breast_batch_correction_UMAP.png")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Paths ─────────────────────────────────────────────────────────────────────
RDS_PATH <- opt$rds_path
OUT_PNG  <- opt$out_png

# ── Load ──────────────────────────────────────────────────────────────────────
log_msg(paste("Loading:", RDS_PATH))
obj <- readRDS(RDS_PATH)
log_msg(sprintf("Loaded: %d cells, %d genes", ncol(obj), nrow(obj)))
log_msg(paste("Reductions present:", paste(names(obj@reductions), collapse = ", ")))
log_msg(paste("Meta cols:", paste(colnames(obj@meta.data), collapse = ", ")))

# ── Determine correct column names ────────────────────────────────────────────
# Study / batch column
study_col <- if ("study" %in% colnames(obj@meta.data)) {
  "study"
} else if ("project_name" %in% colnames(obj@meta.data)) {
  "project_name"
} else {
  stop("Cannot find study/batch column. Available: ", paste(colnames(obj@meta.data), collapse = ", "))
}

# Cell-type label column
celltype_col <- if ("majority_scimilarity" %in% colnames(obj@meta.data)) {
  "majority_scimilarity"
} else {
  stop("Cannot find majority_scimilarity column. Available: ", paste(colnames(obj@meta.data), collapse = ", "))
}

log_msg(paste("Using study col:", study_col, "| cell-type col:", celltype_col))

# ── UMAP key mapping ──────────────────────────────────────────────────────────
# Each entry: list(reduction = <embedding>, umap_name = <umap reduction to use/create>)
method_info <- list(
  "Uncorrected" = list(red = "pca",     umap = "umap.pca"),
  "ComBat"      = list(red = "pca.combat", umap = "umap.combat"),
  "MNN"         = list(red = "mnn",     umap = "umap.mnn"),
  "Harmony"     = list(red = "harmony", umap = "umap_harmony"),
  "scANVI"      = list(red = "scanvi",  umap = "umap.scanvi")
)

all_reds  <- names(obj@reductions)

# ── Compute missing UMAPs ─────────────────────────────────────────────────────
for (label in names(method_info)) {
  info <- method_info[[label]]
  red  <- info$red
  umap <- info$umap

  # Check alternative naming conventions for the embedding
  if (!red %in% all_reds) {
    alternatives <- list(
      "pca.combat" = c("pca.combat", "combat"),
      "mnn"        = c("mnn", "fastmnn"),
      "harmony"    = c("harmony", "integrated.harmony"),
      "scanvi"     = c("scanvi", "SCANVI"),
      "pca"        = c("pca", "PCA")
    )
    for (alt in alternatives[[red]]) {
      if (alt %in% all_reds) { red <- alt; break }
    }
    method_info[[label]]$red <- red
  }

  if (!red %in% all_reds) {
    log_msg(paste("WARNING: reduction", red, "not found for", label, "— skipping"))
    method_info[[label]] <- NULL
    next
  }

  # Check for already-computed UMAP (try multiple naming variants)
  umap_variants <- c(umap,
                     gsub("\\.", "_", umap),
                     gsub("_", ".", umap),
                     paste0("umap.", red),
                     paste0("umap_", red))
  found_umap <- intersect(umap_variants, names(obj@reductions))

  if (length(found_umap) > 0) {
    method_info[[label]]$umap <- found_umap[1]
    log_msg(paste(label, "— using existing UMAP:", found_umap[1]))
  } else {
    log_msg(paste(label, "— computing UMAP from reduction:", red))
    n_dims <- min(30L, ncol(Embeddings(obj, red)))
    obj <- RunUMAP(obj,
                   reduction      = red,
                   dims           = seq_len(n_dims),
                   reduction.name = umap,
                   reduction.key  = paste0(gsub("[^A-Za-z0-9]", "", umap), "_"),
                   verbose        = FALSE)
    method_info[[label]]$umap <- umap
    log_msg(paste(label, "— UMAP done:", umap))
  }
}

# Remove any methods that were dropped
method_info <- Filter(Negate(is.null), method_info)
log_msg(paste("Methods to plot:", paste(names(method_info), collapse = ", ")))

# ── Colour palettes ───────────────────────────────────────────────────────────
n_studies   <- length(unique(obj@meta.data[[study_col]]))
n_celltypes <- length(unique(obj@meta.data[[celltype_col]]))

study_pal    <- setNames(scales::hue_pal()(n_studies),
                         sort(unique(obj@meta.data[[study_col]])))
celltype_pal <- setNames(scales::hue_pal(l = 55, c = 80)(n_celltypes),
                         sort(unique(obj@meta.data[[celltype_col]])))

pt_size <- ifelse(ncol(obj) > 20000, 0.08, ifelse(ncol(obj) > 5000, 0.2, 0.4))

# ── Plot factory ──────────────────────────────────────────────────────────────
make_umap <- function(obj, umap_red, color_by, palette, title = NULL, show_legend = FALSE) {
  emb <- as.data.frame(Embeddings(obj, umap_red)[, 1:2])
  colnames(emb) <- c("U1", "U2")
  emb$fill <- obj@meta.data[[color_by]]
  emb <- emb[sample(nrow(emb)), ]

  p <- ggplot(emb, aes(U1, U2, colour = fill)) +
    geom_point(size = pt_size, alpha = 0.7, stroke = 0) +
    scale_colour_manual(values = palette, name = color_by, na.value = "grey80") +
    theme_classic(base_size = 9) +
    theme(
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      axis.line    = element_line(linewidth = 0.3),
      plot.title   = element_text(size = 9, hjust = 0.5, face = "bold", margin = margin(b = 2)),
      legend.position  = if (show_legend) "right" else "none",
      legend.key.size  = unit(0.3, "cm"),
      legend.text      = element_text(size = 6),
      legend.title     = element_blank(),
      plot.margin      = margin(2, 2, 2, 2)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1))

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

# ── Build rows ────────────────────────────────────────────────────────────────
methods <- names(method_info)
n_methods <- length(methods)

row_label_plot <- function(txt) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = txt,
             angle = 90, fontface = "bold", size = 3.5, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

row1_plots <- lapply(seq_along(methods), function(i) {
  m <- methods[i]
  make_umap(obj,
            umap_red    = method_info[[m]]$umap,
            color_by    = study_col,
            palette     = study_pal,
            title       = m,
            show_legend = (i == n_methods))
})

row2_plots <- lapply(seq_along(methods), function(i) {
  m <- methods[i]
  make_umap(obj,
            umap_red    = method_info[[m]]$umap,
            color_by    = celltype_col,
            palette     = celltype_pal,
            title       = NULL,
            show_legend = (i == n_methods))
})

# ── Assemble panel ────────────────────────────────────────────────────────────
row1 <- (row_label_plot("Study of origin") | wrap_plots(row1_plots, nrow = 1)) +
  plot_layout(widths = c(0.025, 1))

row2 <- (row_label_plot("Cell-type subtype") | wrap_plots(row2_plots, nrow = 1)) +
  plot_layout(widths = c(0.025, 1))

final <- (row1 / row2) +
  plot_annotation(
    title = "Myeloid cells — Breast cancer cohort: batch correction comparison",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                            margin = margin(b = 6)))
  )

# ── Save ──────────────────────────────────────────────────────────────────────
log_msg(paste("Saving:", OUT_PNG))
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PNG, plot = final, width = 20, height = 8, dpi = 300, bg = "white")
log_msg("Done.")
