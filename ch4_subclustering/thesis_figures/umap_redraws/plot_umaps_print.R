#!/usr/bin/env Rscript
# plot_umaps_print.R  (comment 10 redraw, step 2 of 2)
#
# Print-ready redraw of the Ch.4 UMAP figures made by
# create_sc_object_workflow/merge/plot_umaps.R. Data are unchanged: the per-cell
# UMAP coordinates, louvain_res0.50 clusters and scimilarity_cluster_label come
# from the *_umapdata.rds tables written by extract_umap_data.R (Embeddings(obj, "umap")
# and obj@meta.data of <lineage>_clustered.RDS). Only the plotting differs:
#   - two panels (a) cluster, (b) Scimilarity label, at final 5.4 in width
#   - no on-plot labels, no titles, no axes; legends below; text >= 8 pt
#   - colours as in the old figures (Seurat DimPlot default = scales::hue_pal()
#     over the factor levels of the grouping column)
#   - points rasterised in the PDF (ggrastr)
suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(patchwork); library(ggrastr)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--dir", help = "comment10 folder (contains data/)"))))

objects <- data.frame(
  name = c("umap_Tcells_LiverBiliary", "umap_Epithelial_Kidney", "umap_Myeloid_Prostate",
           "umap_Bcells_LiverBiliary", "umap_Fibroblasts_Brain", "umap_Endothelial_Kidney",
           "umap_NKcells_Skin", "umap_HSPC_Breast"),
  expected_clusters = c(16, 43, 19, 7, 16, 26, 4, 7),
  stringsAsFactors = FALSE)

FIG_W  <- 5.4     # inches, final printed width
PANEL  <- 2.55    # inches, UMAP panel side
TXT    <- 8       # legend text pt
TITLE  <- 8.5     # legend title pt
TAG    <- 10      # panel letter pt
LINE_IN <- TXT * 1.25 / 72   # height of one legend text line, inches

log_file <- file.path(opt$dir, "comment10_umap_log.txt")
logl <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = log_file, append = TRUE) }
cat("", file = log_file)
logl("Comment 10 UMAP redraw - ", format(Sys.time(), "%Y-%m-%d %H:%M"))
logl("Script: ", file.path(opt$dir, "plot_umaps_print.R"), " (data via extract_umap_data.R)")
logl("Figure width ", FIG_W, " in; legend text ", TXT, " pt; legend titles ", TITLE,
     " pt; panel letters ", TAG, " pt bold; PNG 400 dpi; PDF with rasterised points (400 dpi)")
logl("")

cap1 <- function(x) sub("^(.)", "\\U\\1", x, perl = TRUE)
# wrap labels longer than ~45 characters onto two lines
wrap45 <- function(x) vapply(x, function(s)
  if (nchar(s) > 45) paste(strwrap(s, width = 45), collapse = "\n") else s, "")

for (i in seq_len(nrow(objects))) {
  nm <- objects$name[i]
  df <- readRDS(file.path(opt$dir, "data", paste0(nm, "_umapdata.rds")))
  info <- attr(df, "info")

  # --- checks -------------------------------------------------------------
  n_cells <- nrow(df)
  n_cl    <- length(unique(df$cluster))
  status  <- if (n_cl == objects$expected_clusters[i]) "OK" else
    sprintf("MISMATCH: expected %d, found %d in column '%s'",
            objects$expected_clusters[i], n_cl, info$cluster_col)

  # --- colours exactly as DimPlot default ---------------------------------
  cl_levels  <- info$cluster_levels                       # DimPlot level order
  cl_cols    <- setNames(scales::hue_pal()(length(cl_levels)), cl_levels)
  cl_order   <- cl_levels[order(suppressWarnings(as.numeric(cl_levels)), cl_levels)]
  lab_levels <- levels(as.factor(df$label))               # DimPlot level order
  lab_cols   <- setNames(scales::hue_pal()(length(lab_levels)), lab_levels)
  lab_show   <- setNames(wrap45(cap1(lab_levels)), lab_levels)

  df$cluster <- factor(as.character(df$cluster), levels = cl_order)
  df$label   <- factor(df$label, levels = lab_levels)

  # --- layout -------------------------------------------------------------
  pt <- if (n_cells > 100000) 0.1 else if (n_cells > 20000) 0.2 else 0.3
  cl_ncol  <- min(length(cl_levels), 8)
  cl_rows  <- ceiling(length(cl_levels) / cl_ncol)
  key_cm   <- if (length(cl_levels) > 30) 0.26 else 0.32
  lab_ncol <- 1
  lab_lines <- sum(lengths(regmatches(lab_show, gregexpr("\n", lab_show))) + 1)
  leg_h <- max(cl_rows * max(LINE_IN, key_cm / 2.54), lab_lines * LINE_IN) + 0.35
  fig_h <- min(max(PANEL + 0.25 + leg_h, 4), 6.5)

  base_theme <- theme_void(base_size = TXT) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.justification = c(0, 1),
          legend.title = element_text(size = TITLE, margin = margin(b = 2)),
          legend.text  = element_text(size = TXT, lineheight = 0.95),
          legend.title.position = "top",
          legend.key.spacing.y = unit(0.5, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(2, 3, 2, 3),
          aspect.ratio = 1)

  pa <- ggplot(df, aes(UMAP_1, UMAP_2, colour = cluster)) +
    rasterise(geom_point(size = pt, shape = 16, stroke = 0), dpi = 400) +
    scale_colour_manual(values = cl_cols, breaks = cl_order, name = "Cluster") +
    guides(colour = guide_legend(ncol = cl_ncol, byrow = TRUE,
                                 override.aes = list(size = 1.8))) +
    base_theme +
    theme(legend.key.size = unit(key_cm, "cm"),
          legend.key.spacing.x = unit(1, "pt"),
          legend.text = element_text(size = TXT, margin = margin(l = 1.5, r = 3)))

  pb <- ggplot(df, aes(UMAP_1, UMAP_2, colour = label)) +
    rasterise(geom_point(size = pt, shape = 16, stroke = 0), dpi = 400) +
    scale_colour_manual(values = lab_cols, breaks = lab_levels, labels = lab_show,
                        name = "Scimilarity label") +
    guides(colour = guide_legend(ncol = lab_ncol, override.aes = list(size = 1.8))) +
    base_theme +
    theme(legend.key.size = unit(0.32, "cm"))

  p <- (pa | pb) + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = TAG, face = "bold"),
          plot.tag.position = c(0, 1))

  png_path <- file.path(opt$dir, paste0(nm, ".png"))
  pdf_path <- file.path(opt$dir, paste0(nm, ".pdf"))
  if (file.exists(png_path) || file.exists(pdf_path)) stop("Refusing to overwrite ", nm)
  ggsave(png_path, p, width = FIG_W, height = fig_h, dpi = 400, units = "in", bg = "white")
  ggsave(pdf_path, p, width = FIG_W, height = fig_h, units = "in", device = cairo_pdf, bg = "white")

  # --- log ----------------------------------------------------------------
  logl("== ", nm, " ==")
  logl("RDS file      : ", info$rds)
  logl("Cluster column: ", info$cluster_col, "   Label column: ", info$scim_col,
       "   Reduction: umap (", paste(info$umap_colnames, collapse = ", "), ")")
  logl("n cells       : ", format(n_cells, big.mark = ","))
  logl("n clusters 0.50: ", n_cl, " (expected ", objects$expected_clusters[i], ") -> ", status)
  cl_tab <- table(df$cluster, useNA = "ifany")
  logl("Cells per cluster: ", paste(sprintf("%s=%d", names(cl_tab), cl_tab), collapse = ", "),
       "   (NA cluster: ", sum(is.na(df$cluster)), ", NA label: ", sum(is.na(df$label)), ")")
  logl("Figure size   : ", FIG_W, " x ", round(fig_h, 2), " in; point size ", pt,
       "; cluster legend ", cl_ncol, " columns")
  logl("Scimilarity labels in legend (", length(lab_levels), "):")
  for (l in lab_levels)
    logl("   - ", cap1(l), "   [", sum(df$label == l), " cells]")
  logl("")
}
logl("Done.")
