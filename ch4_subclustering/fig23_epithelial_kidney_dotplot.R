#!/usr/bin/env Rscript
# fig23_epithelial_kidney_dotplot.R
#
# Ch.4 Fig 23 — Epithelial cells, Kidney: dotplot of marker genes per
# cluster, redesigned for legibility.
#
# The original design put all 44 clusters (louvain_res0.50) x top-5 markers
# each on one x-axis (~220 gene-label columns) -- unreadable at any print
# size. Fix here combines two density reductions:
#   (a) top 3 markers/cluster instead of top 5 (from Epithelial_Kidney_top_markers.csv,
#       already ranked by avg_log2FC via find_markers.R)
#   (b) split into 3 panels of ~15 clusters each, stacked vertically with a
#       shared legend, so each panel's x-axis holds ~45 columns instead of ~220
#
# Run via R container:
#   sbatch --mem=128G ~/rstudio/rstudio.sbatch \
#     --script ch4_subclustering/fig23_epithelial_kidney_dotplot.R

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

rds_path     <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages/Epithelial_Kidney_clustered.RDS"
markers_path <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages/Epithelial_Kidney_top_markers.csv"
out_png      <- "/ofir_grp/users/hadasra/create_sc_object_workflow/plots_and_outputs/Fig23_Epithelial_Kidney_dotplot_top3_3panels.png"
cluster_col  <- "louvain_res0.50"
top_n        <- 3
n_panels     <- 3

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

log_msg("Loading object ...")
obj <- readRDS(rds_path)
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"

# Cluster ids are stored as character/unordered factor, so default sorting
# is alphabetical ("1","10","11",...,"2",...) -- re-level numerically or
# both the axis order and the panel split below come out wrong.
cluster_levels <- sort(as.integer(as.character(unique(obj@meta.data[[cluster_col]]))))
obj@meta.data[[cluster_col]] <- factor(as.character(obj@meta.data[[cluster_col]]),
                                        levels = as.character(cluster_levels))
Idents(obj) <- obj@meta.data[[cluster_col]]
log_msg("Cluster column: ", cluster_col, "  (", length(cluster_levels), " clusters)")

markers <- read.csv(markers_path)
markers$cluster <- as.integer(markers$cluster)

top_genes <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

panel_size    <- ceiling(length(cluster_levels) / n_panels)
cluster_panel <- setNames(ceiling(seq_along(cluster_levels) / panel_size), as.character(cluster_levels))

make_panel <- function(p) {
  cl_ids <- names(cluster_panel)[cluster_panel == p]
  genes  <- top_genes %>% filter(as.character(cluster) %in% cl_ids) %>% pull(gene) %>% unique()
  DotPlot(obj, features = genes, idents = cl_ids, cols = c("lightgrey", "#B2182B"), dot.scale = 5) +
    RotatedAxis() +
    labs(title = paste0("Clusters ", cl_ids[1], "–", cl_ids[length(cl_ids)]), x = NULL, y = "Cluster") +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 8),
          plot.title  = element_text(size = 10, face = "bold"))
}

log_msg("Building ", n_panels, " panels (top ", top_n, " genes/cluster) ...")
panels   <- lapply(seq_len(n_panels), make_panel)
combined <- wrap_plots(panels, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Epithelial cells — Kidney: top marker genes per cluster (louvain res 0.50)")

log_msg("Saving -> ", out_png)
ggsave(out_png, combined, width = 16, height = 4.4 * n_panels, dpi = 400, limitsize = FALSE)
log_msg("Done.")
