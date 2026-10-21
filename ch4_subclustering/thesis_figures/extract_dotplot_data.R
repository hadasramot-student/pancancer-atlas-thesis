#!/usr/bin/env Rscript
# extract_dotplot_data.R  (comment 18)
#
# Exports the Seurat DotPlot values behind the current T cells - Liver/Biliary
# and Epithelial - Kidney marker dot plots (appendix_a_fig51 / fig52), so the
# figures can be redrawn for print. Logic copied unchanged from
#   pancancer-atlas-thesis/ch4_subclustering/fig19_tcells_liver_biliary_dotplot.R
#   pancancer-atlas-thesis/ch4_subclustering/fig23_epithelial_kidney_dotplot.R
# (louvain_res0.50, top 3 markers per cluster from the saved
# <lineage>_top_markers.csv, same panel split, DotPlot per panel with
# idents = the panel's clusters, so avg.exp.scaled is identical to the
# current figures). No markers are recomputed.
#
# Usage (R container):
#   sbatch ~/rstudio/rstudio.sbatch --script extract_dotplot_data.R \
#     --lineage T-cells_Liver_Biliary --n_panels 2 --out_dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--lineage",  type = "character"),
  make_option("--n_panels", type = "integer"),
  make_option("--top_n",    type = "integer", default = 3L),
  make_option("--in_dir",   type = "character",
              default = "/ofir_grp/users/hadasra/temp_objects/clustered_lineages"),
  make_option("--out_dir",  type = "character")
)))

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
cluster_col <- "louvain_res0.50"

rds_path     <- file.path(opt$in_dir, paste0(opt$lineage, "_clustered.RDS"))
markers_path <- file.path(opt$in_dir, paste0(opt$lineage, "_top_markers.csv"))

log_msg("Loading ", rds_path)
obj <- readRDS(rds_path)
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"

cluster_levels <- sort(as.integer(as.character(unique(obj@meta.data[[cluster_col]]))))
obj@meta.data[[cluster_col]] <- factor(as.character(obj@meta.data[[cluster_col]]),
                                        levels = as.character(cluster_levels))
Idents(obj) <- obj@meta.data[[cluster_col]]
log_msg("Cluster column: ", cluster_col, "  (", length(cluster_levels), " clusters, ",
        ncol(obj), " cells)")

markers <- read.csv(markers_path)
markers$cluster <- as.integer(markers$cluster)
top_genes <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = opt$top_n) %>%
  ungroup()

panel_size    <- ceiling(length(cluster_levels) / opt$n_panels)
cluster_panel <- setNames(ceiling(seq_along(cluster_levels) / panel_size), as.character(cluster_levels))

res <- lapply(seq_len(opt$n_panels), function(p) {
  cl_ids <- names(cluster_panel)[cluster_panel == p]
  genes  <- top_genes %>% filter(as.character(cluster) %in% cl_ids) %>% pull(gene) %>% unique()
  d <- DotPlot(obj, features = genes, idents = cl_ids,
               cols = c("lightgrey", "#B2182B"), dot.scale = 5)$data
  d$panel <- p
  d$gene_order <- match(as.character(d$features.plot), genes)
  d
})
out_df <- bind_rows(res)
out_df$features.plot <- as.character(out_df$features.plot)
out_df$id <- as.character(out_df$id)

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(out_df, file.path(opt$out_dir, paste0(opt$lineage, "_dotplot_data.csv")), row.names = FALSE)
writeLines(c(paste0("lineage\t", opt$lineage),
             paste0("n_cells\t", ncol(obj)),
             paste0("n_clusters_object\t", length(cluster_levels)),
             paste0("n_clusters_markers\t", length(unique(markers$cluster))),
             paste0("n_panels\t", opt$n_panels)),
           file.path(opt$out_dir, paste0(opt$lineage, "_dotplot_info.tsv")))
log_msg("Done.")
