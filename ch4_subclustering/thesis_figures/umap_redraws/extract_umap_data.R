#!/usr/bin/env Rscript
# extract_umap_data.R  (comment 10 redraw, step 1 of 2)
# Same data loading as create_sc_object_workflow/merge/plot_umaps.R:
# reads <input_dir>/<lineage>_clustered.RDS, reduction "umap",
# cluster column louvain_res0.50 and scimilarity_cluster_label.
# Writes a small per-cell table (UMAP_1, UMAP_2, cluster, label) so the plotting
# step can be iterated without reloading multi-GB Seurat objects. Nothing is recomputed.
suppressPackageStartupMessages({ library(optparse); library(Seurat) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--input_dir"), make_option("--lineage"), make_option("--out_name"),
  make_option("--out_dir"),
  make_option("--cluster_col", default = "louvain_res0.50"),
  make_option("--scim_col", default = "scimilarity_cluster_label"))))

rds_path <- file.path(opt$input_dir, paste0(opt$lineage, "_clustered.RDS"))
cat("RDS:", rds_path, "\n")
obj <- readRDS(rds_path)
meta <- obj@meta.data
cat("Meta columns:", paste(colnames(meta), collapse = ", "), "\n")
stopifnot("umap" %in% names(obj@reductions))
stopifnot(opt$cluster_col %in% colnames(meta), opt$scim_col %in% colnames(meta))
emb <- Embeddings(obj, "umap")
df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                 cluster = meta[rownames(emb), opt$cluster_col],
                 label   = meta[rownames(emb), opt$scim_col],
                 stringsAsFactors = FALSE)
attr(df, "info") <- list(rds = rds_path, lineage = opt$lineage,
                         cluster_col = opt$cluster_col, scim_col = opt$scim_col,
                         cluster_levels = levels(as.factor(meta[[opt$cluster_col]])),
                         umap_colnames = colnames(emb))
saveRDS(df, file.path(opt$out_dir, paste0(opt$out_name, "_umapdata.rds")))
cat("n cells:", nrow(df), " n clusters:", length(unique(df$cluster)), "\n")
cat("pkgs: ggrastr", requireNamespace("ggrastr", quietly = TRUE),
    " ragg", requireNamespace("ragg", quietly = TRUE),
    " patchwork", requireNamespace("patchwork", quietly = TRUE),
    " cowplot", requireNamespace("cowplot", quietly = TRUE), "\n")
