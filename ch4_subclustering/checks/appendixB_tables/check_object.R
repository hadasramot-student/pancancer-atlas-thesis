# check_object.R — read-only per-object facts for the Appendix B cluster tables.
#
# For one <object>_clustered.RDS this script:
#   1. recomputes, from the object's metadata, n cells / Scimilarity majority label /
#      confidence % per cluster at louvain_res0.50 and at the table resolution
#      (same majority vote as ch4_subclustering/annotate_lineages.R, on scimilarity_celltype);
#   2. checks that the saved <object>_markers.csv belongs to this object's res-0.50 clusters:
#      for the top 3 genes per cluster it recomputes pct.1 / pct.2 (fraction of cells with
#      data > 0 inside / outside the cluster) and compares them with the saved values.
# Nothing is re-clustered, re-annotated or re-tested; nothing is written next to the object.
#
# Usage (inside the R container):
#   Rscript check_object.R <object> <table_res> <clustered_dir> <out_dir>
#   e.g. Rscript check_object.R Epithelial_Kidney 0.30 /ofir_grp/users/hadasra/temp_objects/clustered_lineages out

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
obj_name  <- args[1]
table_res <- args[2]
cl_dir    <- args[3]
out_dir   <- args[4]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

rds <- file.path(cl_dir, paste0(obj_name, "_clustered.RDS"))
log_msg("Reading ", rds)
obj <- readRDS(rds)
md  <- obj@meta.data
log_msg(ncol(obj), " cells")

label_col <- "scimilarity_celltype"
if (!label_col %in% colnames(md)) stop("no ", label_col, " column in ", obj_name)

majority_vote <- function(labels, clusters, res) {
  df <- data.frame(label = as.character(labels), cluster = as.character(clusters))
  do.call(rbind, lapply(sort(unique(df$cluster)), function(cl) {
    sub <- df[df$cluster == cl, ]
    tbl <- sort(table(sub$label), decreasing = TRUE)
    data.frame(resolution = res, cluster = cl, majority_label = names(tbl)[1],
               confidence_pct = round(100 * tbl[[1]] / nrow(sub), 1),
               n_cells = nrow(sub), stringsAsFactors = FALSE)
  }))
}

res_set <- unique(c("0.50", table_res))
votes <- do.call(rbind, lapply(res_set, function(r) {
  col <- paste0("louvain_res", r)
  if (!col %in% colnames(md)) stop("missing ", col, " in ", obj_name)
  majority_vote(md[[label_col]], md[[col]], col)
}))
votes$object <- obj_name
write.table(votes, file.path(out_dir, paste0(obj_name, "_clusters_from_object.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote cluster votes: ", nrow(votes), " rows")

# Marker-file consistency check at res 0.50 (markers were computed only at 0.50)
mk <- read.csv(file.path(cl_dir, paste0(obj_name, "_markers.csv")), check.names = FALSE)
mk <- mk[order(mk$cluster, -mk$avg_log2FC), ]
top3 <- do.call(rbind, lapply(split(mk, mk$cluster), head, 3))
genes <- intersect(unique(top3$gene), rownames(obj))
dat <- LayerData(obj, assay = DefaultAssay(obj), layer = "data")[genes, , drop = FALSE]
pos <- dat > 0
cl50 <- as.character(md[["louvain_res0.50"]])
top3$pct1_obj <- NA_real_
top3$pct2_obj <- NA_real_
for (i in seq_len(nrow(top3))) {
  g <- top3$gene[i]
  if (!g %in% genes) next
  inside <- cl50 == as.character(top3$cluster[i])
  if (!any(inside)) next
  v <- pos[g, ]
  top3$pct1_obj[i] <- round(mean(v[inside]), 3)
  top3$pct2_obj[i] <- round(mean(v[!inside]), 3)
}
top3$object <- obj_name
write.table(top3[, c("object", "cluster", "gene", "avg_log2FC", "pct.1", "pct1_obj", "pct.2", "pct2_obj")],
            file.path(out_dir, paste0(obj_name, "_marker_check.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Marker check: max |pct.1 diff| = ",
        signif(max(abs(top3$pct.1 - top3$pct1_obj), na.rm = TRUE), 3))
log_msg("DONE")
