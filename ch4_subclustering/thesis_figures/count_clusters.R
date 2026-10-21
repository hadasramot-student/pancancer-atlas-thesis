# Comment 18: confirm res0.50 cluster counts in the objects behind captions that disagree with the marker tables
suppressPackageStartupMessages(library(Seurat))
base <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages"
out  <- "/ofir_grp/users/hadasra/thesis_figures/comment18_2026-10-21/data/cluster_counts_check.tsv"
res <- do.call(rbind, lapply(c("Fibroblasts_Brain", "Hematopoietic_Stem_and_Progenitor_cells_Breast"), function(l) {
  md <- readRDS(file.path(base, paste0(l, "_clustered.RDS")))@meta.data
  cl <- md[["louvain_res0.50"]]
  data.frame(lineage = l, n_cells = nrow(md), n_clusters_res0.50 = length(unique(cl)),
             smallest = min(table(cl)))
}))
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
