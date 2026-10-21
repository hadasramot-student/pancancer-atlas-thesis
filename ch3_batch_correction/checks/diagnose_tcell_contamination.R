#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))

path <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/T-cells_Breast_merged_step4_harmony_dev015_nocca.RDS"
cat("Loading:", path, "\n")
obj <- readRDS(path)
meta <- obj@meta.data
cat("N cells:", ncol(obj), "\n")
cat("Colnames of meta.data:\n")
print(colnames(meta))

cat("\n=== table(Consensus_Cell_Type) ===\n")
print(table(meta$Consensus_Cell_Type, useNA = "ifany"))

cat("\n=== table(Meta.scimilarity) ===\n")
print(table(meta$Meta.scimilarity, useNA = "ifany"))

cat("\n=== table(Meta.SingleR) ===\n")
print(table(meta$Meta.SingleR, useNA = "ifany"))

cat("\n=== table(Meta.Sctype) ===\n")
print(table(meta$Meta.Sctype, useNA = "ifany"))

cat("\n=== Cross-tab: Consensus_Cell_Type x Meta.scimilarity ===\n")
print(table(meta$Consensus_Cell_Type, meta$Meta.scimilarity, useNA = "ifany"))

cat("\n=== Doublet-related columns? ===\n")
print(grep("doublet|scrublet|Doublet", colnames(meta), value = TRUE, ignore.case = TRUE))

cat("\n=== Cluster-related columns? ===\n")
print(grep("cluster|seurat_clusters", colnames(meta), value = TRUE, ignore.case = TRUE))

# Focus on cells where Consensus_Cell_Type == T-cells but Meta.scimilarity != T-cells
if ("Consensus_Cell_Type" %in% colnames(meta) && "Meta.scimilarity" %in% colnames(meta)) {
  tcell_labels <- grep("^T", unique(meta$Consensus_Cell_Type), value = TRUE)
  cat("\nDistinct Consensus_Cell_Type values matching T-cell-ish:", paste(tcell_labels, collapse=", "), "\n")

  is_consensus_tcell <- meta$Consensus_Cell_Type == "T-cells"
  cat("\nN cells with Consensus_Cell_Type == 'T-cells':", sum(is_consensus_tcell, na.rm=TRUE), "\n")
  cat("N cells with Consensus_Cell_Type != 'T-cells':", sum(!is_consensus_tcell, na.rm=TRUE), "\n")

  scim_not_tcell <- is_consensus_tcell & !(meta$Meta.scimilarity %in% c("T-cells"))
  cat("\nN cells: Consensus_Cell_Type=='T-cells' BUT Meta.scimilarity != 'T-cells':", sum(scim_not_tcell, na.rm=TRUE), "\n")
  cat("Breakdown of Meta.scimilarity among those cells:\n")
  print(table(meta$Meta.scimilarity[scim_not_tcell], useNA="ifany"))

  cat("\nOf those, how many had SingleR and Sctype agree on T-cells (explaining the 2/3 majority)?\n")
  sub <- meta[scim_not_tcell, c("Meta.SingleR","Meta.Sctype","Meta.scimilarity","maj2_cell","maj2_cluster","Consensus_Cell_Type")]
  print(head(sub, 30))
  cat("\ntable of SingleR==Sctype among these cells:\n")
  print(table(sub$Meta.SingleR == sub$Meta.Sctype))

  cat("\n=== sample project_name / cluster for these contaminant-looking cells ===\n")
  extra_cols <- intersect(c("project_name","seurat_clusters","Harmony_cluster","predicted_doublet","doublet_score"), colnames(meta))
  print(extra_cols)
  if (length(extra_cols) > 0) {
    print(head(meta[scim_not_tcell, extra_cols, drop=FALSE], 30))
  }
}

cat("\nDONE\n")
