#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Probe the 11 epithelial per-cancer-type scANVI objects to confirm
# readiness for the Option-A batch-correction benchmark:
#   - reductions present (need harmony + scanvi already; add mnn + pca.combat)
#   - scimilarity_celltype grading column present & populated
#   - RNA assay layers (counts / data / scale.data) + VariableFeatures
#   - project_name batch column: n batches, smallest batch
# Read-only. Writes a human log to stdout and a tidy CSV summary.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

base <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query"
out_csv <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/epithelial-batch-correction-benchmark-f0c664/epi_benchmark/logs/00_probe_summary.csv"

objs <- c(
  Brain         = file.path(base, "Brain/Epithelial_Brain_merged_full_scanvi_final.RDS"),
  Breast        = file.path(base, "Breast/Epithelial_Breast_subsample450k_merged_full_scanvi_final.RDS"),
  Breast_full   = file.path(base, "Breast/Epithelial_Breast_merged_full_scanvi_final.RDS"),
  Colorectal    = file.path(base, "Colorectal/Epithelial_Colorectal_merged_full_scanvi_final.RDS"),
  Kidney        = file.path(base, "Kidney/Epithelial_Kidney_merged_full_scanvi_final.RDS"),
  Liver_Biliary = file.path(base, "Liver_Biliary/Epithelial_Liver_Biliary_merged_full_scanvi_final.RDS"),
  Lung          = file.path(base, "Lung/Epithelial_Lung_merged_full_scanvi_final.RDS"),
  Ovarian       = file.path(base, "Ovarian/Epithelial_Ovarian_merged_full_scanvi_final.RDS"),
  Pancreas      = file.path(base, "Pancreas/Epithelial_Pancreas_merged_full_scanvi_final.RDS"),
  Prostate      = file.path(base, "Prostate/Epithelial_Prostate_merged_full_scanvi_final.RDS"),
  Sarcoma       = file.path(base, "Sarcoma/Epithelial_Sarcoma_merged_full_scanvi_final.RDS"),
  Skin          = file.path(base, "Skin/Epithelial_Skin_merged_full_scanvi_final.RDS")
)

rows <- list()

for (nm in names(objs)) {
  path <- objs[[nm]]
  cat("\n=================================================================\n")
  cat("### ", nm, "\n### ", path, "\n")
  cat("=================================================================\n")
  if (!file.exists(path)) { cat("  MISSING FILE\n"); next }

  res <- tryCatch({
    obj <- readRDS(path)
    md  <- obj@meta.data
    reds <- SeuratObject::Reductions(obj)
    red_dims <- sapply(reds, function(r) tryCatch(ncol(Embeddings(obj, r)), error = function(e) NA))
    cat("  n_cells:", ncol(obj), "  n_genes:", nrow(obj), "\n")
    cat("  reductions:", paste(sprintf("%s(%s)", reds, red_dims), collapse = ", "), "\n")
    cat("  assays:", paste(Assays(obj), collapse = ", "), "\n")
    for (a in Assays(obj)) {
      ly <- tryCatch(SeuratObject::Layers(obj[[a]]), error = function(e) NA)
      cat("    assay", a, "layers:", paste(ly, collapse = ", "), "\n")
    }
    cat("  n VariableFeatures(RNA):", length(tryCatch(VariableFeatures(obj[["RNA"]]), error = function(e) character(0))), "\n")
    sd_dim <- tryCatch(dim(GetAssayData(obj, assay = "RNA", layer = "scale.data")), error = function(e) NA)
    cat("  RNA scale.data dim:", paste(sd_dim, collapse = " x "), "\n")
    cat("  meta.data columns (", ncol(md), "):\n"); print(colnames(md))

    has_sc <- "scimilarity_celltype" %in% colnames(md)
    cat("  scimilarity_celltype present:", has_sc, "\n")
    sc_na <- NA; sc_ndist <- NA
    if (has_sc) {
      v <- as.character(md$scimilarity_celltype)
      sc_na <- sum(is.na(v) | v == "" | v == "NA")
      sc_ndist <- length(unique(v[!(is.na(v) | v == "" | v == "NA")]))
      cat("  scimilarity_celltype: NA/blank =", sc_na, "/", length(v), "  n_distinct =", sc_ndist, "\n")
      print(head(sort(table(v), decreasing = TRUE), 15))
    }
    for (cand in c("scanvi_pred", "Consensus_Cell_Type", "consensus_meta_type",
                   "Consensus_meta_type", "meta_type", "cell_type")) {
      if (cand %in% colnames(md))
        cat("  [also present] ", cand, " n_distinct =", length(unique(as.character(md[[cand]]))), "\n")
    }

    has_pn <- "project_name" %in% colnames(md)
    cat("  project_name present:", has_pn, "\n")
    nb <- NA; minb <- NA
    if (has_pn) {
      tb <- sort(table(as.character(md$project_name)))
      nb <- length(tb); minb <- as.integer(tb[1])
      cat("  n batches =", nb, "  smallest batch =", minb, " (", names(tb)[1], ")\n")
      cat("  batch size quantiles:\n"); print(quantile(as.integer(tb), c(0, .1, .25, .5, .75, 1)))
    }

    data.frame(
      cancer = nm, file = basename(path), n_cells = ncol(obj), n_genes = nrow(obj),
      reductions = paste(sprintf("%s(%s)", reds, red_dims), collapse = "|"),
      has_harmony = "harmony" %in% reds, has_scanvi = "scanvi" %in% reds,
      has_mnn = "mnn" %in% reds, has_pca_combat = "pca.combat" %in% reds,
      rna_layers = paste(tryCatch(SeuratObject::Layers(obj[["RNA"]]), error = function(e) NA), collapse = "|"),
      n_varfeat = length(tryCatch(VariableFeatures(obj[["RNA"]]), error = function(e) character(0))),
      scale_data_dim = paste(sd_dim, collapse = "x"),
      has_scimilarity_celltype = has_sc, scimilarity_na = sc_na, scimilarity_ndistinct = sc_ndist,
      has_project_name = has_pn, n_batches = nb, smallest_batch = minb,
      stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    data.frame(cancer = nm, file = basename(path), error = conditionMessage(e), stringsAsFactors = FALSE)
  })
  rows[[nm]] <- res
  rm(list = intersect("obj", ls())); gc()
}

all_cols <- unique(unlist(lapply(rows, colnames)))
rows2 <- lapply(rows, function(d) { for (c in setdiff(all_cols, colnames(d))) d[[c]] <- NA; d[all_cols] })
summ <- do.call(rbind, rows2)
write.csv(summ, out_csv, row.names = FALSE)
cat("\n\nWROTE:", out_csv, "\n"); print(summ)
