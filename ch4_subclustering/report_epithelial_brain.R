#!/usr/bin/env Rscript
# report_epithelial_brain.R
#
# Post-hoc report for the Ch.4 subclustering of Epithelial_Brain (run after
# cluster_lineage.R + annotate_lineages.R). Read-only on all inputs.
#
#   1. Per-cluster n_cells + Scimilarity majority label / confidence at
#      louvain_res0.50 and louvain_res0.30; flags clusters < 50% confidence.
#   2. Composition of the Brain object (studies, samples, repositories) and
#      QC metrics, next to the same numbers for the other Epithelial x cancer
#      objects (taken from the atlas metadata table, not the 30+ GB RDSs).
#
# Usage:
#   Rscript report_epithelial_brain.R --clustered <..._clustered.RDS> \
#     --scim_tsv <..._scimilarity_summary.tsv> --atlas_meta <all_cells...RDS> \
#     --out_dir <dir>

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--clustered",  type = "character"),
  make_option("--scim_tsv",   type = "character"),
  make_option("--atlas_meta", type = "character"),
  make_option("--out_dir",    type = "character"),
  make_option("--conf_threshold", type = "double", default = 50)
)))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$out_dir, f)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── 1. Clusters + Scimilarity ─────────────────────────────────────────────────
scim <- read.delim(opt$scim_tsv, stringsAsFactors = FALSE)
for (res in c("louvain_res0.50", "louvain_res0.30")) {
  s <- scim[scim$resolution == res, ]
  s <- s[order(as.integer(s$cluster)), ]
  s$flag_lt50 <- s$confidence_pct < opt$conf_threshold
  write.table(s, out(sprintf("Epithelial_Brain_%s_clusters.tsv", res)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("\n==== %s: %d clusters, %d below %.0f%% confidence ====\n",
              res, nrow(s), sum(s$flag_lt50), opt$conf_threshold))
  print(s, row.names = FALSE)
}

# ── 2. Brain object composition + QC ──────────────────────────────────────────
log_msg("Loading clustered object")
obj <- readRDS(opt$clustered)
md  <- obj@meta.data
rm(obj); invisible(gc())
cat("\nMetadata columns:\n"); print(colnames(md))

tab_top <- function(x, n = 30) {
  t <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
  data.frame(value = names(t), n = as.integer(t),
             pct = round(100 * as.integer(t) / length(x), 2))[seq_len(min(n, length(t))), ]
}
for (col in intersect(c("project_name", "sample", "patient", "batch", "cell_type",
                        "Consensus_Cell_Type", "scimilarity_celltype",
                        "cancer_type", "technology", "platform"), colnames(md))) {
  cat(sprintf("\n-- %s (%d distinct, %.1f%% NA) --\n", col,
              length(unique(md[[col]])), 100 * mean(is.na(md[[col]]))))
  print(tab_top(md[[col]]), row.names = FALSE)
}

# Study mix per cluster: is a cluster essentially one study?
if ("project_name" %in% colnames(md)) {
  for (res in c("louvain_res0.50", "louvain_res0.30")) {
    ct <- table(md[[res]], md$project_name)
    dom <- data.frame(
      cluster   = rownames(ct),
      n_cells   = rowSums(ct),
      n_studies = rowSums(ct > 0),
      top_study = colnames(ct)[apply(ct, 1, which.max)],
      top_study_pct = round(100 * apply(ct, 1, max) / rowSums(ct), 1))
    dom <- dom[order(as.integer(dom$cluster)), ]
    write.table(dom, out(sprintf("Epithelial_Brain_%s_study_mix.tsv", res)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("\n-- %s study mix per cluster --\n", res)); print(dom, row.names = FALSE)
  }
}

qc_cols <- grep("^(nCount_RNA|nFeature_RNA|percent[._]?mt|pct_counts_mt|percent\\.mito)$",
                colnames(md), value = TRUE, ignore.case = TRUE)
qc_summary <- function(d, cols) sapply(cols, function(c) {
  v <- suppressWarnings(as.numeric(d[[c]]))
  c(median = median(v, na.rm = TRUE), q10 = quantile(v, .1, na.rm = TRUE, names = FALSE),
    q90 = quantile(v, .9, na.rm = TRUE, names = FALSE), pct_NA = 100 * mean(is.na(v)))
})
cat("\n-- Brain QC metrics --\n"); print(qc_summary(md, qc_cols))

# ── 3. Comparison against the other Epithelial x cancer objects ───────────────
log_msg("Loading atlas metadata")
am <- as.data.frame(readRDS(opt$atlas_meta))
repo_of <- function(p) ifelse(grepl("_ET_AL_.*_filtered$", p), "Breastscape",
                       ifelse(grepl("^(SRP|ERP|DRP)", p), "ScBaseCount", "3CA"))
epi <- am[am$Consensus_Cell_Type == "Epithelial", ]
am_qc <- intersect(qc_cols, colnames(am))
cmp_row <- function(d, name) {
  sm <- if ("sample" %in% colnames(d)) d$sample else rep(NA, nrow(d))
  pj <- d$project_name
  rp <- table(factor(repo_of(pj), levels = c("3CA", "Breastscape", "ScBaseCount")))
  r <- data.frame(object = name, n_cells = nrow(d),
                  n_studies = length(unique(pj)),
                  n_samples = length(unique(na.omit(sm))),
                  sample_pct_NA = round(100 * mean(is.na(sm)), 1),
                  largest_study_pct = round(100 * max(table(pj)) / length(pj), 1),
                  pct_3CA = round(100 * rp[["3CA"]] / nrow(d), 1),
                  pct_Breastscape = round(100 * rp[["Breastscape"]] / nrow(d), 1),
                  pct_ScBaseCount = round(100 * rp[["ScBaseCount"]] / nrow(d), 1))
  for (c in am_qc) r[[paste0("median_", c)]] <- median(as.numeric(d[[c]]), na.rm = TRUE)
  r
}
cancers <- sort(setdiff(unique(epi$cancer_type), c(NA, "unknown", "Other_Models")))
cmp <- do.call(rbind, lapply(cancers, function(ca)
  cmp_row(epi[epi$cancer_type == ca, ], paste0("Epithelial_", ca))))
brain_cells <- am[am$cancer_type == "Brain" & !is.na(am$cell_type) & am$cell_type == "Malignant", ]
cmp <- rbind(cmp, cmp_row(brain_cells, "Brain Malignant (-> Epithelial_Brain)"))
cat("\n-- Brain Malignant cells: Consensus_Cell_Type in atlas --\n")
print(tab_top(brain_cells$Consensus_Cell_Type))
cmp <- cmp[order(-cmp$n_cells), ]
write.table(cmp, out("Epithelial_objects_composition_qc_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n-- Epithelial objects (atlas meta, Consensus_Cell_Type == Epithelial) --\n")
print(cmp, row.names = FALSE)
log_msg("Done. Outputs in ", opt$out_dir)
