#!/usr/bin/env Rscript
# inspect_immune_annotations.R
#
# Lightweight inspection of the immune lineage Seurat objects to decide:
#   (1) which obs column to use as Spectra's cell_type_key
#   (2) whether features are gene SYMBOLS (required to match Hallmark/TME priors)
#   (3) that 'counts' and 'data' layers exist for h5ad export
#
# For each immune lineage we load the SMALLEST cancer file (cheap) and print
# candidate cell-type columns + their value counts. Column schema is assumed
# consistent across cancers within a lineage.
#
# Usage:
#   Rscript inspect_immune_annotations.R \
#     --input_dir /ofir_grp/users/hadasra/temp_objects/clustered_lineages

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option("--input_dir", type = "character",
              default = "/ofir_grp/users/hadasra/temp_objects/clustered_lineages",
              help = "Directory with *_clustered.RDS files")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Smallest cancer file per immune lineage (chosen for cheap loading)
files <- c(
  "B-cell"        = "B-cell_Sarcoma_clustered.RDS",
  "HSPC"          = "Hematopoietic_Stem_and_Progenitor_cells_Lung_clustered.RDS",
  "NK-cells"      = "NK-cells_Skin_clustered.RDS",
  "T-cells"       = "T-cells_Sarcoma_clustered.RDS",
  "Myeloid"       = "Myeloid_Prostate_clustered.RDS"
)

# Columns most likely to hold immune cell-type annotations
candidate_cols <- c(
  "Consensus_Cell_Type", "cell_type", "cell_name",
  "scimilarity_celltype", "scimilarity_celltype_filtered",
  "SingleR.labels", "sctype", "majority_cell_type",
  "maj2_cell", "maj3_cell", "maj6"
)

cat("==========================================================\n")
cat("Spectra annotation inspection\n")
cat("==========================================================\n")

for (lin in names(files)) {
  f <- file.path(opt$input_dir, files[[lin]])
  cat("\n\n########################################################\n")
  cat(sprintf("# LINEAGE: %s  (%s)\n", lin, files[[lin]]))
  cat("########################################################\n")
  if (!file.exists(f)) { cat("  FILE NOT FOUND, skipping\n"); next }

  obj <- readRDS(f)
  cat(sprintf("dims (genes x cells): %d x %d\n", nrow(obj), ncol(obj)))

  # --- feature naming: symbols vs ensembl ---
  feats <- head(rownames(obj), 15)
  cat("first 15 feature names: ", paste(feats, collapse = ", "), "\n")
  n_ens <- sum(grepl("^ENSG[0-9]+", rownames(obj)))
  cat(sprintf("features matching ENSG\\d+: %d / %d  -> %s\n",
              n_ens, nrow(obj),
              ifelse(n_ens > nrow(obj) * 0.5, "ENSEMBL IDs (need symbol mapping!)",
                     "appear to be SYMBOLS (good)")))

  # --- assay / layers present ---
  da <- DefaultAssay(obj)
  cat("default assay: ", da, " | assays: ", paste(Assays(obj), collapse = ", "), "\n")
  if (inherits(obj[["RNA"]], "Assay5")) {
    cat("RNA is Assay5; layers: ", paste(Layers(obj[["RNA"]]), collapse = ", "), "\n")
  } else {
    cat("RNA is classic Assay (has counts + data slots)\n")
  }

  # --- candidate cell-type columns ---
  meta <- obj@meta.data
  present <- intersect(candidate_cols, colnames(meta))
  cat("\ncandidate cell-type columns present: ", paste(present, collapse = ", "), "\n")
  for (cc in present) {
    v <- as.character(meta[[cc]])
    tab <- sort(table(v), decreasing = TRUE)
    cat(sprintf("\n  --- %s : %d unique values ---\n", cc, length(tab)))
    print(utils::head(tab, 25))
  }
  rm(obj); gc()
}

cat("\n\n==========================================================\n")
cat("DONE. Pick the column with clean, fine-grained immune cell\n")
cat("types and gene SYMBOL features for cell_type_key.\n")
cat("==========================================================\n")
