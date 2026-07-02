#!/usr/bin/env Rscript
# 01_export_counts_for_cnmf.R
#
# Ch.5 (GEP programs) — cNMF step 1/4. Extracts raw integer counts from all
# Epithelial/Cancer _merged.RDS objects (the tumor-cell compartment; see
# ../filter_brain_malignant.R for how Brain's is produced), finds common
# genes, merges into one sparse matrix, and saves:
#   counts.mtx       - raw counts (genes x cells)
#   genes.tsv        - gene names
#   barcodes.tsv     - cell barcodes
#   cell_metadata.csv - cancer_type, project_name, scanvi_pred per cell
#
# Inputs:  ATLAS_DATA_DIR/merged_by_query/<Cancer>/Epithelial_<Cancer>_merged*.RDS
# Outputs: ATLAS_DATA_DIR/cnmf_input/{counts.mtx,genes.tsv,barcodes.tsv,cell_metadata.csv}
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript 01_export_counts_for_cnmf.R

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
.libPaths(c(userlib[nzchar(userlib)], "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

atlas_data_dir <- Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects")
BASE    <- file.path(atlas_data_dir, "merged_by_query")
OUT_DIR <- file.path(atlas_data_dir, "cnmf_input")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Map: cancer_type -> merged RDS (use _scanvi_final for metadata, _merged for raw counts)
OBJECTS <- list(
  Colorectal    = list(counts = file.path(BASE, "Colorectal/Epithelial_Colorectal_merged.RDS"),
                       final  = file.path(BASE, "Colorectal/Epithelial_Colorectal_merged_full_scanvi_final.RDS")),
  Kidney        = list(counts = file.path(BASE, "Kidney/Epithelial_Kidney_merged.RDS"),
                       final  = file.path(BASE, "Kidney/Epithelial_Kidney_merged_full_scanvi_final.RDS")),
  Liver_Biliary = list(counts = file.path(BASE, "Liver_Biliary/Epithelial_Liver_Biliary_merged.RDS"),
                       final  = file.path(BASE, "Liver_Biliary/Epithelial_Liver_Biliary_merged_full_scanvi_final.RDS")),
  Lung          = list(counts = file.path(BASE, "Lung/Epithelial_Lung_merged.RDS"),
                       final  = file.path(BASE, "Lung/Epithelial_Lung_merged_full_scanvi_final.RDS")),
  Ovarian       = list(counts = file.path(BASE, "Ovarian/Epithelial_Ovarian_merged.RDS"),
                       final  = file.path(BASE, "Ovarian/Epithelial_Ovarian_merged_full_scanvi_final.RDS")),
  Pancreas      = list(counts = file.path(BASE, "Pancreas/Epithelial_Pancreas_merged.RDS"),
                       final  = file.path(BASE, "Pancreas/Epithelial_Pancreas_merged_full_scanvi_final.RDS")),
  Prostate      = list(counts = file.path(BASE, "Prostate/Epithelial_Prostate_merged.RDS"),
                       final  = file.path(BASE, "Prostate/Epithelial_Prostate_merged_full_scanvi_final.RDS")),
  Sarcoma       = list(counts = file.path(BASE, "Sarcoma/Epithelial_Sarcoma_merged.RDS"),
                       final  = file.path(BASE, "Sarcoma/Epithelial_Sarcoma_merged_full_scanvi_final.RDS")),
  Skin          = list(counts = file.path(BASE, "Skin/Epithelial_Skin_merged.RDS"),
                       final  = file.path(BASE, "Skin/Epithelial_Skin_merged_full_scanvi_final.RDS")),
  Brain         = list(counts = file.path(BASE, "Brain/Epithelial_Brain_merged_full_scanvi_final.RDS"),
                       final  = file.path(BASE, "Brain/Epithelial_Brain_merged_full_scanvi_final.RDS"))
  # Breast added separately once pipeline finishes
)

all_counts <- list()
all_meta   <- list()

for (cancer in names(OBJECTS)) {
  paths <- OBJECTS[[cancer]]
  log_msg(sprintf("=== %s ===", cancer))

  # Load the counts object
  log_msg(sprintf("  Loading: %s", basename(paths$counts)))
  obj <- readRDS(paths$counts)
  log_msg(sprintf("  Cells: %d  Genes: %d", ncol(obj), nrow(obj)))

  # Extract raw counts
  lyrs <- Layers(obj[["RNA"]])
  count_lyrs <- grep("^counts", lyrs, value = TRUE)

  if (length(count_lyrs) == 0) {
    # Try old-style slot access
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  } else if (length(count_lyrs) == 1) {
    counts <- LayerData(obj, assay = "RNA", layer = count_lyrs[1])
  } else {
    counts <- do.call(cbind, lapply(count_lyrs, function(l)
      LayerData(obj, assay = "RNA", layer = l)))
  }

  # Validate raw integer counts
  sample_vals <- counts@x[seq_len(min(200, length(counts@x)))]
  is_integer  <- all(sample_vals == floor(sample_vals))
  log_msg(sprintf("  Counts are integers: %s  (sample range: %.2f - %.2f)",
                  is_integer, min(sample_vals), max(sample_vals)))
  if (!is_integer) {
    log_msg("  WARNING: counts are log-normalized — back-transforming with expm1() + round()")
    counts@x <- round(expm1(counts@x))
    counts@x[counts@x < 0] <- 0
    log_msg(sprintf("  After back-transform sample: %s", paste(round(counts@x[1:5], 1), collapse=", ")))
  }

  # Pull metadata from final object for scanvi labels
  meta <- data.frame(
    cell_barcode = colnames(counts),
    cancer_type  = cancer,
    row.names    = colnames(counts),
    stringsAsFactors = FALSE
  )

  if (file.exists(paths$final) && paths$final != paths$counts) {
    log_msg(sprintf("  Loading final object for metadata: %s", basename(paths$final)))
    final_obj <- readRDS(paths$final)
    common_cells <- intersect(colnames(counts), colnames(final_obj))
    log_msg(sprintf("  Cells with scanvi metadata: %d / %d", length(common_cells), ncol(counts)))

    for (col in c("project_name", "scanvi_pred", "Consensus_Cell_Type", "orig.ident")) {
      if (col %in% colnames(final_obj@meta.data)) {
        meta[common_cells, col] <- final_obj@meta.data[common_cells, col]
      }
    }
    rm(final_obj); gc()
  } else {
    # Brain: final IS the object, pull metadata directly
    for (col in c("project_name", "scanvi_pred", "Consensus_Cell_Type", "orig.ident")) {
      if (col %in% colnames(obj@meta.data)) {
        meta[, col] <- obj@meta.data[colnames(counts), col]
      }
    }
  }

  all_counts[[cancer]] <- counts
  all_meta[[cancer]]   <- meta
  rm(obj, counts); gc()
  log_msg(sprintf("  Done: %d cells", nrow(meta)))
}

# ── Find common genes ──────────────────────────────────────────────────────
log_msg("Finding common genes across all cancer types...")
common_genes <- Reduce(intersect, lapply(all_counts, rownames))
log_msg(sprintf("Common genes: %d", length(common_genes)))

# ── Merge ──────────────────────────────────────────────────────────────────
log_msg("Merging count matrices...")
merged_counts <- do.call(cbind, lapply(all_counts, function(m) m[common_genes, ]))
merged_meta   <- do.call(rbind, all_meta)
rm(all_counts); gc()

log_msg(sprintf("Final matrix: %d genes x %d cells", nrow(merged_counts), ncol(merged_counts)))

# ── Write outputs ──────────────────────────────────────────────────────────
log_msg("Writing outputs...")
writeMM(merged_counts, file.path(OUT_DIR, "counts.mtx"))
write.table(data.frame(gene = rownames(merged_counts)),
            file.path(OUT_DIR, "genes.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(data.frame(barcode = colnames(merged_counts)),
            file.path(OUT_DIR, "barcodes.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.csv(merged_meta, file.path(OUT_DIR, "cell_metadata.csv"), row.names=FALSE)

log_msg(sprintf("Done! Output: %s", OUT_DIR))
log_msg(sprintf("  counts.mtx      : %d genes x %d cells", nrow(merged_counts), ncol(merged_counts)))
log_msg(sprintf("  genes.tsv       : %d genes", length(common_genes)))
log_msg(sprintf("  barcodes.tsv    : %d cells", ncol(merged_counts)))
log_msg(sprintf("  cell_metadata.csv: %d rows", nrow(merged_meta)))
