#!/usr/bin/env Rscript
# extract_breast_benchmark_composition.R
#
# One-off audit script: reads the 7 Breast-cancer step4 objects used in the
# Ch.3 batch-correction benchmark (B-cell, Endothelial, Fibroblasts, HSPC,
# Myeloid, NK-cells, T-cells) and reports per-object and pooled cell counts,
# unique studies (project_name), and unique patients, so the real numbers
# behind the thesis's "44-study multi-platform dataset" claim can be checked.

suppressPackageStartupMessages(library(Seurat))

BASE <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast"

CELL_TYPES <- c(
  "B-cell", "Endothelial", "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells",
  "Myeloid", "NK-cells", "T-cells"
)

all_studies  <- character(0)
all_patients <- character(0)
total_cells  <- 0L
rows <- list()

for (cel in CELL_TYPES) {
  rds_path <- file.path(BASE, paste0(cel, "_Breast_merged_step4_harmony_dev015_nocca.RDS"))
  cat(sprintf("[READ] %s\n", rds_path))

  if (!file.exists(rds_path)) {
    cat("  MISSING\n")
    rows[[cel]] <- data.frame(cell_type = cel, n_cells = NA, n_studies = NA,
                               n_patients = NA, has_patient_col = NA)
    next
  }

  obj <- readRDS(rds_path)
  meta <- obj@meta.data

  n_cells   <- ncol(obj)
  studies   <- unique(as.character(meta$project_name))
  n_studies <- length(studies)

  has_patient_col <- "patient" %in% colnames(meta)
  patients <- character(0)
  if (has_patient_col) {
    patients <- unique(as.character(meta$patient))
    patients <- patients[!is.na(patients) & patients != ""]
  }

  cat(sprintf("  cells=%d  studies=%d  patient_col=%s  n_patients=%d\n",
              n_cells, n_studies, has_patient_col, length(patients)))

  total_cells  <- total_cells + n_cells
  all_studies  <- union(all_studies, studies)
  all_patients <- union(all_patients, patients)

  rows[[cel]] <- data.frame(cell_type = cel, n_cells = n_cells,
                             n_studies = n_studies,
                             n_patients = length(patients),
                             has_patient_col = has_patient_col)

  rm(obj); gc(verbose = FALSE)
}

summary_df <- do.call(rbind, rows)
print(summary_df)

cat("\n=== POOLED (union across all 7 cell-type objects) ===\n")
cat(sprintf("Total cells (sum across objects, cells may double-count if shared): %d\n", total_cells))
cat(sprintf("Unique studies (project_name) across all 7 objects: %d\n", length(all_studies)))
cat("Studies:\n")
print(sort(all_studies))
cat(sprintf("Unique patients across all 7 objects (only where 'patient' col present/non-empty): %d\n",
            length(all_patients)))

out_csv <- file.path(BASE, "ch3_benchmark_composition_summary.csv")
write.csv(summary_df, out_csv, row.names = FALSE)
cat(sprintf("\nWrote summary to: %s\n", out_csv))

writeLines(sort(all_studies), file.path(BASE, "ch3_benchmark_unique_studies.txt"))
