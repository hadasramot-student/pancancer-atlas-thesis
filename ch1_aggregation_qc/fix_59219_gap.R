#!/usr/bin/env Rscript
# fix_59219_gap.R
#
# Fixes the 59,219-cell gap in all_cells.RDS caused by a bug in
# save_metadata_tsv() (pancancer-atlas-thesis/common/functions.R:81-93):
# for multi-cancer-type studies, `cancer_tp <- unique(obj@meta.data$cancer_type)`
# returns a vector of length >1, and the subsequent
# `mutate(cancer_type = cancer_tp)` throws a hard dplyr/vctrs recycling
# error ("must be size N, not K"). This crashes 05_consensus_annotation.R
# at the tsv-write step (after the .RDS and .html are already written),
# so rebuild_all_cells.R's glob for `_step4.RDS.tsv` never finds these
# studies and they're silently absent from all_cells.RDS.
#
# This script does NOT re-run the buggy function. It reads the 3 affected
# per-project _step4.RDS objects directly, exports their metadata the way
# save_metadata_tsv() *should* have (keeping the existing per-cell
# cancer_type column intact, since these objects already carry correct
# per-cell cancer-type labels - no re-derivation needed), and appends the
# result to a NEW copy of all_cells.RDS. The original all_cells.RDS is
# never opened for writing.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

TO <- "/ofir_grp/users/hadasra/temp_objects"

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

gap_studies <- c("ERP136281", "ERP136992", "SRP288163")

# ── Step 1: load the 3 affected objects, export per-cell metadata correctly ──
new_rows_list <- list()

for (acc in gap_studies) {
  rds_path <- file.path(TO, paste0(acc, "_step4.RDS"))
  log_msg("Loading ", rds_path, " ...")
  obj <- readRDS(rds_path)

  md <- obj@meta.data
  ct_vals <- unique(md$cancer_type)
  log_msg(acc, ": ", nrow(md), " cells, distinct cancer_type values: ",
          paste(sort(ct_vals), collapse = ", "),
          " (", length(ct_vals), " distinct - this is what crashed the original save_metadata_tsv)")

  # Reproduce what save_metadata_tsv() does for every OTHER (working) study,
  # minus the buggy overwrite: rownames -> cell_id, keep cancer_type as-is
  # (it's already correct per-cell, not something to derive or guess).
  meta <- md %>%
    rownames_to_column("cell_id")

  # sanity: per-cancer-type breakdown for this study, for the report
  log_msg("  Per-cancer-type breakdown for ", acc, ":")
  print(sort(table(meta$cancer_type, useNA = "ifany"), decreasing = TRUE))

  new_rows_list[[acc]] <- meta
}

new_rows <- bind_rows(new_rows_list)
log_msg("Total new rows across the 3 studies: ", nrow(new_rows))
stopifnot(nrow(new_rows) == 59219)

# ── Step 2: load the existing (untouched) all_cells.RDS ──────────────────────
orig_path <- file.path(TO, "all_cells.RDS")
log_msg("Loading existing all_cells.RDS from ", orig_path, " (read-only, never opened for writing) ...")
all_cells_orig <- readRDS(orig_path)
log_msg("Existing all_cells.RDS: ", nrow(all_cells_orig), " rows, ", ncol(all_cells_orig), " columns")

# ── Step 3: align columns and combine ─────────────────────────────────────────
orig_cols <- colnames(all_cells_orig)
new_cols  <- colnames(new_rows)
log_msg("Columns in existing table but not in new rows: ",
        paste(setdiff(orig_cols, new_cols), collapse = ", "))
log_msg("Columns in new rows but not in existing table: ",
        paste(setdiff(new_cols, orig_cols), collapse = ", "))

# match the pattern rebuild_all_cells.R uses for freshly-read tsvs:
# patient set to NA (character), predicted_doublet dropped if present
new_rows <- new_rows %>%
  mutate(patient = as.character(NA)) %>%
  select(-any_of("predicted_doublet"))

# ── Harmonize types for every shared column before binding ───────────────────
# all_cells_orig came from read_tsv() (text-inferred types); new_rows comes
# straight from Seurat meta.data (native R types, e.g. RNA_snn_res.* stored
# as double instead of the factor read_tsv produced). Coerce any mismatched
# shared column to character on both sides so bind_rows doesn't choke - this
# only affects representation for comparison/storage, not the underlying
# per-cell values, which are unchanged either way.
shared_cols <- intersect(colnames(all_cells_orig), colnames(new_rows))
mismatched <- character(0)
for (col in shared_cols) {
  cls_orig <- class(all_cells_orig[[col]])[1]
  cls_new  <- class(new_rows[[col]])[1]
  if (cls_orig != cls_new) {
    mismatched <- c(mismatched, col)
    all_cells_orig[[col]] <- as.character(all_cells_orig[[col]])
    new_rows[[col]] <- as.character(new_rows[[col]])
  }
}
log_msg("Type-harmonized ", length(mismatched), " mismatched shared column(s): ",
        paste(mismatched, collapse = ", "))

combined <- bind_rows(all_cells_orig, new_rows)
log_msg("Combined total: ", nrow(combined), " rows (expect 8,270,370)")
stopifnot(nrow(combined) == nrow(all_cells_orig) + 59219)

# ── Step 4: save as a NEW file, original untouched ────────────────────────────
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- file.path(TO, paste0("all_cells.gapfix_", ts, ".RDS"))
saveRDS(combined, out_path)
log_msg("Saved fixed table to: ", out_path)
log_msg("Original all_cells.RDS was never modified: ", orig_path)

# ── Step 5: recompute reports ─────────────────────────────────────────────────
log_msg("\n=== Per-cancer-type counts: BEFORE vs AFTER ===")
before_ct <- table(all_cells_orig$cancer_type, useNA = "ifany")
after_ct  <- table(combined$cancer_type, useNA = "ifany")
all_types <- union(names(before_ct), names(after_ct))
ct_report <- data.frame(
  cancer_type = all_types,
  before = as.integer(before_ct[all_types]),
  after  = as.integer(after_ct[all_types])
)
ct_report$before[is.na(ct_report$before)] <- 0
ct_report$after[is.na(ct_report$after)] <- 0
ct_report$delta <- ct_report$after - ct_report$before
ct_report <- ct_report[order(-ct_report$after), ]
ct_report$pct_after <- round(100 * ct_report$after / sum(ct_report$after), 3)
print(ct_report, row.names = FALSE)
write_csv(ct_report, file.path(TO, paste0("gapfix_cancer_type_report_", ts, ".csv")))

log_msg("\n=== Consensus_Cell_Type (meta-type) counts: BEFORE vs AFTER ===")
before_cct <- table(all_cells_orig$Consensus_Cell_Type, useNA = "ifany")
after_cct  <- table(combined$Consensus_Cell_Type, useNA = "ifany")
all_cct <- union(names(before_cct), names(after_cct))
cct_report <- data.frame(
  meta_type = all_cct,
  before = as.integer(before_cct[all_cct]),
  after  = as.integer(after_cct[all_cct])
)
cct_report$before[is.na(cct_report$before)] <- 0
cct_report$after[is.na(cct_report$after)] <- 0
cct_report$delta <- cct_report$after - cct_report$before
cct_report <- cct_report[order(-cct_report$after), ]
cct_report$pct_before <- round(100 * cct_report$before / sum(cct_report$before), 3)
cct_report$pct_after  <- round(100 * cct_report$after / sum(cct_report$after), 3)
print(cct_report, row.names = FALSE)
write_csv(cct_report, file.path(TO, paste0("gapfix_meta_type_report_", ts, ".csv")))

log_msg("\nDONE. New file: ", out_path)
log_msg("Reports: ",
        file.path(TO, paste0("gapfix_cancer_type_report_", ts, ".csv")), ", ",
        file.path(TO, paste0("gapfix_meta_type_report_", ts, ".csv")))
