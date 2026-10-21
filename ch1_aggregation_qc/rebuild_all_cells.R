#!/usr/bin/env Rscript
# rebuild_all_cells.R
#
# Regenerate all_cells.RDS by reading all *_step4.RDS.tsv files.
# This is the source of truth for Consensus_Cell_Type — the TSVs are written
# directly from each per-project Seurat object and are always current.
#
# Usage:
#   Rscript rebuild_all_cells.R \
#     --tsv_dir   /mnt/.../temp_objects \
#     --out_rds   /mnt/.../temp_objects/all_cells.RDS \
#     [--backup]

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option("--tsv_dir", type = "character",
    default = "/ofir_grp/users/hadasra/temp_objects",
    help = "Directory containing *_step4.RDS.tsv files"),
  make_option("--out_rds", type = "character",
    default = "/ofir_grp/users/hadasra/temp_objects/all_cells.RDS",
    help = "Output path for regenerated all_cells.RDS"),
  make_option("--backup", action = "store_true", default = FALSE,
    help = "Save a backup of existing all_cells.RDS before overwriting")
)

opt <- parse_args(OptionParser(option_list = option_list))

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

# ── Backup existing file ──────────────────────────────────────────────────────
if (isTRUE(opt$backup) && file.exists(opt$out_rds)) {
  bak <- paste0(opt$out_rds, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy(opt$out_rds, bak)
  log_msg("Backed up existing all_cells.RDS → ", bak)
}

# ── Discover TSV files ────────────────────────────────────────────────────────
tsv_files <- list.files(opt$tsv_dir, pattern = "_step4\\.RDS\\.tsv$", full.names = TRUE)
log_msg("Found ", length(tsv_files), " TSV files in: ", opt$tsv_dir)
if (length(tsv_files) == 0) stop("No *_step4.RDS.tsv files found in: ", opt$tsv_dir)

# ── Read and bind ─────────────────────────────────────────────────────────────
log_msg("Reading TSV files ...")

all_cells <- lapply(tsv_files, function(f) {
  tryCatch(
    read_tsv(f, show_col_types = FALSE) %>%
      mutate(patient = as.character(NA)) %>%
      select(-any_of("predicted_doublet")),
    error = function(e) {
      message("  WARNING: failed to read ", basename(f), ": ", conditionMessage(e))
      NULL
    }
  )
})

# Report any failures
failed <- sum(sapply(all_cells, is.null))
if (failed > 0) log_msg("WARNING: ", failed, " TSV(s) failed to read — excluded")

all_cells <- bind_rows(Filter(Negate(is.null), all_cells))
log_msg("Total cells: ", nrow(all_cells), " from ", length(tsv_files) - failed, " projects")

# ── Sanity check ──────────────────────────────────────────────────────────────
if (!"Consensus_Cell_Type" %in% colnames(all_cells))
  stop("'Consensus_Cell_Type' column not found in TSVs — check TSV format")

log_msg("Consensus_Cell_Type distribution:")
cct_tbl <- sort(table(all_cells$Consensus_Cell_Type), decreasing = TRUE)
for (i in seq_along(cct_tbl)) {
  cat(sprintf("  %-50s %d\n", names(cct_tbl)[i], cct_tbl[i]))
}

if ("cancer_type" %in% colnames(all_cells)) {
  log_msg("Cancer types: ", paste(sort(unique(all_cells$cancer_type)), collapse = ", "))
}

# ── Save ──────────────────────────────────────────────────────────────────────
log_msg("Saving → ", opt$out_rds)
saveRDS(all_cells, opt$out_rds)
log_msg("Done. all_cells.RDS regenerated with ", nrow(all_cells), " cells.")
