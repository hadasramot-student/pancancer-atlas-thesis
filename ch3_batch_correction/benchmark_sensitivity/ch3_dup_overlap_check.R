#!/usr/bin/env Rscript
# ch3_dup_overlap_check.R
#
# Chapter 3 benchmark cohort vs. the 5 cross-repository duplicate study pairs.
#
# For each of the 7 Breast compartments used in the Ch.3 batch-correction
# benchmark, load the step4 object (the exact object the composition summary
# and benchmark were built on) and tabulate project_name for every study name
# matching one of the 5 duplicate pairs. Report per-compartment and pooled
# cell counts for the 3CA "kept" copy and the Breastscape "removed" copy of
# each pair.
#
# Read-only. Writes:
#   ch3_dup_overlap_by_compartment.csv   (long: compartment x project_name x n)
#   ch3_dup_overlap_pairsummary.csv      (per pair: kept vs removed, pooled)
#   ch3_dup_overlap_check.log            (full console)

suppressPackageStartupMessages(library(Seurat))

BASE <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast"
OUT  <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/chapter3-benchmark-duplicates-3ba85e"

CELL_TYPES <- c(
  "B-cell", "Endothelial", "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells",
  "Myeloid", "NK-cells", "T-cells"
)

# regex that catches BOTH copies of each of the 5 pairs, in any spelling
DUP_REGEX <- "Wu2021|WU_SZ_ET_AL_2021|Gao2021|GAO_R_ET_AL_2021|Griffiths2021|GRIFFITHS_JI_ET_AL_2021|Gonzalez2022|GONZALEZ_H_ET_AL_2022|Savas2018|SAVAS_P_ET_AL_2018"

# map a project_name to (pair, side)
classify <- function(pn) {
  if (grepl("^Wu2021_Breast",              pn)) return(c("Wu2021",         "3CA_kept"))
  if (grepl("^WU_SZ_ET_AL_2021_P15",       pn)) return(c("Wu2021",         "Breastscape_removed"))
  if (grepl("^Gao2021_Breast",             pn)) return(c("Gao2021",        "3CA_kept"))
  if (grepl("^GAO_R_ET_AL_2021_P17",       pn)) return(c("Gao2021",        "Breastscape_removed"))
  if (grepl("^Griffiths2021_Breast",       pn)) return(c("Griffiths2021",  "3CA_kept"))
  if (grepl("^GRIFFITHS_JI_ET_AL_2021_P18",pn)) return(c("Griffiths2021",  "Breastscape_removed"))
  if (grepl("^Gonzalez2022",               pn)) return(c("Gonzalez2022",   "3CA_kept"))
  if (grepl("^GONZALEZ_H_ET_AL_2022_P31",  pn)) return(c("Gonzalez2022",   "Breastscape_removed"))
  if (grepl("^Savas2018_Breast",           pn)) return(c("Savas2018",      "3CA_kept"))
  if (grepl("^SAVAS_P_ET_AL_2018_P35",     pn)) return(c("Savas2018",      "Breastscape_removed"))
  return(c(NA, NA))
}

long_rows <- list()
grand_total <- 0L

for (cel in CELL_TYPES) {
  rds <- file.path(BASE, paste0(cel, "_Breast_merged_step4_harmony_dev015_nocca.RDS"))
  cat(sprintf("\n[READ] %s\n", rds))
  obj  <- readRDS(rds)
  meta <- obj@meta.data
  n    <- ncol(obj)
  grand_total <- grand_total + n
  cat(sprintf("  cells=%d  reductions={%s}\n", n,
              paste(SeuratObject::Reductions(obj), collapse=",")))

  pn  <- as.character(meta$project_name)
  tab <- sort(table(pn[grepl(DUP_REGEX, pn)]), decreasing = TRUE)
  if (length(tab) == 0) {
    cat("  no duplicate-pair studies in this compartment\n")
  } else {
    for (nm in names(tab)) {
      cc <- classify(nm)
      cat(sprintf("  %-40s n=%-7d  [%s / %s]\n", nm, tab[[nm]], cc[1], cc[2]))
      long_rows[[length(long_rows)+1]] <- data.frame(
        compartment = cel, project_name = nm, n_cells = as.integer(tab[[nm]]),
        pair = cc[1], side = cc[2], stringsAsFactors = FALSE)
    }
  }
  rm(obj); gc(verbose = FALSE)
}

cat(sprintf("\n=== pooled benchmark cells across 7 compartments: %d (expect 138722) ===\n", grand_total))

long <- do.call(rbind, long_rows)
write.csv(long, file.path(OUT, "ch3_dup_overlap_by_compartment.csv"), row.names = FALSE)

# pair-level summary
cat("\n=== PER-PAIR SUMMARY (pooled across compartments) ===\n")
pairs <- c("Wu2021","Gao2021","Griffiths2021","Gonzalez2022","Savas2018")
sum_rows <- list()
for (p in pairs) {
  kept    <- sum(long$n_cells[long$pair == p & long$side == "3CA_kept"])
  removed <- sum(long$n_cells[long$pair == p & long$side == "Breastscape_removed"])
  comps_k <- paste(sort(unique(long$compartment[long$pair==p & long$side=="3CA_kept"])), collapse=";")
  comps_r <- paste(sort(unique(long$compartment[long$pair==p & long$side=="Breastscape_removed"])), collapse=";")
  both <- kept > 0 && removed > 0
  cat(sprintf("  %-14s  3CA_kept=%-7d (%s)   Breastscape_removed=%-7d (%s)   BOTH_IN_BENCHMARK=%s\n",
              p, kept, comps_k, removed, comps_r, both))
  sum_rows[[p]] <- data.frame(pair = p, cells_3CA_kept = kept, cells_BS_removed = removed,
                              compartments_3CA_kept = comps_k,
                              compartments_BS_removed = comps_r,
                              both_copies_in_benchmark = both, stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, sum_rows), file.path(OUT, "ch3_dup_overlap_pairsummary.csv"), row.names = FALSE)
cat("\nWrote ch3_dup_overlap_by_compartment.csv and ch3_dup_overlap_pairsummary.csv\n")
cat("DONE\n")
