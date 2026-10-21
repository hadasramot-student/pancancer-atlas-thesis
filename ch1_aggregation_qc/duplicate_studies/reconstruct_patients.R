suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr)})

worktree <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/cross-repo-duplicate-check-e81250"
project_names <- readLines(file.path(worktree, "work/final_atlas_project_names.txt"))
project_names <- project_names[nchar(trimws(project_names)) > 0]
cat("n final-atlas project_name values (raw, uncollapsed):", length(project_names), "\n")

classify <- function(pn) {
  ifelse(grepl("^(SRP|ERP|DRP)", pn), "scbasecount",
  ifelse(grepl("_ET_AL_|_et_al_", pn), "breastscape", "3ca"))
}
repo <- classify(project_names)
cat("breakdown:\n"); print(table(repo))

targets <- c("WU_SZ_ET_AL_2021_P15_filtered", "GAO_R_ET_AL_2021_P17_filtered",
             "GRIFFITHS_JI_ET_AL_2021_P18_filtered", "GONZALEZ_H_ET_AL_2022_P31_filtered")

# ---------- 3CA: read patient from each project's Meta-data.csv ----------
idx <- read_tsv("/ofir_grp/Dropbox/Cancer_DB/SC.3CA.Tirosh/SC.3CA.Tirosh.project_file_index.tsv",
                 comment = "#", show_col_types = FALSE)
idx <- idx %>% distinct(project_name, metadata_file)

old_prefix <- "/local/ofir/OCohenLab/Dropbox"
new_prefix <- "/ofir_grp/Dropbox"

tca_projects <- project_names[repo == "3ca"]
tca_results <- list()
for (pn in tca_projects) {
  # handle _1/_2/etc split-child suffixes (Dong2020_Neuroendocrine_1 etc.) by
  # stripping to the base project name for the index lookup
  base_pn <- sub("_[0-9]+$", "", pn)
  rows <- idx %>% filter(project_name == base_pn)
  if (nrow(rows) == 0) {
    tca_results[[pn]] <- list(n_patients = NA, note = "not found in project_file_index.tsv")
    next
  }
  mf <- gsub(old_prefix, new_prefix, rows$metadata_file[1], fixed = TRUE)
  if (!file.exists(mf)) {
    tca_results[[pn]] <- list(n_patients = NA, note = paste("Meta-data.csv not found:", mf))
    next
  }
  md <- tryCatch(read_csv(mf, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(md) || !("patient" %in% tolower(names(md)))) {
    tca_results[[pn]] <- list(n_patients = NA, note = "no 'patient' column in Meta-data.csv")
    next
  }
  names(md) <- tolower(names(md))
  npat <- md$patient %>% unique() %>% na.omit() %>% length()
  tca_results[[pn]] <- list(n_patients = npat, note = "ok")
}

cat("\n=== 3CA per-project patient counts ===\n")
for (pn in names(tca_results)) cat(sprintf("  %-45s n_patients=%s  (%s)\n", pn,
    ifelse(is.na(tca_results[[pn]]$n_patients), "NA", tca_results[[pn]]$n_patients), tca_results[[pn]]$note))

tca_ok <- Filter(function(x) !is.na(x$n_patients), tca_results)
tca_sum <- sum(sapply(tca_ok, function(x) x$n_patients))
tca_failed <- length(tca_results) - length(tca_ok)
cat("\n3CA: sum of per-project distinct patients (", length(tca_ok), "of", length(tca_results), "projects readable):", tca_sum, "\n")
cat("3CA projects where patient count could NOT be determined:", tca_failed, "\n")

# ---------- Breastscape: read patient from each project's own _step4.RDS.tsv ----------
bs_projects <- project_names[repo == "breastscape"]
tsv_dir <- "/ofir_grp/users/hadasra/temp_objects"
bs_results <- list()
for (pn in bs_projects) {
  f <- file.path(tsv_dir, paste0(pn, "_step4.RDS.tsv"))
  if (!file.exists(f)) {
    bs_results[[pn]] <- list(n_patients = NA, note = paste("tsv not found:", f))
    next
  }
  # read just the patient column for speed
  hdr <- names(read_tsv(f, n_max = 0, show_col_types = FALSE))
  if (!"patient" %in% hdr) {
    bs_results[[pn]] <- list(n_patients = NA, note = "no 'patient' column in tsv")
    next
  }
  pcol <- read_tsv(f, col_select = "patient", show_col_types = FALSE)
  npat <- pcol$patient %>% unique() %>% na.omit() %>% length()
  bs_results[[pn]] <- list(n_patients = npat, note = "ok")
}

cat("\n=== Breastscape per-project patient counts ===\n")
for (pn in names(bs_results)) cat(sprintf("  %-45s n_patients=%s  (%s)\n", pn,
    ifelse(is.na(bs_results[[pn]]$n_patients), "NA", bs_results[[pn]]$n_patients), bs_results[[pn]]$note))

bs_ok <- Filter(function(x) !is.na(x$n_patients), bs_results)
bs_sum_all <- sum(sapply(bs_ok, function(x) x$n_patients))
bs_failed <- length(bs_results) - length(bs_ok)
cat("\nBreastscape: sum of per-project distinct patients (", length(bs_ok), "of", length(bs_results), "projects readable):", bs_sum_all, "\n")
cat("Breastscape projects where patient count could NOT be determined:", bs_failed, "\n")

bs_sum_excl_targets <- sum(sapply(bs_ok[setdiff(names(bs_ok), targets)], function(x) x$n_patients))
cat("Breastscape sum EXCLUDING the 4 confirmed duplicate studies:", bs_sum_excl_targets, "\n")
cat("(patients removed by dropping the 4 duplicates:", bs_sum_all - bs_sum_excl_targets, ")\n")

# ---------- ScBaseCount: no patient concept ----------
sc_projects <- project_names[repo == "scbasecount"]
cat("\n=== ScBaseCount ===\n")
cat("n projects:", length(sc_projects), "-- no patient-level field exists in SRA metadata (confirmed earlier this session).\n")
cat("NOT counted toward any patient total below.\n")

# ---------- Grand totals ----------
cat("\n=================== GRAND TOTALS (3CA + Breastscape only; ScBaseCount excluded, no patient field) ===================\n")
grand_current <- tca_sum + bs_sum_all
grand_after_dedup <- tca_sum + bs_sum_excl_targets
cat("Current (all 3CA + Breastscape projects, incl. the 4 duplicates counted twice):", grand_current, "\n")
cat("After removing the 4 confirmed Breastscape duplicates:", grand_after_dedup, "\n")
cat("Delta:", grand_current - grand_after_dedup, "\n")
cat("\nTarget from thesis: 11769\n")
cat("Distance (current vs target):", abs(grand_current - 11769), "\n")
cat("Distance (after-dedup vs target):", abs(grand_after_dedup - 11769), "\n")

cat("\nNote: this sum assumes no patient ID string is reused across two DIFFERENT (non-duplicate) studies.",
    "Not independently verified here beyond the 4 already-confirmed duplicate pairs. Also note: 3CA counts are",
    "PER-STUDY unique patient counts, not verified against actual per-cell coverage in the final atlas (i.e. if a",
    "3CA project's Meta-data.csv lists a patient with 0 cells surviving QC, they'd still be counted here) -- a",
    "conservative upper-bound-ish estimate, not a cell-exact figure.\n")

cat("\nDone.\n")
