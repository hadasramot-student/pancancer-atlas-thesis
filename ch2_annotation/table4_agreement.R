## ---------------------------------------------------------------------------
## Regenerate thesis Table 4 (per-meta-type inter-method agreement) on the
## Savas-corrected final atlas.
##
## Correction applied here vs. the stale Table 4: remove the 5th duplicate
## study SAVAS_P_ET_AL_2018_P35_filtered (5,450 cells) in addition to the 4
## Breastscape duplicates already removed for the stale table. NOTHING ELSE
## is changed. Target = the atlas that produced the current Table 3
## (8,107,665 cells).
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

IN_RDS  <- "/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS"
OUT_DIR <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/table-4-savas-correction-6e1e66"

cat("Loading:", IN_RDS, "\n")
md <- readRDS(IN_RDS)
if (!is.data.frame(md)) md <- as.data.frame(md)
cat("Raw rows:", nrow(md), " | distinct project_name:", dplyr::n_distinct(md$project_name), "\n\n")

## --- the 5 cross-repository duplicate studies to exclude -------------------
dup_studies <- c(
  "WU_SZ_ET_AL_2021_P15_filtered",
  "GAO_R_ET_AL_2021_P17_filtered",
  "GRIFFITHS_JI_ET_AL_2021_P18_filtered",
  "GONZALEZ_H_ET_AL_2022_P31_filtered",
  "SAVAS_P_ET_AL_2018_P35_filtered"     # <-- the 5th duplicate, newly removed
)

cat("Cells per duplicate study (to be excluded):\n")
print(md %>% filter(project_name %in% dup_studies) %>% count(project_name))
missing_dups <- setdiff(dup_studies, unique(md$project_name))
if (length(missing_dups)) cat("\n!! WARNING - dup study id(s) not found verbatim:", paste(missing_dups, collapse=", "), "\n")

md2 <- md %>% filter(!(project_name %in% dup_studies))
cat("\nRows after excluding the 5 duplicate studies:", nrow(md2),
    "  (target 8107665; delta =", nrow(md2) - 8107665, ")\n")
cat("distinct project_name after removal:", dplyr::n_distinct(md2$project_name), "\n\n")

## --- per-cell inter-method agreement status -------------------------------
##  Full    (3/3): maj3_cell is a real meta-type (not "equivocal")
##  Partial (2/3): maj3_cell == "equivocal" but maj2_cell is a real meta-type
##  None   (<=1/3): maj2_cell == "equivocal"
md2 <- md2 %>%
  mutate(
    cell_status = case_when(
      maj3_cell != "equivocal"                          ~ "Full",
      maj3_cell == "equivocal" & maj2_cell != "equivocal" ~ "Partial",
      maj2_cell == "equivocal"                          ~ "None",
      TRUE ~ NA_character_
    )
  )
cat("NA cell_status (expect 0):", sum(is.na(md2$cell_status)), "\n\n")

## --- counts by final Consensus_Cell_Type x status ------------------------
raw_tbl <- md2 %>%
  count(Consensus_Cell_Type, cell_status) %>%
  pivot_wider(names_from = cell_status, values_from = n, values_fill = 0)
for (c in c("Full","Partial","None")) if (!c %in% names(raw_tbl)) raw_tbl[[c]] <- 0L

raw_tbl <- raw_tbl %>%
  mutate(total = Full + Partial + None) %>%
  select(Consensus_Cell_Type, Full, Partial, None, total)

cat("=== Raw table (internal Consensus_Cell_Type labels) ===\n")
print(as.data.frame(raw_tbl))
cat("\nGrand total:", sum(raw_tbl$total), " (target 8107665)\n\n")

## --- map to thesis meta-type display names + order ----------------------
disp <- c(
  "Epithelial"                             = "Epithelial/Cancer",
  "T-cells"                                = "T cells",
  "Myeloid"                                = "Myeloid",
  "B-cell"                                 = "B cells",
  "Fibroblasts"                            = "Fibroblasts/Stroma",
  "Endothelial"                            = "Endothelial/Vascular",
  "NK-cells"                               = "NK cells",
  "Hematopoietic_Stem_and_Progenitor_cells"= "HSC/Stem-like",
  "other"                                  = "Other/Unknown"
)
row_order <- c("Epithelial/Cancer","T cells","Myeloid","B cells",
               "Fibroblasts/Stroma","Endothelial/Vascular","NK cells",
               "HSC/Stem-like","Other/Unknown")

tbl <- raw_tbl %>%
  mutate(meta_type = dplyr::recode(Consensus_Cell_Type, !!!disp)) %>%
  group_by(meta_type) %>%
  summarise(Full = sum(Full), Partial = sum(Partial), None = sum(None),
            total = sum(total), .groups = "drop") %>%
  mutate(meta_type = factor(meta_type, levels = row_order)) %>%
  arrange(meta_type) %>%
  mutate(pct_full = round(100 * Full / total, 1))

## --- reconcile against current Table 3 row totals ----------------------
expected_tot <- c(
  "Epithelial/Cancer"    = 2443133,
  "T cells"              = 1510735,
  "Myeloid"              = 979202,
  "B cells"              = 371455,
  "Fibroblasts/Stroma"   = 264397,
  "Endothelial/Vascular" = 257365,
  "NK cells"             = 90334,
  "HSC/Stem-like"        = 32492,
  "Other/Unknown"        = 2158552
)
tbl$table3_total <- expected_tot[as.character(tbl$meta_type)]
tbl$reconciles   <- tbl$total == tbl$table3_total

cat("=== CORRECTED TABLE 4 (Savas-corrected atlas) ===\n")
print(as.data.frame(tbl))

cat("\n--- Reconciliation vs Table 3 ---\n")
if (all(tbl$reconciles)) {
  cat("ALL 9 rows reconcile exactly against Table 3.\n")
} else {
  bad <- tbl[!tbl$reconciles, c("meta_type","total","table3_total")]
  bad$diff <- bad$total - bad$table3_total
  cat("!! MISMATCH in", nrow(bad), "row(s):\n"); print(as.data.frame(bad))
}
gt  <- sum(tbl$total)
cat("Grand total:", gt, " vs Table 3 total 8107665  -> ",
    ifelse(gt == 8107665, "OK", paste0("MISMATCH (diff ", gt - 8107665, ")")), "\n")

## --- overall percentages (for abstract / Conclusions) ------------------
tot_full    <- sum(tbl$Full)
tot_partial <- sum(tbl$Partial)
tot_none    <- sum(tbl$None)
grand       <- tot_full + tot_partial + tot_none
other_row   <- tbl[tbl$meta_type == "Other/Unknown", ]

cat("\n=== Overall agreement (all meta-types) ===\n")
cat(sprintf("Full  (3/3) : %10d  (%.1f%%)\n", tot_full,    100*tot_full/grand))
cat(sprintf("Partial(2/3): %10d  (%.1f%%)\n", tot_partial, 100*tot_partial/grand))
cat(sprintf("None (<=1/3): %10d  (%.1f%%)\n", tot_none,    100*tot_none/grand))
cat(sprintf("Grand total : %10d\n", grand))
cat(sprintf("\nOther/Unknown share of atlas: %d  (%.1f%%)\n",
            other_row$total, 100*other_row$total/grand))
cat(sprintf("Full+Partial (some agreement): %.1f%%\n", 100*(tot_full+tot_partial)/grand))

## --- markdown rendering ----------------------------------------------------
fmt <- function(x) formatC(x, big.mark = ",", format = "d")
md_lines <- c(
  "| Meta-type | Full (3/3) | Partial (2/3) | None (<=1/3) | Row total | % Full agreement |",
  "|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(tbl))) {
  r <- tbl[i, ]
  none_cell <- if (r$None == 0) "--" else fmt(r$None)
  md_lines <- c(md_lines, sprintf("| %s | %s | %s | %s | %s | %.1f%% |",
    r$meta_type, fmt(r$Full), fmt(r$Partial), none_cell, fmt(r$total), r$pct_full))
}
md_lines <- c(md_lines, sprintf("| **Total** | **%s** | **%s** | **%s** | **%s** | **%.1f%%** |",
  fmt(tot_full), fmt(tot_partial), fmt(tot_none), fmt(grand), 100*tot_full/grand))

cat("\n=== Markdown ===\n"); cat(md_lines, sep = "\n"); cat("\n")

## --- write outputs -------------------------------------------------------
out_csv <- file.path(OUT_DIR, "table4_savas_corrected.csv")
out_md  <- file.path(OUT_DIR, "table4_savas_corrected.md")
write.csv(tbl %>% mutate(meta_type = as.character(meta_type)), out_csv, row.names = FALSE)
writeLines(md_lines, out_md)
cat("\nWrote:\n ", out_csv, "\n ", out_md, "\n\nDONE\n")
