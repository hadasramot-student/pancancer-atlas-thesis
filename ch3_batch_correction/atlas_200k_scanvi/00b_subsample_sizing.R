# 01_subsample_sizing.R
# Metadata-only sizing pass for the atlas-scale scANVI Figure 1/2 subsample.
# NO modelling. Single group-by over the canonical atlas metadata table.
#
# Questions (from Prompt 15):
#   1. cells per cancer_type, smallest -> largest
#   2. cells per lineage (Consensus_Cell_Type, the Ch.2 consensus meta-type)
#   3. full cancer_type x lineage cross-tab
#   4. expected cell count per cross-tab cell under a uniform (proportional)
#      stratified draw at 200k; flag cells < ~200-300
#   5. same at 400k and 800k
#   6. explicit flag for Head_and_Neck and Neuroendocrine

options(width = 220, max.print = 200000, stringsAsFactors = FALSE)

TO  <- "/ofir_grp/users/hadasra/temp_objects"
f   <- file.path(TO, "all_cells.gapfix_20260727_164657.RDS")
OUT <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/atlas-subsample-sizing-ab0533/work"

cat("Reading:", f, "\n")
md <- as.data.frame(readRDS(f))
cat("raw dim:", paste(dim(md), collapse = " x "), "\n\n")

## --- apply the 5-study cross-repo dedup (Breastscape side), matching the
##     canonical corrected atlas count of 8,107,665 -----------------------------
dups <- c("WU_SZ_ET_AL_2021_P15_filtered",
          "GRIFFITHS_JI_ET_AL_2021_P18_filtered",
          "GONZALEZ_H_ET_AL_2022_P31_filtered",
          "SAVAS_P_ET_AL_2018_P35_filtered",
          "GAO_R_ET_AL_2021_P17_filtered")
n_before <- nrow(md)
md <- md[!(md$project_name %in% dups), , drop = FALSE]
cat(sprintf("dropped %d cross-repo duplicate rows (5 studies); atlas N = %d\n\n",
            n_before - nrow(md), nrow(md)))
N_ATLAS <- nrow(md)

lin_col <- "Consensus_Cell_Type"
can_col <- "cancer_type"
stopifnot(all(c(lin_col, can_col) %in% colnames(md)))

## =====================================================================
## 1. cells per cancer type (smallest -> largest)
## =====================================================================
can <- as.data.frame(table(md[[can_col]]), responseName = "n")
colnames(can)[1] <- "cancer_type"
can <- can[order(can$n), ]
can$pct_atlas <- round(100 * can$n / N_ATLAS, 3)
cat("=== 1. CELLS PER CANCER TYPE (asc) ===\n"); print(can, row.names = FALSE)
write.csv(can, file.path(OUT, "counts_by_cancer_type.csv"), row.names = FALSE)

## =====================================================================
## 2. cells per lineage / consensus meta-type
## =====================================================================
lin <- as.data.frame(table(md[[lin_col]]), responseName = "n")
colnames(lin)[1] <- "lineage"
lin <- lin[order(lin$n), ]
lin$pct_atlas <- round(100 * lin$n / N_ATLAS, 3)
cat("\n=== 2. CELLS PER LINEAGE (Consensus_Cell_Type, asc) ===\n"); print(lin, row.names = FALSE)
write.csv(lin, file.path(OUT, "counts_by_lineage.csv"), row.names = FALSE)

## =====================================================================
## 3. full cancer_type x lineage cross-tab (long form)
## =====================================================================
ct <- as.data.frame(table(cancer_type = md[[can_col]], lineage = md[[lin_col]]),
                     responseName = "n")
ct <- ct[order(ct$n), ]
ct$pct_atlas <- round(100 * ct$n / N_ATLAS, 4)

## expected counts under a uniform / proportional stratified draw:
##   f = N_sub / N_atlas ;  E[cell] = n_ij * f
for (S in c(2e5, 4e5, 8e5)) {
  frac <- S / N_ATLAS
  ct[[sprintf("exp_%dk", S / 1e3)]] <- round(ct$n * frac, 1)
}
ct$flag_200k <- ifelse(ct$`exp_200k` < 200, "BELOW_200",
                ifelse(ct$`exp_200k` < 300, "below_300", ""))

cat(sprintf("\n=== 3+4+5. CROSS-TAB (cancer_type x lineage), asc; sampling fracs: 200k=%.3f%%  400k=%.3f%%  800k=%.3f%% ===\n",
            100 * 2e5 / N_ATLAS, 100 * 4e5 / N_ATLAS, 100 * 8e5 / N_ATLAS))
print(ct, row.names = FALSE)
write.csv(ct, file.path(OUT, "crosstab_cancer_x_lineage_with_expected.csv"), row.names = FALSE)

## wide count matrix too, for eyeballing
wide <- xtabs(n ~ cancer_type + lineage, data = ct)
cat("\n--- wide count matrix (cancer_type x lineage) ---\n"); print(addmargins(wide))
write.csv(as.data.frame.matrix(addmargins(wide)),
          file.path(OUT, "crosstab_cancer_x_lineage_counts_wide.csv"))

## =====================================================================
## 4. the smallest 15 cross-tab cells, explicitly
## =====================================================================
cat("\n=== 4. SMALLEST 15 (cancer_type x lineage) CELLS + expected draw ===\n")
print(head(ct[, c("cancer_type","lineage","n","pct_atlas","exp_200k","exp_400k","exp_800k","flag_200k")], 15),
      row.names = FALSE)

cat("\n--- how many cross-tab cells fall below thresholds, by subsample size ---\n")
thr_tab <- data.frame(
  size      = c("200k","400k","800k"),
  below_100 = c(sum(ct$`exp_200k` < 100), sum(ct$`exp_400k` < 100), sum(ct$`exp_800k` < 100)),
  below_200 = c(sum(ct$`exp_200k` < 200), sum(ct$`exp_400k` < 200), sum(ct$`exp_800k` < 200)),
  below_300 = c(sum(ct$`exp_200k` < 300), sum(ct$`exp_400k` < 300), sum(ct$`exp_800k` < 300)),
  nonzero_cells = sum(ct$n > 0)
)
print(thr_tab, row.names = FALSE)
write.csv(thr_tab, file.path(OUT, "threshold_summary.csv"), row.names = FALSE)

## =====================================================================
## 6. Head_and_Neck and Neuroendocrine, explicit
## =====================================================================
cat("\n=== 6. THIN COHORTS: Head_and_Neck & Neuroendocrine ===\n")
for (cc in c("Head_and_Neck","Neuroendocrine")) {
  sub <- ct[ct$cancer_type == cc, c("lineage","n","exp_200k","exp_400k","exp_800k")]
  sub <- sub[order(-sub$n), ]
  cat(sprintf("\n-- %s : total = %d cells (%.3f%% of atlas); contributing projects: %d\n",
              cc, sum(sub$n), 100 * sum(sub$n) / N_ATLAS,
              length(unique(md$project_name[md[[can_col]] == cc]))))
  print(sub, row.names = FALSE)
}

cat("\n=== done ===\n")
