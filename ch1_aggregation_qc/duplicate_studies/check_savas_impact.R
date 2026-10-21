suppressPackageStartupMessages({library(dplyr)})

obj <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")

cat("class:", class(obj), "\n")
cat("total rows (current, before removal):", nrow(obj), "\n")
cat("total distinct project_name (current, before removal):", n_distinct(obj$project_name), "\n\n")

target <- "SAVAS_P_ET_AL_2018_P35_filtered"

meta <- as.data.frame(obj)

cat("Sanity check - does this exact string exist in project_name?\n")
print(target %in% unique(meta$project_name))

type_col_candidates <- c("Consensus_Cell_Type", "meta_type", "consensus_cell_type", "cell_type")
type_col <- type_col_candidates[type_col_candidates %in% colnames(meta)][1]
cat("meta-type column used:", ifelse(is.na(type_col), "NONE FOUND", type_col), "\n\n")

sub <- meta[meta$project_name == target, , drop = FALSE]
n <- nrow(sub)
cat("=== Target:", target, "===\n")
cat("n cells:", n, "\n")

cat("\ncancer_type breakdown:\n")
print(sort(table(sub$cancer_type), decreasing = TRUE))

if (!is.na(type_col)) {
  cat("\nmeta-type breakdown:\n")
  print(sort(table(sub[[type_col]]), decreasing = TRUE))
}

# Also print a few other candidate columns for context (patient, sample) if present
for (cand in c("patient", "sample", "orig.ident")) {
  if (cand %in% colnames(sub)) {
    cat("\n", cand, "breakdown:\n")
    print(sort(table(sub[[cand]]), decreasing = TRUE))
  }
}

summary_df <- data.frame(project_name = target, n_cells = n,
                          cancer_type = paste(names(sort(table(sub$cancer_type), decreasing=TRUE)), collapse=";"))

cat("\n=== Corrected final-atlas totals if removed ===\n")
cat("Corrected total cells (current - target):", nrow(obj) - n, "\n")
cat("Corrected total distinct project_name (current - 1):", n_distinct(obj$project_name) - 1, "\n")

cat("\n=== Corrected per-cancer-type counts (Breast) ===\n")
before_after <- meta %>%
  filter(cancer_type == "Breast") %>%
  mutate(is_target = project_name == target) %>%
  group_by(cancer_type) %>%
  summarise(before = n(), removed = sum(is_target), after = before - removed, .groups = "drop")
print(as.data.frame(before_after), row.names = FALSE)

write.csv(summary_df, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/savas2018-duplicate-removal-8a3b6b/work/savas_impact_summary.csv", row.names = FALSE)
write.csv(before_after, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/savas2018-duplicate-removal-8a3b6b/work/savas_impact_cancer_type_before_after.csv", row.names = FALSE)

if (!is.na(type_col)) {
  mt_df <- as.data.frame(sort(table(sub[[type_col]]), decreasing = TRUE))
  colnames(mt_df) <- c("meta_type", "n_cells")
  write.csv(mt_df, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/savas2018-duplicate-removal-8a3b6b/work/savas_meta_type_breakdown.csv", row.names = FALSE)
}

cat("\nDone. Wrote savas_impact_summary.csv, savas_impact_cancer_type_before_after.csv, savas_meta_type_breakdown.csv\n")
