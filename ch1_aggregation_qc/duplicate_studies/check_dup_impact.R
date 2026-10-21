suppressPackageStartupMessages({library(dplyr)})

obj <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")

cat("class:", class(obj), "\n")
cat("total rows (current):", nrow(obj), "\n")
cat("total distinct project_name (current):", n_distinct(obj$project_name), "\n\n")

targets <- c("WU_SZ_ET_AL_2021_P15_filtered", "GAO_R_ET_AL_2021_P17_filtered",
             "GRIFFITHS_JI_ET_AL_2021_P18_filtered", "GONZALEZ_H_ET_AL_2022_P31_filtered")

meta <- as.data.frame(obj)

cat("Sanity check - do these exact strings exist in project_name?\n")
print(targets %in% unique(meta$project_name))

# figure out a meta-type / consensus cell type column name if present
type_col_candidates <- c("Consensus_Cell_Type", "meta_type", "consensus_cell_type", "cell_type")
type_col <- type_col_candidates[type_col_candidates %in% colnames(meta)][1]
cat("meta-type column used:", ifelse(is.na(type_col), "NONE FOUND", type_col), "\n\n")

cat("=== Per-target: exact final-atlas cell count + cancer_type + (meta-type if available) ===\n")
summary_rows <- list()
for (t in targets) {
  sub <- meta[meta$project_name == t, , drop = FALSE]
  n <- nrow(sub)
  ct <- if (n > 0) sort(table(sub$cancer_type), decreasing = TRUE) else NA
  cat("\n--", t, "--\n")
  cat("  n cells:", n, "\n")
  cat("  cancer_type breakdown:\n"); print(ct)
  if (!is.na(type_col)) {
    mt <- sort(table(sub[[type_col]]), decreasing = TRUE)
    cat("  meta-type breakdown:\n"); print(mt)
  }
  summary_rows[[t]] <- data.frame(project_name = t, n_cells = n,
                                   cancer_type = paste(names(ct), collapse=";"))
}

summary_df <- do.call(rbind, summary_rows)
cat("\n=== Summary table ===\n")
print(summary_df, row.names = FALSE)

total_removed <- sum(summary_df$n_cells)
cat("\nTotal cells across the 4 targets:", total_removed, "\n")
cat("Corrected total cells (current - removed):", nrow(obj) - total_removed, "\n")
cat("Corrected total distinct project_name (current - 4):", n_distinct(obj$project_name) - 4, "\n")

cat("\n=== Corrected per-cancer-type counts (Breast and any other touched type) ===\n")
touched_types <- unique(meta$cancer_type[meta$project_name %in% targets])
cat("Cancer types touched by these 4 studies:", paste(touched_types, collapse=", "), "\n\n")

before_after <- meta %>%
  filter(cancer_type %in% touched_types) %>%
  mutate(is_target = project_name %in% targets) %>%
  group_by(cancer_type) %>%
  summarise(before = n(), removed = sum(is_target), after = before - removed, .groups = "drop")
print(as.data.frame(before_after), row.names = FALSE)

write.csv(summary_df, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/cross-repo-duplicate-check-e81250/work/dup_impact_summary.csv", row.names = FALSE)
write.csv(before_after, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/cross-repo-duplicate-check-e81250/work/dup_impact_cancer_type_before_after.csv", row.names = FALSE)

cat("\nDone. Wrote dup_impact_summary.csv and dup_impact_cancer_type_before_after.csv\n")
