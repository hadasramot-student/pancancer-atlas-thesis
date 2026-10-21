#!/usr/bin/env Rscript
# query_project_status.R
# Pulls distinct project_name x n_cells from the (gap-fixed) atlas table.
suppressPackageStartupMessages({ library(dplyr); library(readr) })

TO <- "/ofir_grp/users/hadasra/temp_objects"
df <- readRDS(file.path(TO, "all_cells.gapfix_20260727_164657.RDS"))
cat("Total rows:", nrow(df), "\n")

per_proj <- df %>%
  count(project_name, name = "n_cells_in_atlas") %>%
  arrange(desc(n_cells_in_atlas))

write_tsv(per_proj, file.path(TO, "gapfix_project_name_counts.tsv"))
cat("Distinct project_name values:", nrow(per_proj), "\n")
cat("Saved to gapfix_project_name_counts.tsv\n")
