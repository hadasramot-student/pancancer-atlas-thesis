suppressPackageStartupMessages({
  library(dplyr)
})

path <- "/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS"
cat("Loading:", path, "\n")
md <- readRDS(path)
cat("Raw rows:", nrow(md), "\n")

dup_studies <- c(
  "WU_SZ_ET_AL_2021_P15_filtered",
  "GAO_R_ET_AL_2021_P17_filtered",
  "GRIFFITHS_JI_ET_AL_2021_P18_filtered",
  "GONZALEZ_H_ET_AL_2022_P31_filtered"
)

cat("\nCells per duplicate study (to be excluded):\n")
print(md %>% filter(project_name %in% dup_studies) %>% count(project_name))

md2 <- md %>% filter(!(project_name %in% dup_studies))
cat("\nRows after excluding duplicate studies:", nrow(md2), "(expect 8113115)\n")

# --- cell-level agreement status ---
md2 <- md2 %>%
  mutate(
    cell_status = case_when(
      maj3_cell != "equivocal" ~ "Full",
      maj3_cell == "equivocal" & maj2_cell != "equivocal" ~ "Partial",
      maj2_cell == "equivocal" ~ "None",
      TRUE ~ NA_character_
    )
  )

cat("\nSanity: any NA cell_status?", sum(is.na(md2$cell_status)), "\n")

# --- main table: Full / Partial / None counts per FINAL Consensus_Cell_Type ---
main_tbl <- md2 %>%
  count(Consensus_Cell_Type, cell_status) %>%
  tidyr::pivot_wider(names_from = cell_status, values_from = n, values_fill = 0)

cat("\n=== Main table: Full/Partial/None by Consensus_Cell_Type ===\n")
print(as.data.frame(main_tbl))

main_tbl <- main_tbl %>%
  mutate(total = Full + Partial + ifelse("None" %in% names(main_tbl), None, 0))

cat("\n=== With totals ===\n")
print(as.data.frame(main_tbl))

cat("\nGrand total check:", sum(main_tbl$total), "(expect 8113115)\n")

# --- Other/Unknown override sub-counts cited in prose ---
# "full/partial agreement on a REAL (non-'other', non-equivocal) type that was
#  nonetheless overridden to Other/Unknown by the cluster-level check"
other_full_overridden <- md2 %>%
  filter(Consensus_Cell_Type == "other",
         maj3_cell != "equivocal", maj3_cell != "other") %>%
  nrow()

other_partial_overridden <- md2 %>%
  filter(Consensus_Cell_Type == "other",
         maj3_cell == "equivocal",
         maj2_cell != "equivocal", maj2_cell != "other") %>%
  nrow()

cat("\n=== Other/Unknown override sub-counts ===\n")
cat("Full (3/3) cell-level agreement, overridden by cluster check:", other_full_overridden, "\n")
cat("Partial (2/3) cell-level agreement, overridden by cluster check:", other_partial_overridden, "\n")

# For context: how many Other/Unknown Full/Partial cells directly agreed on "other" itself
other_full_direct <- md2 %>% filter(Consensus_Cell_Type == "other", maj3_cell == "other") %>% nrow()
other_partial_direct <- md2 %>% filter(Consensus_Cell_Type == "other", maj3_cell == "equivocal", maj2_cell == "other") %>% nrow()
cat("Full cells that directly agreed label = 'other' (not an override):", other_full_direct, "\n")
cat("Partial cells that directly agreed label = 'other' (not an override):", other_partial_direct, "\n")
cat("Check Full sums:", other_full_overridden + other_full_direct,
    "vs main_tbl Full for other:",
    main_tbl$Full[main_tbl$Consensus_Cell_Type == "other"], "\n")
cat("Check Partial sums:", other_partial_overridden + other_partial_direct,
    "vs main_tbl Partial for other:",
    main_tbl$Partial[main_tbl$Consensus_Cell_Type == "other"], "\n")

# --- % full agreement per meta-type ---
main_tbl <- main_tbl %>%
  mutate(pct_full = 100 * Full / total) %>%
  arrange(desc(pct_full))

cat("\n=== Ranked by %% Full agreement ===\n")
print(as.data.frame(main_tbl %>% select(Consensus_Cell_Type, Full, Partial, any_of("None"), total, pct_full)))

saveRDS(main_tbl, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/table-4-agreement-counts-69fa2f/table4_main_tbl.RDS")

cat("\nDONE\n")
