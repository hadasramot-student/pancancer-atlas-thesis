#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(dplyr) })

base <- "/ofir_grp/users/hadasra/temp_objects"
all_cells <- readRDS(file.path(base, "all_cells.RDS"))

cat("=== all_cells.RDS columns ===\n")
print(colnames(all_cells))

cat("\n=== Unique values of ORIGIN (if exists) ===\n")
if ("ORIGIN" %in% colnames(all_cells)) print(table(all_cells$ORIGIN, useNA="ifany"))

cat("\n=== Unique values of C.TYPE (if exists) ===\n")
if ("C.TYPE" %in% colnames(all_cells)) print(table(all_cells$C.TYPE, useNA="ifany"))

cat("\n=== Unique values of paper (sample) ===\n")
if ("paper" %in% colnames(all_cells)) print(head(unique(all_cells$paper), 20))

lineages <- c("T-cells","Myeloid","Epithelial","B-cell","Fibroblasts",
              "Endothelial","NK-cells","Hematopoietic_Stem_and_Progenitor_cells")

breast <- all_cells %>% filter(cancer_type == "Breast", Consensus_Cell_Type %in% lineages)

cat("\n=== Breast (lineage cells only) total: ", nrow(breast), " ===\n")

proj_summary <- breast %>%
  group_by(project_name) %>%
  summarise(n_cells = n(),
            ORIGIN = paste(unique(ORIGIN), collapse="|"),
            C.TYPE = paste(unique(C.TYPE), collapse="|"),
            paper  = paste(unique(paper), collapse="|"),
            n_patients = n_distinct(patient)) %>%
  arrange(desc(n_cells))

cat("\n=== Per-project breakdown (project_name, n_cells, ORIGIN, C.TYPE, paper, n_patients) ===\n")
print(as.data.frame(proj_summary), row.names = FALSE)

cat("\n=== Total distinct projects: ", nrow(proj_summary), " ===\n")
cat("=== Total distinct patients: ", n_distinct(breast$patient), " ===\n")

# Classify by naming convention as a cross-check
proj_summary$name_pattern <- case_when(
  grepl("^ERP|^SRP|^ERX|^SRX", proj_summary$project_name) ~ "ScBaseCount-style (ERP/SRP accession)",
  grepl("^[A-Z][a-z]+[0-9]{4}_", proj_summary$project_name) ~ "3CA-style (Author+Year_CancerType)",
  TRUE ~ "OTHER/UNRECOGNIZED"
)
cat("\n=== Classification by project_name pattern ===\n")
print(table(proj_summary$name_pattern))
cat("\n=== Projects classified as OTHER/UNRECOGNIZED (potential Breastscape or unknown) ===\n")
print(as.data.frame(proj_summary[proj_summary$name_pattern == "OTHER/UNRECOGNIZED", c("project_name","n_cells","ORIGIN","C.TYPE","paper")]), row.names = FALSE)

write.csv(proj_summary, file.path(base, "breast_cohort_project_origins.csv"), row.names = FALSE)
cat("\nSaved: ", file.path(base, "breast_cohort_project_origins.csv"), "\n")
cat("DONE\n")
