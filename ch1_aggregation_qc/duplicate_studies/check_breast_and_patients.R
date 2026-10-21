suppressPackageStartupMessages({library(dplyr)})

obj <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")
meta <- as.data.frame(obj)
cat("nrow:", nrow(meta), " n_distinct(project_name):", n_distinct(meta$project_name), "\n\n")

targets <- c("WU_SZ_ET_AL_2021_P15_filtered", "GAO_R_ET_AL_2021_P17_filtered",
             "GRIFFITHS_JI_ET_AL_2021_P18_filtered", "GONZALEZ_H_ET_AL_2022_P31_filtered")
stopifnot(all(targets %in% unique(meta$project_name)))

cat("=================== PATIENT COUNT INVESTIGATION ===================\n")
cat("Target: figure out which method reproduces the thesis's cited 11,769 patients\n\n")

# raw patient field stats
cat("class(meta$patient):", class(meta$patient), "\n")
n_na_patient <- sum(is.na(meta$patient) | meta$patient == "" | meta$patient == "NA")
cat("rows with NA/empty patient:", n_na_patient, "out of", nrow(meta), "\n")

# which project_names have NA/empty patient (sanity: expect ScBaseCount SRP/ERP rows)
na_projects <- meta %>% filter(is.na(patient) | patient == "" | patient == "NA") %>%
  distinct(project_name) %>% pull(project_name)
cat("n distinct project_name with NA/empty patient:", length(na_projects), "\n")
cat("sample of those project_names:\n"); print(head(sort(na_projects), 15))
cat("Do all of these look like SRA accessions (SRP/ERP/DRP prefix)? ",
    all(grepl("^(SRP|ERP|DRP)", na_projects)), "\n\n")

# Method A: global raw distinct patient string (ignoring NA)
methodA <- meta %>% filter(!is.na(patient), patient != "", patient != "NA") %>%
  summarise(n = n_distinct(patient)) %>% pull(n)
cat("METHOD A - global n_distinct(patient) [raw string, NAs dropped]:", methodA, "\n")

# Method B: compound key project_name + patient (per-study patient IDs, avoids cross-study collisions)
methodB <- meta %>% filter(!is.na(patient), patient != "", patient != "NA") %>%
  summarise(n = n_distinct(paste(project_name, patient, sep = "___"))) %>% pull(n)
cat("METHOD B - n_distinct(project_name + patient) [compound key]:", methodB, "\n")

# Method C: same as A/B but including NA patients each as their own singleton via project_name+sample fallback
# (only run if neither A nor B is close to 11769)
methodC <- meta %>%
  mutate(patient_key = ifelse(!is.na(patient) & patient != "" & patient != "NA",
                               paste(project_name, patient, sep = "___"),
                               paste(project_name, sample, sep = "___SAMPLE___"))) %>%
  summarise(n = n_distinct(patient_key)) %>% pull(n)
cat("METHOD C - compound key with sample-level fallback for NA-patient rows:", methodC, "\n\n")

cat("Target value from thesis text: 11769\n")
cat("Closest method to target -> A:", abs(methodA-11769), " B:", abs(methodB-11769), " C:", abs(methodC-11769), "\n\n")

cat("=================== RECOMPUTE AFTER REMOVING THE 4 BREASTSCAPE ENTRIES ===================\n")
meta_after <- meta %>% filter(!project_name %in% targets)
cat("nrow after removal:", nrow(meta_after), " (expect", nrow(meta) - sum(meta$project_name %in% targets), ")\n\n")

methodA_after <- meta_after %>% filter(!is.na(patient), patient != "", patient != "NA") %>%
  summarise(n = n_distinct(patient)) %>% pull(n)
methodB_after <- meta_after %>% filter(!is.na(patient), patient != "", patient != "NA") %>%
  summarise(n = n_distinct(paste(project_name, patient, sep = "___"))) %>% pull(n)
methodC_after <- meta_after %>%
  mutate(patient_key = ifelse(!is.na(patient) & patient != "" & patient != "NA",
                               paste(project_name, patient, sep = "___"),
                               paste(project_name, sample, sep = "___SAMPLE___"))) %>%
  summarise(n = n_distinct(patient_key)) %>% pull(n)

cat("METHOD A after removal:", methodA_after, " delta:", methodA - methodA_after, "\n")
cat("METHOD B after removal:", methodB_after, " delta:", methodB - methodB_after, "\n")
cat("METHOD C after removal:", methodC_after, " delta:", methodC - methodC_after, "\n\n")

# Explicit check requested in task point 4: do the 4 removed Breastscape studies'
# patient ID strings collide with any 3CA patient ID strings (raw string, method A logic)?
cat("=== Explicit patient-ID collision check for the 4 removed studies (raw string overlap with rest of atlas) ===\n")
for (t in targets) {
  removed_patients <- unique(meta$patient[meta$project_name == t])
  removed_patients <- removed_patients[!is.na(removed_patients) & removed_patients != "" & removed_patients != "NA"]
  rest_patients <- unique(meta$patient[meta$project_name != t])
  collide <- intersect(removed_patients, rest_patients)
  cat("\n--", t, "--\n")
  cat("  n distinct patient IDs in this study:", length(removed_patients), "\n")
  cat("  of which also appear (as raw strings) elsewhere in the atlas:", length(collide), "\n")
  if (length(collide) > 0 && length(collide) <= 20) {
    cat("  colliding IDs:", paste(collide, collapse=", "), "\n")
    # show which OTHER project_names those colliding IDs belong to
    other_hits <- meta %>% filter(patient %in% collide, project_name != t) %>%
      distinct(patient, project_name)
    print(as.data.frame(other_hits))
  }
}

cat("\n=================== BREAST RAW-BASELINE CROSS-CHECK (raw h5ad vs final-atlas cell counts) ===================\n")
raw_h5ad_counts <- c(WU_SZ_ET_AL_2021_P15_filtered = 46305,
                     GAO_R_ET_AL_2021_P17_filtered = 4617,
                     GRIFFITHS_JI_ET_AL_2021_P18_filtered = 86334,
                     GONZALEZ_H_ET_AL_2022_P31_filtered = 19999)
final_atlas_counts <- sapply(targets, function(t) sum(meta$project_name == t))
names(final_atlas_counts) <- targets
cmp <- data.frame(project_name = targets,
                   raw_h5ad = raw_h5ad_counts[targets],
                   final_atlas = final_atlas_counts[targets],
                   identical = raw_h5ad_counts[targets] == final_atlas_counts[targets])
print(cmp, row.names = FALSE)
cat("\nAll identical (raw == final-atlas, no QC loss for these 4)?", all(cmp$identical), "\n")

write.csv(cmp, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/cross-repo-duplicate-check-e81250/work/raw_vs_final_atlas_4studies.csv", row.names = FALSE)

cat("\nDone.\n")
