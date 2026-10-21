suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
})

M <- "/gpfs0/bgu-ofircohen/users/hadasra/Raw_Data/scBaseCount/meta_data"
meta <- read_csv(file.path(M, "ScBaseCount_metadata.csv"), show_col_types = FALSE)
merged_meta <- read_csv(file.path(M, "merged_meta_data.csv"), show_col_types = FALSE)

# same function as merge_scbasecount_per_proj.R
infer_cancer_type <- function(disease, tissue) {
  txt <- tolower(paste(disease, tissue))
  if (grepl("breast", txt))                      return("Breast")
  if (grepl("colorectal|colon|rectal", txt))      return("Colorectal")
  if (grepl("ovary|ovarian", txt))                return("Ovarian")
  if (grepl("pancreas|pancreatic", txt))          return("Pancreas")
  if (grepl("sarcoma", txt))                      return("Sarcoma")
  if (grepl("liver|hepatocellular|bile duct|cholangiocarcinoma", txt)) return("Liver_Biliary")
  if (grepl("kidney|renal", txt))                 return("Kidney")
  if (grepl("skin|melanoma|cutaneous", txt))      return("Skin")
  if (grepl("lung|bronch", txt))                  return("Lung")
  if (grepl("head and neck|hypopharyn|oropharyn|larynx|oral cavity", txt)) return("Head_and_Neck")
  if (grepl("prostate", txt))                     return("Prostate")
  if (grepl("neuroendocrine", txt))                return("Neuroendocrine")
  if (grepl("brain|glioblastoma|glioma|astrocytoma", txt)) return("Brain")
  return("unknown")
}

meta$cancer_type <- mapply(infer_cancer_type, meta$disease, meta$tissue)

# join SRX -> Study.Accession (the SRP/ERP-level project id, i.e. our raw filenames)
j <- merged_meta %>%
  select(srx_accession = `Experiment.Accession`, study_accession = `Study.Accession`) %>%
  distinct() %>%
  inner_join(meta %>% select(srx_accession, cancer_type), by = "srx_accession")

cat("SRX rows joined:", nrow(j), "of", nrow(merged_meta), "merged_meta rows,", nrow(meta), "metadata rows\n")

per_study <- j %>%
  group_by(study_accession) %>%
  summarise(
    n_srx = n(),
    n_distinct_types = n_distinct(cancer_type),
    types = paste(sort(unique(cancer_type)), collapse = ";"),
    .groups = "drop"
  )

cat("\nTotal distinct study_accession found via metadata join:", nrow(per_study), "\n")
cat("Studies with >1 distinct cancer_type (incl. 'unknown' as its own type):", sum(per_study$n_distinct_types > 1), "\n")
cat("Studies with exactly 1 distinct cancer_type:", sum(per_study$n_distinct_types == 1), "\n")
cat("  ...of which 'unknown' only:", sum(per_study$n_distinct_types == 1 & per_study$types == "unknown"), "\n")
cat("  ...of which a real single cancer type:", sum(per_study$n_distinct_types == 1 & per_study$types != "unknown"), "\n")

write_tsv(per_study, "~/claude_scbasecount_reconcile.tsv")
cat("\nWrote ~/claude_scbasecount_reconcile.tsv\n")
