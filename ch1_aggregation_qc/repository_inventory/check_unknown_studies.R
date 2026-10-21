suppressPackageStartupMessages({ library(readr); library(dplyr) })

M <- "/gpfs0/bgu-ofircohen/users/hadasra/Raw_Data/scBaseCount/meta_data"
meta <- read_csv(file.path(M, "ScBaseCount_metadata.csv"), show_col_types = FALSE)
merged_meta <- read_csv(file.path(M, "merged_meta_data.csv"), show_col_types = FALSE)

infer_cancer_type <- function(disease, tissue) {
  txt <- tolower(paste(disease, tissue))
  if (grepl("breast", txt)) return("Breast")
  if (grepl("colorectal|colon|rectal", txt)) return("Colorectal")
  if (grepl("ovary|ovarian", txt)) return("Ovarian")
  if (grepl("pancreas|pancreatic", txt)) return("Pancreas")
  if (grepl("sarcoma", txt)) return("Sarcoma")
  if (grepl("liver|hepatocellular|bile duct|cholangiocarcinoma", txt)) return("Liver_Biliary")
  if (grepl("kidney|renal", txt)) return("Kidney")
  if (grepl("skin|melanoma|cutaneous", txt)) return("Skin")
  if (grepl("lung|bronch", txt)) return("Lung")
  if (grepl("head and neck|hypopharyn|oropharyn|larynx|oral cavity", txt)) return("Head_and_Neck")
  if (grepl("prostate", txt)) return("Prostate")
  if (grepl("neuroendocrine", txt)) return("Neuroendocrine")
  if (grepl("brain|glioblastoma|glioma|astrocytoma", txt)) return("Brain")
  return("unknown")
}
meta$cancer_type <- mapply(infer_cancer_type, meta$disease, meta$tissue)

j <- merged_meta %>%
  select(srx_accession = `Experiment.Accession`, study_accession = `Study.Accession`) %>%
  distinct() %>%
  inner_join(meta %>% select(srx_accession, disease, tissue, cancer_type), by = "srx_accession")

per_study <- j %>% group_by(study_accession) %>%
  summarise(n_distinct_types = n_distinct(cancer_type),
            types = paste(sort(unique(cancer_type)), collapse=";"), .groups="drop")

unknown_only <- per_study %>% filter(n_distinct_types==1, types=="unknown") %>% pull(study_accession)
cat("N unknown-only studies:", length(unknown_only), "\n\n")

detail <- j %>% filter(study_accession %in% unknown_only) %>%
  distinct(study_accession, disease, tissue) %>%
  arrange(study_accession)

write_tsv(detail, "~/claude_unknown_studies_detail.tsv")
print(detail, n=100)
