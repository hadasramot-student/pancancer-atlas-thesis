suppressPackageStartupMessages({library(dplyr)})
obj <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")
meta <- as.data.frame(obj)

classify_repo <- function(pn) {
  ifelse(grepl("^(SRP|ERP|DRP)", pn), "scbasecount",
  ifelse(grepl("_ET_AL_|_et_al_", pn), "breastscape", "3ca"))
}
meta$repo <- classify_repo(meta$project_name)

cat("=== population rate of 'sample' and 'SRX_accession' by repo ===\n")
meta %>% group_by(repo) %>%
  summarise(n_rows = n(),
            n_distinct_project = n_distinct(project_name),
            sample_nonNA = sum(!is.na(sample) & sample != "" & sample != "NA"),
            srx_nonNA = sum(!is.na(SRX_accession) & SRX_accession != "" & SRX_accession != "NA"),
            .groups = "drop") %>% as.data.frame() %>% print()

targets <- c("WU_SZ_ET_AL_2021_P15_filtered", "GAO_R_ET_AL_2021_P17_filtered",
             "GRIFFITHS_JI_ET_AL_2021_P18_filtered", "GONZALEZ_H_ET_AL_2022_P31_filtered")

cat("\n=== candidate proxy counts (whole atlas, current) ===\n")
d1 <- meta %>% filter(!is.na(SRX_accession), SRX_accession != "", SRX_accession != "NA") %>%
  summarise(n = n_distinct(SRX_accession)) %>% pull(n)
cat("distinct SRX_accession (global, non-NA only):", d1, "\n")

d2 <- meta %>% filter(!is.na(sample), sample != "", sample != "NA") %>%
  summarise(n = n_distinct(paste(project_name, sample, sep="___"))) %>% pull(n)
cat("distinct (project_name, sample) [non-NA sample rows only]:", d2, "\n")

# combined: use SRX_accession where sample is NA (scbasecount), else project_name+sample
meta2 <- meta %>% mutate(
  unit = case_when(
    !is.na(sample) & sample != "" & sample != "NA" ~ paste(project_name, sample, sep="___S___"),
    !is.na(SRX_accession) & SRX_accession != "" & SRX_accession != "NA" ~ paste(project_name, SRX_accession, sep="___X___"),
    TRUE ~ paste(project_name, "UNKNOWN_UNIT", sep="___U___")
  )
)
d3 <- n_distinct(meta2$unit)
cat("combined unit (sample if present else SRX_accession else fallback-per-project):", d3, "\n")

# same three, but also broken down by repo to see which repos actually contribute
cat("\n=== combined-unit distinct count broken down by repo ===\n")
meta2 %>% group_by(repo) %>% summarise(n_units = n_distinct(unit), .groups="drop") %>% as.data.frame() %>% print()

cat("\nTarget: 11769\n")
cat("Distances -> d1:", abs(d1-11769), " d2:", abs(d2-11769), " d3:", abs(d3-11769), "\n")

cat("\n=== same three, AFTER excluding the 4 confirmed Breastscape duplicates ===\n")
meta2_after <- meta2 %>% filter(!project_name %in% targets)
d1a <- meta2_after %>% filter(!is.na(SRX_accession), SRX_accession != "", SRX_accession != "NA") %>% summarise(n=n_distinct(SRX_accession)) %>% pull(n)
d2a <- meta2_after %>% filter(!is.na(sample), sample != "", sample != "NA") %>% summarise(n=n_distinct(paste(project_name, sample, sep="___"))) %>% pull(n)
d3a <- n_distinct(meta2_after$unit)
cat("d1 after:", d1a, " delta:", d1-d1a, "\n")
cat("d2 after:", d2a, " delta:", d2-d2a, "\n")
cat("d3 after:", d3a, " delta:", d3-d3a, "\n")

cat("\nDone.\n")
