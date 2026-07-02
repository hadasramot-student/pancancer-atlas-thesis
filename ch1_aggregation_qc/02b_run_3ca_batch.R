# 02b_run_3ca_batch.R
#
# Ch.1 (Aggregation & QC) — batch driver for the 3CA source.
#
# Reads the 3CA project file index, restricts to 10x technology (non-
# hematologic cancers handled elsewhere), and invokes 00_main_dispatch.R
# (--ptype 3CA) for a selected set of row numbers in that index.
#
# Inputs:
#   SC3CA_PROJECT_INDEX - TSV index of all 3CA projects (path, technology,
#                          Cancer_type, project_name columns)
# Outputs:
#   Invokes 00_main_dispatch.R --ptype 3CA --row_number_3CA <i> once per
#   selected project row (each call writes its own ingested/QC'd RDS via
#   02_ingest_3ca.R).
#
# Usage:
#   SC3CA_PROJECT_INDEX=/path/to/SC.3CA.Tirosh.project_file_index.tsv \
#     Rscript 02b_run_3ca_batch.R

library(readr)
library(dplyr)

projects_df <- read_tsv(
  Sys.getenv("SC3CA_PROJECT_INDEX", unset = "config/SC.3CA.Tirosh.project_file_index.tsv"),
  comment = "#"
)

sub_df <- projects_df %>%
  filter(technology == "10x", Cancer_type != "Hematologic") %>%
  group_by(project_name) %>%
  mutate(rn = row_number(), n = n()) %>%
  mutate(project_name_i = glue::glue("{project_name}_{rn}")) %>%
  mutate(project_name_i = if_else(n == 1, project_name, project_name_i))

# Row numbers (into sub_df) to (re-)process this run — edit per batch.
number_row_to_proccess <- c(7, 10, 19, 25, 26, 40)

for (i in seq_len(nrow(sub_df))) {
  if (i %in% number_row_to_proccess) {
    cmd <- sprintf(
      "Rscript 00_main_dispatch.R --ptype 3CA --row_number_3CA %s",
      shQuote(i)
    )
    message("Running: ", cmd)
    system(cmd)
  }
}
