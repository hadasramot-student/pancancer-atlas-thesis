suppressPackageStartupMessages({library(dplyr)})
obj <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")
meta <- as.data.frame(obj)

cat("=== raw table(patient) top 20, no filtering at all ===\n")
print(head(sort(table(meta$patient, useNA="always"), decreasing=TRUE), 20))

cat("\n=== str() of patient column ===\n")
str(meta$patient)

cat("\n=== GAO_R_ET_AL_2021_P17_filtered rows: patient values (expect TNBC1..TNBC5) ===\n")
sub <- meta[meta$project_name == "GAO_R_ET_AL_2021_P17_filtered", c("project_name","sample","patient","paper")]
print(table(sub$patient, useNA="always"))
cat("\nfirst 5 raw rows:\n")
print(head(sub, 5))

cat("\n=== Azizi2018_Breast rows (3CA, known patient col in source Meta-data.csv) ===\n")
sub2 <- meta[meta$project_name == "Azizi2018_Breast", c("project_name","sample","patient")]
print(table(sub2$patient, useNA="always"))
print(head(sub2, 5))
