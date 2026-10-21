suppressPackageStartupMessages(library(Seurat))

files <- c(
  "ERP136281_step4.RDS",
  "ERP136992_step4.RDS",
  "SRP288163_step4.RDS"
)

total <- 0
for (f in files) {
  path <- paste0("/ofir_grp/users/hadasra/temp_objects/", f)
  cat("=== ", f, " ===\n")
  if (!file.exists(path)) {
    cat("FILE DOES NOT EXIST\n")
    next
  }
  obj <- readRDS(path)
  md <- obj@meta.data
  n <- nrow(md)
  cat("nrow meta.data:", n, "\n")
  total <- total + n
  if ("Consensus_Cell_Type" %in% colnames(md)) {
    print(table(md$Consensus_Cell_Type, useNA = "ifany"))
  } else {
    cat("No Consensus_Cell_Type column. Columns:\n")
    print(colnames(md))
  }
  cat("\n")
}
cat("TOTAL across these 3 projects:", total, "\n")

k_path <- "/ofir_grp/users/hadasra/temp_objects/KARAAYVAZ_M_ET_AL_2018_P20_filtered_step2.h5ad"
cat("\nKARAAYVAZ step2.h5ad exists:", file.exists(k_path), "\n")
