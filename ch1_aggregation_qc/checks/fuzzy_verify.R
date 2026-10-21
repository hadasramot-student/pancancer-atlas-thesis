#!/usr/bin/env Rscript
# fuzzy_verify.R - independent fuzzy-match verification of "not yet processed"
# projects against all_cells.RDS, searching multiple metadata columns with
# loose matching (not exact project_name equality).

suppressPackageStartupMessages({ library(dplyr) })

TO <- "/ofir_grp/users/hadasra/temp_objects"
df <- readRDS(file.path(TO, "all_cells.gapfix_20260727_164657.RDS"))
cat("Loaded all_cells.gapfix:", nrow(df), "rows,", ncol(df), "cols\n")
cat("Columns:", paste(colnames(df), collapse=", "), "\n\n")

# candidate columns to search - anything that could carry a source identifier
search_cols <- intersect(
  c("project_name", "orig.ident", "sample", "paper", "SRX_accession",
    "patient", "ORIGIN", "cell_name", "BIOPSY", "BIOPSY_GROUP"),
  colnames(df)
)
cat("Searching columns:", paste(search_cols, collapse=", "), "\n\n")

# pre-extract unique values per column (fast lookup, avoids scanning 8.27M cells per term)
uniq_vals <- lapply(search_cols, function(cn) unique(as.character(df[[cn]])))
names(uniq_vals) <- search_cols

samples <- list(
  list(id="ZHANG_Y_ET_AL_2021_P39", terms=c("ZHANG_Y", "P39", "Zhang.*2021")),
  list(id="KARAAYVAZ_M_ET_AL_2018_P20", terms=c("KARAAYVAZ", "P20")),
  list(id="Pal2021_Breast", terms=c("Pal2021", "^Pal$")),
  list(id="Cillo2020_Head and Neck", terms=c("Cillo2020", "Cillo")),
  list(id="Liu2021_Hematologic", terms=c("Liu2021")),
  list(id="Kürten2021_Head and Neck", terms=c("Kurten2021", "K.rten2021", "Kurten")),
  list(id="Ferrari de Andrade2019_Skin", terms=c("Ferrari", "Andrade2019")),
  list(id="Caron2020_Hematologic", terms=c("Caron2020")),
  list(id="Roider2020_Hematologic", terms=c("Roider2020")),
  list(id="Nam2019_Hematologic", terms=c("Nam2019")),
  list(id="SRP329496", terms=c("SRP329496")),
  list(id="SRP445031", terms=c("SRP445031")),
  list(id="SRP421766", terms=c("SRP421766")),
  list(id="SRP519697", terms=c("SRP519697")),
  list(id="SRP457645", terms=c("SRP457645")),
  list(id="SRP484945", terms=c("SRP484945")),
  list(id="SRP498406", terms=c("SRP498406")),
  list(id="SRP519065", terms=c("SRP519065")),
  list(id="SRP484621", terms=c("SRP484621")),
  list(id="SRP485164", terms=c("SRP485164")),
  list(id="SRP506243", terms=c("SRP506243")),
  list(id="SRP515136", terms=c("SRP515136"))
)

results <- list()
for (s in samples) {
  hits <- list()
  for (cn in search_cols) {
    vals <- uniq_vals[[cn]]
    for (term in s$terms) {
      m <- grep(term, vals, ignore.case = TRUE, value = TRUE, perl = TRUE)
      if (length(m) > 0) {
        hits[[paste0(cn, ":", term)]] <- m
      }
    }
  }
  cat("===", s$id, "===\n")
  if (length(hits) == 0) {
    cat("  NO MATCHES found in any searched column for any search term.\n")
  } else {
    for (h in names(hits)) {
      cat("  MATCH in", h, "->", paste(head(hits[[h]], 5), collapse=" | "),
          if (length(hits[[h]]) > 5) paste0(" ...(+", length(hits[[h]])-5, " more)") else "", "\n")
    }
  }
  cat("\n")
  results[[s$id]] <- hits
}

n_with_hits <- sum(sapply(results, length) > 0)
cat("\n=== SUMMARY: ", n_with_hits, "of", length(samples), "sampled 'not yet processed' projects had ANY fuzzy match in all_cells.RDS ===\n")
