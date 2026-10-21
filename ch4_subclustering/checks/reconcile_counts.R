#!/usr/bin/env Rscript
# reconcile_counts.R
# Trace ncol() through each pipeline stage for every (lineage x cancer_type)
# combo, and compare against the expected counts implied by all_cells.RDS
# (the Ch.2 Consensus_Cell_Type / Table 3 source of truth).

suppressPackageStartupMessages({
  library(dplyr)
})

base <- "/ofir_grp/users/hadasra/temp_objects"
mbq  <- file.path(base, "merged_by_query")
clus <- file.path(base, "clustered_lineages")
ann  <- file.path(base, "annotated_lineages")

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── Expected counts from all_cells.RDS ──────────────────────────────────────
all_cells <- readRDS(file.path(base, "all_cells.RDS"))
lineages  <- c("T-cells","Myeloid","Epithelial","B-cell","Fibroblasts",
               "Endothelial","NK-cells","Hematopoietic_Stem_and_Progenitor_cells")

expected <- all_cells %>%
  filter(Consensus_Cell_Type %in% lineages) %>%
  count(Consensus_Cell_Type, cancer_type, name = "expected_n")

log_msg("Expected atlas-wide totals per lineage (Table 3 / Ch2 source of truth):")
print(all_cells %>% count(Consensus_Cell_Type, sort = TRUE))

cancer_types <- sort(unique(all_cells$cancer_type))
safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

results <- list()

for (lineage in lineages) {
  for (ct in cancer_types) {
    lin_s <- safe(lineage)
    ct_s  <- safe(ct)

    exp_n <- expected$expected_n[expected$Consensus_Cell_Type == lineage & expected$cancer_type == ct]
    exp_n <- if (length(exp_n) == 0) 0 else exp_n
    if (exp_n == 0) next

    raw_path   <- file.path(mbq, ct, paste0(lin_s, "_", ct_s, "_merged.RDS"))
    scanvi_path<- file.path(mbq, ct, paste0(lin_s, "_", ct_s, "_merged_full_scanvi_final.RDS"))
    clus_path  <- file.path(clus, paste0(lin_s, "_", ct_s, "_clustered.RDS"))
    ann_path   <- file.path(ann, paste0(lin_s, "_", ct_s, "_annotated.RDS"))

    get_n <- function(p) {
      if (!file.exists(p)) return(NA_integer_)
      obj <- tryCatch(readRDS(p), error = function(e) NULL)
      if (is.null(obj)) return(NA_integer_)
      n <- tryCatch(ncol(obj), error = function(e) NA_integer_)
      rm(obj); gc(FALSE)
      n
    }

    n_raw    <- get_n(raw_path)
    n_scanvi <- get_n(scanvi_path)
    n_clus   <- get_n(clus_path)
    n_ann    <- get_n(ann_path)

    results[[length(results) + 1]] <- data.frame(
      lineage = lineage, cancer_type = ct, expected_n = exp_n,
      raw_merge_n = n_raw, scanvi_final_n = n_scanvi,
      clustered_n = n_clus, annotated_n = n_ann
    )
    log_msg(sprintf("%-42s %-15s expected=%-8d raw=%-8s scanvi=%-8s clustered=%-8s annotated=%-8s",
                    lineage, ct, exp_n,
                    ifelse(is.na(n_raw),"NA",n_raw),
                    ifelse(is.na(n_scanvi),"NA",n_scanvi),
                    ifelse(is.na(n_clus),"NA",n_clus),
                    ifelse(is.na(n_ann),"NA",n_ann)))
  }
}

res_df <- do.call(rbind, results)
out_csv <- file.path(base, "reconcile_counts_result.csv")
write.csv(res_df, out_csv, row.names = FALSE)
log_msg("Saved: ", out_csv)

log_msg("=== Per-lineage totals (summed across cancer types, using best-available stage) ===")
res_df$best_n <- ifelse(!is.na(res_df$clustered_n), res_df$clustered_n,
                  ifelse(!is.na(res_df$annotated_n), res_df$annotated_n,
                  ifelse(!is.na(res_df$scanvi_final_n), res_df$scanvi_final_n,
                  ifelse(!is.na(res_df$raw_merge_n), res_df$raw_merge_n, 0))))

summary_df <- res_df %>%
  group_by(lineage) %>%
  summarise(expected_total = sum(expected_n),
            raw_total      = sum(raw_merge_n, na.rm = TRUE),
            scanvi_total   = sum(scanvi_final_n, na.rm = TRUE),
            clustered_total= sum(clustered_n, na.rm = TRUE),
            best_total     = sum(best_n))
print(as.data.frame(summary_df))
write.csv(summary_df, file.path(base, "reconcile_counts_summary.csv"), row.names = FALSE)
log_msg("Saved summary: ", file.path(base, "reconcile_counts_summary.csv"))
log_msg("DONE")
