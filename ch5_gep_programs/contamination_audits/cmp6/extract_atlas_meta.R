#!/usr/bin/env Rscript
# Extract a slim per-cell metadata table from the Chapter-2 consensus atlas
# object for the cMP6 / myeloid-contamination check.

RDS <- "/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS"
OUT <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/cmp6-myeloid-contamination-83948d/cmp6_check/out/atlas_meta_slim.csv.gz"

cat("[", format(Sys.time()), "] loading", RDS, "\n")
df <- readRDS(RDS)
cat("class:", paste(class(df), collapse=","), " dim:", paste(dim(df), collapse=" x "), "\n")

want <- c("cell_id","cell_name","orig.ident","cancer_type","project_name",
          "sample","patient","paper","ORIGIN","Consensus_Cell_Type",
          "maj6","nCount_RNA","nFeature_RNA","doublet_score","pct_counts_mt")
have <- intersect(want, colnames(df))
cat("keeping columns:", paste(have, collapse=", "), "\n")
missing <- setdiff(want, colnames(df))
if (length(missing)) cat("MISSING:", paste(missing, collapse=", "), "\n")

slim <- as.data.frame(df)[, have, drop=FALSE]

# quick sanity summaries
cat("\n--- Consensus_Cell_Type ---\n"); print(table(slim$Consensus_Cell_Type, useNA="ifany"))
for (cc in c("cell_id","cell_name","sample","patient")) {
  if (cc %in% names(slim)) {
    v <- slim[[cc]]
    cat(sprintf("\n%s: n_unique=%d  n_NA=%d  examples: %s\n", cc, length(unique(v)),
                sum(is.na(v)), paste(utils::head(unique(v[!is.na(v)]),4), collapse=" | ")))
  }
}

cat("\n[", format(Sys.time()), "] writing", OUT, "\n")
con <- gzfile(OUT, "w")
write.csv(slim, con, row.names = FALSE)
close(con)
cat("done. nrow =", nrow(slim), "\n")
