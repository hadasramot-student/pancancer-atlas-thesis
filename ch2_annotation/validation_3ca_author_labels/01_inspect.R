# 01_inspect.R -- diagnostic pass over the current merged atlas metadata table.
# Goal: identify exact column names for (a) consensus meta-type, (b) repository
# of origin, (c) original 3CA author cell-type label, plus cancer-type column;
# confirm total row count and the 5-study dedup arithmetic.

options(width = 200, max.print = 100000)
TO <- "/ofir_grp/users/hadasra/temp_objects"
f  <- file.path(TO, "all_cells.gapfix_20260727_164657.RDS")
cat("Reading:", f, "\n")
cat("file mtime:", format(file.info(f)$mtime), " size:", file.info(f)$size, "\n\n")

obj <- readRDS(f)
cat("class:", paste(class(obj), collapse=","), "\n")
cat("dim  :", paste(dim(obj), collapse=" x "), "\n")
cat("nrow :", nrow(obj), "\n\n")

md <- if (inherits(obj, "Seurat")) obj@meta.data else as.data.frame(obj)
cat("=== metadata dim:", paste(dim(md), collapse=" x "), "===\n\n")

cn <- colnames(md)
cat("=== ALL COLUMN NAMES (", length(cn), ") ===\n", sep="")
print(cn)
cat("\n")

cat("=== PER-COLUMN SUMMARY ===\n")
for (c in cn) {
  v <- md[[c]]
  cl <- paste(class(v), collapse="/")
  nna <- sum(is.na(v))
  if (is.numeric(v)) {
    cat(sprintf("\n--- %s [%s]  NA=%d  range=[%s, %s]\n", c, cl, nna,
                format(min(v, na.rm=TRUE)), format(max(v, na.rm=TRUE))))
  } else {
    vv <- as.character(v)
    nd <- length(unique(vv[!is.na(vv)]))
    cat(sprintf("\n--- %s [%s]  NA=%d  n_distinct=%d\n", c, cl, nna, nd))
    tb <- sort(table(vv, useNA="no"), decreasing=TRUE)
    print(utils::head(tb, 40))
  }
}

cat("\n\n=== CANDIDATE REPOSITORY COLUMNS: cross-tab with any 3CA-ish label col ===\n")
repo_cands  <- grep("repo|source|dataset|cohort|ptype|origin|db$|platform", cn, ignore.case=TRUE, value=TRUE)
label_cands <- grep("^cell_type$|m\\.cell_type|author|orig.*type|original|celltype|cell_type", cn, ignore.case=TRUE, value=TRUE)
meta_cands  <- grep("meta_type|meta.type|consensus|maj2|maj3|maj6|Meta\\.", cn, ignore.case=TRUE, value=TRUE)
cat("repo candidates :", paste(repo_cands, collapse=", "), "\n")
cat("label candidates:", paste(label_cands, collapse=", "), "\n")
cat("meta  candidates:", paste(meta_cands, collapse=", "), "\n\n")

for (rc in repo_cands) {
  cat("\n#### table(", rc, ") ####\n", sep="")
  print(sort(table(md[[rc]], useNA="ifany"), decreasing=TRUE))
  for (lc in label_cands) {
    nonNA <- tapply(md[[lc]], md[[rc]], function(x) sum(!is.na(x) & x != "" & x != "NA"))
    cat("  non-missing '", lc, "' by ", rc, ":\n", sep="")
    print(nonNA)
  }
}

cat("\n\n=== project_name: the 5 flagged Breastscape duplicate studies ===\n")
dups <- c("WU_SZ_ET_AL_2021_P15","GAO_R_ET_AL_2021_P17","GRIFFITHS_JI_ET_AL_2021_P18",
          "GONZALEZ_H_ET_AL_2022_P31","SAVAS_P_ET_AL_2018_P35")
if ("project_name" %in% cn) {
  pn <- md$project_name
  for (d in dups) cat(sprintf("  %-32s rows=%d\n", d, sum(pn == d, na.rm=TRUE)))
  cat("  total rows in these 5 :", sum(pn %in% dups, na.rm=TRUE), "\n")
  cat("  nrow after removing 5 :", sum(!(pn %in% dups)), "\n")
  cat("  (target corrected count: 8,107,665)\n")
  cat("  n distinct project_name:", length(unique(pn)), "\n")
} else cat("  NO project_name column!\n")

cat("\n=== done ===\n")
