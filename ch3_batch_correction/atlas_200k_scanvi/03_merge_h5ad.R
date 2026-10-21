# 03_merge_h5ad.R -- cbind the per-project count submatrices into one
# 28,231 x ~200,000 matrix, attach the drawn metadata, LogNormalize, write the
# h5ad the Python scVI/scANVI + PCA step consumes.
#
# "counts" fed downstream = the step4 counts layer AS-IS (SCimilarity log-norm,
# 28,231-gene panel) per the match-the-atlas decision (user, 2026-08-31).
# We still LogNormalize it (scale.factor 1e4) exactly as Seurat's NormalizeData
# does -- that mirrors the standard pipeline step and gives Figure 2's raw-PCA
# branch its input. No Seurat dependency (avoids cell-name munging of the
# "project|cell_id" row_uid).
suppressPackageStartupMessages({
  library(optparse); library(Matrix)
  library(SingleCellExperiment); library(zellkonverter); library(S4Vectors)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--root", type = "character",
              default = "/ofir_grp/users/hadasra/temp_objects/figure1_2_scanvi")
)))
lg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
ROOT  <- opt$root
CDIR  <- file.path(ROOT, "per_project_counts")
SCALE <- 10000

drawn <- read.csv(file.path(ROOT, "drawn_cells.csv"), stringsAsFactors = FALSE)
lg("drawn cells: ", nrow(drawn))
stopifnot(!anyDuplicated(drawn$row_uid))

## ---- assemble -------------------------------------------------------------
proj <- sort(unique(drawn$project_name))
lg("projects to gather: ", length(proj))
mats <- vector("list", length(proj)); names(mats) <- proj
g0 <- NULL; tot <- 0L; miss_tot <- 0L; rep_rows <- list()
for (p in proj) {
  f <- file.path(CDIR, paste0(p, ".rds"))
  if (!file.exists(f)) stop("MISSING per-project counts: ", f)
  m <- readRDS(f)
  if (is.null(g0)) g0 <- rownames(m) else if (!identical(rownames(m), g0))
    stop("gene order mismatch in ", p)
  want_p <- drawn$cell_id[drawn$project_name == p]
  if (!setequal(colnames(m), want_p)) {
    miss <- setdiff(want_p, colnames(m)); miss_tot <- miss_tot + length(miss)
    lg("  ", p, ": have ", ncol(m), " / ", length(want_p), " (", length(miss), " missing)")
  }
  colnames(m) <- paste(p, colnames(m), sep = "|")
  mats[[p]] <- m; tot <- tot + ncol(m)
  rt <- file.path(CDIR, paste0(p, ".report.tsv"))
  if (file.exists(rt)) rep_rows[[p]] <- read.delim(rt, stringsAsFactors = FALSE)
}
lg("cells gathered: ", tot, " ; missing barcodes total: ", miss_tot)
if (length(rep_rows)) {
  allrep <- do.call(rbind, rep_rows)
  write.csv(allrep, file.path(ROOT, "assembly_report.csv"), row.names = FALSE)
  lg("assembly_report.csv written (", sum(allrep$n_missing), " missing across ",
     sum(allrep$n_missing > 0), " projects)")
}

counts <- do.call(cbind, mats); rm(mats); gc(verbose = FALSE)
counts <- as(counts, "CsparseMatrix")
lg("merged counts: ", nrow(counts), " x ", ncol(counts))

mi <- match(colnames(counts), drawn$row_uid)
stopifnot(!anyNA(mi))
meta <- drawn[mi, ]; rownames(meta) <- meta$row_uid
if (miss_tot > 0)
  lg("NOTE: ", ncol(counts), " of ", nrow(drawn), " drawn cells assembled (",
     nrow(drawn) - ncol(counts), " dropped as missing barcodes)")

## ---- LogNormalize (== Seurat NormalizeData, LogNormalize, scale.factor 1e4) --
lg("LogNormalize(scale.factor = ", SCALE, ")")
cs <- Matrix::colSums(counts)
stopifnot(all(cs > 0))
logn <- counts %*% Matrix::Diagonal(x = SCALE / cs)
logn <- as(logn, "CsparseMatrix")
logn@x <- log1p(logn@x)
dimnames(logn) <- dimnames(counts)

## ---- h5ad via zellkonverter (X = lognorm, layer 'counts' = step4 counts) ----
for (cc in colnames(meta)) if (is.factor(meta[[cc]])) meta[[cc]] <- as.character(meta[[cc]])
meta[] <- lapply(meta, function(x) if (is.character(x)) ifelse(is.na(x), "NA", x) else x)

sce <- SingleCellExperiment(
  assays  = list(logcounts = logn, counts = counts),
  colData = DataFrame(meta, row.names = rownames(meta))
)
saveRDS(sce, file.path(ROOT, "fig12_200k_sce.rds"))
out_h5 <- file.path(ROOT, "fig12_200k.h5ad")
writeH5AD(sce, file = out_h5, X_name = "logcounts")
lg("wrote ", out_h5, "  (X=logcounts, layers$counts=step4 counts, ", ncol(sce), " cells)")
lg("DONE 03_merge_h5ad")
