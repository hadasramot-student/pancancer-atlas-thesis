# 02_assemble_counts.R -- pull the drawn barcodes' expression columns out of one
# per-project _step4.RDS and save a slim genes x cells sparse matrix.
#
# DECISION (user, 2026-08-31): match the existing atlas -> use the step4 "counts"
# layer AS-IS. It holds SCimilarity log-normalised values (28,231-gene panel),
# which is exactly what the Ch.3 §3.4.5 benchmarked/selected atlas scANVI and every
# per-lineage x cancer-type embedding were trained on. Disclosed in methods.
#
# One SLURM array task = one project. Restartable: skips if output exists.
#
# args: --project <name> --root <figure1_2_scanvi dir> --to <temp_objects dir>
suppressPackageStartupMessages({ library(optparse); library(Seurat); library(Matrix) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--project", type = "character"),
  make_option("--root",    type = "character",
              default = "/ofir_grp/users/hadasra/temp_objects/figure1_2_scanvi"),
  make_option("--to",      type = "character",
              default = "/ofir_grp/users/hadasra/temp_objects")
)))
stopifnot(!is.null(opt$project))
lg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

bc_file  <- file.path(opt$root, "per_project_barcodes", paste0(opt$project, ".txt"))
out_rds  <- file.path(opt$root, "per_project_counts",   paste0(opt$project, ".rds"))
rep_file <- file.path(opt$root, "per_project_counts",   paste0(opt$project, ".report.tsv"))
step4    <- file.path(opt$to, paste0(opt$project, "_step4.RDS"))
stopifnot(file.exists(bc_file), file.exists(step4))
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

want <- readLines(bc_file)
want <- want[nzchar(want)]
lg("project ", opt$project, " : ", length(want), " barcodes requested")

if (file.exists(out_rds)) {
  m <- readRDS(out_rds)
  if (ncol(m) == length(want) && setequal(colnames(m), want)) {
    lg("output already complete -> skip"); quit(save = "no", status = 0)
  }
  lg("stale output present, rebuilding")
}

lg("readRDS ", step4, "  (", round(file.info(step4)$size / 1e9, 1), " GB)")
o <- readRDS(step4)
DefaultAssay(o) <- "RNA"
cn <- colnames(o)
found   <- want[want %in% cn]
missing <- want[!want %in% cn]
lg("found ", length(found), " / ", length(want), " ; missing ", length(missing))

cm <- tryCatch(
  SeuratObject::GetAssayData(o, assay = "RNA", layer = "counts"),
  error = function(e) SeuratObject::LayerData(o, assay = "RNA", layer = "counts"))
stopifnot(!is.null(cm), length(cm@x) > 0)
sub <- cm[, found, drop = FALSE]
sub <- as(sub, "CsparseMatrix")
lg("submatrix: ", nrow(sub), " genes x ", ncol(sub), " cells ; max=", signif(max(sub@x), 4),
   " ; integer-valued=", isTRUE(all(sub@x == floor(sub@x))))

saveRDS(sub, out_rds, compress = FALSE)
lg("wrote ", out_rds)

rep <- data.frame(
  project      = opt$project,
  n_requested  = length(want),
  n_found      = length(found),
  n_missing    = length(missing),
  n_genes      = nrow(sub),
  counts_max   = signif(max(sub@x), 6),
  counts_integer = isTRUE(all(sub@x == floor(sub@x))),
  missing_examples = paste(utils::head(missing, 5), collapse = ";")
)
write.table(rep, rep_file, sep = "\t", quote = FALSE, row.names = FALSE)
if (length(missing) > 0)
  lg("WARNING: ", length(missing), " requested barcodes not in step4 object for ", opt$project)
lg("DONE ", opt$project)
