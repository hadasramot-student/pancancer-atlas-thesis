#!/usr/bin/env Rscript
# subsample_epithelial_breast.R
#
# Chapter 4 gap-fix for Epithelial_Breast (1,164,224 cells) — the only
# lineage x cancer-type object that never clustered, because JoinLayers on the
# full merged object overflows the 32-bit dgCMatrix non-zero index limit
# (`p[length(p)] cannot exceed 2^31-1`, see
#  temp_objects/clustered_lineages/logs/cl_Epithelial_Breast_19998318.err).
#
# Supervisor-approved option 3: subsample to 450,000 cells FOR CLUSTERING ONLY,
# as a stated design choice, then run the standard Ch.4 pipeline unchanged.
#
# What this script does (and nothing else):
#   1. Load the existing scANVI-corrected object (embedding already computed —
#      batch correction is NOT redone).
#   2. Proportional STRATIFIED random subsample to exactly 450,000 cells,
#      stratified by `project_name` (study), seed 42. Study is the only usable
#      stratum key: `patient` is 77% NA and `sample` 35% NA in this object.
#   3. JoinLayers on the 450k subset (now well under the 2^31 ceiling) and
#      verify the non-zero count.
#   4. Save the subsampled object in the exact format the Ch.4 clustering
#      script expects (`*_merged_full_scanvi_final.RDS`, scANVI reduction kept).
#   5. Write before/after composition reports (per study / patient / sample /
#      Scimilarity cell type) for an accurate figure/table legend.
#
# Usage:
#   Rscript subsample_epithelial_breast.R \
#     --input  /ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/Epithelial_Breast_merged_full_scanvi_final.RDS \
#     --output /ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/Epithelial_Breast_subsample450k_merged_full_scanvi_final.RDS \
#     --report_dir /ofir_grp/users/hadasra/temp_objects/clustered_lineages \
#     --n 450000 --seed 42

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
})

option_list <- list(
  make_option("--input",      type = "character"),
  make_option("--output",     type = "character"),
  make_option("--report_dir", type = "character"),
  make_option("--n",          type = "double",  default = 450000),
  make_option("--seed",       type = "integer", default = 42L),
  make_option("--strat_col",  type = "character", default = "project_name")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), !is.null(opt$output), !is.null(opt$report_dir))
target <- as.integer(opt$n)
dir.create(opt$report_dir, recursive = TRUE, showWarnings = FALSE)

lg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...))); flush.console() }

# ── 1. Load ──────────────────────────────────────────────────────────────────
lg("Loading ", opt$input)
obj <- readRDS(opt$input)
lg("Loaded: ", ncol(obj), " cells x ", nrow(obj), " genes")
lg("Assay class : ", class(obj[["RNA"]])[1])
lg("RNA layers  : ", paste(SeuratObject::Layers(obj[["RNA"]]), collapse = ", "))
lg("Reductions  : ", paste(Reductions(obj), collapse = ", "))
scanvi_name <- intersect(c("scanvi","scANVI","SCANVI","latent.scanvi"), Reductions(obj))
if (length(scanvi_name) == 0) stop("No scANVI reduction found in object — cannot proceed.")
scanvi_name <- scanvi_name[1]
lg("scANVI reduction: '", scanvi_name, "'  dims = ", ncol(Embeddings(obj, scanvi_name)))
if (!opt$strat_col %in% colnames(obj@meta.data))
  stop("Stratification column '", opt$strat_col, "' not in metadata.")

md <- obj@meta.data
n_full <- nrow(md)

# ── 2. Proportional stratified allocation (largest-remainder rounding) ────────
strat <- as.character(md[[opt$strat_col]]); strat[is.na(strat)] <- "__NA__"
tab   <- table(strat)
raw   <- target * as.numeric(tab) / sum(as.numeric(tab))
base  <- pmin(floor(raw), as.numeric(tab))          # never ask for more than a stratum has
rem   <- target - sum(base)
if (rem > 0) {
  headroom <- as.numeric(tab) - base                # strata that can still take more
  ord <- order(raw - floor(raw), decreasing = TRUE)
  ord <- ord[headroom[ord] > 0]
  take <- ord[seq_len(min(rem, length(ord)))]
  base[take] <- base[take] + 1
}
alloc <- setNames(as.integer(base), names(tab))
lg("Target n = ", target, "  | allocated sum = ", sum(alloc),
   "  | sampling fraction = ", sprintf("%.4f", target / n_full))

set.seed(opt$seed)
cells_all <- colnames(obj)
sel_idx <- unlist(lapply(names(alloc), function(s) {
  pool <- which(strat == s)
  k <- min(alloc[[s]], length(pool))
  if (k <= 0) integer(0) else sample(pool, k)
}), use.names = FALSE)

# top up (or trim) to land on exactly `target`, respecting seed
if (length(sel_idx) < target) {
  extra_pool <- setdiff(seq_len(n_full), sel_idx)
  sel_idx <- c(sel_idx, sample(extra_pool, target - length(sel_idx)))
} else if (length(sel_idx) > target) {
  sel_idx <- sample(sel_idx, target)
}
sel_idx  <- sort(unique(sel_idx))
stopifnot(length(sel_idx) == target)
sel_cells <- cells_all[sel_idx]
lg("Selected ", length(sel_cells), " unique cells")

# ── 3. Subset + JoinLayers (the step that overflowed on the full object) ─────
lg("subset() to the 450k cells ...")
obj_sub <- subset(obj, cells = sel_cells)
rm(obj); gc(verbose = FALSE)
lg("Subset object: ", ncol(obj_sub), " cells x ", nrow(obj_sub), " genes")
lg("Subset reductions: ", paste(Reductions(obj_sub), collapse = ", "))

if (inherits(obj_sub[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(obj_sub[["RNA"]])) > 1) {
  lg("JoinLayers on the subset ...")
  ok <- tryCatch({ obj_sub[["RNA"]] <- JoinLayers(obj_sub[["RNA"]]); TRUE },
                 error = function(e) { lg("JoinLayers FAILED: ", conditionMessage(e)); FALSE })
  if (!ok) stop("JoinLayers still failed at 450k — stop here and report (do not retry blindly).")
}
DefaultAssay(obj_sub) <- "RNA"
cnts <- SeuratObject::GetAssayData(obj_sub, assay = "RNA", layer = "counts")
nnz <- length(cnts@x)
lg("Joined counts matrix: ", nrow(cnts), " x ", ncol(cnts),
   "  | non-zero entries = ", format(nnz, big.mark = ","),
   "  | 2^31-1 = ", format(2^31 - 1, big.mark = ","))
if (nnz >= 2^31 - 1) stop("Non-zero count still exceeds 2^31-1 after subsampling — unexpected.")
lg("OK: under the 32-bit ceiling with ", sprintf("%.1f%%", 100 * nnz / (2^31 - 1)), " headroom used.")

# ── 4. Save in the standard pipeline input format ───────────────────────────
lg("Saving subsampled object -> ", opt$output)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_sub, opt$output)
lg("Saved (", round(file.info(opt$output)$size / 1e9, 2), " GB)")

# ── 5. Composition reports (before vs after) ────────────────────────────────
comp_report <- function(col) {
  if (!col %in% colnames(md)) return(NULL)
  bf <- as.character(md[[col]]);            bf[is.na(bf)] <- "NA"
  af <- as.character(obj_sub@meta.data[[col]]); af[is.na(af)] <- "NA"
  keys <- sort(union(unique(bf), unique(af)))
  bt <- table(factor(bf, levels = keys)); at <- table(factor(af, levels = keys))
  data.frame(
    variable         = col,
    level            = keys,
    n_full           = as.integer(bt),
    pct_full         = round(100 * as.integer(bt) / length(bf), 4),
    n_subsample      = as.integer(at),
    pct_subsample    = round(100 * as.integer(at) / length(af), 4),
    pct_point_diff   = round(100 * as.integer(at) / length(af) -
                             100 * as.integer(bt) / length(bf), 4),
    row.names = NULL
  )
}
reports <- do.call(rbind, lapply(
  c(opt$strat_col, "patient", "sample", "paper",
    "scimilarity_celltype", "Consensus_Cell_Type"), comp_report))
out_csv <- file.path(opt$report_dir, "Epithelial_Breast_subsample_composition.csv")
write.csv(reports, out_csv, row.names = FALSE)
lg("Composition report -> ", out_csv)

# max drift on the stratification variable (should be ~0 for proportional strat)
st <- reports[reports$variable == opt$strat_col, ]
lg(sprintf("Max %s drift: %.4f pct-points (over %d studies)",
           opt$strat_col, max(abs(st$pct_point_diff)), nrow(st)))

# human-readable note that travels with the outputs
note <- c(
  "Epithelial_Breast — Chapter 4 clustering was run on a SUBSAMPLE, not the full object.",
  "",
  sprintf("Full object            : %s cells (largest lineage x cancer-type object in the atlas)", format(n_full, big.mark = ",")),
  sprintf("Subsample size         : %s cells", format(target, big.mark = ",")),
  sprintf("Sampling method        : proportional stratified random sampling by '%s' (study)", opt$strat_col),
  sprintf("                         within-stratum: uniform random, largest-remainder rounding to hit exactly %s", format(target, big.mark = ",")),
  sprintf("Random seed            : %d (set.seed, R %s)", opt$seed, paste(R.version$major, R.version$minor, sep = ".")),
  sprintf("Sampling fraction      : %.4f", target / n_full),
  sprintf("Strata (studies)       : %d; max per-study composition drift %.3f pct-points", nrow(st), max(abs(st$pct_point_diff))),
  "Reason                 : JoinLayers on the full object overflows the 32-bit dgCMatrix",
  "                         non-zero index limit (p[length(p)] cannot exceed 2^31-1).",
  "                         Subsampling keeps the merged sparse matrix under that ceiling.",
  "Batch correction       : NOT redone — the existing scANVI embedding for these cells is reused.",
  "Pipeline               : identical to the other lineage x cancer-type objects",
  "                         (cluster_lineage.R -> annotate_lineages.R -> find_markers.R -> gsea_hallmarks.R).",
  sprintf("Generated              : %s", format(Sys.time())),
  sprintf("Composition detail     : %s", basename(out_csv))
)
writeLines(note, file.path(opt$report_dir, "Epithelial_Breast_SUBSAMPLE_NOTE.txt"))
lg("Wrote Epithelial_Breast_SUBSAMPLE_NOTE.txt")
lg("DONE.")
