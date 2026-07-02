#!/usr/bin/env Rscript
# run_genenmf_per_cancer.R
#
# Ch.5 (GEP programs) — GeneNMF (Carmona lab 2023, carmonalab/GeneNMF)
# meta-program discovery for one cancer type (epithelial/tumor compartment),
# plus a Jaccard comparison of the resulting meta-programs (MPs) against the
# cNMF consensus programs (K=20) already computed for that cancer.
#
# Pipeline per cancer type:
#   1. Load the epithelial / tumour Seurat object (Brain = malignant-filtered other_Brain,
#      already saved under the Epithelial_Brain name).
#   2. JoinLayers (Seurat v5) + NormalizeData.
#   3. Split into per-sample objects, drop samples with too few cells.
#   4. multiNMF(): run NMF per sample over a range of k.
#   5. getMetaPrograms(): cluster per-sample programs into cross-sample meta-programs.
#   6. Extract top-N (default 50) genes per meta-program.
#   7. Jaccard(top50 MP, top50 cNMF program); MPs with best Jaccard > thr are
#      "recovered by both methods" = high-confidence programs.
#
# All outputs go to <output_dir>/<Cancer>/.
#
# Container note: this typically runs inside an apptainer Rstudio image with
# a bind mount — set --base_dir/--cnmf_base/--output_dir (or ATLAS_DATA_DIR)
# to whatever path the container sees, which may differ from the host path.
#
# Usage:
#   Rscript run_genenmf_per_cancer.R --cancer Lung
#   Rscript run_genenmf_per_cancer.R --cancer Brain --sample_col orig.ident
#
# (Submitted on SLURM via submit_genenmf_all_cancers.sh + ~/rstudio/rstudio.sbatch.)

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
if (nzchar(userlib)) .libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--cancer",      type = "character",
              help = "Cancer type, e.g. Lung (matches merged_by_query/<Cancer>/)"),
  make_option("--input",       type = "character", default = NULL,
              help = "Input Seurat RDS. Default: <base>/<Cancer>/Epithelial_<Cancer>_merged_full_scanvi_final.RDS"),
  make_option("--base_dir",    type = "character",
              default = file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "merged_by_query"),
              help = "Base dir holding per-cancer merged objects"),
  make_option("--cnmf_base",   type = "character",
              default = file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "cnmf_output"),
              help = "Base dir holding per-cancer cNMF output"),
  make_option("--output_dir",  type = "character",
              default = file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "genenmf_output"),
              help = "Base output dir; results land in <output_dir>/<Cancer>/"),
  make_option("--filter_celltype", type = "character", default = NULL,
              help = "If set, subset to cells where cell_type == this value (e.g. Malignant for Brain)"),
  make_option("--sample_col",  type = "character", default = NULL,
              help = "Metadata column defining a sample. Default: auto (sample > orig.ident > project_name)"),
  make_option("--k_range",     type = "character", default = "4:9",
              help = "Per-sample NMF k values, R expression, e.g. '4:9' or 'c(5,6,7)' [4:9]"),
  make_option("--n_mp",        type = "integer",   default = 20,
              help = "Number of meta-programs to extract (matches cNMF K=20) [20]"),
  make_option("--nfeatures",   type = "integer",   default = 2000,
              help = "HVGs used per-sample by multiNMF [2000]"),
  make_option("--top_n",       type = "integer",   default = 50,
              help = "Top genes per meta-program / cNMF program for Jaccard [50]"),
  make_option("--min_cells",   type = "integer",   default = 50,
              help = "Drop samples with fewer than this many cells [50]"),
  make_option("--min_samples", type = "integer",   default = 5,
              help = "Minimum #samples required to run cross-sample meta-programs [5]"),
  make_option("--jaccard_thr", type = "double",    default = 0.1,
              help = "Jaccard threshold for 'recovered by both' [0.1]"),
  make_option("--cnmf_k",      type = "integer",   default = 20,
              help = "cNMF K used in the consensus program filenames [20]"),
  make_option("--cnmf_dt",     type = "character", default = "0_1",
              help = "cNMF density-threshold tag in filenames (dt_<tag>) [0_1]"),
  make_option("--seed",        type = "integer",   default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$cancer)) stop("--cancer is required")

set.seed(opt$seed)
log_msg <- function(...) cat(sprintf("[%s] %s\n", Sys.time(), sprintf(...)))

if (!requireNamespace("GeneNMF", quietly = TRUE)) {
  stop("GeneNMF is not installed. Run install_genenmf.R first ",
       "(devtools::install_github('carmonalab/GeneNMF')).")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(GeneNMF)
  library(Matrix)
})

cancer  <- opt$cancer
out_dir <- file.path(opt$output_dir, cancer)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input <- opt$input
if (is.null(input)) {
  input <- file.path(opt$base_dir, cancer,
                     sprintf("Epithelial_%s_merged_full_scanvi_final.RDS", cancer))
}
if (!file.exists(input)) stop("Input object not found: ", input)

k_values <- eval(parse(text = opt$k_range))

log_msg("==== GeneNMF: %s ====", cancer)
log_msg("Input : %s", input)
log_msg("Output: %s", out_dir)
log_msg("k = %s | n_mp = %d | nfeatures = %d | top_n = %d",
        paste(range(k_values), collapse=":"), opt$n_mp, opt$nfeatures, opt$top_n)

# ── 1. Load object ────────────────────────────────────────────────────────────
log_msg("Loading Seurat object...")
obj <- readRDS(input)
DefaultAssay(obj) <- "RNA"
log_msg("Loaded: %d cells x %d genes", ncol(obj), nrow(obj))

# Optional cell-type filter (Brain: other_Brain -> cell_type == "Malignant").
if (!is.null(opt$filter_celltype)) {
  if (!"cell_type" %in% colnames(obj@meta.data))
    stop("--filter_celltype given but 'cell_type' column not in metadata. Available: ",
         paste(colnames(obj@meta.data), collapse = ", "))
  log_msg("cell_type distribution:")
  print(sort(table(obj@meta.data$cell_type), decreasing = TRUE))
  keep_ct <- obj@meta.data$cell_type == opt$filter_celltype
  if (sum(keep_ct, na.rm = TRUE) == 0)
    stop("No cells with cell_type == '", opt$filter_celltype, "'")
  obj <- obj[, which(keep_ct)]
  log_msg("Filtered to cell_type == '%s': %d cells", opt$filter_celltype, ncol(obj))
}

# Seurat v5: collapse split count layers, then (re)normalise from counts.
if (inherits(obj[["RNA"]], "Assay5")) {
  log_msg("JoinLayers (Seurat v5 multi-layer object)")
  obj <- JoinLayers(obj)
}
log_msg("NormalizeData (log-normalize)")
obj <- NormalizeData(obj, verbose = FALSE)

# ── 2. Choose the per-sample split column ─────────────────────────────────────
# The 'sample' column is unreliable here (1 level for some cancers, <5 for others);
# 'orig.ident' is the true per-patient/per-sample id. Auto-pick the candidate with
# the most groups passing the cell-count floor, preferring finer granularity.
n_usable <- function(meta, col, min_cells) {
  if (!col %in% colnames(meta)) return(-1L)
  v <- meta[[col]]
  v <- v[!is.na(v) & !v %in% c("", "NA", "NaN")]
  if (length(v) == 0) return(-1L)
  sum(table(v) >= min_cells)
}
pick_sample_col <- function(obj, requested, min_cells) {
  meta <- obj@meta.data
  if (!is.null(requested)) {
    if (!requested %in% colnames(meta))
      stop("--sample_col '", requested, "' not in metadata. Available: ",
           paste(colnames(meta), collapse = ", "))
    return(requested)
  }
  cands  <- c("orig.ident", "sample", "project_name")
  scores <- vapply(cands, function(c) n_usable(meta, c, min_cells), integer(1))
  log_msg("Usable samples (>= %d cells) per candidate column: %s",
          min_cells, paste(sprintf("%s=%d", cands, scores), collapse = ", "))
  if (max(scores) < 1)
    stop("No usable sample column found (tried ", paste(cands, collapse = ", "), ").")
  cands[which.max(scores)]   # ties -> first (orig.ident preferred)
}
sample_col <- pick_sample_col(obj, opt$sample_col, opt$min_cells)
log_msg("Sample column: '%s'", sample_col)

# Drop cells with missing sample label.
keep <- !is.na(obj@meta.data[[sample_col]]) &
        !obj@meta.data[[sample_col]] %in% c("", "NA", "NaN")
if (sum(!keep) > 0) {
  log_msg("Dropping %d cells with missing '%s'", sum(!keep), sample_col)
  obj <- obj[, keep]
}

# Drop samples below the cell-count floor.
tab <- table(obj@meta.data[[sample_col]])
good_samples <- names(tab)[tab >= opt$min_cells]
log_msg("Samples: %d total, %d with >= %d cells",
        length(tab), length(good_samples), opt$min_cells)
if (length(good_samples) < opt$min_samples) {
  stop(sprintf("Only %d samples have >= %d cells (need >= %d). Lower --min_cells or --min_samples.",
               length(good_samples), opt$min_cells, opt$min_samples))
}
obj <- obj[, obj@meta.data[[sample_col]] %in% good_samples]
obj@meta.data[[sample_col]] <- droplevels(factor(obj@meta.data[[sample_col]]))
log_msg("Cells after sample filtering: %d", ncol(obj))

# ── 3. Split per sample ───────────────────────────────────────────────────────
log_msg("Splitting into %d per-sample objects", length(good_samples))
seu.list <- SplitObject(obj, split.by = sample_col)
rm(obj); gc()

# ── 4. Per-sample NMF across k ────────────────────────────────────────────────
log_msg("multiNMF() over k = {%s} ...", paste(k_values, collapse = ","))
geneNMF.programs <- multiNMF(
  seu.list,
  assay     = "RNA",
  slot      = "data",
  k         = k_values,
  nfeatures = opt$nfeatures,
  hvg       = NULL,
  min.cells.per.sample = opt$min_cells
)
saveRDS(geneNMF.programs, file.path(out_dir, sprintf("%s_genenmf_programs.rds", cancer)))
log_msg("Per-sample programs saved.")

# ── 5. Cluster programs into meta-programs ────────────────────────────────────
log_msg("getMetaPrograms(nMP = %d, max.genes = %d) ...", opt$n_mp, opt$top_n)
geneNMF.mp <- getMetaPrograms(
  geneNMF.programs,
  nMP       = opt$n_mp,
  max.genes = opt$top_n,
  hclust.method = "ward.D2"
)
saveRDS(geneNMF.mp, file.path(out_dir, sprintf("%s_genenmf_metaprograms.rds", cancer)))

# metaprograms.genes: named list (MP1..MPn) of gene vectors (<= top_n)
mp_genes <- geneNMF.mp$metaprograms.genes
mp_genes <- lapply(mp_genes, function(g) head(g, opt$top_n))
log_msg("Extracted %d meta-programs.", length(mp_genes))

# Long table: cancer, MP, rank, gene
mp_long <- do.call(rbind, lapply(names(mp_genes), function(mp) {
  g <- mp_genes[[mp]]
  if (length(g) == 0) return(NULL)
  data.frame(cancer = cancer, meta_program = mp,
             rank = seq_along(g), gene = g, stringsAsFactors = FALSE)
}))
write.csv(mp_long, file.path(out_dir, sprintf("%s_genenmf_MP_top%d_genes.csv", cancer, opt$top_n)),
          row.names = FALSE)

# Per-MP QC metrics, if available
if (!is.null(geneNMF.mp$metaprograms.metrics)) {
  m <- as.data.frame(geneNMF.mp$metaprograms.metrics)
  m$cancer <- cancer
  write.csv(m, file.path(out_dir, sprintf("%s_genenmf_MP_metrics.csv", cancer)), row.names = TRUE)
}

# ── 6. Load cNMF consensus programs (top genes) ───────────────────────────────
cnmf_score <- file.path(
  opt$cnmf_base, cancer, sprintf("cnmf_%s", cancer),
  sprintf("cnmf_%s.gene_spectra_score.k_%d.dt_%s.txt", cancer, opt$cnmf_k, opt$cnmf_dt))

compare_done <- FALSE
if (!file.exists(cnmf_score)) {
  log_msg("WARNING: cNMF spectra file not found, skipping comparison: %s", cnmf_score)
} else {
  log_msg("Loading cNMF spectra: %s", basename(cnmf_score))
  # rows = programs (1..K), columns = genes, values = z-scored gene loadings
  spectra <- read.table(cnmf_score, header = TRUE, row.names = 1,
                        sep = "\t", check.names = FALSE)
  cnmf_genes <- apply(spectra, 1, function(r) {
    names(sort(r, decreasing = TRUE))[seq_len(opt$top_n)]
  })  # matrix: top_n x K
  cnmf_list <- lapply(seq_len(ncol(cnmf_genes)), function(j) cnmf_genes[, j])
  names(cnmf_list) <- paste0("cNMF_", rownames(spectra))

  # cNMF top genes table
  cnmf_long <- do.call(rbind, lapply(names(cnmf_list), function(p) {
    data.frame(cancer = cancer, cnmf_program = p,
               rank = seq_along(cnmf_list[[p]]), gene = cnmf_list[[p]],
               stringsAsFactors = FALSE)
  }))
  write.csv(cnmf_long, file.path(out_dir, sprintf("%s_cnmf_top%d_genes.csv", cancer, opt$top_n)),
            row.names = FALSE)

  # ── 7. Jaccard MP vs cNMF ───────────────────────────────────────────────────
  jaccard <- function(a, b) {
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    length(intersect(a, b)) / length(union(a, b))
  }
  J <- outer(seq_along(mp_genes), seq_along(cnmf_list),
             Vectorize(function(i, j) jaccard(mp_genes[[i]], cnmf_list[[j]])))
  rownames(J) <- names(mp_genes)
  colnames(J) <- names(cnmf_list)
  write.csv(round(J, 4),
            file.path(out_dir, sprintf("%s_jaccard_MP_vs_cNMF.csv", cancer)))

  # Best cNMF match per MP + high-confidence flag
  best_j   <- apply(J, 1, max)
  best_hit <- colnames(J)[apply(J, 1, which.max)]
  summary_df <- data.frame(
    cancer         = cancer,
    meta_program   = rownames(J),
    best_cnmf_match = best_hit,
    best_jaccard   = round(best_j, 4),
    recovered_by_both = best_j > opt$jaccard_thr,
    n_mp_genes     = vapply(mp_genes, length, integer(1)),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df,
            file.path(out_dir, sprintf("%s_highconfidence_programs.csv", cancer)),
            row.names = FALSE)
  compare_done <- TRUE
  log_msg("High-confidence (Jaccard > %.2f): %d / %d meta-programs.",
          opt$jaccard_thr, sum(summary_df$recovered_by_both), nrow(summary_df))
}

log_msg("==== Done: %s ====", cancer)
log_msg("Outputs in %s :", out_dir)
log_msg("  %s_genenmf_metaprograms.rds", cancer)
log_msg("  %s_genenmf_MP_top%d_genes.csv", cancer, opt$top_n)
if (compare_done) {
  log_msg("  %s_jaccard_MP_vs_cNMF.csv", cancer)
  log_msg("  %s_highconfidence_programs.csv", cancer)
}
