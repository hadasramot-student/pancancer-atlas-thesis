#!/usr/bin/env Rscript
# run_genenmf_per_lineage.R
#
# Ch.5 (GEP programs) — GeneNMF (carmonalab/GeneNMF) meta-program discovery
# for ONE immune lineage, on the SAME merged immune h5ad that Spectra uses,
# so the cross-method comparison is a valid 2-method comparison (Spectra vs
# GeneNMF) on the SAME immune cells.
#
# (cNMF/GeneNMF were originally run on epithelial cells per cancer type — a
#  different compartment — so we re-run GeneNMF here on the immune populations.)
#
# Pipeline:
#   1. Read <lineage>_immune.h5ad (anndata) -> Seurat (counts), NormalizeData.
#   2. Pick a per-sample split column; drop tiny samples.
#   3. multiNMF() per sample over k; getMetaPrograms() -> MPs.
#   4. Top-N genes per MP.
#   5. Jaccard(MP top genes, Spectra factor top genes) for the same lineage.
#
# Usage (inside apptainer Rstudio container via ~/rstudio/rstudio.sbatch):
#   Rscript run_genenmf_per_lineage.R --lineage B-cell \
#     --h5ad   /path/to/clustered_lineages/spectra_out/B-cell_immune.h5ad \
#     --spectra_factors /path/to/clustered_lineages/spectra_out/B-cell_spectra_gene_factors.csv \
#     --output_dir /path/to/temp_objects/genenmf_immune_output

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
if (nzchar(userlib)) .libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--lineage",  type = "character", help = "Lineage name, e.g. B-cell"),
  make_option("--h5ad",     type = "character", help = "Merged immune h5ad (Spectra input)"),
  make_option("--spectra_factors", type = "character", default = NULL,
              help = "Spectra <lineage>_spectra_gene_factors.csv (factors x genes) for Jaccard"),
  make_option("--output_dir", type = "character",
              default = file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "genenmf_immune_output")),
  make_option("--sample_col", type = "character", default = NULL,
              help = "Per-sample column. Default auto: project_name > sample > orig.ident"),
  make_option("--k_range",  type = "character", default = "4:9"),
  make_option("--n_mp",     type = "integer",   default = 20),
  make_option("--nfeatures",type = "integer",   default = 2000),
  make_option("--top_n",    type = "integer",   default = 50),
  make_option("--min_cells",type = "integer",   default = 50),
  make_option("--min_samples", type = "integer", default = 5)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$lineage)) stop("--lineage is required")
if (is.null(opt$h5ad))    stop("--h5ad is required")

log_msg <- function(fmt, ...) cat(sprintf(paste0("[%s] ", fmt, "\n"), Sys.time(), ...))

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(anndata); library(GeneNMF)
})

k_values <- eval(parse(text = opt$k_range))
lineage  <- opt$lineage
out_dir  <- file.path(opt$output_dir, lineage)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_msg("==== GeneNMF (immune): %s ====", lineage)
log_msg("k = %s | n_mp = %d | nfeatures = %d | top_n = %d",
        paste(range(k_values), collapse = ":"), opt$n_mp, opt$nfeatures, opt$top_n)

# ── 1. h5ad -> Seurat ─────────────────────────────────────────────────────────
log_msg("Reading h5ad: %s", opt$h5ad)
ad <- read_h5ad(opt$h5ad)
counts <- ad$layers[["counts"]]
if (is.null(counts)) stop("h5ad has no layers['counts']")
counts <- Matrix::t(counts)                      # anndata is cells x genes -> genes x cells
rownames(counts) <- ad$var_names
colnames(counts) <- ad$obs_names
meta <- as.data.frame(ad$obs)
rownames(meta) <- ad$obs_names
log_msg("Building Seurat: %d genes x %d cells", nrow(counts), ncol(counts))
obj <- CreateSeuratObject(counts = counts, meta.data = meta)
rm(counts, ad); gc()
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) { log_msg("JoinLayers"); obj <- JoinLayers(obj) }
log_msg("NormalizeData")
obj <- NormalizeData(obj, verbose = FALSE)

# ── 2. Per-sample split column ────────────────────────────────────────────────
n_usable <- function(meta, col, min_cells) {
  if (!col %in% colnames(meta)) return(-1L)
  v <- meta[[col]]; v <- v[!is.na(v) & !v %in% c("", "NA", "NaN", "Unknown")]
  if (length(v) == 0) return(-1L)
  sum(table(v) >= min_cells)
}
pick_sample_col <- function(obj, requested, min_cells) {
  meta <- obj@meta.data
  if (!is.null(requested)) {
    if (!requested %in% colnames(meta)) stop("--sample_col '", requested, "' not in metadata")
    return(requested)
  }
  cands  <- c("project_name", "sample", "orig.ident")
  scores <- vapply(cands, function(c) n_usable(meta, c, min_cells), integer(1))
  log_msg("Usable samples (>= %d cells): %s", min_cells,
          paste(sprintf("%s=%d", cands, scores), collapse = ", "))
  if (max(scores) < 1) stop("No usable sample column found.")
  cands[which.max(scores)]
}
sample_col <- pick_sample_col(obj, opt$sample_col, opt$min_cells)
log_msg("Sample column: '%s'", sample_col)

keep <- !is.na(obj@meta.data[[sample_col]]) &
        !obj@meta.data[[sample_col]] %in% c("", "NA", "NaN", "Unknown")
if (sum(!keep) > 0) { log_msg("Dropping %d cells with missing sample", sum(!keep)); obj <- obj[, keep] }

tab <- table(obj@meta.data[[sample_col]])
good_samples <- names(tab)[tab >= opt$min_cells]
log_msg("Samples: %d total, %d with >= %d cells", length(tab), length(good_samples), opt$min_cells)
if (length(good_samples) < opt$min_samples)
  stop(sprintf("Only %d samples >= %d cells (need %d).", length(good_samples), opt$min_cells, opt$min_samples))
obj <- obj[, obj@meta.data[[sample_col]] %in% good_samples]
obj@meta.data[[sample_col]] <- droplevels(factor(obj@meta.data[[sample_col]]))
log_msg("Cells after sample filtering: %d", ncol(obj))

# ── 3. Per-sample NMF + meta-programs ─────────────────────────────────────────
seu.list <- SplitObject(obj, split.by = sample_col)
rm(obj); gc()
log_msg("multiNMF() over k = {%s} on %d samples", paste(k_values, collapse = ","), length(seu.list))
geneNMF.programs <- multiNMF(seu.list, assay = "RNA", slot = "data", k = k_values,
                             nfeatures = opt$nfeatures, hvg = NULL,
                             min.cells.per.sample = opt$min_cells)
saveRDS(geneNMF.programs, file.path(out_dir, sprintf("%s_genenmf_programs.rds", lineage)))

log_msg("getMetaPrograms(nMP = %d) ...", opt$n_mp)
geneNMF.mp <- getMetaPrograms(geneNMF.programs, nMP = opt$n_mp, max.genes = opt$top_n,
                              hclust.method = "ward.D2")
saveRDS(geneNMF.mp, file.path(out_dir, sprintf("%s_genenmf_metaprograms.rds", lineage)))

mp_genes <- lapply(geneNMF.mp$metaprograms.genes, function(g) head(g, opt$top_n))
log_msg("Extracted %d meta-programs.", length(mp_genes))
mp_long <- do.call(rbind, lapply(names(mp_genes), function(mp) {
  g <- mp_genes[[mp]]; if (length(g) == 0) return(NULL)
  data.frame(lineage = lineage, meta_program = mp, rank = seq_along(g), gene = g)
}))
write.csv(mp_long, file.path(out_dir, sprintf("%s_genenmf_MP_top%d_genes.csv", lineage, opt$top_n)),
          row.names = FALSE)
if (!is.null(geneNMF.mp$metaprograms.metrics)) {
  m <- as.data.frame(geneNMF.mp$metaprograms.metrics); m$lineage <- lineage
  write.csv(m, file.path(out_dir, sprintf("%s_genenmf_MP_metrics.csv", lineage)), row.names = TRUE)
}

# ── 4. Jaccard MP vs Spectra factors (same lineage, same cells) ───────────────
if (!is.null(opt$spectra_factors) && file.exists(opt$spectra_factors)) {
  log_msg("Loading Spectra factors: %s", basename(opt$spectra_factors))
  sf <- read.csv(opt$spectra_factors, row.names = 1, check.names = FALSE)  # factors x genes
  spectra_list <- lapply(seq_len(nrow(sf)), function(i) {
    r <- unlist(sf[i, ]); names(sort(r, decreasing = TRUE))[seq_len(opt$top_n)]
  })
  names(spectra_list) <- rownames(sf)
  jaccard <- function(a, b) { a <- a[!is.na(a)]; b <- b[!is.na(b)]
    length(intersect(a, b)) / length(union(a, b)) }
  J <- outer(seq_along(mp_genes), seq_along(spectra_list),
             Vectorize(function(i, j) jaccard(mp_genes[[i]], spectra_list[[j]])))
  rownames(J) <- names(mp_genes); colnames(J) <- names(spectra_list)
  write.csv(round(J, 4), file.path(out_dir, sprintf("%s_jaccard_MP_vs_Spectra.csv", lineage)))
  best_j <- apply(J, 1, max); best_hit <- colnames(J)[apply(J, 1, which.max)]
  summary_df <- data.frame(lineage = lineage, meta_program = names(mp_genes),
                           best_spectra_factor = best_hit, best_jaccard = round(best_j, 4))
  write.csv(summary_df, file.path(out_dir, sprintf("%s_MP_vs_Spectra_summary.csv", lineage)),
            row.names = FALSE)
  log_msg("Jaccard comparison written (%d MPs vs %d Spectra factors).",
          length(mp_genes), length(spectra_list))
} else {
  log_msg("No Spectra factors file yet — skipping comparison (run after Spectra finishes).")
}

log_msg("Done. Outputs in %s", out_dir)
