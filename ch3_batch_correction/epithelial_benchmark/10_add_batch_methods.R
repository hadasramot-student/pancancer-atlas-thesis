#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Epithelial batch-correction benchmark (Option A) - per cancer type.
#
# Takes the existing per-cancer scANVI object
#   Epithelial_<Cancer>_merged_full_scanvi_final.RDS
# (which already carries the atlas `scanvi` embedding, and for most
#  cancers also `harmony`), and ADDS the two missing method embeddings
#   - mnn        (batchelor::fastMNN)
#   - pca.combat (sva::ComBat on scaled HVG -> PCA)
# rebuilding harmony / normalization / scaling only when absent
# (Breast-450k-subsample and Kidney objects ship with scanvi only).
#
# Output object carries reductions: pca, harmony, scanvi, mnn, pca.combat
# plus meta cols project_name + scimilarity_celltype, ready for
# Batch_Correction_Validation_Metrics.R (2000-cell / seed-42 F1 etc.).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(SeuratWrappers)
  library(batchelor)
  library(SingleCellExperiment)
  library(sva)
  library(harmony)
})

option_list <- list(
  make_option("--input_rds",  type = "character", default = NULL),
  make_option("--output_rds", type = "character", default = NULL),
  make_option("--batch_col",  type = "character", default = "project_name"),
  make_option("--label_col",  type = "character", default = "scimilarity_celltype"),
  make_option("--n_hvg",      type = "integer",   default = 3000L),
  make_option("--n_pcs",      type = "integer",   default = 30L),
  make_option("--min_batch",  type = "integer",   default = 20L,
              help = "Batches with fewer cells than this are dropped from fastMNN [default 20]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_rds) || is.null(opt$output_rds)) stop("need --input_rds and --output_rds")

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time()), paste0(...)))
t0 <- Sys.time()

log_msg("Loading ", opt$input_rds)
obj <- readRDS(opt$input_rds)
DefaultAssay(obj) <- "RNA"
log_msg(sprintf("Loaded: %d cells x %d genes", ncol(obj), nrow(obj)))

## --- 0. join split layers -------------------------------------------------
obj[["RNA"]] <- tryCatch(JoinLayers(obj[["RNA"]]), error = function(e) { log_msg("JoinLayers skipped: ", conditionMessage(e)); obj[["RNA"]] })
layers_now <- SeuratObject::Layers(obj[["RNA"]])
log_msg("RNA layers after join: ", paste(layers_now, collapse = ", "))

batch <- as.character(obj[[opt$batch_col]][, 1])
stopifnot(!anyNA(batch))
log_msg(sprintf("batch '%s': %d levels; label '%s': %d NA",
                opt$batch_col, length(unique(batch)),
                opt$label_col, sum(is.na(obj[[opt$label_col]][, 1]))))

## --- 1. ensure data (log-norm) layer ------------------------------------
if (!("data" %in% SeuratObject::Layers(obj[["RNA"]]))) {
  log_msg("No 'data' layer -> NormalizeData()")
  obj <- NormalizeData(obj, verbose = FALSE)
}

## --- 2. ensure HVG + scale.data ---------------------------------------
sd_dim <- tryCatch(dim(GetAssayData(obj, assay = "RNA", layer = "scale.data")),
                   error = function(e) c(0, 0))
if (is.null(sd_dim) || sd_dim[1] == 0) {
  log_msg(sprintf("No scale.data -> FindVariableFeatures(%d) + ScaleData()", opt$n_hvg))
  obj <- FindVariableFeatures(obj, nfeatures = opt$n_hvg, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
} else {
  log_msg(sprintf("Reusing existing scale.data (%d x %d)", sd_dim[1], sd_dim[2]))
  if (length(VariableFeatures(obj[["RNA"]])) == 0)
    VariableFeatures(obj[["RNA"]]) <- rownames(GetAssayData(obj, assay = "RNA", layer = "scale.data"))
}

## --- 3. ensure pca ------------------------------------------------------
if (!("pca" %in% Reductions(obj))) {
  log_msg("No pca -> RunPCA()")
  obj <- RunPCA(obj, npcs = max(opt$n_pcs, 30L), verbose = FALSE)
}

## --- 4. ensure harmony ------------------------------------------------
if (!("harmony" %in% Reductions(obj))) {
  log_msg("No harmony -> RunHarmony(group.by.vars = ", opt$batch_col, ")")
  obj <- RunHarmony(obj, group.by.vars = opt$batch_col,
                    reduction = "pca", reduction.save = "harmony", verbose = TRUE)
} else {
  log_msg(sprintf("Reusing existing harmony (%d dims)", ncol(Embeddings(obj, "harmony"))))
}

## --- 5. fastMNN -> 'mnn' --------------------------------------------------
log_msg("fastMNN ...")
tb <- table(batch)
keep_batches <- names(tb)[tb >= opt$min_batch]
dropped <- setdiff(names(tb), keep_batches)
if (length(dropped))
  log_msg(sprintf("  dropping %d batch(es) < %d cells from MNN: %s",
                  length(dropped), opt$min_batch, paste(dropped, collapse = ", ")))
cells_mnn <- colnames(obj)[batch %in% keep_batches]

obj_mnn <- obj[, cells_mnn]
object.list <- SplitObject(obj_mnn, split.by = opt$batch_col)
common <- Reduce(intersect, lapply(object.list, rownames))
object.list <- lapply(object.list, function(x) x[common, ])
sce.list <- lapply(object.list, as.SingleCellExperiment)
min_cells <- min(vapply(sce.list, ncol, integer(1)))
k_mnn <- max(5L, min(20L, min_cells - 1L))
d_mnn <- min(opt$n_pcs, min_cells - 1L, 30L)
log_msg(sprintf("  fastMNN k=%d d=%d on %d batches / %d cells", k_mnn, d_mnn,
                length(sce.list), ncol(obj_mnn)))
mnn_out <- batchelor::fastMNN(sce.list, k = k_mnn, d = d_mnn)
mnn_mat <- reducedDim(mnn_out, "corrected")

full_mnn <- matrix(NA_real_, nrow = ncol(obj), ncol = ncol(mnn_mat),
                   dimnames = list(colnames(obj), paste0("mnn_", seq_len(ncol(mnn_mat)))))
full_mnn[rownames(mnn_mat), ] <- as.matrix(mnn_mat)
obj[["mnn"]] <- CreateDimReducObject(embeddings = full_mnn, key = "mnn_",
                                     assay = DefaultAssay(obj))
rm(obj_mnn, object.list, sce.list, mnn_out, mnn_mat, full_mnn); gc()

## --- 6. ComBat -> 'pca.combat' ----------------------------------------
log_msg("ComBat ...")
mat <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "scale.data"))
log_msg(sprintf("  scale.data for ComBat: %d genes x %d cells", nrow(mat), ncol(mat)))
combat_mat <- sva::ComBat(dat = mat, batch = batch)
obj[["RNA.combat"]] <- CreateAssayObject(data = combat_mat)
rm(mat, combat_mat); gc()
DefaultAssay(obj) <- "RNA.combat"
VariableFeatures(obj) <- intersect(rownames(obj[["RNA.combat"]]),
                                   VariableFeatures(obj[["RNA"]]))
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = max(opt$n_pcs, 30L), reduction.name = "pca.combat",
              verbose = FALSE)
DefaultAssay(obj) <- "RNA"

## --- 7. slim + save --------------------------------------------------
obj[["RNA.combat"]] <- NULL
keep_red <- intersect(c("pca", "harmony", "scanvi", "mnn", "pca.combat"), Reductions(obj))
for (r in setdiff(Reductions(obj), keep_red)) obj[[r]] <- NULL
log_msg("Final reductions: ", paste(Reductions(obj), collapse = ", "))

dir.create(dirname(opt$output_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(obj, opt$output_rds)
log_msg(sprintf("Saved %s  (%.1f min elapsed)", opt$output_rds,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
