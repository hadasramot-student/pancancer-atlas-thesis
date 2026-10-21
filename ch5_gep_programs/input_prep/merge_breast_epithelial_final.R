#!/usr/bin/env Rscript
# merge_breast_epithelial_final.R
#
# Completes the Breast Epithelial merge from existing chunks.
# Chunks already exist in tmp_merge/Epithelial/ — skips Pass 1.
# Uses gene INTERSECTION across chunks to avoid sparse matrix overflow.
# Then runs normalization → PCA → Harmony → clustering → UMAP.

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(Matrix)
  source("/ofir_grp/users/hadasra/create_sc_object_workflow/functions.R", chdir = FALSE)
})

CHUNK_DIR  <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/tmp_merge/Epithelial"
OUTPUT_RDS <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/Epithelial_Breast_merged.RDS"

chunk_files <- sort(list.files(CHUNK_DIR, pattern = "^chunk_.*\\.RDS$", full.names = TRUE))
log_msg(sprintf("Found %d chunk files", length(chunk_files)))

# ── Step 1: find common genes across all chunks ────────────────────────────
log_msg("Pass 0: scanning gene names across all chunks...")
gene_lists <- lapply(chunk_files, function(f) {
  obj <- readRDS(f)
  g <- rownames(obj)
  rm(obj); gc()
  g
})
common_genes <- Reduce(intersect, gene_lists)
log_msg(sprintf("Common genes: %d", length(common_genes)))
rm(gene_lists); gc()

# ── Step 2: load chunks, subset to common genes, collect as list ───────────
log_msg("Pass 1: loading chunks and subsetting to common genes...")
obj_list <- list()
for (i in seq_along(chunk_files)) {
  f <- chunk_files[[i]]
  log_msg(sprintf("  Loading chunk %d/%d: %s", i, length(chunk_files), basename(f)))
  obj <- readRDS(f)
  # ensure only counts layer remains (DietSeurat)
  obj <- DietSeurat(obj, assays = "RNA", dimreducs = NULL, graphs = NULL, misc = FALSE)
  obj <- obj[common_genes, ]
  log_msg(sprintf("  Chunk %d: %d cells x %d genes", i, ncol(obj), nrow(obj)))
  obj_list[[i]] <- obj
  rm(obj); gc()
}

log_msg(sprintf("Merging %d objects...", length(obj_list)))
merged <- merge(obj_list[[1]], y = obj_list[-1], merge.data = FALSE)
rm(obj_list); gc()
log_msg(sprintf("Final merged object: %d cells x %d genes", ncol(merged), nrow(merged)))

# ── Step 3: normalization → HVG → PCA → Harmony → clustering → UMAP ───────
log_msg("Joining layers...")
merged[["RNA"]] <- JoinLayers(merged[["RNA"]], layers = "counts")
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$project_name)

batch_sizes <- table(merged$project_name)
log_msg(sprintf("Projects: %d", length(batch_sizes)))
print(batch_sizes)

# drop tiny batches
small <- names(batch_sizes[batch_sizes < 30])
if (length(small) > 0) {
  log_msg(sprintf("Dropping small batches: %s", paste(small, collapse = ", ")))
  merged <- subset(merged, subset = !(project_name %in% small))
}

single_batch <- length(unique(merged$project_name)) < 2
max_dims <- min(20L, floor(min(table(merged$project_name)) / 2L))
max_dims <- max(max_dims, 5L)
log_msg(sprintf("max_dims = %d, single_batch = %s", max_dims, single_batch))

log_msg("NormalizeData...")
merged <- NormalizeData(merged)
log_msg("FindVariableFeatures...")
merged <- FindVariableFeatures(merged, nfeatures = 3000)
log_msg("ScaleData...")
merged <- ScaleData(merged)
log_msg("RunPCA...")
merged <- RunPCA(merged, npcs = max_dims)

if (!single_batch) {
  log_msg("Harmony integration...")
  merged <- IntegrateLayers(
    object         = merged,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    assay          = "RNA",
    dims           = 1:max_dims,
    batch          = "project_name",
    verbose        = TRUE
  )
  reduction_to_use <- "harmony"
} else {
  reduction_to_use <- "pca"
}

log_msg("FindNeighbors...")
n_cells <- ncol(merged)
k_param <- min(20L, max(5L, n_cells - 1L))
merged <- FindNeighbors(merged, reduction = reduction_to_use, dims = 1:max_dims,
                        graph.name = "Harmony_nn", k.param = k_param, nn.method = "rann")

log_msg("FindClusters...")
merged <- FindClusters(merged, graph.name = "Harmony_nn", resolution = 0.6,
                       cluster.name = "Harmony_cluster")

log_msg("RunUMAP...")
merged <- RunUMAP(merged, reduction = reduction_to_use, dims = 1:max_dims,
                  reduction.name = "umap_harmony")

log_msg(sprintf("Saving to: %s", OUTPUT_RDS))
saveRDS(merged, OUTPUT_RDS)
log_msg("Done!")
