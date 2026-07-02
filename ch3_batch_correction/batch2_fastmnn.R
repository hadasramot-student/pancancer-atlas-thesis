###_______________Script 2: fastMNN_______________


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(SeuratWrappers)
  library(batchelor)
  library(SingleCellExperiment)
})

option_list <- list(
  make_option(c("--input_rds"),  type="character", default=NULL),
  make_option(c("--output_rds"), type="character", default=NULL),
  make_option(c("--batch_col"),  type="character", default="project_name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_rds)) stop("Please provide --input_rds")

if (is.null(opt$output_rds)) {
  opt$output_rds <- sub("\\.RDS$", "_mnn.RDS", opt$input_rds)
}

batch_col <- opt$batch_col

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

log_msg(paste("Loading object:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Loaded object with %d cells and %d genes", ncol(obj), nrow(obj)))

# --------------------------------------------------
# fastMNN
# --------------------------------------------------
log_msg("Starting fastMNN")

object.list <- SplitObject(obj, split.by = batch_col)

feature_sets    <- lapply(object.list, rownames)
common_features <- Reduce(intersect, feature_sets)
object.list     <- lapply(object.list, function(x) x[common_features, ])

sce.list <- lapply(object.list, as.SingleCellExperiment)

cell_counts <- sapply(sce.list, ncol)
min_cells   <- min(cell_counts)
k_mnn       <- max(5L,  min(30L, min_cells - 1L))
ndims_mnn   <- min(20L, min_cells - 1L)

log_msg(sprintf("fastMNN parameters: k=%d, d=%d", k_mnn, ndims_mnn))

mnn_out <- batchelor::fastMNN(sce.list, k = k_mnn, d = ndims_mnn)

mnn_mat <- reducedDim(mnn_out, "corrected")
mnn_mat <- mnn_mat[colnames(obj), , drop = FALSE]

obj[["mnn"]] <- CreateDimReducObject(
  embeddings = as.matrix(mnn_mat),
  key        = "mnn_",
  assay      = DefaultAssay(obj)
)

n_mnn_pcs <- ncol(Embeddings(obj, "mnn"))
ndims_mnn <- min(ndims_mnn, n_mnn_pcs)

obj <- FindNeighbors(obj, reduction = "mnn", dims = 1:ndims_mnn)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "mnn_clusters")
obj <- RunUMAP(obj, reduction = "mnn", dims = 1:ndims_mnn,
               reduction.name = "umap.mnn")

log_msg("Finished fastMNN")

# --------------------------------------------------
# Save checkpoint
# --------------------------------------------------
log_msg(paste("Saving MNN checkpoint:", opt$output_rds))
saveRDS(obj, opt$output_rds)
log_msg("Script 2 (fastMNN) done")

# sbatch rstudio.sbatch --script  merge/batch2_fastMNN.R \
# --input_rds $ATLAS_DATA_DIR/merged_by_query/Myeloid_Breast_merged_cca.RDS \
# --output_rds $ATLAS_DATA_DIR/merged_by_query/Myeloid_Breast_merged_fastMNN.RDS \
# --batch_col project_name