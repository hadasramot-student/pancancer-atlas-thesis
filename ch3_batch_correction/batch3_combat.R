###_______________Script 3: ComBat_______________


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(sva)
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
  opt$output_rds <- sub("\\.RDS$", "_combat.RDS", opt$input_rds)
}

batch_col <- opt$batch_col

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

choose_max_dims <- function(obj, max_cap = 30L) {
  max(min(max_cap, floor(ncol(obj) / 2L)), 5L)
}

log_msg(paste("Loading object:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Loaded object with %d cells and %d genes", ncol(obj), nrow(obj)))

# --------------------------------------------------
# ComBat
# --------------------------------------------------
log_msg("Starting ComBat")

mat <- GetAssayData(object = obj, assay = "RNA", layer = "scale.data")
batch <- obj[[batch_col]][,1]

combat_mat <- ComBat(dat = as.matrix(mat), batch = batch)

obj[["RNA.combat"]] <- CreateAssayObject(data = combat_mat)
DefaultAssay(obj) <- "RNA.combat"

VariableFeatures(obj) <- intersect(
  rownames(obj[["RNA.combat"]]),
  VariableFeatures(obj[["RNA"]])
)

obj <- ScaleData(obj, verbose = FALSE)

max_dims <- choose_max_dims(obj, 50L)
obj <- RunPCA(obj, npcs = max_dims, reduction.name = "pca.combat")

use_dims <- min(30L, max_dims)
obj <- RunUMAP(obj, reduction = "pca.combat", dims = 1:use_dims,
               reduction.name = "umap.combat")
obj <- FindNeighbors(obj, reduction = "pca.combat", dims = 1:use_dims)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "combat_clusters")

log_msg("Finished ComBat")

# --------------------------------------------------
# Save final object
# --------------------------------------------------
log_msg(paste("Saving final object:", opt$output_rds))
saveRDS(obj, opt$output_rds)
log_msg("Script 3 (ComBat) done")

# sbatch rstudio.sbatch --script  merge/batch2_fastMNN.R \
# --input_rds $ATLAS_DATA_DIR/merged_by_query/Myeloid_Breast_merged_fastMNN.RDS \
# --output_rds $ATLAS_DATA_DIR/merged_by_query/Myeloid_Breast_merged_ComBat.RDS \
# --batch_col project_name