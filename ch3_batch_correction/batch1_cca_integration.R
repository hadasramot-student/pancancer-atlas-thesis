#!/usr/bin/env Rscript
# batch1_cca_integration.R
#
# Ch.3 (Batch correction) — Script 1: Seurat CCA integration.
# Inputs:  --input_rds  (merged, un-integrated Seurat object)
# Outputs: --output_rds (adds pca.cca/umap.cca reductions + cca_clusters)
#
# Usage:
#   Rscript batch1_cca_integration.R --input_rds obj_merged.RDS \
#     --output_rds obj_merged_cca.RDS --batch_col project_name

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(SeuratWrappers)
})

# --------------------------------------------------
# Arguments
# --------------------------------------------------
option_list <- list(
  make_option(c("--input_rds"), type="character", default=NULL),
  make_option(c("--output_rds"), type="character", default=NULL),
  make_option(c("--batch_col"), type="character", default="project_name"),
  make_option(c("--subsample"), type="logical", default=FALSE,
              help="Subsample 90%% of cells for development [default=FALSE]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_rds)) stop("Please provide --input_rds")

if (is.null(opt$output_rds)) {
  opt$output_rds <- sub("\\.RDS$", "_cca.RDS", opt$input_rds)
}

batch_col <- opt$batch_col

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

# --------------------------------------------------
# Load object
# --------------------------------------------------
log_msg(paste("Loading object:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Loaded object with %d cells and %d genes", ncol(obj), nrow(obj)))

# --------------------------------------------------
# Optional subsample for development
# --------------------------------------------------
if (opt$subsample) {
  set.seed(42)
  n_keep <- floor(0.9 * ncol(obj))
  cells_keep <- sample(colnames(obj), n_keep)
  obj <- obj[, cells_keep]
  log_msg(sprintf("Subsampled to %d cells (90%% of original)", ncol(obj)))
}

# --------------------------------------------------
# CCA integration
# --------------------------------------------------
log_msg("Starting CCA integration")

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)

object.list <- SplitObject(obj, split.by = batch_col)
cell_counts <- sapply(object.list, ncol)
log_msg("Cell counts per batch:")
print(cell_counts[order(cell_counts)])

min_cells <- min(cell_counts)

if (min_cells < 20L) {
  log_msg("Minimum cells per batch < 20; skipping CCA integration.")
} else {
  
  max_dims_cca <- min(20L, min_cells - 1L)
  k.weight     <- max(5L,  min(50L, min_cells - 5L))
  k.filter     <- max(10L, min(50L, min_cells - 1L))
  k.score      <- max(5L,  min(30L, floor(min_cells / 2L)))
  
  log_msg(sprintf(
    "CCA parameters: max_dims_cca=%d, k.weight=%d, k.filter=%d, k.score=%d",
    max_dims_cca, k.weight, k.filter, k.score
  ))
  
  object.list <- lapply(object.list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE)
    x
  })
  
  features <- SelectIntegrationFeatures(object.list, nfeatures = 2000)
  
  log_msg("Finding CCA anchors")
  anchors <- tryCatch({
    FindIntegrationAnchors(
      object.list     = object.list,
      anchor.features = features,
      reduction       = "cca",
      dims            = 1:max_dims_cca,
      k.filter        = k.filter,
      k.anchor        = 5L,
      k.score         = k.score,
      verbose         = TRUE
    )
  }, error = function(e) {
    log_msg(paste("FindIntegrationAnchors (CCA) failed:", e$message))
    NULL
  })
  
  if (is.null(anchors)) {
    log_msg("Anchors NULL; skipping CCA.")
  } else {
    n_anchors <- nrow(anchors@anchors)
    log_msg(sprintf("Total anchors after filtering: %d", n_anchors))
    
    if (!is.na(n_anchors) && n_anchors > 0L) {
      
      log_msg("Running IntegrateLayers (CCA)")
      obj <- tryCatch({
        IntegrateLayers(
          object         = obj,
          method         = CCAIntegration,
          anchors        = anchors,
          orig.reduction = "pca",
          new.reduction  = "integrated.cca",
          assay          = "RNA",
          dims           = 1:max_dims_cca,
          k.weight       = k.weight,
          verbose        = TRUE
        )
      }, error = function(e) {
        log_msg(paste("IntegrateLayers (CCA) failed:", e$message))
        obj
      })
      
      if ("integrated.cca" %in% Reductions(obj)) {
        cca_emb    <- Embeddings(obj, "integrated.cca")
        n_cca_dims <- ncol(cca_emb)
        log_msg(sprintf("integrated.cca has %d dimensions", n_cca_dims))
        
        if (n_cca_dims >= 2L) {
          used_dims     <- 1:min(max_dims_cca, n_cca_dims)
          n_cells_total <- ncol(obj)
          k_param       <- min(20L, max(5L, n_cells_total - 1L))
          
          obj <- FindNeighbors(obj, reduction = "integrated.cca",
                               dims = used_dims, k.param = k_param,
                               nn.method = "rann")
          obj <- FindClusters(obj, resolution = 0.6,
                              cluster.name = "cca_clusters")
          obj <- RunUMAP(obj, reduction = "integrated.cca",
                         dims = used_dims, reduction.name = "umap.cca")
          
          log_msg("Finished CCA integration")
        }
      }
    }
  }
}

# --------------------------------------------------
# Save checkpoint
# --------------------------------------------------
log_msg(paste("Saving CCA checkpoint:", opt$output_rds))
saveRDS(obj, opt$output_rds)
log_msg("Script 1 (CCA) done")

# sbatch rstudio.sbatch --script  merge/batch1_CCA_integration.R \
# --input_rds $ATLAS_DATA_DIR/merged_by_query/Myeloid_Breast_merged.RDS \
# --output_rds $ATLAS_DATA_DIR/merged_by_query/Myeloid_Breast_merged_cca.RDS \
# --batch_col project_name
