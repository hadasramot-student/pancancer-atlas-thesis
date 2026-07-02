# 03_ingest_breastscape.R
#
# Ch.1 (Aggregation & QC) — ingest one Breastscape (Technion) project from
# its pre-processed .h5ad into a Seurat object.
#
# Inputs:
#   --h5ad_path        : path to one Breastscape project's .h5ad file
#   --temp_objects_dir : output directory
#   --output            : output RDS filename
# Outputs:
#   <temp_objects_dir>/<output> — Seurat object with project_name set
#
# Usage:
#   Rscript 03_ingest_breastscape.R \
#     --h5ad_path path/to/breastscape/DING_K_ET_AL_2022_P23_filtered.h5ad \
#     --temp_objects_dir data/temp_objects --output DING_K_ET_AL_2022_P23.RDS

library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(readr)
library(optparse)

option_list <- list(
  make_option("--h5ad_path", type = "character", default = NULL,
              help = "Path to a Breastscape project's .h5ad file [required]"),
  make_option("--temp_objects_dir", type = "character", default = "data/temp_objects/",
              help = "Output directory [default %default]"),
  make_option("--output", type = "character", default = NULL,
              help = "Output RDS filename [required]")
)

opt_parser <- optparse::OptionParser(usage = "usage: %prog [options]",
                                      description = "Ingest one Breastscape project into a Seurat object.",
                                      option_list = option_list)
opt <- optparse::parse_args(opt_parser)


suppressPackageStartupMessages({
  library(Seurat)
  library(zellkonverter)
})

source("common/functions.R")
obj_path <- opt$h5ad_path

save_path <- paste0(opt$temp_objects_dir, "/", opt$output)

print(obj_path)
### Read the .h5ad object:
sce <- readH5AD(obj_path)       # returns SingleCellExperiment
### Make it Seurat object, and add the columns:
### project_name, cancer_type = "Breast" :
expr_mat <- assay(sce, "X")
meta_df  <- as.data.frame(colData(sce))
seurat_object <- CreateSeuratObject(
  counts = expr_mat,
  meta.data = meta_df
)
project_name <- get_project_name(obj_path)
seurat_object@project.name <- project_name

saveRDS(seurat_object, file = save_path)
