# 02_ingest_3ca.R
#
# Ch.1 (Aggregation & QC) — ingest one 3CA (Tirosh lab "Curated Cancer Cell
# Atlas") project from its raw 10x-style mtx/cells/genes triplet into a QC'd
# Seurat object, and cross-links the original per-cell "cell_type" label
# from the 3CA metadata (used later for the original-metadata cross-
# validation step in Ch.2 — 3CA is the only one of the 3 sources whose raw
# metadata carries a usable original cell-type column).
#
# QC: min.cells = 3, min.features = 200 (Seurat::CreateSeuratObject).
#
# Requires: a lab-internal `Utils.R` providing
# `my.str_split_fixed.and.paste()` (string helper for de-duplicating cell
# barcodes) — point DROPBOX_ROOT at wherever your copy lives, or replace
# that one call with an equivalent barcode-dedup helper.
#
# Inputs:
#   --mtx_path, --cells_path, --genes_path : one project's raw 3CA files
#   --project_name, --cancer_type          : labels for this project
#   --temp_objects_dir, --output           : where to write the QC'd RDS
# Outputs:
#   <temp_objects_dir>/<output> — QC'd Seurat object with cancer_type,
#   cell_type (original 3CA label) in meta.data
#
# Usage:
#   Rscript 02_ingest_3ca.R \
#     --mtx_path Exp_data_UMIcounts.mtx --cells_path Cells.csv --genes_path Genes.txt \
#     --project_name Azizi2018_Breast --cancer_type Breast \
#     --temp_objects_dir data/temp_objects/ --output Azizi2018_Breast.RDS

suppressPackageStartupMessages({
  library(Seurat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(readr)
  library(optparse)
})


option_list <- list(
  make_option("--mtx_path",        type = "character", default = NULL,
              help = "Path to the 3CA project's UMI-count .mtx file [required]"),
  make_option("--cells_path",      type = "character", default = NULL,
              help = "Path to the 3CA project's Cells.csv (barcodes + original annotation) [required]"),
  make_option("--project_name",    type = "character", default = NULL,
              help = "Project identifier, e.g. Azizi2018_Breast [required]"),
  make_option("--cancer_type",     type = "character", default = NULL,
              help = "Cancer type label for this project [required]"),
  make_option("--genes_path",      type = "character", default = NULL,
              help = "Path to the 3CA project's Genes.txt [required]"),
  make_option("--output",          type = "character", default = NULL,
              help = "Output RDS filename [required]"),
  make_option("--temp_objects_dir", type = "character", default = "data/temp_objects/",
              help = "Output directory [default %default]"),
  make_option("--failed_log",      type = "character",
              default = "data/temp_objects/3ca_failed_projects.tsv",
              help = "Log of skipped/failed projects [default %default]")
)

opt_parser <- optparse::OptionParser(usage = "usage: %prog [options]",
                                      description = "Ingest one 3CA project into a QC'd Seurat object.",
                                      option_list = option_list)
opt <- optparse::parse_args(opt_parser)

source(file.path(Sys.getenv("DROPBOX_ROOT", unset = "config/lab_utils"), "Shared_Code/R_Generic_Code/Utils.R"))

# Example call (see Usage above for the CLI form):
# opt <- list()
# opt$mtx_path      <- "path/to/3CA/Breast/Data_Azizi2018_Breast/10x/Exp_data_UMIcounts.mtx"
# opt$cells_path    <- "path/to/3CA/Breast/Data_Azizi2018_Breast/10x/Cells.csv"
# opt$project_name  <- "Azizi2018_Breast"
# opt$cancer_type   <- "Breast"
# opt$genes_path    <- "path/to/3CA/Breast/Data_Azizi2018_Breast/10x/Genes.txt"

mtx_path     <- opt$mtx_path
cells_path   <- opt$cells_path
project_name <- opt$project_name
cancer_type  <- opt$cancer_type
genes_path   <- opt$genes_path

expression_matrix <- ReadMtx(
  mtx = mtx_path,
  cell.sep = ",",
  features = genes_path,
  cells = cells_path, feature.column = 1, skip.cell = 1
)

min.cells <- 3
min.features <- 200

seurat_object <- CreateSeuratObject(counts = expression_matrix,
                                  project = project_name,
                                  min.cells = min.cells,
                                  min.features = min.features)

cells_org_annot <- data.table::fread(cells_path)

log_path <- opt$failed_log
log_skip <- function(reason) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("%s\t%s\t%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              project_name,
              reason),
      file = log_path, append = TRUE)
  return(NULL)  # <-- skip this object cleanly (works in lapply/map/inside fn)
}

has_cell_type <- "cell_type" %in% names(cells_org_annot)
has_cell_name <- any(c("cell_name", "Cell_names") %in% names(cells_org_annot))

if (!has_cell_type) { log_skip("missing 'cell_type'"); quit(save = "no", status = 0) }
if (!has_cell_name) { log_skip("missing 'cell_name'/'Cell_names'"); quit(save = "no", status = 0) }

seurat_object@meta.data$cancer_type <- cancer_type
cells_org_annot <- cells_org_annot[cells_org_annot$cell_name %in% rownames(seurat_object@meta.data), ]
meta.data <- seurat_object@meta.data
meta.data <- meta.data[order(meta.data$nCount_RNA), ]
meta.data$cell_name <- my.str_split_fixed.and.paste(row.names(meta.data), ",", str_count(row.names(meta.data)[1], ",") + 1, 1)
meta.data.merged <- merge(meta.data, cells_org_annot, by = "cell_name", all = TRUE, sort = FALSE)
row.names(seurat_object@meta.data) <- meta.data.merged$cell_name
m.cell_type <- meta.data.merged$cell_type
names(m.cell_type) <- row.names(seurat_object@meta.data)
seurat_object <- AddMetaData(seurat_object, m.cell_type, "m.cell_type")
seurat_object@meta.data <- meta.data
cell_type <- meta.data.merged$cell_type
names(cell_type) <- row.names(seurat_object@meta.data)
seurat_object <- AddMetaData(seurat_object, cell_type, "cell_type")
Idents(seurat_object) <- "cell_type"

save_path <- paste0(opt$temp_objects_dir, opt$output)
saveRDS(seurat_object, file = save_path)
