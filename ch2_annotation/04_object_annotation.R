# 04_object_annotation.R
#
# Ch.2 (Annotation) — first annotation stage for one object, run after Ch.1
# ingestion. Runs SingleR (reference: HumanPrimaryCellAtlasData), the 3CA
# original-metadata cross-validation (3CA only — the only source with a
# usable original per-cell label), and SCType, then hands off to
# my_scimilarity.py (Scimilarity) by writing an .h5ad.
#
# Inputs:
#   --input : ingested, QC'd Seurat RDS (Ch.1 step1 output)
#   --db_for_sctype : SCType marker database (ScTypeDB_full.xlsx)
#   --dedup_path : 3CA original-vs-deduplicated label lookup (3CA only)
#   --lookup_table : Ensembl ID -> gene symbol table (ScBaseCount only)
# Outputs:
#   <temp_objects_dir>/<output> — h5ad with SingleR.labels, sctype columns
#   (and cell_type cross-validated against 3CA original labels, if ptype=="3CA")
#
# Usage:
#   Rscript 04_object_annotation.R \
#     --project_name Azizi2018_Breast --ptype 3CA \
#     --input Azizi2018_Breast_step1.RDS --output Azizi2018_Breast_step2.h5ad

suppressPackageStartupMessages({
  library(Seurat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(readr)
  library(optparse)
})

option_list <- list(

  make_option("--project_name",
              type = "character",
              default = NULL,
              help = "[default %default]",
              metavar = "string"),

  make_option("--dedup_path",
              type = "character",
              default = "config/original_vs_dedup.xlsx",
              help = "3CA original-vs-deduplicated cell-type label lookup [default %default]",
              metavar = "path"),

  make_option("--db_for_sctype",
              type = "character",
              default = "ch2_annotation/ScTypeDB_full.xlsx",
              help = "SCType marker database [default %default]",
              metavar = "path"),

  make_option("--ptype",
              type = "character",
              default = NULL,
              help = "breastscape OR scbasecount OR 3CA [default %default]",
              metavar = "string"),

  make_option("--lookup_table",
              type = "character",
              default = "config/lookup_table_full.tsv",
              help = "Ensembl ID -> gene symbol lookup table (ScBaseCount only) [default %default]",
              metavar = "path"),

  make_option("--temp_objects_dir",
              type = "character",
              default = "data/temp_objects/",
              help = "[default %default]",
              metavar = "path"),

  make_option("--input",
              type = "character",
              default = NULL,
              help = "[default %default]",
              metavar = "path"),

  make_option("--output",
              type = "character",
              default = NULL,
              help = "[default %default]",
              metavar = "path")
)

opt_parser <- optparse::OptionParser(usage = "usage: %prog [options]",
                                    description = "SingleR + 3CA cross-validation + SCType annotation (Ch.2).",
                                    option_list = option_list)
opt <- optparse::parse_args(opt_parser)

source("common/functions.R")

project_path <- paste0(opt$temp_objects_dir, opt$input)

seurat_object <- readRDS(project_path)
project_name <- seurat_object@project.name
ptype <- opt$ptype

if (ptype == "scbasecount") {
  lookup_table <- read_tsv(opt$lookup_table)
  seurat_object <- rownames_as_gene_symbols(seurat_object, lookup_table = lookup_table)
  print("finished rownames_as_gene_symbols")
}
seurat_object <- creating_object(mSeurat.obj = seurat_object, project_name = project_name)
print("finished creating_object")

seurat_object <- singleR_to_metadata(seurat_object)
print("finished singleR_to_metadata")

if (ptype == "3CA") {
  dedup_path <- opt$dedup_path
  seurat_object <- original_de_dup(seurat_object, dedup_path)
}
db_for_sctype <- opt$db_for_sctype
tissue    <- "Immune_plus"
seurat_object <- sctype_cell_annotation(seurat_object, db_for_sctype, tissue)
print("finished sctype_cell_annotation")

# save as h5ad, handoff to my_scimilarity.py:
save_path <- paste0(opt$temp_objects_dir, opt$output)

DefaultAssay(seurat_object) <- "RNA"
sce <- as.SingleCellExperiment(seurat_object)
writeH5AD(sce, file = save_path, X_name = "counts")

# delete the previous-stage object:
file.remove(project_path)

print("finished and saved")
