# 05_consensus_annotation.R
#
# Ch.2 (Annotation) — final annotation stage for one object, run after
# my_scimilarity.py. Builds the Scimilarity meta-type mapping, majority-vote
# aggregation for SingleR/SCType/Scimilarity, sample/patient ID assignment,
# and the final maj2/maj3/maj6 consensus vote (Consensus_Cell_Type).
#
# Requires an Rmd report template (see --rmd_template) for the per-object
# QC/annotation HTML summary — optional, only used for the reporting step.
#
# Inputs:
#   --input  : h5ad written by my_scimilarity.py (Ch.2 step 3)
#   --ontology : meta-type ontology mapping (Ontology_4_methods.xlsx)
#   --table_3CA + --row_number_3CA : only for ptype == "3CA", to recover the
#     original Cells.csv path for sample/patient ID assignment
# Outputs:
#   <temp_objects_dir>/<output>       - final annotated Seurat RDS
#   <temp_objects_dir>/<output>.tsv   - flat metadata dump
#   <temp_objects_dir>/<output>.html  - per-object QC/annotation report (if
#                                        --rmd_template is available)
#
# Usage:
#   Rscript 05_consensus_annotation.R \
#     --project_name Azizi2018_Breast --ptype 3CA \
#     --table_3CA config/SC.3CA.Tirosh.project_file_index.tsv --row_number_3CA 26 \
#     --input Azizi2018_Breast_step3.h5ad --output Azizi2018_Breast_step4.RDS

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(zellkonverter)
  library(SummarizedExperiment)
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

  make_option("--ontology",
              type = "character",
              default = "ch2_annotation/Ontology_4_methods.xlsx",
              help = "Meta-type ontology mapping (SingleR/SCType/Scimilarity -> meta-type) [default %default]",
              metavar = "path"),

  make_option("--temp_objects_dir",
              type = "character",
              default = "data/temp_objects/",
              help = "[default %default]",
              metavar = "path"),

  make_option("--table_3CA",
              type = "character",
              default = Sys.getenv("SC3CA_PROJECT_INDEX", unset = "config/SC.3CA.Tirosh.project_file_index.tsv"),
              help = "The path to the table of all the 3CA objects [default %default]",
              metavar = "path"),

  make_option("--row_number_3CA",
              type = "integer",
              default = NULL,
              help = "The number of the row in the 3CA table [default %default]",
              metavar = "number"),

  make_option("--rmd_template",
              type = "character",
              default = "common/Rmd_report.Rmd",
              help = "Per-object QC/annotation report template; skipped if missing [default %default]",
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
                                    description = "Final consensus annotation stage (Ch.2).",
                                    option_list = option_list)
opt <- optparse::parse_args(opt_parser)

source("common/functions.R")


project_name <- opt$project_name
h5ad_after_scimilarity <- paste0(opt$temp_objects_dir, opt$input)
sce <- readH5AD(h5ad_after_scimilarity)
expr_mat <- assay(sce, "X")
meta_df  <- as.data.frame(colData(sce))
seurat_object <- CreateSeuratObject(
  counts = expr_mat,
  project = project_name,
  meta.data = meta_df
)
seurat_object <- creating_object(seurat_object)
seurat_object@project.name <- project_name
print("finished create seurat from h5ad")

seurat_object <- into_Meta(seurat_object, opt$ontology)

seurat_object <- add_meta_from_scimilarity(seurat_object, opt$ontology)
print("finished add_meta_from_scimilarity")

if (opt$ptype == "3CA") {
  seurat_object <- add_precentage_and_project_name(seurat_object)
  print("finished add_precentage_and_project_name")
}
seurat_object <- annotate_majority_singler(seurat_object)
print("finished annotate_majority_singler")

seurat_object <- add_precentage_sctype(seurat_object)
print("finished add_precentage_sctype")

seurat_object <- add_precentage_scimilarity(seurat_object)
print("finished add_precentage_scimilarity")

if (opt$ptype == "3CA") {
  projects_df <- read_tsv(opt$table_3CA, comment = "#")
  project_details <- extract_3CA_data_from_row(projects_df = projects_df,
                                               rownum = opt$row_number_3CA)
  cells_path <- project_details$cells_path
  seurat_object <- add_sample_for_project(seurat_object, cells_path)
} else {
  seurat_object <- add_sample_for_project_h5ad(seurat_object)
}
print("finished add sample")

seurat_object <- add_consensuses(seurat_object, opt$ontology)
print("finished add_consensuses")

if (file.exists(opt$rmd_template)) {
  save_rmd_path <- paste0(opt$output, ".html")
  rmarkdown::render(input = opt$rmd_template,
                    output_file = save_rmd_path,
                    output_dir = opt$temp_objects_dir,
                    params = list(seurat_obj = seurat_object,
                                 project_name = unique(opt$project_name),
                                 data_source = opt$ptype,
                                 author = Sys.getenv("USER")))
} else {
  message("rmd_template not found (", opt$rmd_template, "), skipping HTML report")
}

save_path <- paste0(opt$temp_objects_dir, opt$output)
saveRDS(seurat_object, file = save_path)

save_tsv_path <- paste0(opt$temp_objects_dir, opt$output, ".tsv")
save_metadata_tsv(seurat_object, save_tsv_path)
print("finished saved")

file.remove(h5ad_after_scimilarity)
print(paste0("finished delete: ", h5ad_after_scimilarity))
