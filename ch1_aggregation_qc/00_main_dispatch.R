# 00_main_dispatch.R
#
# Ch.1->Ch.2 orchestrator — takes one raw object from any of the 3 sources
# (Breastscape / 3CA / ScBaseCount) through the full per-object pipeline:
#   ingest + QC (Ch.1)
#     -> SingleR + SCType annotation (Ch.2, 04_object_annotation.R)
#     -> Scimilarity annotation (Ch.2, my_scimilarity.py)
#     -> original-metadata cross-validation + consensus (Ch.2, 05_object_annotion.R)
#
# Each stage is shelled out to as its own Rscript/python call so it can be
# resubmitted independently on the cluster; intermediate files are named
# <project_name>_step{1,2,3,4}.{RDS,h5ad}.
#
# Inputs (one of, depending on --ptype):
#   --h5ad_path  : Breastscape .h5ad
#   --rds        : ScBaseCount .RDS
#   --table_3CA + --row_number_3CA : 3CA project index + row to ingest
# Outputs:
#   <temp_objects_dir>/<project_name>_step4.RDS — fully annotated object
#
# Usage (run from the repo root):
#   Rscript ch1_aggregation_qc/00_main_dispatch.R --ptype breastscape --h5ad_path <path>
#   Rscript ch1_aggregation_qc/00_main_dispatch.R --ptype 3CA --row_number_3CA 10
#   Rscript ch1_aggregation_qc/00_main_dispatch.R --ptype scbasecount --rds <path>

library(optparse)

option_list <- list(

  make_option("--ptype",
              type = "character",
              default = "scbasecount",
              help = "breastscape OR scbasecount OR 3CA",
              metavar = "character"),

  make_option("--h5ad_path",
              type = "character",
              default = NULL,
              help = "Path to the h5ad object (for Breastscape) [default %default]",
              metavar = "path"),

  make_option("--rds",
              type = "character",
              default = NULL,
              help = "Path to the RDS object (for ScBaseCount) [default %default]",
              metavar = "path"),

  make_option("--table_3CA",
              type = "character",
              default = Sys.getenv("SC3CA_PROJECT_INDEX", unset = "config/SC.3CA.Tirosh.project_file_index.tsv"),
              help = "The path to the table of all the 3CA objects [default %default]",
              metavar = "path"),

  make_option("--row_number_3CA",
              type = "integer",
              default = 1,
              help = "The number of the row in the 3CA table [default %default]",
              metavar = "number"),

  make_option("--temp_objects_dir",
              type = "character",
              default = "data/temp_objects/",
              help = "Path to the directory to save the temp objects [default %default]",
              metavar = "path")
)

opt_parser <- optparse::OptionParser(usage = "usage: %prog [options]",
                                    description = "Ingest one raw object (any source) through QC + annotation + consensus.",
                                    option_list = option_list)
opt <- optparse::parse_args(opt_parser)


suppressPackageStartupMessages({
  library(tidyverse)
  library(logger)
})

source("common/functions.R")

temp_objects_dir <- opt$temp_objects_dir
ptype <- opt$ptype

if (opt$ptype == "breastscape") {
  stopifnot("file h5ad_path does not exist" = file.exists(opt$h5ad_path))
  output1 <- basename(file_no_ext <- tools::file_path_sans_ext(basename(opt$h5ad_path)))
  output1 <- paste0(output1, "_step1.RDS")
  CMD1 <- paste0("Rscript ch1_aggregation_qc/03_ingest_breastscape.R --h5ad_path ", opt$h5ad_path, " --output ", output1)
  project_name <- get_project_name(opt$h5ad_path)

  log_info(paste0("Running ", CMD1))
  system(CMD1)
} else if (opt$ptype == "3CA") {
  projects_df <- read_tsv(opt$table_3CA, comment = "#")

  project_details <- extract_3CA_data_from_row(projects_df = projects_df,
                                               rownum = opt$row_number_3CA)

  project_name <- project_details$project_name
  output1 <- paste0(project_details$project_name, "_step1.RDS")

  CMD1 <- glue::glue("
  Rscript ch1_aggregation_qc/02_ingest_3ca.R \\
    --mtx_path {project_details$mtx_path} \\
    --cells_path {project_details$cells_path} \\
    --project_name {project_details$project_name} \\
    --cancer_type {project_details$cancer_type} \\
   --genes_path {project_details$genes_path} \\
   --output  {output1}")


  log_info(paste0("Running ", CMD1))
  system(CMD1)

} else if (opt$ptype == "scbasecount") {
  stopifnot("file h5ad_path does not exist" = file.exists(opt$rds))
  project_name <- get_project_name(opt$rds)
  output1 <- paste0(project_name, ".RDS")
  CMD1 <- paste0("cp  ", opt$rds, " ", temp_objects_dir, project_name, ".RDS")
  log_info(paste0("Running ", CMD1))
  system(CMD1)
} else {
  stop("invalid --ptype value")
}

output2 <- paste0(project_name, "_step2.h5ad")
CMD2 <- glue::glue("
  Rscript ch2_annotation/04_object_annotation.R \\
    --project_name {project_name} \\
    --ptype {ptype} \\
   --input {output1} \\
   --output {output2}")
log_info(paste0("Running: ", CMD2))
system(CMD2)

output3 <- paste0(project_name, "_step3.h5ad")
py_script <- file.path(Sys.getenv("SC_OBJECT_WORKFLOW", unset = "ch2_annotation"), "my_scimilarity.py")
CMD3 <- paste(
  "conda run -n scimilarity",
  "python ", py_script,
  " --input", shQuote(output2),
  " --project-name", shQuote(project_name),
  " --output", shQuote(output3)
)

log_info(paste0("Running ", CMD3))
if (system(CMD3) != 0) {
  stop("Error in CMD3")
}

table_3CA <- opt$table_3CA
row_number_3CA <- opt$row_number_3CA

output4 <- paste0(project_name, "_step4.RDS")
CMD4 <- glue::glue("
  Rscript ch2_annotation/05_consensus_annotation.R \\
    --project_name {project_name} \\
    --ptype {ptype} \\
    --table_3CA {table_3CA} \\
    --row_number_3CA {row_number_3CA} \\
   --input {output3} \\
   --output {output4}
                   ")


log_info(paste0("Running ", CMD4))
system(CMD4)
