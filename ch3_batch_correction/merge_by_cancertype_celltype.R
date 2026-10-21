#!/usr/bin/env Rscript
# print(.libPaths())
# quit()
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}
log_msg("Script started")

on.exit({
  log_msg("Script finished")
})

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
cat(".libPaths() =\n")
print(.libPaths())
cat("fs package path:\n")
print(find.package("fs", lib.loc = .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(SeuratObject)
  library(scales)
  library(harmony)
  library(SeuratDisk)
  library(SeuratWrappers)   # for some wrappers if needed
  library(batchelor)        # fastMNN, mnnCorrect
  library(sva)              # ComBat
  library(Matrix)

  source("/ofir_grp/users/hadasra/create_sc_object_workflow/functions.R")
  
})


#--------------------------------------------------
# Helper: streaming merge by meta table
#--------------------------------------------------
merge_subset_streaming <- function(meta_df,
                                   input_dir,
                                   output_file,
                                   project_col = "project_name",
                                   cellid_col  = "cell_id",
                                   file_suffix = "_step4.RDS",
                                   diet        = TRUE,
                                   keep_counts_only = TRUE,
                                   verbose     = TRUE) {
  meta_df <- meta_df %>%
    filter(!is.na(.data[[project_col]]), .data[[project_col]] != "",
           !is.na(.data[[cellid_col]]),  .data[[cellid_col]]  != "")
  
  proj_list <- split(as.character(meta_df[[cellid_col]]),
                     as.character(meta_df[[project_col]]))
  
  if (verbose) print(paste0("Projects to process: ", length(proj_list)))
  
  read_obj <- function(path) {
    tryCatch(
      readRDS(path),
      error = function(e) {
        print("[RDS READ ERROR] ", path, " | ", conditionMessage(e))
        NULL
      }
    )
  }
 
  
  # 1) Build merged object from RDS files (similar to your current code)
  merged <- merge_cancer_cell_types(proj_list, input_dir, file_suffix,   max_cells_per_project = 2000L,   # tune: 2000–5000 for Epithelial
                                    chunk_size            = 4L,      # tune: fewer = less peak RAM per chunk
                                    tmp_dir               = file.path(opt$output_dir, "tmp_merge", gsub("[^A-Za-z0-9._-]", "_", opt$cell_type)))
  
  # 1.1) ---- drop tiny batches BEFORE split/integration ----
  min_cells_threshold <- 30L
  batch_sizes <- table(merged$project_name)
  small_batches <- names(batch_sizes[batch_sizes < min_cells_threshold])
  
  if (length(small_batches) > 0) {
    log_msg(sprintf(
      "Dropping %d small batches: %s",
      length(small_batches),
      paste(small_batches, collapse = ", ")
    ))
    merged <- subset(merged, subset = !(project_name %in% small_batches))
  }
  
  batch_sizes <- table(merged$project_name)
  print(batch_sizes)
  
  single_batch <- (length(batch_sizes) < 2)
  if (single_batch) {
    log_msg("Only 1 batch after filtering — skipping Harmony, using PCA directly.")
  }

  # ---- prepare proper v5 layered assay ----
  log_msg("Joining layers before split...")
  merged[["RNA"]] <- JoinLayers(merged[["RNA"]], layers = "counts")
  log_msg(sprintf("Layers after join: %s",
                  paste(Layers(merged[["RNA"]]), collapse = ", ")))
  merged[["RNA"]] <- split(merged[["RNA"]], f = merged$project_name)
  log_msg(sprintf("Layers after split: %s",
                  paste(Layers(merged[["RNA"]]), collapse = ", ")))

  # ---- recompute sizes from actual split layers used for integration ----
  layer_names <- Layers(merged[["RNA"]])
  count_layers <- grep("^counts", layer_names, value = TRUE)
  data_layers  <- grep("^data",   layer_names, value = TRUE)

  if (length(data_layers) > 0) {
    layer_set <- data_layers
  } else if (length(count_layers) > 0) {
    layer_set <- count_layers
  } else {
    stop("No split RNA layers found after split(...).")
  }

  layer_sizes <- sapply(layer_set, function(ly) {
    mat <- LayerData(merged, assay = "RNA", layer = ly)
    ncol(mat)
  })

  names(layer_sizes) <- sub("^(counts|data)\\.", "", layer_set)

  log_msg("Layer sizes after split:")
  print(layer_sizes)

  small_layers <- names(layer_sizes[layer_sizes < min_cells_threshold])
  if (length(small_layers) > 0 && !single_batch) {
    stop(
      "Some split RNA layers still have < ", min_cells_threshold, " cells: ",
      paste0(small_layers, " (", layer_sizes[small_layers], ")", collapse = ", "),
      ". Lower dims or drop these batches before integration."
    )
  }

  min_cells <- min(layer_sizes)
  if (min_cells < 30L) {
    stop("Some layers have < 30 cells; too small for reliable processing. Layers: ",
         paste(names(layer_sizes), layer_sizes, collapse = ", "))
  }

  max_dims <- min(20L, floor(min_cells / 2L))
  if (max_dims < 5L) {
    stop("max_dims < 5 after checking layer sizes; not enough cells.")
  }

  log_msg(sprintf("single_batch=%s, min_cells=%d, max_dims=%d", single_batch, min_cells, max_dims))

  # ---- normalization + PCA on layered object ----
  merged <- basic_normalizations(merged, max_dims)

  if (!single_batch) {
    # ---- Harmony integration on PCA ----
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
    reduction_to_use  <- "harmony"
    graph_name        <- "Harmony_nn"
    cluster_name      <- "Harmony_cluster"
    umap_name         <- "umap_harmony"
  } else {
    # single batch: use PCA directly (no batch correction needed)
    reduction_to_use  <- "pca"
    graph_name        <- "Harmony_nn"
    cluster_name      <- "Harmony_cluster"
    umap_name         <- "umap_harmony"
  }

  # ---- neighbors / clusters / UMAP ----
  n_cells_total <- ncol(merged)
  k_param <- min(20L, max(5L, n_cells_total - 1L))

  log_msg(sprintf(
    "Neighbors: reduction=%s, n_cells_total=%d, max_dims=%d, k.param=%d",
    reduction_to_use, n_cells_total, max_dims, k_param
  ))

  merged <- FindNeighbors(
    merged,
    reduction  = reduction_to_use,
    dims       = 1:max_dims,
    graph.name = graph_name,
    k.param    = k_param,
    nn.method  = "rann"
  )

  merged <- FindClusters(
    merged,
    graph.name   = graph_name,
    resolution   = 0.6,
    cluster.name = cluster_name
  )

  merged <- RunUMAP(
    merged,
    reduction      = reduction_to_use,
    dims           = 1:max_dims,
    reduction.name = umap_name
  )

  saveRDS(merged, output_file)
  if (verbose) print(paste0("Saved merged object: ", output_file,
                            " | total cells: ", ncol(merged)))
  
  invisible(merged)
}

#--------------------------------------------------
# CLI (optparse)
#--------------------------------------------------
option_list <- list(
  make_option(c("-m", 
                "--meta_rds"), 
              type = "character",
              help = "Path to all_cells metadata RDS (data.frame). Required."),
  make_option(c("-i", "--input_dir"), 
              type = "character",
              help = "Directory containing per-project Seurat RDS files. Required."),
  make_option(c("-o", "--output_dir"), 
              type = "character",
              help = "Directory to write merged object. Required."),
  make_option(c("-c", "--cancer_type"), 
              type = "character",
              help = "Cancer type value to filter (exact match). Required."),
  make_option(c("-t", "--cell_type"), 
              type = "character",
              help = "Cell type label to filter (exact match). Required."),
  
  make_option(c("--meta_cancer_col"), 
              type = "character", 
              default = "cancer_type",
              help = "Column in metadata for cancer type [default: %default]"),
  make_option(c("--meta_celltype_col"), 
              type = "character", 
              default = "Consensus_Cell_Type",
              help = "Column in metadata for cell type [default: %default]"),
  make_option(c("--project_col"), 
              type = "character", 
              default = "project_name",
              help = "Column in metadata for project name [default: %default]"),
  make_option(c("--cellid_col"), 
              type = "character", 
              default = "cell_id",
              help = "Column in metadata for cell id/barcode [default: %default]"),
  
  make_option(c("--file_suffix"), 
              type = "character", 
              default = "_step4.RDS",
              help = "Suffix of per-project RDS files [default: %default]"),
  
  make_option(c("--diet"), 
              action = "store_true", 
              default = TRUE,
              help = "Run DietSeurat on subsets [default: TRUE]"),
  make_option(c("--no_diet"), 
              action = "store_false", 
              dest = "diet",
              help = "Disable DietSeurat"),
  make_option(c("--keep_counts_only"), 
              action = "store_true", 
              default = TRUE,
              help = "When dieting, keep counts layer only [default: TRUE]"),
  make_option(c("--keep_all_layers"), 
              action = "store_false", 
              dest = "keep_counts_only",
              help = "When dieting, do not restrict to counts only"),
  
  make_option(c("--out_name"), 
              type = "character", default = NA,
              help = "Optional output filename (without dir). If not given, auto-generated."),
  make_option(c("-q", "--quiet"), 
              action = "store_true", 
              default = FALSE,
              help = "Less verbose output.")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
# save(opt, file = "opt.R")
# quit()
# load("opt.R")

print("started")
print(opt$cancer_type)
print(opt$cell_type)
# validate required
req <- c("meta_rds", "input_dir", "output_dir", "cancer_type", "cell_type")
missing_req <- req[!nzchar(sapply(req, function(x) ifelse(is.null(opt[[x]]), "", opt[[x]])))]
if (length(missing_req) > 0) {
  print_help(parser)
  stop("Missing required arguments: ", paste(missing_req, collapse = ", "))
}

verbose <- !isTRUE(opt$quiet)

# load metadata
all_cells <- readRDS(opt$meta_rds)
if (!is.data.frame(all_cells)) {
  stop("meta_rds did not load a data.frame. Got: ", class(all_cells)[1])
}
if (verbose) print(paste0("Loaded metadata: ", nrow(all_cells), " rows"))

# filter by cancer type + cell type
meta_sub <- all_cells %>%
  filter(.data[[opt$meta_cancer_col]] == opt$cancer_type,
         .data[[opt$meta_celltype_col]] == opt$cell_type)

if (nrow(meta_sub) == 0) {
  stop("No rows matched: ",
       opt$meta_cancer_col, " == '", opt$cancer_type, "' AND ",
       opt$meta_celltype_col, " == '", opt$cell_type, "'.")
}

# output path
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
default_name <- paste0(safe(opt$cell_type), "_", safe(opt$cancer_type), "_merged.RDS")
out_name <- if (!is.na(opt$out_name) && nzchar(opt$out_name)) opt$out_name else default_name
out_file <- file.path(opt$output_dir, out_name)

if (verbose) {
  print(paste0("Filtering rows: ", nrow(meta_sub)))
  print(paste0("Writing to: ", out_file))
}

# merge and save
merge_subset_streaming(
  meta_df          = meta_sub,
  input_dir        = opt$input_dir,
  output_file      = out_file,
  project_col      = opt$project_col,
  cellid_col       = opt$cellid_col,
  file_suffix      = opt$file_suffix,
  diet             = opt$diet,
  keep_counts_only = opt$keep_counts_only,
  verbose          = verbose
)

if (verbose) print("Done.")

##usage:
# Rscript merge_by_cancertype_celltype.R \
# --meta_rds /gpfs0/bgu-ofircohen/users/hadasra/thesis/create_sc_object_workflow/temp_objects/all_cells.RDS \
# --input_dir /gpfs0/bgu-ofircohen/users/hadasra/thesis/create_sc_object_workflow/temp_objects \
# --output_dir /gpfs0/bgu-ofircohen/users/hadasra/thesis/create_sc_object_workflow/temp_objects/merged_by_query \
# --cancer_type Breast \
# --cell_type "Epithelial"

#usage in sise:
# Rscript merge_by_cancertype_celltype.R \
# --meta_rds /ofir_grp/users/hadasra/temp_objects/all_cells.RDS \
# --input_dir /ofir_grp/users/hadasra/temp_objects \
# --output_dir /ofir_grp/users/hadasra/temp_objects/merged_by_query \
# --cancer_type Breast \
# --cell_type "Epithelial"