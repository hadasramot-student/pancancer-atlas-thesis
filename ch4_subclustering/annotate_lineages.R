#!/usr/bin/env Rscript
# annotate_lineages.R
#
# Job 1 — SCimilarity cell-type characterisation for pan-cancer per-lineage objects.
#
# Adds new columns directly to each existing *_clustered.RDS (in-place overwrite).
# No duplicate files created.
#
# For each lineage RDS file the script:
#   1. Loads the Seurat object (expects UMAP + Louvain/Leiden cluster columns)
#   2. Obtains per-cell SCimilarity labels — two modes:
#        a) REUSE : column already present (--scimilarity_col exists in metadata)
#        b) RERUN : export → h5ad → Python run_scimilarity.py → import CSV
#   3. Majority-vote aggregation per cluster across ALL louvain_res* columns found
#   4. Overwrites the input RDS with the 4 new columns added
#   5. Saves cluster summary TSV alongside input → <input_dir>/<lineage>_scimilarity_summary.tsv
#      columns: resolution | cluster | majority_label | confidence_pct | n_cells
#
# HPC job-array usage (one lineage per task):
#   Rscript annotate_lineages.R --input_dir <dir> --lineage T-cells_Breast
#
# Batch usage (all lineages sequentially):
#   Rscript annotate_lineages.R --input_dir <dir>
#
# Run via R container:
#   sbatch [--mem=Xg] ~/rstudio/rstudio.sbatch \
#     --script ch4_subclustering/annotate_lineages.R \
#     --input_dir /path/to/temp_objects/clustered_lineages \
#     [--lineage T-cells_Breast]

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
})

# ── CLI ──────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input_dir",     type = "character", default = NULL,
    help = "Directory containing *_clustered.RDS files [required]"),
  make_option("--lineage",       type = "character", default = NULL,
    help = "Process only this lineage (e.g. T-cells_Breast). Default: all"),
  make_option("--cluster_col",   type = "character", default = NULL,
    help = "Primary cluster column for summary. Default: auto-detect louvain_res0.50 or first louvain_res*"),
  make_option("--scimilarity_col", type = "character", default = "scimilarity_celltype",
    help = "Metadata column with existing SCimilarity labels. If present, skip Python re-run [scimilarity_celltype]"),
  make_option("--py_script",     type = "character",
    default = "ch4_subclustering/run_scimilarity.py",
    help = "Path to run_scimilarity.py (only needed when re-running SCimilarity)"),
  make_option("--model_path",    type = "character",
    default = Sys.getenv("SCIMILARITY_MODEL_PATH", unset = "models/scimilarity/model_v1.1"),
    help = "SCimilarity model directory (only needed when re-running) [$SCIMILARITY_MODEL_PATH]"),
  make_option("--conda_env",     type = "character", default = "scimilarity",
    help = "Conda env that has scimilarity installed [scimilarity]"),
  make_option("--threshold",     type = "double",    default = 0.60,
    help = "vs2nd confidence threshold for filtering [0.60]"),
  make_option("--k",             type = "integer",   default = 50L,
    help = "KNN neighbours for SCimilarity [50]"),
  make_option("--use_gpu",       action = "store_true", default = FALSE,
    help = "Pass --use_gpu to run_scimilarity.py"),
  make_option("--skip_if_done",  action = "store_true", default = FALSE,
    help = "Skip lineage if annotated RDS already exists in output_dir"),
  make_option("--overwrite",     action = "store_true", default = FALSE,
    help = "Overwrite existing annotated RDS (default: skip)")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_dir)) stop("--input_dir is required")

# ── Logging ──────────────────────────────────────────────────────────────────
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

# ── Discover RDS files ────────────────────────────────────────────────────────
all_rds <- list.files(opt$input_dir, pattern = "\\.RDS$", full.names = TRUE)
if (length(all_rds) == 0) stop("No .RDS files found in: ", opt$input_dir)

if (!is.null(opt$lineage)) {
  # Job-array mode: match files whose basename contains the lineage name
  all_rds <- all_rds[grepl(opt$lineage, basename(all_rds), fixed = TRUE)]
  if (length(all_rds) == 0) stop("No RDS found matching lineage: ", opt$lineage)
}

log_msg("Found ", length(all_rds), " RDS file(s) to process")

# ── Helper: majority vote per cluster column ──────────────────────────────────
majority_vote <- function(labels, clusters, resolution_label) {
  df <- data.frame(
    label   = as.character(labels),
    cluster = as.character(clusters),
    stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(sort(unique(df$cluster)), function(cl) {
    sub <- df[df$cluster == cl, ]
    tbl <- sort(table(sub$label), decreasing = TRUE)
    maj_label  <- names(tbl)[1]
    confidence <- round(100 * tbl[1] / nrow(sub), 1)
    data.frame(
      resolution    = resolution_label,
      cluster       = cl,
      majority_label = maj_label,
      confidence_pct = as.numeric(confidence),
      n_cells       = nrow(sub),
      stringsAsFactors = FALSE
    )
  }))
}

# ── Export minimal h5ad for Python scimilarity ────────────────────────────────
export_h5ad <- function(obj, h5ad_path) {
  if (!requireNamespace("anndata", quietly = TRUE))
    stop("R package 'anndata' is required to export h5ad. Install via: install.packages('anndata')")
  library(anndata)

  if (inherits(obj[["RNA"]], "Assay5"))
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

  lognorm <- t(GetAssayData(obj, assay = "RNA", layer = "data"))
  counts  <- t(GetAssayData(obj, assay = "RNA", layer = "counts"))

  # Fallback: compute lognorm from counts if data layer is empty
  if (sum(lognorm) == 0) {
    log_msg("  WARNING: data layer empty — computing lognorm from counts")
    raw <- t(counts)
    cs  <- Matrix::colSums(raw); cs[cs == 0] <- 1
    lognorm <- t(log1p(Matrix::t(raw) / cs * 1e4))
  }

  ad <- AnnData(
    X      = lognorm,
    layers = list(counts = counts),
    obs    = obj@meta.data
  )
  ad$obs_names <- colnames(obj)
  ad$var_names <- rownames(obj)
  ad$write_h5ad(h5ad_path)
  log_msg("  Exported h5ad: ", h5ad_path,
          "  (", nrow(lognorm), " cells × ", ncol(lognorm), " genes)")
}

# ── Call Python scimilarity, return label data.frame ─────────────────────────
run_scimilarity_py <- function(h5ad_path, csv_path, opt) {
  gpu_flag <- if (opt$use_gpu) "--use_gpu" else ""
  cmd <- paste(
    "source /storage/modules/packages/anaconda/etc/profile.d/conda.sh &&",
    "conda run -n", shQuote(opt$conda_env),
    "python", shQuote(opt$py_script),
    "--h5ad",       shQuote(h5ad_path),
    "--out_csv",    shQuote(csv_path),
    "--model_path", shQuote(opt$model_path),
    "--k",          opt$k,
    "--threshold",  opt$threshold,
    gpu_flag
  )
  log_msg("  Running SCimilarity (conda env: ", opt$conda_env, ") ...")
  ret <- system(cmd, intern = FALSE)
  if (ret != 0) stop("run_scimilarity.py exited with code: ", ret)
  read.csv(csv_path, stringsAsFactors = FALSE)
}

# ── Main loop ─────────────────────────────────────────────────────────────────
for (rds_path in all_rds) {

  fname    <- basename(rds_path)
  lineage  <- sub("_clustered\\.RDS$", "", fname)   # e.g. T-cells_Breast
  lineage  <- sub("\\.RDS$", "", lineage)            # fallback for non-standard names

  log_msg("═══════════════════════════════════════════")
  log_msg("Processing: ", lineage)
  log_msg("  Input: ", rds_path)

  # Save back to the same file (in-place); TSV sits alongside it
  out_rds <- rds_path
  out_tsv <- file.path(dirname(rds_path), paste0(lineage, "_scimilarity_summary.tsv"))

  if (file.exists(out_tsv) && !isTRUE(opt$overwrite)) {
    log_msg("  Summary TSV already exists, skipping (use --overwrite to force): ", out_tsv)
    next
  }

  # ── Load ────────────────────────────────────────────────────────────────────
  log_msg("  Loading RDS ...")
  obj <- readRDS(rds_path)
  log_msg("  Loaded: ", ncol(obj), " cells × ", nrow(obj), " genes")
  log_msg("  Reductions: ", paste(names(obj@reductions), collapse = ", "))
  log_msg("  Metadata cols: ", paste(colnames(obj@meta.data), collapse = ", "))

  # ── SCimilarity labels ───────────────────────────────────────────────────────
  scim_col_in_obj <- opt$scimilarity_col %in% colnames(obj@meta.data)

  if (scim_col_in_obj) {
    log_msg("  Using existing SCimilarity column: ", opt$scimilarity_col)
    scim_labels <- as.character(obj@meta.data[[opt$scimilarity_col]])
    obj$scimilarity_label          <- scim_labels
    obj$scimilarity_label_filtered <- ifelse(
      !is.na(obj$scimilarity_vs2nd) & obj$scimilarity_vs2nd >= opt$threshold,
      scim_labels, "Unknown"
    )
    if (!"scimilarity_vs2nd" %in% colnames(obj@meta.data)) {
      log_msg("  NOTE: 'scimilarity_vs2nd' not found — filtered column = raw labels")
      obj$scimilarity_label_filtered <- scim_labels
    }

  } else {
    # Re-run Python scimilarity
    log_msg("  Column '", opt$scimilarity_col, "' not found — running SCimilarity via Python")
    if (!file.exists(opt$py_script))
      stop("--py_script not found: ", opt$py_script,
           "\nEither provide --scimilarity_col (if labels already exist) or a valid --py_script")
    if (!file.exists(opt$model_path))
      stop("--model_path not found: ", opt$model_path)

    h5ad_tmp <- file.path(dirname(rds_path), paste0(lineage, "_tmp.h5ad"))
    csv_tmp  <- file.path(dirname(rds_path), paste0(lineage, "_scim_pred.csv"))

    export_h5ad(obj, h5ad_tmp)
    pred_df <- run_scimilarity_py(h5ad_tmp, csv_tmp, opt)
    file.remove(h5ad_tmp)

    # Match cells by name
    common <- intersect(colnames(obj), pred_df$cell)
    if (length(common) < 0.9 * ncol(obj))
      stop("Only ", length(common), "/", ncol(obj), " cells matched between Seurat and CSV")
    rownames(pred_df) <- pred_df$cell

    obj$scimilarity_label          <- pred_df[colnames(obj), "scimilarity_label"]
    obj$scimilarity_label_filtered <- pred_df[colnames(obj), "scimilarity_label_filtered"]
    if ("confidence" %in% colnames(pred_df))
      obj$scimilarity_vs2nd        <- pred_df[colnames(obj), "confidence"]
    log_msg("  Labels added (", sum(obj$scimilarity_label_filtered != "Unknown"),
            " above threshold=", opt$threshold, ")")
  }

  # ── Identify cluster columns ──────────────────────────────────────────────
  louvain_cols <- grep("^louvain_res", colnames(obj@meta.data), value = TRUE)
  leiden_cols  <- grep("^leiden_res",  colnames(obj@meta.data), value = TRUE)
  cluster_cols <- sort(c(louvain_cols, leiden_cols))

  if (length(cluster_cols) == 0) {
    # Fall back to Seurat_clusters or user-specified col
    if (!is.null(opt$cluster_col) && opt$cluster_col %in% colnames(obj@meta.data)) {
      cluster_cols <- opt$cluster_col
    } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
      cluster_cols <- "seurat_clusters"
    } else {
      log_msg("  WARNING: No cluster columns found — skipping majority vote")
      cluster_cols <- character(0)
    }
  }

  # Primary column for quick reference
  primary_col <- if (!is.null(opt$cluster_col) && opt$cluster_col %in% cluster_cols) {
    opt$cluster_col
  } else if ("louvain_res0.50" %in% cluster_cols) {
    "louvain_res0.50"
  } else {
    cluster_cols[1]
  }

  log_msg("  Cluster columns: ", paste(cluster_cols, collapse = ", "))
  log_msg("  Primary column:  ", primary_col)

  # ── Majority vote across all resolutions ─────────────────────────────────
  all_summaries <- lapply(cluster_cols, function(col) {
    majority_vote(obj$scimilarity_label, obj@meta.data[[col]], col)
  })
  summary_df <- do.call(rbind, all_summaries)

  # Add majority label back to object for primary resolution
  if (length(cluster_cols) > 0) {
    vote_primary <- subset(summary_df, resolution == primary_col)[, c("cluster", "majority_label", "confidence_pct")]
    names(vote_primary) <- c("cluster_id", "scimilarity_cluster_label", "scimilarity_cluster_conf")
    rownames(vote_primary) <- vote_primary$cluster_id
    cell_clusters  <- as.character(obj@meta.data[[primary_col]])
    obj$scimilarity_cluster_label <- vote_primary[cell_clusters, "scimilarity_cluster_label"]
    obj$scimilarity_cluster_conf  <- vote_primary[cell_clusters, "scimilarity_cluster_conf"]
  }

  # ── Save summary TSV ──────────────────────────────────────────────────────
  write.table(summary_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("  Summary TSV: ", out_tsv)

  # Print summary for primary resolution to console
  primary_summary <- subset(summary_df, resolution == primary_col)
  cat("\n  ── Cluster summary (", primary_col, ") ──\n", sep = "")
  cat(sprintf("  %-10s  %-35s  %8s  %8s\n", "Cluster", "Majority label", "Conf %", "N cells"))
  for (i in seq_len(nrow(primary_summary))) {
    r <- primary_summary[i, ]
    cat(sprintf("  %-10s  %-35s  %7.1f%%  %8d\n",
                r$cluster, r$majority_label, r$confidence_pct, r$n_cells))
  }
  cat("\n")

  # ── Overwrite input RDS with new columns added ────────────────────────────
  log_msg("  Saving updated RDS (in-place) → ", out_rds)
  saveRDS(obj, out_rds)
  log_msg("  Done: ", lineage)
  cat("\n")
}

log_msg("═══════════════════════════════════════════")
log_msg("All done. Columns added in-place to: ", opt$input_dir)
