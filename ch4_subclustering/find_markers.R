#!/usr/bin/env Rscript
# find_markers.R
#
# Ch.4 (Subclustering) — FindAllMarkers (Wilcoxon) for pan-cancer per-lineage
# x cancer-type objects.
#
# For each *_clustered.RDS:
#   1. Runs FindAllMarkers() at a chosen cluster resolution
#   2. Filters: avg_log2FC > 0.25, p_val_adj < 0.05, min.pct >= 0.1
#   3. Saves full marker table  -> <input_dir>/<lineage>_markers.csv
#   4. Saves top-10 per cluster -> <input_dir>/<lineage>_top_markers.csv
#      columns: cluster | gene | avg_log2FC | pct.1 | pct.2 | p_val_adj
#
# Run via R container (one lineage):
#   sbatch --mem=128G ~/rstudio/rstudio.sbatch \
#     --script ch4_subclustering/find_markers.R \
#     --input_dir /path/to/temp_objects/clustered_lineages \
#     --lineage T-cells_Breast

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(future)
})

option_list <- list(
  make_option("--input_dir",   type = "character", default = NULL,
    help = "Directory containing *_clustered.RDS files [required]"),
  make_option("--lineage",     type = "character", default = NULL,
    help = "Process only files matching this string (e.g. T-cells_Breast). Default: all"),
  make_option("--cluster_col", type = "character", default = "louvain_res0.50",
    help = "Cluster column to use for FindAllMarkers [louvain_res0.50]"),
  make_option("--logfc",       type = "double",    default = 0.25,
    help = "log2FC threshold [0.25]"),
  make_option("--min_pct",     type = "double",    default = 0.1,
    help = "min.pct threshold [0.1]"),
  make_option("--padj",        type = "double",    default = 0.05,
    help = "Adjusted p-value cutoff [0.05]"),
  make_option("--top_n",       type = "integer",   default = 10L,
    help = "Top N markers per cluster [10]"),
  make_option("--n_cores",     type = "integer",   default = 8L,
    help = "Parallel workers for FindAllMarkers [8]"),
  make_option("--overwrite",   action = "store_true", default = FALSE,
    help = "Overwrite existing marker CSVs")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_dir)) stop("--input_dir is required")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

# Parallelise FindAllMarkers across cores
plan("multicore", workers = opt$n_cores)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB per worker

# ── Discover files ────────────────────────────────────────────────────────────
all_rds <- list.files(opt$input_dir, pattern = "_clustered\\.RDS$", full.names = TRUE)
if (!is.null(opt$lineage))
  all_rds <- all_rds[grepl(opt$lineage, basename(all_rds), fixed = TRUE)]
if (length(all_rds) == 0) stop("No matching *_clustered.RDS files found in: ", opt$input_dir)

log_msg("Found ", length(all_rds), " file(s) to process")

# ── Main loop ─────────────────────────────────────────────────────────────────
for (rds_path in all_rds) {

  lineage <- sub("_clustered\\.RDS$", "", basename(rds_path))
  out_markers <- file.path(opt$input_dir, paste0(lineage, "_markers.csv"))
  out_top     <- file.path(opt$input_dir, paste0(lineage, "_top_markers.csv"))

  log_msg("═══════════════════════════════════════════")
  log_msg("Processing: ", lineage)

  if (file.exists(out_markers) && !isTRUE(opt$overwrite)) {
    log_msg("  Markers CSV already exists, skipping (use --overwrite to force)")
    next
  }

  # ── Load ──────────────────────────────────────────────────────────────────
  log_msg("  Loading RDS ...")
  obj <- readRDS(rds_path)
  log_msg("  Loaded: ", ncol(obj), " cells × ", nrow(obj), " genes")

  # Join Seurat v5 split layers
  if (inherits(obj[["RNA"]], "Assay5")) {
    log_msg("  JoinLayers (Seurat v5) ...")
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  }
  DefaultAssay(obj) <- "RNA"

  # ── Set cluster identity ──────────────────────────────────────────────────
  cluster_col <- opt$cluster_col

  # Fall back if requested column not present
  if (!cluster_col %in% colnames(obj@meta.data)) {
    avail <- grep("^louvain_res|^leiden_res", colnames(obj@meta.data), value = TRUE)
    if (length(avail) == 0) stop("No cluster columns found in object: ", lineage)
    cluster_col <- sort(avail)[ceiling(length(avail) / 2)]  # pick middle resolution
    log_msg("  '", opt$cluster_col, "' not found — using '", cluster_col, "' instead")
  }

  Idents(obj) <- obj@meta.data[[cluster_col]]
  n_clusters  <- nlevels(Idents(obj))
  log_msg("  Cluster column: ", cluster_col, "  (", n_clusters, " clusters)")

  # ── FindAllMarkers ────────────────────────────────────────────────────────
  log_msg("  Running FindAllMarkers (Wilcoxon, logfc=", opt$logfc,
          ", min.pct=", opt$min_pct, ", workers=", opt$n_cores, ") ...")

  markers <- FindAllMarkers(
    obj,
    test.use        = "wilcox",
    logfc.threshold = opt$logfc,
    min.pct         = opt$min_pct,
    only.pos        = FALSE,
    verbose         = FALSE
  )

  log_msg("  Raw markers: ", nrow(markers), " rows")

  # ── Filter ────────────────────────────────────────────────────────────────
  markers_filt <- subset(markers,
    avg_log2FC > opt$logfc & p_val_adj < opt$padj
  )
  log_msg("  After filter (log2FC>", opt$logfc, ", padj<", opt$padj, "): ",
          nrow(markers_filt), " rows")

  # Add lineage column for easy stacking later
  markers_filt$lineage <- lineage
  markers_filt         <- markers_filt[order(markers_filt$cluster,
                                              -markers_filt$avg_log2FC), ]

  # ── Save full markers ─────────────────────────────────────────────────────
  write.csv(markers_filt, out_markers, row.names = FALSE, quote = FALSE)
  log_msg("  Full markers → ", out_markers)

  # ── Top N per cluster ─────────────────────────────────────────────────────
  top_markers <- do.call(rbind, lapply(split(markers_filt, markers_filt$cluster), function(df) {
    df <- df[order(-df$avg_log2FC), ]
    head(df[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")],
         opt$top_n)
  }))
  rownames(top_markers) <- NULL

  write.csv(top_markers, out_top, row.names = FALSE, quote = FALSE)
  log_msg("  Top ", opt$top_n, " markers → ", out_top)

  # ── Console summary ───────────────────────────────────────────────────────
  cat("\n  ── Top 3 markers per cluster ──\n")
  for (cl in sort(unique(top_markers$cluster))) {
    sub <- head(top_markers[top_markers$cluster == cl, ], 3)
    genes <- paste(sprintf("%s(%.2f)", sub$gene, sub$avg_log2FC), collapse = ", ")
    cat(sprintf("  Cluster %-4s: %s\n", cl, genes))
  }
  cat("\n")

  rm(obj, markers, markers_filt, top_markers); gc(verbose = FALSE)
}

log_msg("═══════════════════════════════════════════")
log_msg("All done.")
