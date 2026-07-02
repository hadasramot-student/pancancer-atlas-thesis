#!/usr/bin/env Rscript

#' Batch-Correction Validation Metrics (improved version)
#'
#' Evaluates all batch correction methods (CCA, MNN, ComBat, Harmony, scANVI)
#' Computes per-method metrics and saves CSV for cross-object comparison
#'
#' Metrics:
#'   - batch_lisi (iLISI)  : higher = better batch mixing
#'   - label_lisi (cLISI)  : lower = better cell-type separation
#'   - kbet_accept         : higher = better batch mixing
#'   - graph_connectivity  : higher = better cell-type preservation
#'   - f1_label           : higher = better cluster-label alignment
#'
#' Install once if needed:
#'   remotes::install_github("immunogenomics/lisi")
#'   remotes::install_github("theislab/kBET")

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(lisi)
  library(caret)
})

option_list <- list(
  make_option(c("--input_rds"),  type="character", default=NULL,
              help="Input Seurat object RDS file"),
  make_option(c("--output_dir"), type="character", default="./validation",
              help="Output directory for results [default=./validation]"),
  make_option(c("--batch_col"),  type="character", default="project_name",
              help="Metadata column for batch [default=project_name]"),
  make_option(c("--label_col"),  type="character", default="scimilarity_celltype",
              help="Metadata column for cell type [default=scimilarity_celltype]"),
  make_option(c("--n_cells"),    type="integer",   default=NULL,
              help="Optional: subsample to N cells for faster computation"),
  make_option(c("--seed"),       type="integer",   default=42)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_rds)) stop("Please provide --input_rds")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

log_msg(paste("Loading:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Loaded: %d cells, %d genes", ncol(obj), nrow(obj)))

# Optional subsampling for speed
if (!is.null(opt$n_cells) && ncol(obj) > opt$n_cells) {
  set.seed(opt$seed)
  cells_keep <- sample(colnames(obj), opt$n_cells)
  obj <- obj[, cells_keep]
  log_msg(sprintf("Subsampled to %d cells for validation", ncol(obj)))
}

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# kBET helper
# ============================================================
compute_kbet <- function(emb, batch, k0 = 30) {
  if (!requireNamespace("kBET", quietly = TRUE)) {
    message("kBET not installed. Run: remotes::install_github('theislab/kBET')")
    return(NA_real_)
  }
  tryCatch({
    res <- kBET::kBET(
      df      = as.matrix(emb),
      batch   = batch,
      k0      = k0,
      plot    = FALSE,
      do.pca  = FALSE,
      verbose = FALSE
    )
    1 - mean(res$stats$kBET.observed, na.rm = TRUE)
  }, error = function(e) {
    warning("kBET failed: ", e$message)
    NA_real_
  })
}

# ============================================================
# Graph connectivity helper
# ============================================================
compute_graph_connectivity <- function(emb, label, k = 10) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("igraph not installed. Install with install.packages('igraph').")
    return(NA_real_)
  }

  if (nrow(emb) < 2L || length(unique(label)) < 2L) {
    return(NA_real_)
  }

  d <- as.matrix(dist(emb))
  diag(d) <- Inf

  gc_scores <- c()

  for (lab in unique(label)) {
    idx <- which(label == lab)
    if (length(idx) <= 1L) next

    d_sub <- d[idx, idx, drop = FALSE]
    # Use lapply (not apply) so the result is always a list, regardless of
    # whether groups return equal numbers of neighbors. apply() collapses to
    # an integer/numeric vector when each row yields a single neighbor (e.g.
    # groups with only 2 cells), which breaks matrix-style `knn_idx[, i]` with
    # "incorrect number of dimensions".
    knn_list <- lapply(seq_len(nrow(d_sub)), function(j) {
      row <- d_sub[j, ]
      ord <- order(row, decreasing = FALSE)
      head(ord, n = min(k, sum(is.finite(row))))
    })

    edges <- c()
    for (i in seq_along(idx)) {
      neighbors <- knn_list[[i]]
      neighbors <- neighbors[is.finite(d_sub[i, neighbors])]
      if (length(neighbors) > 0L) {
        edges <- c(edges, rbind(i, neighbors))
      }
    }
    if (length(edges) == 0L) next

    g <- igraph::make_empty_graph(n = length(idx), directed = FALSE)
    g <- igraph::add_edges(g, edges)

    lcc_size <- max(igraph::components(g)$csize)
    gc_scores <- c(gc_scores, lcc_size / length(idx))
  }

  if (length(gc_scores) == 0L) {
    return(NA_real_)
  }
  mean(gc_scores, na.rm = TRUE)
}

# ============================================================
# Label-cluster F1 helper
# ============================================================
compute_f1_label <- function(label, cluster) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    message("caret not installed. Install with install.packages('caret').")
    return(NA_real_)
  }

  keep <- !(is.na(label) | is.na(cluster))
  label   <- factor(label[keep])
  cluster <- factor(cluster[keep])

  if (length(unique(label)) < 2L || length(unique(cluster)) < 2L) {
    return(NA_real_)
  }

  # Cluster IDs (0,1,2,...) and cell-type labels are different label sets, so
  # confusionMatrix(cluster, label) errors because their factor levels do not
  # overlap. Map each cluster to its majority cell type to obtain predicted
  # labels in the same space as the reference, then compute macro F1.
  cluster_to_label <- tapply(as.character(label), cluster, function(x) {
    names(sort(table(x), decreasing = TRUE))[1]
  })
  predicted <- cluster_to_label[as.character(cluster)]

  predicted <- factor(predicted, levels = levels(label))
  reference <- factor(as.character(label), levels = levels(label))

  cm <- caret::confusionMatrix(predicted, reference)
  f1_per_class <- 2 * cm$byClass[, "Sensitivity"] * cm$byClass[, "Pos Pred Value"] /
    (cm$byClass[, "Sensitivity"] + cm$byClass[, "Pos Pred Value"])
  mean(f1_per_class, na.rm = TRUE)
}

# ============================================================
# Main metrics function
# ============================================================
compute_metrics_for_reduction <- function(seu, reduction, batch_var, label_var,
                                          k_neighbors = 30) {
  emb   <- Embeddings(seu, reduction = reduction)
  batch <- seu[[batch_var]][, 1]
  label <- seu[[label_var]][, 1]

  finite_cells <- apply(emb, 1, function(x) all(is.finite(x)))
  emb   <- emb[finite_cells, , drop = FALSE]
  batch <- batch[finite_cells]
  label <- label[finite_cells]

  # Each metric is wrapped independently so that one failing metric does not
  # blank out the entire row (returns NA for the failed metric instead).

  # LISI
  lisi_vals <- tryCatch({
    lisi_input        <- as.data.frame(emb)
    lisi_input$batch  <- as.factor(batch)
    lisi_input$label  <- as.factor(label)

    lisi_res <- compute_lisi(
      X              = as.matrix(lisi_input[, colnames(emb), drop = FALSE]),
      meta_data      = lisi_input[, c("batch", "label")],
      label_colnames = c("batch", "label"),
      perplexity     = k_neighbors
    )
    c(batch = mean(lisi_res$batch, na.rm = TRUE),
      label = mean(lisi_res$label, na.rm = TRUE))
  }, error = function(e) {
    log_msg(sprintf("    LISI failed: %s", e$message))
    c(batch = NA_real_, label = NA_real_)
  })
  batch_lisi <- lisi_vals[["batch"]]
  label_lisi <- lisi_vals[["label"]]

  # kBET
  kbet_accept <- tryCatch(compute_kbet(emb, batch, k0 = k_neighbors),
                          error = function(e) {
                            log_msg(sprintf("    kBET failed: %s", e$message))
                            NA_real_
                          })

  # Graph connectivity
  graph_conn <- tryCatch(compute_graph_connectivity(emb, label, k = 10),
                         error = function(e) {
                           log_msg(sprintf("    graph_connectivity failed: %s", e$message))
                           NA_real_
                         })

  # Clustering + F1
  f1_label <- tryCatch({
    seu_tmp <- seu
    seu_tmp[[reduction]] <- CreateDimReducObject(
      embeddings = emb,
      key       = paste0(toupper(reduction), "_"),
      assay     = DefaultAssay(seu_tmp)
    )

    seu_tmp <- FindNeighbors(
      seu_tmp,
      reduction = reduction,
      dims      = 1:ncol(emb),
      k.param   = k_neighbors,
      verbose   = FALSE
    )
    seu_tmp <- FindClusters(
      seu_tmp,
      resolution = 0.5,
      verbose    = FALSE
    )

    cluster_vec <- seu_tmp$seurat_clusters[finite_cells]
    compute_f1_label(label, cluster_vec)
  }, error = function(e) {
    log_msg(sprintf("    F1 failed: %s", e$message))
    NA_real_
  })

  data.frame(
    reduction         = reduction,
    batch_lisi        = batch_lisi,
    label_lisi        = label_lisi,
    kbet_accept       = kbet_accept,
    graph_connectivity = graph_conn,
    f1_label          = f1_label
  )
}

# ============================================================
# Run across all reductions
# ============================================================
reductions <- c(
  cca    = "integrated.cca",
  mnn    = "mnn",
  combat = "pca.combat",
  harmony = "harmony",
  scanvi  = "scanvi"
)

log_msg("Computing metrics for all methods...")
all_metrics <- do.call(rbind, lapply(names(reductions), function(m) {
  red <- reductions[[m]]
  log_msg(sprintf("  %s (%s)...", m, red))
  tryCatch({
    df <- compute_metrics_for_reduction(
      seu         = obj,
      reduction   = red,
      batch_var   = opt$batch_col,
      label_var   = opt$label_col,
      k_neighbors = 30
    )
    df$method <- m
    df[, c("method", setdiff(colnames(df), "method"))]
  }, error = function(e) {
    log_msg(sprintf("  ✗ Skipping %s: %s", m, e$message))
    NULL
  })
}))

# ============================================================
# Results
# ============================================================
log_msg("")
log_msg("VALIDATION RESULTS")
log_msg("==================")
print(all_metrics)

# Save to CSV
output_file <- file.path(opt$output_dir, "batch_correction_metrics.csv")
write.csv(all_metrics, output_file, row.names = FALSE)
log_msg(paste("Saved:", output_file))

# Summary interpretation
log_msg("")
log_msg("INTERPRETATION GUIDE:")
log_msg("  batch_lisi (iLISI): higher = better batch mixing")
log_msg("  label_lisi (cLISI): lower = better cell-type separation")
log_msg("  kbet_accept: higher = better batch mixing")
log_msg("  graph_connectivity: higher = better cell-type preservation")
log_msg("  f1_label: higher = better cluster-label alignment")
log_msg("")
