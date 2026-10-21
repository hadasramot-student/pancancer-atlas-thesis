#!/usr/bin/env Rscript
# ch3_dedup_rerun.R
#
# Ch.3 batch-correction benchmark: sensitivity re-run of the Step-8 validation
# metrics on a duplicate-free cohort.
#
# The published metrics (all_validation_results.csv) were produced by
# merge/Batch_Correction_Validation_Metrics.R on
#   <CellType>_Breast_merged_step7_final_dev015_nocca.RDS
# with --n_cells 2000 --seed 42 --batch_col project_name
#      --label_col scimilarity_celltype.
#
# This script reproduces that exactly ("before") and re-runs it after dropping
# the Breastscape-side copy of the three cross-repository duplicate pairs whose
# BOTH copies sit in the benchmark cohort ("after"):
#     WU_SZ_ET_AL_2021_P15_filtered      (twin of 3CA Wu2021_Breast)
#     GAO_R_ET_AL_2021_P17_filtered      (twin of 3CA Gao2021_Breast_1/_2)
#     SAVAS_P_ET_AL_2018_P35_filtered    (twin of 3CA Savas2018_Breast)
#
# NOTE: the integration embeddings themselves are NOT recomputed. Metrics are
# recomputed on the existing mnn / pca.combat / harmony / scanvi embeddings
# after removing the duplicate-labelled cells from the 2,000-cell validation
# subsample. This isolates the batch-mixing-metric inflation caused by the same
# cells carrying two project_name labels; it is a sensitivity check, not a full
# re-integration.

suppressPackageStartupMessages({
  library(Seurat)
  library(lisi)
  library(caret)
})

BASE <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast"
OUT  <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/chapter3-benchmark-duplicates-3ba85e"

CELL_TYPES <- c(
  "B-cell", "Endothelial", "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells",
  "Myeloid", "NK-cells", "T-cells"
)

REMOVED_COPIES <- c(
  "WU_SZ_ET_AL_2021_P15_filtered",
  "GAO_R_ET_AL_2021_P17_filtered",
  "SAVAS_P_ET_AL_2018_P35_filtered"
)

BATCH_VAR <- "project_name"
LABEL_VAR <- "scimilarity_celltype"
N_CELLS   <- 2000L
SEED      <- 42L

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

# ---- helpers copied verbatim from merge/Batch_Correction_Validation_Metrics.R ----
compute_kbet <- function(emb, batch, k0 = 30) {
  if (!requireNamespace("kBET", quietly = TRUE)) return(NA_real_)
  tryCatch({
    res <- kBET::kBET(df = as.matrix(emb), batch = batch, k0 = k0,
                      plot = FALSE, do.pca = FALSE, verbose = FALSE)
    1 - mean(res$stats$kBET.observed, na.rm = TRUE)
  }, error = function(e) { warning("kBET failed: ", e$message); NA_real_ })
}

compute_graph_connectivity <- function(emb, label, k = 10) {
  if (!requireNamespace("igraph", quietly = TRUE)) return(NA_real_)
  if (nrow(emb) < 2L || length(unique(label)) < 2L) return(NA_real_)
  d <- as.matrix(dist(emb)); diag(d) <- Inf
  gc_scores <- c()
  for (lab in unique(label)) {
    idx <- which(label == lab)
    if (length(idx) <= 1L) next
    d_sub <- d[idx, idx, drop = FALSE]
    knn_list <- lapply(seq_len(nrow(d_sub)), function(j) {
      row <- d_sub[j, ]; ord <- order(row, decreasing = FALSE)
      head(ord, n = min(k, sum(is.finite(row))))
    })
    edges <- c()
    for (i in seq_along(idx)) {
      neighbors <- knn_list[[i]]
      neighbors <- neighbors[is.finite(d_sub[i, neighbors])]
      if (length(neighbors) > 0L) edges <- c(edges, rbind(i, neighbors))
    }
    if (length(edges) == 0L) next
    g <- igraph::make_empty_graph(n = length(idx), directed = FALSE)
    g <- igraph::add_edges(g, edges)
    lcc_size <- max(igraph::components(g)$csize)
    gc_scores <- c(gc_scores, lcc_size / length(idx))
  }
  if (length(gc_scores) == 0L) return(NA_real_)
  mean(gc_scores, na.rm = TRUE)
}

compute_f1_label <- function(label, cluster) {
  if (!requireNamespace("caret", quietly = TRUE)) return(NA_real_)
  keep <- !(is.na(label) | is.na(cluster))
  label <- factor(label[keep]); cluster <- factor(cluster[keep])
  if (length(unique(label)) < 2L || length(unique(cluster)) < 2L) return(NA_real_)
  cluster_to_label <- tapply(as.character(label), cluster, function(x)
    names(sort(table(x), decreasing = TRUE))[1])
  predicted <- cluster_to_label[as.character(cluster)]
  predicted <- factor(predicted, levels = levels(label))
  reference <- factor(as.character(label), levels = levels(label))
  cm <- caret::confusionMatrix(predicted, reference)
  f1_per_class <- 2 * cm$byClass[, "Sensitivity"] * cm$byClass[, "Pos Pred Value"] /
    (cm$byClass[, "Sensitivity"] + cm$byClass[, "Pos Pred Value"])
  mean(f1_per_class, na.rm = TRUE)
}

compute_metrics_for_reduction <- function(seu, reduction, batch_var, label_var,
                                          k_neighbors = 30) {
  emb   <- Embeddings(seu, reduction = reduction)
  batch <- seu[[batch_var]][, 1]
  label <- seu[[label_var]][, 1]
  finite_cells <- apply(emb, 1, function(x) all(is.finite(x)))
  emb <- emb[finite_cells, , drop = FALSE]
  batch <- batch[finite_cells]; label <- label[finite_cells]

  lisi_vals <- tryCatch({
    lisi_input <- as.data.frame(emb)
    lisi_input$batch <- as.factor(batch); lisi_input$label <- as.factor(label)
    lisi_res <- compute_lisi(
      X = as.matrix(lisi_input[, colnames(emb), drop = FALSE]),
      meta_data = lisi_input[, c("batch", "label")],
      label_colnames = c("batch", "label"), perplexity = k_neighbors)
    c(batch = mean(lisi_res$batch, na.rm = TRUE),
      label = mean(lisi_res$label, na.rm = TRUE))
  }, error = function(e) { log_msg(sprintf("    LISI failed: %s", e$message));
    c(batch = NA_real_, label = NA_real_) })

  kbet_accept <- tryCatch(compute_kbet(emb, batch, k0 = k_neighbors),
                          error = function(e) NA_real_)
  graph_conn  <- tryCatch(compute_graph_connectivity(emb, label, k = 10),
                          error = function(e) NA_real_)

  f1_label <- tryCatch({
    seu_tmp <- seu
    seu_tmp[[reduction]] <- CreateDimReducObject(
      embeddings = emb, key = paste0(toupper(reduction), "_"),
      assay = DefaultAssay(seu_tmp))
    seu_tmp <- FindNeighbors(seu_tmp, reduction = reduction, dims = 1:ncol(emb),
                             k.param = k_neighbors, verbose = FALSE)
    seu_tmp <- FindClusters(seu_tmp, resolution = 0.5, verbose = FALSE)
    cluster_vec <- seu_tmp$seurat_clusters[finite_cells]
    compute_f1_label(label, cluster_vec)
  }, error = function(e) { log_msg(sprintf("    F1 failed: %s", e$message)); NA_real_ })

  data.frame(reduction = reduction, batch_lisi = lisi_vals[["batch"]],
             label_lisi = lisi_vals[["label"]], kbet_accept = kbet_accept,
             graph_connectivity = graph_conn, f1_label = f1_label)
}
# --------------------------------------------------------------------------------

reductions <- c(mnn = "mnn", combat = "pca.combat", harmony = "harmony",
                scanvi = "scanvi")

run_one <- function(obj, mode) {
  do.call(rbind, lapply(names(reductions), function(m) {
    red <- reductions[[m]]
    if (!red %in% SeuratObject::Reductions(obj)) {
      log_msg(sprintf("  [%s] %s: reduction '%s' absent, skipping", mode, m, red))
      return(NULL)
    }
    log_msg(sprintf("  [%s] %s (%s)...", mode, m, red))
    df <- compute_metrics_for_reduction(obj, red, BATCH_VAR, LABEL_VAR, 30)
    df$method <- m; df$mode <- mode
    df
  }))
}

all_rows <- list()
for (cel in CELL_TYPES) {
  rds <- file.path(BASE, paste0(cel, "_Breast_merged_step7_final_dev015_nocca.RDS"))
  log_msg(sprintf("==== %s ====", cel))
  log_msg(paste("Loading:", rds))
  obj <- readRDS(rds)
  n_total <- ncol(obj)
  n_dup   <- sum(obj[[BATCH_VAR]][, 1] %in% REMOVED_COPIES)
  log_msg(sprintf("  cells=%d  duplicate-copy cells to remove=%d  reductions={%s}",
                  n_total, n_dup, paste(SeuratObject::Reductions(obj), collapse = ",")))

  # ---- BEFORE: exact reproduction of the published run ----
  obj_before <- obj
  if (n_total > N_CELLS) {
    set.seed(SEED)
    obj_before <- obj[, sample(colnames(obj), N_CELLS)]
  }
  rb <- run_one(obj_before, "before")
  if (!is.null(rb)) { rb$cell_type <- cel; rb$n_subsample <- ncol(obj_before)
                      rb$n_dup_removed <- 0L; all_rows[[length(all_rows)+1]] <- rb }

  # ---- AFTER: drop the Breastscape duplicate copy, then subsample ----
  keep_cells <- colnames(obj)[!(obj[[BATCH_VAR]][, 1] %in% REMOVED_COPIES)]
  obj_dd <- obj[, keep_cells]
  n_dd <- ncol(obj_dd)
  if (n_dd > N_CELLS) {
    set.seed(SEED)
    obj_dd <- obj_dd[, sample(colnames(obj_dd), N_CELLS)]
  }
  ra <- run_one(obj_dd, "after")
  if (!is.null(ra)) { ra$cell_type <- cel; ra$n_subsample <- ncol(obj_dd)
                      ra$n_dup_removed <- n_dup; all_rows[[length(all_rows)+1]] <- ra }

  rm(obj, obj_before, obj_dd); gc(verbose = FALSE)
}

res <- do.call(rbind, all_rows)
res <- res[, c("cell_type", "method", "mode", "reduction",
               "batch_lisi", "label_lisi", "kbet_accept",
               "graph_connectivity", "f1_label", "n_subsample", "n_dup_removed")]
write.csv(res, file.path(OUT, "ch3_dedup_rerun_results.csv"), row.names = FALSE)

log_msg("")
log_msg("=== BEFORE vs AFTER (batch_lisi / kbet_accept are the batch-mixing metrics) ===")
w <- reshape(res[, c("cell_type","method","mode","batch_lisi","kbet_accept",
                     "label_lisi","graph_connectivity","f1_label")],
             idvar = c("cell_type","method"), timevar = "mode", direction = "wide")
print(w, row.names = FALSE)
log_msg("")
log_msg("Wrote ch3_dedup_rerun_results.csv")
log_msg("DONE")
