#!/usr/bin/env Rscript
# ch3_seed_control.R
# Control: how much do the Ch.3 benchmark metrics move purely from the choice of
# the random 2,000-cell validation subsample? Re-runs merge/
# Batch_Correction_Validation_Metrics.R logic on the UNCHANGED (non-deduplicated)
# step7 objects for several seeds, so the before/after shifts in ch3_dedup_rerun
# can be read against subsample noise.

suppressPackageStartupMessages({ library(Seurat); library(lisi); library(caret) })

BASE <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast"
OUT  <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/chapter3-benchmark-duplicates-3ba85e"
CELL_TYPES <- c("B-cell", "Myeloid", "T-cells")   # compartments where our F1 harness matches the published CSV exactly
SEEDS <- c(42, 43, 44, 101)
BATCH_VAR <- "project_name"; LABEL_VAR <- "scimilarity_celltype"; N_CELLS <- 2000L

log_msg <- function(m) cat(sprintf("[%s] %s\n", Sys.time(), m))

compute_kbet <- function(emb, batch, k0 = 30) {
  if (!requireNamespace("kBET", quietly = TRUE)) return(NA_real_)
  tryCatch({ res <- kBET::kBET(df = as.matrix(emb), batch = batch, k0 = k0,
      plot = FALSE, do.pca = FALSE, verbose = FALSE)
    1 - mean(res$stats$kBET.observed, na.rm = TRUE) },
    error = function(e) NA_real_)
}
compute_graph_connectivity <- function(emb, label, k = 10) {
  if (!requireNamespace("igraph", quietly = TRUE)) return(NA_real_)
  if (nrow(emb) < 2L || length(unique(label)) < 2L) return(NA_real_)
  d <- as.matrix(dist(emb)); diag(d) <- Inf; gc_scores <- c()
  for (lab in unique(label)) {
    idx <- which(label == lab); if (length(idx) <= 1L) next
    d_sub <- d[idx, idx, drop = FALSE]
    knn_list <- lapply(seq_len(nrow(d_sub)), function(j) {
      row <- d_sub[j, ]; head(order(row), n = min(k, sum(is.finite(row)))) })
    edges <- c()
    for (i in seq_along(idx)) {
      nb <- knn_list[[i]]; nb <- nb[is.finite(d_sub[i, nb])]
      if (length(nb) > 0L) edges <- c(edges, rbind(i, nb)) }
    if (length(edges) == 0L) next
    g <- igraph::make_empty_graph(n = length(idx), directed = FALSE)
    g <- igraph::add_edges(g, edges)
    gc_scores <- c(gc_scores, max(igraph::components(g)$csize) / length(idx)) }
  if (length(gc_scores) == 0L) return(NA_real_); mean(gc_scores, na.rm = TRUE)
}
compute_f1_label <- function(label, cluster) {
  keep <- !(is.na(label) | is.na(cluster))
  label <- factor(label[keep]); cluster <- factor(cluster[keep])
  if (length(unique(label)) < 2L || length(unique(cluster)) < 2L) return(NA_real_)
  c2l <- tapply(as.character(label), cluster, function(x)
    names(sort(table(x), decreasing = TRUE))[1])
  predicted <- factor(c2l[as.character(cluster)], levels = levels(label))
  reference <- factor(as.character(label), levels = levels(label))
  cm <- caret::confusionMatrix(predicted, reference)
  f1 <- 2 * cm$byClass[, "Sensitivity"] * cm$byClass[, "Pos Pred Value"] /
    (cm$byClass[, "Sensitivity"] + cm$byClass[, "Pos Pred Value"])
  mean(f1, na.rm = TRUE)
}
metrics_for <- function(seu, reduction, k = 30) {
  emb <- Embeddings(seu, reduction = reduction)
  batch <- seu[[BATCH_VAR]][, 1]; label <- seu[[LABEL_VAR]][, 1]
  fin <- apply(emb, 1, function(x) all(is.finite(x)))
  emb <- emb[fin, , drop = FALSE]; batch <- batch[fin]; label <- label[fin]
  lv <- tryCatch({ li <- as.data.frame(emb); li$batch <- factor(batch); li$label <- factor(label)
    r <- compute_lisi(as.matrix(li[, colnames(emb)]), li[, c("batch","label")],
                      c("batch","label"), perplexity = k)
    c(batch = mean(r$batch, na.rm = TRUE), label = mean(r$label, na.rm = TRUE)) },
    error = function(e) c(batch = NA, label = NA))
  kb <- tryCatch(compute_kbet(emb, batch, k0 = k), error = function(e) NA_real_)
  gc <- tryCatch(compute_graph_connectivity(emb, label, 10), error = function(e) NA_real_)
  f1 <- tryCatch({ st <- seu
    st[[reduction]] <- CreateDimReducObject(embeddings = emb,
        key = paste0(toupper(reduction), "_"), assay = DefaultAssay(st))
    st <- FindNeighbors(st, reduction = reduction, dims = 1:ncol(emb), k.param = k, verbose = FALSE)
    st <- FindClusters(st, resolution = 0.5, verbose = FALSE)
    compute_f1_label(label, st$seurat_clusters[fin]) }, error = function(e) NA_real_)
  data.frame(batch_lisi = lv[["batch"]], label_lisi = lv[["label"]],
             kbet_accept = kb, graph_connectivity = gc, f1_label = f1)
}
reductions <- c(mnn = "mnn", combat = "pca.combat", harmony = "harmony", scanvi = "scanvi")

out <- list()
for (cel in CELL_TYPES) {
  rds <- file.path(BASE, paste0(cel, "_Breast_merged_step7_final_dev015_nocca.RDS"))
  log_msg(sprintf("==== %s ====", cel)); obj <- readRDS(rds)
  for (s in SEEDS) {
    o <- obj
    if (ncol(obj) > N_CELLS) { set.seed(s); o <- obj[, sample(colnames(obj), N_CELLS)] }
    for (m in names(reductions)) {
      log_msg(sprintf("  seed %d  %s", s, m))
      df <- metrics_for(o, reductions[[m]]); df$cell_type <- cel; df$seed <- s; df$method <- m
      out[[length(out)+1]] <- df
    }
  }
  rm(obj); gc(verbose = FALSE)
}
res <- do.call(rbind, out)
res <- res[, c("cell_type","method","seed","batch_lisi","label_lisi","kbet_accept","graph_connectivity","f1_label")]
write.csv(res, file.path(OUT, "ch3_seed_control_results.csv"), row.names = FALSE)

log_msg("=== spread across seeds (min / max / sd) per cell_type x method ===")
for (cel in CELL_TYPES) for (m in names(reductions)) {
  sub <- res[res$cell_type == cel & res$method == m, ]
  cat(sprintf("%-9s %-7s  iLISI[%.2f-%.2f sd=%.3f]  kBET[%.3f-%.3f sd=%.3f]  F1[%.3f-%.3f sd=%.3f]\n",
    cel, m, min(sub$batch_lisi), max(sub$batch_lisi), sd(sub$batch_lisi),
    min(sub$kbet_accept), max(sub$kbet_accept), sd(sub$kbet_accept),
    min(sub$f1_label), max(sub$f1_label), sd(sub$f1_label)))
}
log_msg("DONE")
