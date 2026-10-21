#!/usr/bin/env Rscript
# extract_batch_umap_data.R  (comment 18)
#
# Exports, for one Breast benchmark lineage, the 2-D UMAP coordinates of the
# five batch-correction views (Uncorrected, MNN, ComBat, Harmony, scANVI)
# plus the colouring metadata, so the figures can be redrawn without
# touching the analysis. Data loading, scANVI attachment and UMAP settings
# are copied unchanged from
#   pancancer-atlas-thesis/ch3_batch_correction/plot_publication_figures.R
# (step4 RDS + step6 scANVI latent CSV; RunUMAP n.neighbors 30, min.dist 0.3,
# Seurat default seed.use = 42). A UMAP is recomputed only if it is not
# stored in the object; every recomputation is logged.
#
# Usage (R container):
#   sbatch ~/rstudio/rstudio.sbatch --script extract_batch_umap_data.R \
#     --lineage T-cells --out_dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--lineage", type = "character"),
  make_option("--base",    type = "character",
              default = "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast"),
  make_option("--out_dir", type = "character")
)))

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

recomputed <- character(0)
ensure_umap <- function(obj, from_reduction, umap_name,
                        n_neighbors = 30L, min_dist = 0.3) {
  if (umap_name %in% names(obj@reductions)) {
    log_msg("  UMAP '", umap_name, "' stored in object - reused")
    return(obj)
  }
  log_msg("  Computing UMAP '", umap_name, "' from '", from_reduction, "'...")
  recomputed <<- c(recomputed, umap_name)
  ndim <- ncol(Embeddings(obj, from_reduction))
  RunUMAP(obj,
          reduction      = from_reduction,
          dims           = seq_len(min(50L, ndim)),
          reduction.name = umap_name,
          reduction.key  = paste0(toupper(gsub("[^A-Za-z0-9]", "", umap_name)), "_"),
          n.neighbors    = n_neighbors,
          min.dist       = min_dist,
          verbose        = FALSE)
}

cel <- opt$lineage
step4_path <- file.path(opt$base, paste0(cel, "_Breast_merged_step4_harmony_dev015_nocca.RDS"))
lat_csv    <- file.path(opt$base, paste0(cel, "_Breast_merged_step6_scanvi_dev015_nocca_scanvi_latent.csv"))

log_msg("Loading ", step4_path)
obj <- readRDS(step4_path)
log_msg("  Cells: ", ncol(obj), " | Reductions: ", paste(names(obj@reductions), collapse = ", "))

obj <- ensure_umap(obj, "pca",        "umap.pca")
obj <- ensure_umap(obj, "mnn",        "umap.mnn")
obj <- ensure_umap(obj, "pca.combat", "umap.combat")
obj <- ensure_umap(obj, "harmony",    "umap.harmony")

lat_df <- read.csv(lat_csv, row.names = 1, stringsAsFactors = FALSE)
common <- intersect(colnames(obj), rownames(lat_df))
log_msg("  scANVI cells in common: ", length(common), " / ", ncol(obj))
stopifnot(length(common) >= 0.85 * ncol(obj))
obj     <- obj[, common]
lat_mat <- as.matrix(lat_df[common, , drop = FALSE])
colnames(lat_mat) <- paste0("SCANVI_", seq_len(ncol(lat_mat)))
obj[["scanvi"]] <- CreateDimReducObject(embeddings = lat_mat, key = "SCANVI_",
                                        assay = DefaultAssay(obj))
obj <- ensure_umap(obj, "scanvi", "umap.scanvi")

umaps <- c(Uncorrected = "umap.pca", MNN = "umap.mnn", ComBat = "umap.combat",
           Harmony = "umap.harmony", scANVI = "umap.scanvi")
df <- data.frame(cell                 = colnames(obj),
                 project_name         = as.character(obj$project_name),
                 scimilarity_celltype = as.character(obj$scimilarity_celltype),
                 Meta.scimilarity     = as.character(obj$Meta.scimilarity),
                 stringsAsFactors = FALSE)
for (m in names(umaps)) {
  e <- Embeddings(obj, umaps[[m]])
  df[[paste0(m, "_1")]] <- e[, 1]
  df[[paste0(m, "_2")]] <- e[, 2]
}

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(opt$out_dir, paste0(cel, "_Breast_umapdata.csv.gz"))
con <- gzfile(out, "w"); write.csv(df, con, row.names = FALSE); close(con)
writeLines(c(paste0("lineage\t", cel),
             paste0("n_cells\t", nrow(df)),
             paste0("n_projects\t", length(unique(df$project_name))),
             paste0("recomputed_umaps\t", if (length(recomputed)) paste(recomputed, collapse = ",") else "none")),
           file.path(opt$out_dir, paste0(cel, "_Breast_umapdata_info.tsv")))
log_msg("Saved ", out, " (recomputed: ", paste(recomputed, collapse = ","), ")")
