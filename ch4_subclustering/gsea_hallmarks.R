#!/usr/bin/env Rscript
# gsea_hallmarks.R
#
# Run fgsea Hallmarks enrichment on FindMarkers results for one lineage.
#
# For each cluster in the lineage's _markers.csv:
#   - Build a gene ranking by avg_log2FC (all genes passing thresholds)
#   - Run fgsea against MSigDB H collection (Hallmarks)
#   - Collect results across all clusters into one table
#
# Output: <output_dir>/<lineage>_hallmarks_gsea.csv
#
# Usage:
#   Rscript gsea_hallmarks.R --lineage T-cells_Breast \
#     --input_dir /path/to/temp_objects/clustered_lineages \
#     --output_dir /path/to/temp_objects/clustered_lineages/FindMarkers_results/enrichment

if (nzchar(Sys.getenv("R_USER_LIBRARY"))) .libPaths(c(Sys.getenv("R_USER_LIBRARY"), .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--lineage",    type="character", help="Lineage name, e.g. T-cells_Breast"),
  make_option("--input_dir",  type="character",
              default=file.path(Sys.getenv("ATLAS_DATA_DIR", unset="data/temp_objects"), "clustered_lineages")),
  make_option("--output_dir", type="character",
              default=file.path(Sys.getenv("ATLAS_DATA_DIR", unset="data/temp_objects"), "clustered_lineages/FindMarkers_results/enrichment"))
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$lineage)) stop("--lineage is required")

# msigdbr is in the user library; fgsea is in the container's system library

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
})

lineage    <- opt$lineage
input_dir  <- opt$input_dir
output_dir <- opt$output_dir

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

markers_file <- file.path(input_dir, paste0(lineage, "_markers.csv"))
if (!file.exists(markers_file)) stop("Markers file not found: ", markers_file)

message("[", lineage, "] Loading markers: ", markers_file)
markers <- read.csv(markers_file, stringsAsFactors=FALSE)

# Hallmarks gene sets (human)
h_sets <- msigdbr(species="Homo sapiens", collection="H")
pathways <- split(h_sets$gene_symbol, h_sets$gs_name)
message("[", lineage, "] Loaded ", length(pathways), " Hallmark gene sets")

clusters <- sort(unique(markers$cluster))
message("[", lineage, "] Running fgsea for ", length(clusters), " clusters")

results_list <- lapply(clusters, function(cl) {
  sub <- markers[markers$cluster == cl, ]
  if (nrow(sub) < 5) {
    message("  Cluster ", cl, ": too few genes (", nrow(sub), "), skipping")
    return(NULL)
  }
  # Named vector: gene -> avg_log2FC, deduplicated (keep max abs)
  sub <- sub[order(abs(sub$avg_log2FC), decreasing=TRUE), ]
  sub <- sub[!duplicated(sub$gene), ]
  ranks <- setNames(sub$avg_log2FC, sub$gene)
  ranks <- sort(ranks, decreasing=TRUE)

  score_type <- if (all(ranks >= 0)) "pos" else if (all(ranks <= 0)) "neg" else "std"
  set.seed(42)
  res <- tryCatch(
    fgsea(pathways=pathways, stats=ranks, minSize=10, maxSize=500, nPermSimple=10000, scoreType=score_type),
    error=function(e) { message("  Cluster ", cl, " fgsea error: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)

  res$cluster  <- cl
  res$lineage  <- lineage
  res$n_genes_ranked <- length(ranks)
  res
})

results <- bind_rows(Filter(Negate(is.null), results_list))

if (nrow(results) == 0) {
  message("[", lineage, "] No enrichment results produced.")
  quit(status=0)
}

# Flatten leadingEdge list column to semicolon-separated string
results$leadingEdge <- sapply(results$leadingEdge, paste, collapse=";")

# Reorder columns
results <- results[, c("lineage","cluster","pathway","pval","padj","ES","NES",
                        "size","n_genes_ranked","leadingEdge")]

out_file <- file.path(output_dir, paste0(lineage, "_hallmarks_gsea.csv"))
write.csv(results, out_file, row.names=FALSE)
message("[", lineage, "] Saved ", nrow(results), " rows -> ", out_file)
