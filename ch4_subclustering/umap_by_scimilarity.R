# umap_by_scimilarity.R
#
# Ch.4 (Subclustering) — cluster-purity check against Scimilarity labels
# across the 5 Louvain/Leiden resolutions, for one lineage object.
#
# Interactive/sourced script: expects `seurat_obj` (a clustered lineage
# object, e.g. output of clustering_jobs/run_clustering.R +
# annotate_lineages.R) to already be loaded.
#
# Usage:
#   seurat_obj <- readRDS("data/temp_objects/clustered_lineages/Fibroblasts_clustered.RDS")
#   source("umap_by_scimilarity.R")

library(Seurat)
library(dplyr)
library(ggplot2)

# Your Seurat object
# seurat_obj <- fibro

# Exact columns
cluster_cols <- c("RNA_snn_res.0.25", "RNA_snn_res.0.1", "RNA_snn_res.0.5", 
                  "RNA_snn_res.0.75", "RNA_snn_res.1", "RNA_snn_res.2")

valid_cols <- cluster_cols[cluster_cols %in% colnames(seurat_obj@meta.data)]
# similarity_col <- "scimilarity_celltype"
similarity_col <- "Meta.scimilarity"

# SIMPLIFIED purity function - NO filter(total > 10) causing issues
calc_purity <- function(res_col) {
  clusters <- as.character(seurat_obj[[res_col]][[1]])
  cell_types <- seurat_obj[[similarity_col]][[1]]
  
  # Valid data only
  valid <- !is.na(clusters) & !is.na(cell_types) & clusters != ""
  clusters <- clusters[valid]
  cell_types <- cell_types[valid]
  
  # Table
  tbl <- table(clusters, cell_types)
  
  # Per cluster majority
  cluster_maj <- tbl %>%
    as.data.frame() %>%
    group_by(clusters) %>%
    summarise(
      total = sum(Freq),
      maj_pct = max(Freq / total * 100),
      .groups = "drop"
    )
  
  # Metrics
  n_clust <- nrow(cluster_maj)
  mean_maj <- mean(cluster_maj$maj_pct)
  pure50 <- sum(cluster_maj$maj_pct > 50)
  
  list(
    res = gsub("RNA_snn_res\\.", "", res_col),
    col = res_col,
    clusters = n_clust,
    mean_maj_pct = round(mean_maj, 1),
    pure50 = pure50,
    pure75 = sum(cluster_maj$maj_pct > 75)
  )
}

# Process all
purity_results <- do.call(rbind, lapply(valid_cols, function(col) {
  as.data.frame(calc_purity(col))
}))

print(purity_results)

# Sort and select optimal
purity_results <- purity_results[order(-purity_results$mean_maj_pct), ]
optimal_col <- purity_results$col[1]
optimal_res <- purity_results$res[1]

cat("\nOPTIMAL:", optimal_col, "(res", optimal_res, ") Mean maj:", purity_results$mean_maj_pct[1], "%\n")

# Set optimal
Idents(seurat_obj) <- seurat_obj[[optimal_col]][[1]]
seurat_obj$optimal_clusters <- Idents(seurat_obj)

# Use existing umap.mnn (your MNN UMAP)
# DefaultDimReduc(seurat_obj) <- "umap.mnn"

# Plots using umap.mnn
p1 <- DimPlot(seurat_obj, reduction = "umap.mnn", group.by = "project_name", 
              label = TRUE, label.size = 3, repel = TRUE)
p2 <- DimPlot(seurat_obj, reduction = "umap.mnn", group.by = similarity_col, 
              label = TRUE, label.size = 3, repel = TRUE)

print(p1)
print(p2)

saveRDS(seurat_obj, "fibro_mnn_optimal_final.rds")
ggsave("final_umap.png", p1, width = 10, height = 8)

cat("✅ Using umap.mnn - DONE!\n")


###___ Basic statistics____
similarity_col <- "scimilarity_celltype"

# 1. TOTAL cells and projects
total_cells <- ncol(seurat_obj)
cat("Total cells:", total_cells, "\n")

# Check for project/batch column (common names)
project_cols <- c( "project_name","sample")
project_col <- project_cols[project_cols %in% colnames(seurat_obj@meta.data)][1]
if(length(project_col) > 0) {
  n_projects <- length(unique(seurat_obj[[project_col]][[1]]))
  cat("Projects:", n_projects, " (column:", project_col, ")\n")
  print(table(seurat_obj[[project_col]][[1]]))
} else {
  cat("No project column found. Available: ", paste(head(colnames(seurat_obj@meta.data), 20), collapse = ", "), "\n")
}

# 2. Cells per Scimilarity type
table(seurat_obj@meta.data$Meta.scimilarity) %>% View()

# 3. Detailed table
scim_table <- data.frame(
  Cell_Type = names(scim_summary),
  Num_Cells = as.numeric(scim_summary),
  Pct_Total = round(as.numeric(scim_summary) / total_cells * 100, 1)
)

print("Scimilarity Summary:")
print(scim_table)

# 4. If project column exists, cells per type per project
if(exists("project_col") && length(project_col) > 0) {
  cat("\nCells per type per project:\n")
  cross_tab <- table(seurat_obj[[similarity_col]][[1]], seurat_obj[[project_col]][[1]])
  print(cross_tab)
}


#another try 

library(Seurat)
library(ggplot2)
library(dplyr)

# Use fewer dimensions first
mnn_dims <- ncol(Embeddings(seurat_obj, reduction = "mnn"))
cat("MNN dims available:", mnn_dims, "\n")
dims_use <- 1:min(30, mnn_dims)

# Recompute neighbors/clusters if needed
seurat_obj <- FindNeighbors(seurat_obj, reduction = "mnn", dims = dims_use)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Re-run UMAP with gentler settings
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "mnn",
  dims = dims_use,
  n.neighbors = 15,
  min.dist = 0.1,
  spread = 1
)

# Better colors for many cell types
cell_levels <- unique(seurat_obj[[similarity_col]][[1]])
cell_colors <- setNames(scales::hue_pal(l = 55, c = 100)(length(cell_levels)), cell_levels)

# Clean Scimilarity plot
p2 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = similarity_col,
  pt.size = 0.25,
  cols = cell_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +
  ggtitle("Scimilarity cell types") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

print(p2)
