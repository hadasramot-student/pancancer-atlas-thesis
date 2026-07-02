# batch_correction_methods_summary.R
#
# Ch.3 (Batch correction) — cluster-purity summary across batch methods.
# For each reduction in `batch_methods` present on `seurat_obj`, computes,
# per clustering resolution, the majority-cell-type purity of each cluster
# (against Meta.scimilarity), then ranks methods by mean purity.
#
# Interactive/sourced script: expects `seurat_obj` (a Seurat object with
# harmony/mnn/pca.combat/integrated.cca/scanvi reductions and RNA_snn_res.*
# cluster columns already computed) to be loaded in the R session.
#
# Output: batch_methods_comparison.png (in the working directory)
#
# Usage:
#   seurat_obj <- readRDS("obj_with_all_batch_methods.RDS")
#   source("batch_correction_methods_summary.R")

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

similarity_col <- "Meta.scimilarity"

# YOUR batch methods (from reductions)
batch_methods <- c("harmony", "mnn", "pca.combat", "integrated.cca", "scanvi")
cluster_cols <- grep("RNA_snn_res\\.", colnames(seurat_obj@meta.data), value = TRUE)

cat("Testing methods:", paste(batch_methods, collapse = ", "), "\n")
cat("Resolutions:", length(cluster_cols), "\n")

# Purity per method/cluster
all_purities <- data.frame()

for(method in batch_methods) {
  if(method %in% names(seurat_obj@reductions)) {
    cat("Processing", method, "...\n")
    for(clust_col in cluster_cols) {
      # Cluster majority calculation
      clusters <- seurat_obj[[clust_col]][[1]]
      cell_types <- seurat_obj[[similarity_col]][[1]]
      
      valid <- !is.na(clusters) & !is.na(cell_types)
      clusters <- as.character(clusters[valid])
      cell_types <- cell_types[valid]
      
      purity_tbl <- table(clusters, cell_types) %>%
        as.data.frame() %>%
        group_by(clusters) %>%
        mutate(
          total = sum(Freq),
          pct = Freq / total * 100
        ) %>%
        slice_max(pct, n = 1) %>%
        ungroup() %>%
        select(cluster = clusters, majority_pct = pct, total_cells = total) %>%
        mutate(
          reduction = method,
          resolution = gsub("RNA_snn_res\\.", "", clust_col)
        )
      
      all_purities <- rbind(all_purities, purity_tbl)
    }
  } else {
    cat("⚠️", method, "not found\n")
  }
}

# Summary table
summary_stats <- all_purities %>%
  group_by(reduction, resolution) %>%
  summarise(
    mean_purity = round(mean(majority_pct), 1),
    median_purity = round(median(majority_pct), 1),
    pure50_clusters = sum(majority_pct > 50),
    total_clusters = n(),
    pure50_pct = round(pure50_clusters / total_clusters * 100, 1),
    .groups = "drop"
  ) %>%
  ungroup()

print("Method Performance Summary:")
print(summary_stats)

# BEST METHOD OVERALL
overall_best <- summary_stats %>%
  group_by(reduction) %>%
  summarise(avg_purity = mean(mean_purity), .groups = "drop") %>%
  arrange(desc(avg_purity))

cat("\n🏆 RANKING (avg mean purity %):\n")
print(overall_best)

# BOXPLOT 1: Purity distribution by method
p1 <- ggplot(all_purities, aes(x = reorder(reduction, majority_pct, median), y = majority_pct, fill = reduction)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +
  labs(title = "Cluster Purity by Batch Method", y = "Majority Cell Type %", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# BOXPLOT 2: % Pure clusters (>50%) by method
p2 <- ggplot(summary_stats, aes(x = reorder(reduction, pure50_pct, median), y = pure50_pct, fill = reduction)) +
  geom_col(alpha = 0.7, show.legend = FALSE) +
  labs(title = "% Clusters >50% Pure", y = "% Pure Clusters", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# BOXPLOT 3: Mean purity by resolution
p3 <- ggplot(summary_stats, aes(x = resolution, y = mean_purity, color = reduction, group = reduction)) +
  geom_point(size = 3) + geom_line() +
  labs(title = "Mean Purity vs Resolution", y = "Mean Purity %", x = "Resolution") +
  theme_minimal()

# Combined
combined <- (p1 | p2) / p3
print(combined)

ggsave("batch_methods_comparison.png", combined, width = 16, height = 12, dpi = 300)

cat("\n✅ Plots saved! Higher boxes = better batch correction.\n")