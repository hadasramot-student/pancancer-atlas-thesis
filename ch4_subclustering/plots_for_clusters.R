# plots_for_clusters.R
#
# Ch.4 (Subclustering) — cluster x Scimilarity-celltype composition heatmap
# and z-scored profile plot, for one lineage object at resolution 2.0.
#
# Interactive/sourced script: expects `seurat_obj` (a clustered lineage
# object) to already be loaded.
#
# Usage:
#   seurat_obj <- readRDS("data/temp_objects/clustered_lineages/<lineage>_clustered.RDS")
#   source("plots_for_clusters.R")

library(pheatmap)
library(tidyverse)
table(seurat_obj@meta.data$scimilarity_celltype, seurat_obj@meta.data$RNA_snn_res.2) %>%
  as.data.frame() %>%
  pivot_wider(id_cols = Var2,names_from = Var1,values_from = Freq) %>% 
  dplyr::select(-Var2) %>% 
  pheatmap()
  
table(seurat_obj@meta.data$scimilarity_celltype, seurat_obj@meta.data$Meta.scimilarity) 


table(seurat_obj@meta.data$scimilarity_celltype, seurat_obj@meta.data$RNA_snn_res.2) %>% 
  as.data.frame() %>% 
  pivot_wider(id_cols = Var2,names_from = Var1,values_from = Freq) %>% 
  column_to_rownames("Var2") %>% 
  apply(MARGIN = 2,FUN = scale) %>%
  as.data.frame() %>% 
  mutate(cluster_new_num=row_number()) %>% 
  pivot_longer(cols = -cluster_new_num,names_to = "cell_type",values_to = "zscore") %>% 
  ggplot(aes(x=cluster_new_num,y = zscore,color = cell_type)) +
  geom_line() +
  geom_point()
