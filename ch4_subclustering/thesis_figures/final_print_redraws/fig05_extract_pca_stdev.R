# Extract the stored PCA standard deviations (first 30 PCs) of the
# T-cells Liver/Biliary clustered object for the elbow plot (Figure A.2).
# Same object and reduction as ch4_subclustering/elbow_check.R; read-only.
obj <- readRDS("/ofir_grp/users/hadasra/temp_objects/clustered_lineages/T-cells_Liver_Biliary_clustered.RDS")
pca_name <- if ("pca" %in% names(obj@reductions)) "pca" else grep("pca", names(obj@reductions), value = TRUE)[1]
sd <- obj[[pca_name]]@stdev
n <- min(30, length(sd))
write.csv(data.frame(PC = seq_len(n), stdev = sd[seq_len(n)]),
          "/ofir_grp/users/hadasra/thesis_figures/final_figure_fixes_2026-10-21/data/Tcells_LiverBiliary_pca_stdev.csv",
          row.names = FALSE)
cat("pca reduction:", pca_name, " n stdev:", length(sd), "\n")
