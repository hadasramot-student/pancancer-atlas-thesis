#!/usr/bin/env Rscript
# elbow_check.R
#
# Load 8 lineage Seurat objects, extract PCA stdev, save ElbowPlots + summary CSV.
# Output: temp_objects/clustered_lineages/elbow_plots/

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

IN_DIR  <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages"
OUT_DIR <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages/elbow_plots"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

lineages <- c(
  "T-cells_Liver_Biliary",
  "Myeloid_Prostate",
  "Epithelial_Kidney",
  "B-cell_Liver_Biliary",
  "Fibroblasts_Brain",
  "Endothelial_Kidney",
  "NK-cells_Skin",
  "Hematopoietic_Stem_and_Progenitor_cells_Breast"
)

summary_rows <- list()

for (lin in lineages) {
  rds <- file.path(IN_DIR, paste0(lin, "_clustered.RDS"))
  message("\n[", lin, "] Loading...")
  obj <- readRDS(rds)

  # Find PCA reduction name
  pca_name <- grep("^pca", names(obj@reductions), value=TRUE, ignore.case=TRUE)[1]
  if (is.na(pca_name)) {
    message("  No PCA reduction found, skipping")
    next
  }

  stdev <- obj[[pca_name]]@stdev
  ndims <- min(30, length(stdev))
  stdev30 <- stdev[seq_len(ndims)]
  pcs     <- seq_len(ndims)

  # Detect elbow: largest drop in stdev
  diffs      <- -diff(stdev30)
  elbow_pc   <- which.max(diffs) + 1
  # Also check if stdev is still declining at PC30
  last_drop  <- stdev30[ndims-1] - stdev30[ndims]
  first_drop <- stdev30[1] - stdev30[2]
  pct_remaining <- round(last_drop / first_drop * 100, 1)
  stabilized <- pct_remaining < 5   # <5% of initial drop = essentially flat

  message(sprintf("  PCA: %s | ndims available: %d | elbow ~PC%d | drop at PC30: %.1f%% of initial",
                  pca_name, length(stdev), elbow_pc, pct_remaining))
  message("  Stabilized before PC30: ", stabilized)

  # ElbowPlot
  df <- data.frame(PC=pcs, stdev=stdev30)
  p <- ggplot(df, aes(x=PC, y=stdev)) +
    geom_line(color="steelblue") +
    geom_point(color="steelblue", size=2) +
    geom_vline(xintercept=elbow_pc, linetype="dashed", color="firebrick", linewidth=0.8) +
    annotate("text", x=elbow_pc+0.5, y=max(stdev30)*0.95,
             label=paste0("elbow ~PC", elbow_pc), hjust=0, color="firebrick", size=3.5) +
    labs(title=lin, x="PC", y="Std Dev") +
    theme_classic(base_size=13)

  png(file.path(OUT_DIR, paste0(lin, "_elbow.png")), width=700, height=450)
  print(p)
  dev.off()

  summary_rows[[lin]] <- data.frame(
    lineage          = lin,
    pca_reduction    = pca_name,
    n_pcs_available  = length(stdev),
    elbow_pc         = elbow_pc,
    stdev_pc1        = round(stdev30[1], 3),
    stdev_pc30       = round(stdev30[ndims], 3),
    pct_drop_at_pc30 = pct_remaining,
    stabilized_by_pc30 = stabilized
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR, "elbow_summary.csv"), row.names=FALSE)
message("\nDone. Summary saved to: ", file.path(OUT_DIR, "elbow_summary.csv"))
print(summary_df[, c("lineage","elbow_pc","pct_drop_at_pc30","stabilized_by_pc30")])
