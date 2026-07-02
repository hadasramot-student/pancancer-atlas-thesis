#!/usr/bin/env Rscript

#' compare_batch_methods.R
#'
#' Ch.3 (Batch correction) — cross-object method comparison. Reads the
#' per-object validation CSVs written by batch_correction_validation_metrics.R
#' and creates:
#' 1. Summary comparison table
#' 2. Ranking of methods
#' 3. Heatmap visualization
#'
#' Usage:
#'   Rscript compare_batch_methods.R \
#'     --base_dir /path/to/temp_objects/merged_by_query \
#'     --pattern "*/validation/batch_correction_metrics.csv" \
#'     --output_dir ./comparison

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--base_dir"), type="character",
              default=Sys.getenv("ATLAS_DATA_DIR_MERGED", unset="data/temp_objects/merged_by_query"),
              help="Base directory containing validation folders"),
  make_option(c("--output_dir"), type="character", default="./comparison",
              help="Output directory for results [default=./comparison]"),
  make_option(c("--tag"), type="character", default="dev30",
              help="Tag for validation folders (e.g., dev30, dev50) [default=dev30]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Find and read all validation CSVs
# ============================================================
log_msg("Searching for validation results...")

validation_pattern <- sprintf("*_validation_%s/batch_correction_metrics.csv", opt$tag)
csv_files <- Sys.glob(file.path(opt$base_dir, validation_pattern))

if (length(csv_files) == 0) {
  log_msg("ERROR: No validation CSV files found matching pattern:")
  log_msg(sprintf("  %s/%s", opt$base_dir, validation_pattern))
  stop("No validation files found")
}

log_msg(sprintf("Found %d validation files:", length(csv_files)))
for (f in csv_files) log_msg(sprintf("  %s", basename(dirname(f))))

# Read and combine
all_data <- do.call(rbind, lapply(csv_files, function(f) {
  object_name <- sub(sprintf("_validation_%s", opt$tag), "",
                     basename(dirname(f)))
  df <- read.csv(f)
  df$object <- object_name
  df
}))

log_msg("")
log_msg(sprintf("Combined data: %d methods × %d objects",
                length(unique(all_data$method)),
                length(unique(all_data$object))))

# ============================================================
# Summary statistics per method (across all objects)
# ============================================================
log_msg("")
log_msg("METHOD RANKINGS (across all objects)")
log_msg("====================================")

summary_stats <- data.frame()

for (metric in c("batch_lisi", "label_lisi", "kbet_accept",
                 "graph_connectivity", "f1_label")) {
  if (!metric %in% colnames(all_data)) next

  # Determine if higher or lower is better
  better_higher <- !grepl("lisi", metric) || metric == "batch_lisi"

  method_means <- tapply(all_data[[metric]], all_data$method, mean, na.rm = TRUE)
  method_sds   <- tapply(all_data[[metric]], all_data$method, sd, na.rm = TRUE)

  if (better_higher) {
    ranking <- rank(-method_means)  # higher is better
  } else {
    ranking <- rank(method_means)   # lower is better
  }

  temp <- data.frame(
    method    = names(method_means),
    metric    = metric,
    mean      = as.numeric(method_means),
    sd        = as.numeric(method_sds),
    rank      = as.numeric(ranking),
    direction = if(better_higher) "↑ higher" else "↓ lower"
  )
  summary_stats <- rbind(summary_stats, temp)
}

# Pretty print
for (metric in unique(summary_stats$metric)) {
  sub_data <- summary_stats[summary_stats$metric == metric, ]
  sub_data <- sub_data[order(sub_data$rank), ]

  cat(sprintf("\n%s (%s is better):\n", metric, sub_data$direction[1]))
  cat(sprintf("  %-10s %8s %8s\n", "Method", "Mean", "SD"))
  cat(strrep("-", 35), "\n")
  for (i in seq_len(nrow(sub_data))) {
    cat(sprintf("  %-10s %8.4f %8.4f  ← Rank %d\n",
                sub_data$method[i], sub_data$mean[i], sub_data$sd[i],
                sub_data$rank[i]))
  }
}

# ============================================================
# Overall ranking (average rank across all metrics)
# ============================================================
log_msg("")
log_msg("OVERALL METHOD RANKING")
log_msg("======================")

overall_rank <- aggregate(rank ~ method, summary_stats, mean)
overall_rank <- overall_rank[order(overall_rank$rank), ]
overall_rank$overall_rank <- seq_len(nrow(overall_rank))

cat("\n")
for (i in seq_len(nrow(overall_rank))) {
  cat(sprintf("%d. %s (avg rank: %.2f)\n",
              overall_rank$overall_rank[i],
              overall_rank$method[i],
              overall_rank$rank[i]))
}

# ============================================================
# Save detailed results
# ============================================================
output_file <- file.path(opt$output_dir, "method_comparison.csv")
write.csv(all_data, output_file, row.names = FALSE)
log_msg(paste("\nDetailed results saved:", output_file))

output_file2 <- file.path(opt$output_dir, "method_summary.csv")
write.csv(summary_stats, output_file2, row.names = FALSE)
log_msg(paste("Summary table saved:", output_file2))

output_file3 <- file.path(opt$output_dir, "overall_ranking.csv")
write.csv(overall_rank, output_file3, row.names = FALSE)
log_msg(paste("Rankings saved:", output_file3))

# ============================================================
# Heatmap visualization
# ============================================================
log_msg("")
log_msg("Creating heatmap...")

# Normalize metrics for visualization
viz_data <- all_data
for (metric in c("batch_lisi", "kbet_accept", "graph_connectivity", "f1_label")) {
  if (metric %in% colnames(viz_data)) {
    min_val <- min(viz_data[[metric]], na.rm = TRUE)
    max_val <- max(viz_data[[metric]], na.rm = TRUE)
    viz_data[[metric]] <- (viz_data[[metric]] - min_val) / (max_val - min_val)
  }
}

# For label_lisi, invert so higher is better
if ("label_lisi" %in% colnames(viz_data)) {
  min_val <- min(viz_data$label_lisi, na.rm = TRUE)
  max_val <- max(viz_data$label_lisi, na.rm = TRUE)
  viz_data$label_lisi <- 1 - ((viz_data$label_lisi - min_val) / (max_val - min_val))
}

# Average across metrics and objects for heatmap
heatmap_data <- aggregate(cbind(batch_lisi, label_lisi, kbet_accept,
                                graph_connectivity, f1_label) ~ method,
                          data = viz_data,
                          mean,
                          na.rm = TRUE)

rownames(heatmap_data) <- heatmap_data$method
heatmap_data$method <- NULL

# Create heatmap using base R
pdf(file.path(opt$output_dir, "batch_methods_heatmap.pdf"), width = 8, height = 6)
heatmap(as.matrix(heatmap_data),
        main = "Batch Correction Methods Comparison\n(Normalized scores, higher = better)",
        scale = "none", margins = c(10, 10))
dev.off()
log_msg(paste("Heatmap saved:", file.path(opt$output_dir, "batch_methods_heatmap.pdf")))

log_msg("")
log_msg("COMPARISON COMPLETE")
log_msg("==================")
log_msg("")
log_msg("Recommendation:")
best_method <- overall_rank$method[1]
log_msg(sprintf("  Based on aggregated metrics, '%s' ranks highest.", best_method))
log_msg("")
log_msg("Next steps:")
log_msg(sprintf("  1. Review detailed metrics in: %s", output_file2))
log_msg(sprintf("  2. Check per-object results in: %s", output_file))
log_msg(sprintf("  3. View visualization in: %s", file.path(opt$output_dir, "batch_methods_heatmap.pdf")))
log_msg("")
