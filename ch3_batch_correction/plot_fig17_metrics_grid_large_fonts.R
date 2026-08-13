#!/usr/bin/env Rscript
# plot_fig17_metrics_grid_large_fonts.R
#
# Recreates thesis Fig. 17 ("Batch correction method comparison across five
# metrics and seven cell types in the breast cancer benchmark cohort") directly
# from the combined validation CSV, with large typography baked in from the
# start for legibility at 6in print display width (no smaller-font baseline
# version exists to derive this from — the original rendering script was not
# found in the repo, so this rebuilds the figure from the underlying data,
# which is unchanged).
#
# Layout: 5 metric rows x 7 cell-type columns, each panel a 4-bar chart
# (MNN | ComBat | Harmony | scANVI), bars ordered within each panel from
# best- to worst-performing (accounting for each metric's higher/lower-is-
# better direction), with the best-performing bar per panel given a bold
# border. Same data/plot type as described; only typography (axis/tick/
# legend/strip text) is sized up ~1.5-2x relative to panel size vs. a
# "standard" ggplot facet grid at this size.
#
# Input:  ATLAS_DATA_DIR/merged_by_query/all_validation_results.csv
#         (columns: cell_type, method, reduction, batch_lisi, label_lisi,
#          kbet_accept, graph_connectivity, f1_label)
# Output: ATLAS_DATA_DIR/merged_by_query/batch_correction_comparison_plots/
#           all_celltypes_Breast_fig17_metrics_grid_large_fonts.{pdf,png}
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript plot_fig17_metrics_grid_large_fonts.R

suppressPackageStartupMessages({
  library(ggplot2)
})

BASE    <- file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "merged_by_query")
OUT_DIR <- file.path(BASE, "batch_correction_comparison_plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

in_csv <- file.path(BASE, "all_validation_results.csv")
log_msg("Reading: ", in_csv)
df <- read.csv(in_csv, stringsAsFactors = FALSE)

# ---------------------------------------------------------------
# Ordering / display labels (kept consistent with the other Ch.3
# method-comparison figures: MNN, ComBat, Harmony, scANVI column order;
# same 7-cell-type set as plot_batch_correction_comparison.R)
# ---------------------------------------------------------------
method_levels <- c("mnn", "combat", "harmony", "scanvi")
method_labels <- c(mnn = "MNN", combat = "ComBat", harmony = "Harmony", scanvi = "scANVI")
method_colors <- c(MNN = "#377EB8", ComBat = "#4DAF4A", Harmony = "#FF7F00", scANVI = "#E41A1C")

celltype_levels <- c(
  "B-cell", "Endothelial", "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells",
  "Myeloid", "NK-cells", "T-cells"
)
celltype_labels <- c(
  "B-cell" = "B-cell", "Endothelial" = "Endothelial", "Fibroblasts" = "Fibroblasts",
  "Hematopoietic_Stem_and_Progenitor_cells" = "HSPCs",
  "Myeloid" = "Myeloid", "NK-cells" = "NK-cells", "T-cells" = "T-cells"
)

metric_levels <- c("batch_lisi", "label_lisi", "kbet_accept", "graph_connectivity", "f1_label")
metric_labels <- c(
  batch_lisi         = "Batch LISI\n(higher = better)",
  label_lisi         = "Label LISI\n(lower = better)",
  kbet_accept        = "kBET Acceptance\n(higher = better)",
  graph_connectivity = "Graph Connectivity\n(higher = better)",
  f1_label           = "F1 Label Score\n(higher = better)"
)
# Direction of "better" per metric, used both to rank bars best->worst
# within a panel and to flag which bar gets the bold "winner" border.
metric_direction <- c(
  batch_lisi = "higher", label_lisi = "lower", kbet_accept = "higher",
  graph_connectivity = "higher", f1_label = "higher"
)

df$method    <- factor(df$method, levels = method_levels, labels = method_labels[method_levels])
df$cell_type <- factor(df$cell_type, levels = celltype_levels, labels = celltype_labels[celltype_levels])

long <- do.call(rbind, lapply(metric_levels, function(m) {
  sub <- data.frame(
    cell_type = df$cell_type,
    method    = df$method,
    metric    = factor(metric_labels[[m]], levels = metric_labels[metric_levels]),
    value     = df[[m]]
  )

  direction   <- metric_direction[[m]]
  sub$rank    <- NA_integer_
  sub$is_best <- FALSE

  # Rank + best-method flag computed independently within each cell-type
  # panel (4 methods per panel), so bars can be sorted best -> worst per
  # panel regardless of the fixed method order used for fill/legend.
  for (ct in levels(sub$cell_type)) {
    idx  <- which(sub$cell_type == ct)
    vals <- sub$value[idx]
    ord  <- if (direction == "higher") {
      order(-vals, is.na(vals))
    } else {
      order(vals, is.na(vals))
    }
    sub$rank[idx[ord]] <- seq_along(idx)

    if (!all(is.na(vals))) {
      best_val <- if (direction == "higher") max(vals, na.rm = TRUE) else min(vals, na.rm = TRUE)
      sub$is_best[idx] <- !is.na(vals) & vals == best_val
    }
  }
  sub
}))

log_msg(sprintf("Panels: %d metrics x %d cell types x %d methods (%d rows, %d NA)",
                length(metric_levels), length(celltype_levels), length(method_levels),
                nrow(long), sum(is.na(long$value))))

# ---------------------------------------------------------------
# Build the 5 x 7 bar-chart grid
# ---------------------------------------------------------------
fig <- ggplot(long, aes(x = factor(rank), y = value, fill = method, linewidth = is_best)) +
  geom_col(colour = "black", width = 0.75, na.rm = TRUE) +
  facet_grid(metric ~ cell_type, scales = "free_y", switch = "y") +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_linewidth_manual(values = c(`TRUE` = 2.0, `FALSE` = 0.35), guide = "none") +
  labs(
    title    = "Batch correction method comparison across five metrics and seven cell types",
    subtitle = "Breast cancer benchmark cohort  |  bars ranked best -> worst per panel; bold border = best method",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 26) +
  theme(
    plot.title      = element_text(face = "bold", size = 32, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle   = element_text(size = 22, hjust = 0.5, colour = "grey35", margin = margin(b = 10)),
    strip.text.x    = element_text(size = 22, face = "bold", margin = margin(t = 6, b = 6)),
    strip.text.y.left = element_text(size = 20, face = "bold", angle = 0, hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", colour = NA),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.text.y     = element_text(size = 17),
    axis.title      = element_blank(),
    panel.spacing   = unit(0.6, "lines"),
    legend.position = "bottom",
    legend.title    = element_text(size = 24, face = "bold"),
    legend.text     = element_text(size = 22),
    legend.key.size = unit(1.1, "cm"),
    plot.margin     = margin(12, 14, 12, 14)
  )

# ---------------------------------------------------------------
# Save — physical canvas sized so 300dpi export meets/exceeds the
# ~5972x4115px native resolution of the original figure.
# ---------------------------------------------------------------
out_png <- file.path(OUT_DIR, "all_celltypes_Breast_fig17_metrics_grid_large_fonts.png")
out_pdf <- file.path(OUT_DIR, "all_celltypes_Breast_fig17_metrics_grid_large_fonts.pdf")

ggsave(out_png, fig, width = 20.5, height = 14.1, dpi = 300, bg = "white")
ggsave(out_pdf, fig, width = 20.5, height = 14.1, device = pdf, useDingbats = FALSE)

log_msg("Saved: ", out_png)
log_msg("Saved: ", out_pdf)
