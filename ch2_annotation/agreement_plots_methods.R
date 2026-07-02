# agreement_plots_methods.R
#
# Ch.2 (Annotation) — pairwise agreement/contradiction analysis between the
# 3 annotation methods (SingleR, SCType, Scimilarity) at the meta-type
# level, ahead of the maj2/maj3/maj6 consensus vote (see
# add_consensuses() in common/functions.R). For each method, computes the
# % of cells where the other two methods agree vs. contradict it, per
# meta-cell-type and per project.
#
# Interactive/sourced script: expects a metadata data.frame (e.g. from
# rbind_tsvs.R's all_metadata_merged.tsv) with Meta.SingleR, Meta.Sctype,
# Meta.scimilarity columns, assigned to `df` before use.
#
# Usage:
#   df <- readr::read_tsv("data/temp_objects/all_metadata_merged.tsv")
#   source("agreement_plots_methods.R")

library(dplyr)
library(tidyr)
library(ggplot2)

methods <- c("Meta.SingleR", "Meta.Sctype", "Meta.scimilarity")

calc_agree_contradict <- function(df, focal_method, study_col = "project_name") {
  
  other_methods <- setdiff(methods, focal_method)
  m1 <- focal_method
  m2 <- other_methods[1]
  m3 <- other_methods[2]
  
  df %>%
    dplyr::mutate(
      focal_label = .data[[m1]],
      other1 = .data[[m2]],
      other2 = .data[[m3]],
      status = dplyr::case_when(
        is.na(focal_label) | is.na(other1) | is.na(other2) ~ NA_character_,
        focal_label == other1 | focal_label == other2 ~ "agreement",
        other1 == other2 & focal_label != other1 ~ "contradicted",
        focal_label != other1 & focal_label != other2 & other1 != other2 ~ "three_way",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::group_by(method = m1, project = .data[[study_col]]) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_agreement = sum(status == "agreement", na.rm = TRUE),
      n_contradicted = sum(status == "contradicted", na.rm = TRUE),
      n_three_way = sum(status == "three_way", na.rm = TRUE),
      pct_agreement = 100 * n_agreement / n_total,
      pct_contradicted = 100 * n_contradicted / n_total,
      pct_three_way = 100 * n_three_way / n_total,
      .groups = "drop"
    )
}

agree_contradict_df <- dplyr::bind_rows(lapply(methods, function(m) {
  calc_agree_contradict(all_cells, focal_method = m, study_col = "project_name")
}))

plot_counts_df <- agree_contradict_df %>%
  dplyr::select(method, project, n_agreement, n_contradicted) %>%
  tidyr::pivot_longer(
    cols = c(n_agreement, n_contradicted),
    names_to = "type",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    type = dplyr::recode(
      type,
      "n_agreement" = "Agreement",
      "n_contradicted" = "Contradicted"
    )
  )

p_counts <- ggplot(plot_counts_df, aes(x = project, y = count, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(~ method, ncol = 1, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Project",
    y = "Number of cells",
    fill = "Status",
    title = "Agreement and contradiction counts per project for each annotation method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  )

p_counts

plot_pct_df <- agree_contradict_df %>%
  dplyr::select(method, project, pct_agreement, pct_contradicted) %>%
  tidyr::pivot_longer(
    cols = c(pct_agreement, pct_contradicted),
    names_to = "type",
    values_to = "percent"
  ) %>%
  dplyr::mutate(
    type = dplyr::recode(
      type,
      "pct_agreement" = "Agreement",
      "pct_contradicted" = "Contradicted"
    )
  )
p_pct <- ggplot(plot_pct_df, aes(x = project, y = percent, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(~ method, ncol = 1, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Project",
    y = "Percent of cells",
    fill = "Status",
    title = "Agreement and contradiction percentages per project for each annotation method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  )

p_pct


project_key <- agree_contradict_df %>%
  dplyr::distinct(project) %>%
  dplyr::arrange(project) %>%
  dplyr::mutate(project_num = row_number())

plot_pct_df_num <- plot_pct_df %>%
  dplyr::left_join(project_key, by = "project") %>%
  dplyr::mutate(project_num = factor(project_num, levels = project_key$project_num))


p_pct_num <- ggplot(plot_pct_df_num, aes(x = project_num, y = percent, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(~ method, ncol = 1) +
  theme_classic() +
  labs(
    x = "Project number",
    y = "Percent of cells",
    fill = "Status",
    title = "Agreement and contradiction percentages per project for each annotation method"
  ) +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1)
  )

p_pct_num

celltype_cols <- c(
  "Cancer/Epithelial cells" = "#D62728",
  "T-cells" = "#2CA02C",
  "B-cell" = "#1F77B4",
  "NK-cells" = "#17BECF",
  "Myeloid" = "#9467BD",
  "Fibroblasts/Stroma" = "#8C564B",
  "Endothelial/Vascular" = "#FF7F0E",
  "HSPC/Stem" = "#E377C2",
  "other" = "#B0B0B0"
)
celltype_map <- c(
  "Epithelial" = "Cancer/Epithelial cells",
  "Cancer/Epithelial cells" = "Cancer/Epithelial cells",
  "T-cells" = "T-cells",
  "B-cell" = "B-cell",
  "NK-cells" = "NK-cells",
  "Myeloid" = "Myeloid",
  "Fibroblasts" = "Fibroblasts/Stroma",
  "Fibroblasts/Stroma" = "Fibroblasts/Stroma",
  "Endothelial" = "Endothelial/Vascular",
  "Endothelial/Vascular" = "Endothelial/Vascular",
  "Hematopoietic_Stem_and_Progenitor_cells" = "HSPC/Stem",
  "HSPC/Stem" = "HSPC/Stem",
  "other" = "other"
)
calc_one_vs_two_by_celltype <- function(df, focal_method,
                                        study_col = "project_name",
                                        celltype_col = "Consensus_Cell_Type") {
  
  other_methods <- setdiff(methods, focal_method)
  m1 <- focal_method
  m2 <- other_methods[1]
  m3 <- other_methods[2]
  
  df %>%
    dplyr::mutate(
      focal_label = .data[[m1]],
      other1 = .data[[m2]],
      other2 = .data[[m3]],
      meta_cell_type = .data[[celltype_col]],
      status = dplyr::case_when(
        is.na(focal_label) | is.na(other1) | is.na(other2) ~ NA_character_,
        focal_label == other1 | focal_label == other2 ~ "agreement",
        other1 == other2 & focal_label != other1 ~ "contradicted",
        focal_label != other1 & focal_label != other2 & other1 != other2 ~ "three_way",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::group_by(
      method = m1,
      project = .data[[study_col]],
      meta_cell_type
    ) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_agreement = sum(status == "agreement", na.rm = TRUE),
      n_contradicted = sum(status == "contradicted", na.rm = TRUE),
      n_three_way = sum(status == "three_way", na.rm = TRUE),
      denom = n_agreement + n_contradicted,
      pct_not_contradicted = ifelse(denom > 0, 100 * n_agreement / denom, NA_real_),
      .groups = "drop"
    )
}
agree_df <- dplyr::bind_rows(lapply(methods, function(m) {
  calc_one_vs_two_by_celltype(
    df = all_cells,
    focal_method = m,
    study_col = "project_name",
    celltype_col = "Consensus_Cell_Type"
  )
}))
agree_df$meta_cell_type_std <- dplyr::recode(
  agree_df$meta_cell_type,
  !!!celltype_map
)

agree_df$meta_cell_type_std[is.na(agree_df$meta_cell_type_std)] <- "other"

agree_df <- agree_df %>%
  dplyr::filter(n_total >= 50)

median_rank <- agree_df %>%
  dplyr::group_by(method, meta_cell_type_std) %>%
  dplyr::summarise(
    median_value = median(pct_not_contradicted, na.rm = TRUE),
    .groups = "drop"
  )
agree_df_ranked <- agree_df %>%
  dplyr::left_join(median_rank, by = c("method", "meta_cell_type_std"))

p_box <- ggplot(
  agree_df_ranked,
  aes(
    x = reorder(meta_cell_type_std, median_value),
    y = pct_not_contradicted,
    fill = meta_cell_type_std
  )
) +
  geom_boxplot(outlier.size = 0.4, width = 0.7) +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  theme_classic() +
  labs(
    x = "Meta-cell-type",
    y = "% not contradicted by the other two methods",
    title = "Agreement distribution by meta-cell-type across projects"
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 11)
  )

p_box
