#!/usr/bin/env Rscript

#_________Subsetting script___________
# Strategy: joint batch x cell-subtype stratification
#   - Proportional sampling (--fraction) per stratum
#   - Soft floor (--min_cells): strata below the floor are kept whole + warned
#   - Optional hard cap (--max_cells): total budget allocated proportionally
#     across batches before subtype floor is applied

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option(c("--input_rds"),  type="character", default=NULL),
  make_option(c("--output_rds"), type="character", default=NULL),
  make_option(c("--fraction"),   type="double",    default=0.8,
              help="Fraction of cells to keep per stratum [default=0.8]"),
  make_option(c("--batch_col"),  type="character", default="project_name",
              help="Metadata column for batch [default=project_name]"),
  make_option(c("--label_col"),  type="character", default=NULL,
              help="Optional metadata column for cell subtype (joint stratification)"),
  make_option(c("--min_cells"),  type="integer",   default=10,
              help="Min cells to keep per stratum; strata below this are kept whole [default=10]"),
  make_option(c("--max_cells"),  type="integer",   default=NULL,
              help="Optional hard cap on total output cells (overrides --fraction if smaller)"),
  make_option(c("--seed"),       type="integer",   default=42)
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$input_rds))  stop("Please provide --input_rds")
if (is.null(opt$output_rds)) stop("Please provide --output_rds")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

log_msg(paste("Loading:", opt$input_rds))
obj <- readRDS(opt$input_rds)
log_msg(sprintf("Original: %d cells, %d genes", ncol(obj), nrow(obj)))

set.seed(opt$seed)

batch_vec <- obj[[opt$batch_col]][, 1]

# Build strata: batch alone or batch x cell subtype
if (!is.null(opt$label_col)) {
  label_vec <- obj[[opt$label_col]][, 1]
  strata_vec <- paste(batch_vec, label_vec, sep = "___")
  log_msg(sprintf("Stratifying by '%s' x '%s' (%d strata)",
                  opt$batch_col, opt$label_col, length(unique(strata_vec))))
} else {
  strata_vec <- batch_vec
  log_msg(sprintf("Stratifying by '%s' only (%d strata)",
                  opt$batch_col, length(unique(strata_vec))))
}

# If --max_cells is set, compute per-batch fraction to hit the budget,
# then apply whichever fraction is smaller (--fraction or budget fraction)
if (!is.null(opt$max_cells) && opt$max_cells < ncol(obj)) {
  batch_sizes   <- table(batch_vec)
  batch_budgets <- round(opt$max_cells * batch_sizes / ncol(obj))
  budget_fracs  <- setNames(
    as.numeric(batch_budgets) / as.numeric(batch_sizes),
    names(batch_sizes)
  )
  log_msg(sprintf("--max_cells=%d active; effective per-batch fractions range %.2f-%.2f",
                  opt$max_cells, min(budget_fracs), max(budget_fracs)))
} else {
  budget_fracs <- NULL
}

# Sample each stratum
sparse_strata <- character(0)

cells_keep <- unlist(lapply(unique(strata_vec), function(s) {
  s_cells <- colnames(obj)[strata_vec == s]
  b        <- batch_vec[strata_vec == s][1]

  frac <- if (!is.null(budget_fracs)) {
    min(opt$fraction, budget_fracs[b])
  } else {
    opt$fraction
  }

  n_target <- floor(length(s_cells) * frac)

  if (n_target < opt$min_cells) {
    sparse_strata <<- c(sparse_strata, s)
    n_target <- length(s_cells)   # keep whole stratum, never inflate
  }

  sample(s_cells, min(n_target, length(s_cells)))
}))

if (length(sparse_strata) > 0) {
  warning(sprintf(
    "%d sparse strata kept whole (below min_cells=%d): %s",
    length(sparse_strata), opt$min_cells,
    paste(sparse_strata, collapse=", ")
  ))
}

obj_sub <- obj[, cells_keep]
log_msg(sprintf("Subset: %d cells (%.0f%% of original)",
                ncol(obj_sub), 100 * ncol(obj_sub) / ncol(obj)))

log_msg("Cells per batch after subsetting:")
print(table(obj_sub[[opt$batch_col]][, 1]))

if (!is.null(opt$label_col)) {
  log_msg("Cells per batch x subtype after subsetting:")
  print(table(obj_sub[[opt$batch_col]][, 1], obj_sub[[opt$label_col]][, 1]))
}

log_msg(paste("Saving:", opt$output_rds))
saveRDS(obj_sub, opt$output_rds)
log_msg("Done")

# Example usage:
#
# Batch-only (original behaviour):
# sbatch rstudio.sbatch \
#   --script .../create_dev_subset.R \
#   --input_rds  .../Myeloid_Breast_merged.RDS \
#   --output_rds .../Myeloid_Breast_merged_dev.RDS \
#   --batch_col project_name \
#   --fraction 0.7
#
# Joint batch x subtype, hard cap at 50k cells:
# sbatch rstudio.sbatch \
#   --script .../create_dev_subset.R \
#   --input_rds  .../Myeloid_Breast_merged.RDS \
#   --output_rds .../Myeloid_Breast_merged_dev.RDS \
#   --batch_col project_name \
#   --label_col cell_type \
#   --max_cells 50000 \
#   --min_cells 10
