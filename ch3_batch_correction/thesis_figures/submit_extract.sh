#!/bin/bash
# Submit the comment-18 data-extraction jobs (one per object).
set -euo pipefail
D_HOST=/groups/ofircohen-group/users/hadasra/thesis_figures/comment18_2026-10-21
D_R=/ofir_grp/users/hadasra/thesis_figures/comment18_2026-10-21

while read -r lin mem; do
  sbatch --parsable --job-name="c18_umap_${lin}" --time=0-06:00:00 --mem=$mem --cpus-per-task=4 \
    --output="$D_HOST/logs/umap_${lin}_%j.out" --error="$D_HOST/logs/umap_${lin}_%j.err" \
    ~/rstudio/rstudio.sbatch --script "$D_R/scripts/extract_batch_umap_data.R" \
    --lineage "$lin" --out_dir "$D_R/data"
done <<LIST
B-cell 64G
Endothelial 64G
Fibroblasts 96G
Hematopoietic_Stem_and_Progenitor_cells 32G
Myeloid 128G
NK-cells 32G
T-cells 192G
LIST

while read -r lin np mem; do
  sbatch --parsable --job-name="c18_dot_${lin}" --time=0-03:00:00 --mem=$mem --cpus-per-task=2 \
    --output="$D_HOST/logs/dot_${lin}_%j.out" --error="$D_HOST/logs/dot_${lin}_%j.err" \
    ~/rstudio/rstudio.sbatch --script "$D_R/scripts/extract_dotplot_data.R" \
    --lineage "$lin" --n_panels "$np" --out_dir "$D_R/data"
done <<LIST
T-cells_Liver_Biliary 2 96G
Epithelial_Kidney 3 256G
LIST
