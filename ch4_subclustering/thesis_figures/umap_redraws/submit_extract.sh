#!/bin/bash
# Submit one extraction job per object (memory sized from RDS file size).
set -euo pipefail
D_HOST=/groups/ofircohen-group/users/hadasra/thesis_figures/comment10_2026-10-19
D_R=/ofir_grp/users/hadasra/thesis_figures/comment10_2026-10-19
IN_HOST=/groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages
IN_R=/ofir_grp/users/hadasra/temp_objects/clustered_lineages
while read -r lin name; do
  bytes=$(stat -c%s "$IN_HOST/${lin}_clustered.RDS")
  if (( bytes > 5000000000 )); then mem=256G; elif (( bytes > 1000000000 )); then mem=96G; else mem=48G; fi
  job=$(sbatch --parsable --job-name="ex_${name#umap_}" --time=0-03:00:00 --mem=$mem --cpus-per-task=2 \
    --output="$D_HOST/logs/extract_${name}_%j.out" --error="$D_HOST/logs/extract_${name}_%j.err" \
    ~/rstudio/rstudio.sbatch --script "$D_R/extract_umap_data.R" \
    --input_dir "$IN_R" --lineage "$lin" --out_name "$name" --out_dir "$D_R/data")
  echo "$name [$mem] job=$job"
done <<LIST
T-cells_Liver_Biliary umap_Tcells_LiverBiliary
Epithelial_Kidney umap_Epithelial_Kidney
Myeloid_Prostate umap_Myeloid_Prostate
B-cell_Liver_Biliary umap_Bcells_LiverBiliary
Fibroblasts_Brain umap_Fibroblasts_Brain
Endothelial_Kidney umap_Endothelial_Kidney
NK-cells_Skin umap_NKcells_Skin
Hematopoietic_Stem_and_Progenitor_cells_Breast umap_HSPC_Breast
LIST
