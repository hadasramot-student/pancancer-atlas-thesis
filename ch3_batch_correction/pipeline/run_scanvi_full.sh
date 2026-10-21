#!/bin/bash

#============================================================
# scANVI-only pipeline on FULL (non-subsampled) objects
# Uses GPU (RTX 3090) for fast training
# Runs all 7 objects in parallel on separate GPUs
#
# Usage:
#   ./run_scanvi_full.sh
#============================================================

set -euo pipefail

BASE=/ofir_grp/users/hadasra/temp_objects/merged_by_query
SCRIPT_DIR=/ofir_grp/users/hadasra/create_sc_object_workflow/merge
SEURAT_TO_H5AD=/ofir_grp/users/hadasra/temp_objects/seurat_to_h5ad.R
SBATCH_R=~/rstudio/rstudio.sbatch
SBATCH_PY=~/rstudio/python.sbatch
TAG="full"

OBJECTS=(
  "Myeloid_Breast_merged"
  "B-cell_Breast_merged"
  "T-cells_Breast_merged"
  "NK-cells_Breast_merged"
  "Fibroblasts_Breast_merged"
  "Endothelial_Breast_merged"
  "Hematopoietic_Stem_and_Progenitor_cells_Breast_merged"
)

echo "=============================================="
echo "scANVI Full Pipeline - All 7 Objects"
echo "GPU: RTX 3090 (24GB VRAM)"
echo "Running in parallel (1 GPU per object)"
echo "=============================================="
echo ""

for obj in "${OBJECTS[@]}"; do
  INPUT="${BASE}/${obj}.RDS"
  echo "--- Submitting: $obj ---"

  # Step 5: Seurat → h5ad (CPU, fast)
  JOB5=$(sbatch --parsable \
    --job-name="f5_h5ad_${obj}" \
    --time=0-04:00:00 \
    ${SBATCH_R} \
      --script ${SEURAT_TO_H5AD} \
      --input_rds   "${INPUT}" \
      --output_h5ad "${BASE}/${obj}_${TAG}.h5ad")
  echo "  Step 5 (h5ad):   $JOB5"

  # Step 6: scANVI with GPU - no mem/cpu specified, let system assign
  JOB6=$(sbatch --parsable \
    --job-name="f6_scanvi_${obj}" \
    --dependency=afterok:${JOB5} \
    --time=1-00:00:00 \
    --gres=gpu:rtx_3090:1 \
    ${SBATCH_PY} \
      --script /mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/export_scanvi_to_csv.py \
      --h5ad              "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/${obj}_${TAG}.h5ad" \
      --output_prefix     "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/${obj}_${TAG}_scanvi" \
      --batch_key         project_name \
      --labels_key        scimilarity_celltype \
      --unlabeled_category Unknown \
      --layer_counts      counts \
      --n_latent          30 \
      --max_epochs_scvi   400 \
      --max_epochs_scanvi 200 \
      --label_key_out     scanvi_pred)
  echo "  Step 6 (scANVI): $JOB6"

  # Step 7: Import scANVI back to Seurat (CPU)
  JOB7=$(sbatch --parsable \
    --job-name="f7_seurat_${obj}" \
    --dependency=afterok:${JOB6} \
    --time=0-04:00:00 \
    ${SBATCH_R} \
      --script ${SCRIPT_DIR}/scanvi_into_seurat.R \
      --seurat      "${INPUT}" \
      --latent_csv  "/ofir_grp/users/hadasra/temp_objects/merged_by_query/${obj}_${TAG}_scanvi_scanvi_latent.csv" \
      --labels_csv  "/ofir_grp/users/hadasra/temp_objects/merged_by_query/${obj}_${TAG}_scanvi_scanvi_labels.csv" \
      --output      "${BASE}/${obj}_${TAG}_scanvi_final.RDS" \
      --label_col   scanvi_pred \
      --meta_col    scanvi_pred \
      --reduction   scanvi \
      --overwrite)
  echo "  Step 7 (Seurat): $JOB7"

  echo ""
done

echo "=============================================="
echo "All jobs submitted!"
echo "Monitor: squeue -u hadasra | grep -E 'f[567]_'"
echo ""
echo "Expected timeline:"
echo "  Step 5 (h5ad):   ~10-30 min per object (parallel)"
echo "  Step 6 (scANVI): ~3-5 hrs per object on RTX 3090 (parallel)"
echo "  Step 7 (Seurat): ~10-20 min per object (parallel)"
echo "  TOTAL:           ~4-6 hours for ALL 7 objects"
echo "=============================================="
