#!/bin/bash

#============================================================
# Fast Batch Correction Pipeline (NO CCA)
# Runs: fastMNN → ComBat → Harmony → scANVI → validation
# Skips CCA for speed (6-8 hours instead of 30+)
#
# Usage:
#   ./run_batch_correction_no_cca.sh "Myeloid_Breast_merged" 0.15
#
# Or submit multiple in parallel:
#   for CELL_TYPE in "Myeloid_Breast_merged" "B_cell_Breast_merged"; do
#     ./run_batch_correction_no_cca.sh "$CELL_TYPE" 0.15 &
#   done
#   wait
#============================================================

set -euo pipefail

# ============================================================
# Configuration
# ============================================================
# Note: Use /ofir_grp/ for paths in SLURM/docker container
BASE=/ofir_grp/users/hadasra/temp_objects/merged_by_query
SCRIPT_DIR=/ofir_grp/users/hadasra/create_sc_object_workflow/merge
SEURAT_TO_H5AD=/ofir_grp/users/hadasra/temp_objects/seurat_to_h5ad.R
SBATCH_R=~/rstudio/rstudio.sbatch
SBATCH_PY=~/rstudio/python.sbatch

# Parse arguments
OBJECT_NAME="${1:?Usage: $0 OBJECT_NAME [FRACTION]}"
FRACTION="${2:-0.15}"

# Derive file paths
INPUT_FULL="${BASE}/${OBJECT_NAME}.RDS"
TAG="dev${FRACTION/./}"  # e.g., "dev015" for fraction 0.15

# Verify input exists (check from login node before submitting)
if [[ ! -f "$INPUT_FULL" ]]; then
  # Try alternate path on login node
  ALT_INPUT="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/${OBJECT_NAME}.RDS"
  if [[ -f "$ALT_INPUT" ]]; then
    echo "Note: Found input at alternate path, but will use /ofir_grp in container"
  else
    echo "ERROR: Input file not found: $OBJECT_NAME"
    exit 1
  fi
fi

echo "=========================================="
echo "NO-CCA Pipeline: ${OBJECT_NAME}"
echo "Fraction: ${FRACTION}"
echo "Tag:      ${TAG}"
echo "Methods:  fastMNN → ComBat → Harmony → scANVI"
echo "=========================================="

# ============================================================
# Step 0 — Create fraction dev subset
# ============================================================
JOB0=$(sbatch --parsable \
  --job-name="s0_subset_${OBJECT_NAME}" \
  --time=0-02:00:00 \
  --mem=128G \
  --cpus-per-task=8 \
  ${SBATCH_R} \
    --script ${SCRIPT_DIR}/create_dev_subset.R \
    --input_rds  "${INPUT_FULL}" \
    --output_rds "${BASE}/${OBJECT_NAME}_step0_${TAG}.RDS" \
    --batch_col  project_name \
    --label_col  scimilarity_celltype \
    --fraction   ${FRACTION} \
    --min_cells  10)
echo "Step 0 (subset):   $JOB0"

# ============================================================
# Step 2 — fastMNN (skipping CCA)
# ============================================================
JOB2=$(sbatch --parsable \
  --job-name="s2_mnn_${OBJECT_NAME}" \
  --dependency=afterok:${JOB0} \
  --time=0-12:00:00 \
  --mem=128G \
  --cpus-per-task=15 \
  ${SBATCH_R} \
    --script ${SCRIPT_DIR}/batch2_fastmnn.R \
    --input_rds  "${BASE}/${OBJECT_NAME}_step0_${TAG}.RDS" \
    --output_rds "${BASE}/${OBJECT_NAME}_step2_mnn_${TAG}_nocca.RDS" \
    --batch_col  project_name)
echo "Step 2 (MNN):      $JOB2"

# ============================================================
# Step 3 — ComBat
# ============================================================
JOB3=$(sbatch --parsable \
  --job-name="s3_combat_${OBJECT_NAME}" \
  --dependency=afterok:${JOB2} \
  --time=0-06:00:00 \
  --mem=128G \
  --cpus-per-task=15 \
  ${SBATCH_R} \
    --script ${SCRIPT_DIR}/batch3_combat.R \
    --input_rds  "${BASE}/${OBJECT_NAME}_step2_mnn_${TAG}_nocca.RDS" \
    --output_rds "${BASE}/${OBJECT_NAME}_step3_combat_${TAG}_nocca.RDS" \
    --batch_col  project_name)
echo "Step 3 (ComBat):   $JOB3"

# ============================================================
# Step 4 — Harmony
# ============================================================
JOB4=$(sbatch --parsable \
  --job-name="s4_harmony_${OBJECT_NAME}" \
  --dependency=afterok:${JOB3} \
  --time=0-06:00:00 \
  --mem=128G \
  --cpus-per-task=15 \
  ${SBATCH_R} \
    --script ${SCRIPT_DIR}/batch4_harmony.R \
    --input_rds  "${BASE}/${OBJECT_NAME}_step3_combat_${TAG}_nocca.RDS" \
    --output_rds "${BASE}/${OBJECT_NAME}_step4_harmony_${TAG}_nocca.RDS" \
    --batch_col  project_name)
echo "Step 4 (Harmony):  $JOB4"

# ============================================================
# Step 5 — Convert to h5ad
# ============================================================
JOB5=$(sbatch --parsable \
  --job-name="s5_h5ad_${OBJECT_NAME}" \
  --dependency=afterok:${JOB4} \
  --time=0-04:00:00 \
  --mem=128G \
  --cpus-per-task=8 \
  ${SBATCH_R} \
    --script ${SEURAT_TO_H5AD} \
    --input_rds   "${BASE}/${OBJECT_NAME}_step4_harmony_${TAG}_nocca.RDS" \
    --output_h5ad "${BASE}/${OBJECT_NAME}_step5_h5ad_${TAG}_nocca.h5ad")
echo "Step 5 (h5ad):     $JOB5"

# ============================================================
# Step 6 — scANVI (Python)
# ============================================================
JOB6=$(sbatch --parsable \
  --job-name="s6_scanvi_${OBJECT_NAME}" \
  --dependency=afterok:${JOB5} \
  --time=1-00:00:00 \
  --mem=64G \
  --cpus-per-task=8 \
  ${SBATCH_PY} \
    --script /mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/export_scanvi_to_csv.py \
    --h5ad              "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/${OBJECT_NAME}_step5_h5ad_${TAG}_nocca.h5ad" \
    --output_prefix     "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/${OBJECT_NAME}_step6_scanvi_${TAG}_nocca" \
    --batch_key         project_name \
    --labels_key        scimilarity_celltype \
    --unlabeled_category Unknown \
    --layer_counts      counts \
    --n_latent          30 \
    --max_epochs_scvi   200 \
    --max_epochs_scanvi 100 \
    --label_key_out     scanvi_pred)
echo "Step 6 (scANVI):   $JOB6"

# ============================================================
# Step 7 — Add scANVI back to Seurat
# ============================================================
JOB7=$(sbatch --parsable \
  --job-name="s7_scanvi_seurat_${OBJECT_NAME}" \
  --dependency=afterok:${JOB6} \
  --time=0-04:00:00 \
  --mem=128G \
  --cpus-per-task=8 \
  ${SBATCH_R} \
    --script ${SCRIPT_DIR}/scanvi_into_seurat.R \
    --seurat      "${BASE}/${OBJECT_NAME}_step4_harmony_${TAG}_nocca.RDS" \
    --latent_csv  "${BASE}/${OBJECT_NAME}_step6_scanvi_${TAG}_nocca_scanvi_latent.csv" \
    --labels_csv  "${BASE}/${OBJECT_NAME}_step6_scanvi_${TAG}_nocca_scanvi_labels.csv" \
    --output      "${BASE}/${OBJECT_NAME}_step7_final_${TAG}_nocca.RDS" \
    --label_col   scanvi_pred \
    --meta_col    scanvi_pred \
    --reduction   scanvi \
    --overwrite)
echo "Step 7 (->Seurat): $JOB7"

# ============================================================
# Step 8 — Validate
# ============================================================
JOB8=$(sbatch --parsable \
  --job-name="s8_validate_${OBJECT_NAME}" \
  --dependency=afterok:${JOB7} \
  --time=0-06:00:00 \
  --mem=128G \
  --cpus-per-task=8 \
  ${SBATCH_R} \
    --script ${SCRIPT_DIR}/Batch_Correction_Validation_Metrics.R \
    --input_rds  "${BASE}/${OBJECT_NAME}_step7_final_${TAG}_nocca.RDS" \
    --output_dir "${BASE}/${OBJECT_NAME}_validation_${TAG}_nocca" \
    --batch_col  project_name \
    --label_col  scimilarity_celltype \
    --n_cells    2000)
echo "Step 8 (validate): $JOB8"

# ============================================================
# Summary
# ============================================================
echo ""
echo "======================================"
echo "NO-CCA Pipeline submitted: ${OBJECT_NAME}"
echo "Fraction: ${FRACTION}, Tag: ${TAG}"
echo "======================================"
echo "Input:  ${INPUT_FULL}"
echo "Output: ${BASE}/${OBJECT_NAME}_step7_final_${TAG}_nocca.RDS"
echo ""
echo "Timeline: ~6-8 hours (vs 30+ hours with CCA)"
echo ""
echo "Monitor: squeue -u hadasra"
echo "         squeue -j $JOB8  (final job)"
echo ""
