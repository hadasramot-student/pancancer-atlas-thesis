#!/bin/bash
# ------------------------------------------------------------------
# Submit the epithelial batch-correction benchmark (Option A) for the
# 10 epithelial cancer types (Brain excluded: glioma, non-epithelial
# scimilarity labels; Head&Neck excluded: no scANVI object).
#
# Per cancer:  10_add_batch_methods.R  ->  Batch_Correction_Validation_Metrics.R
# ------------------------------------------------------------------
set -euo pipefail

WT=/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/epithelial-batch-correction-benchmark-f0c664
WT_HOST=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/epithelial-batch-correction-benchmark-f0c664
MB=/ofir_grp/users/hadasra/temp_objects/merged_by_query
OUT=/ofir_grp/users/hadasra/temp_objects/merged_by_query/epi_benchmark
OUT_HOST=/groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/epi_benchmark
LOGS_HOST=$WT_HOST/epi_benchmark/logs

ADD=$WT/epi_benchmark/scripts/10_add_batch_methods.R
VAL=/ofir_grp/users/hadasra/create_sc_object_workflow/merge/Batch_Correction_Validation_Metrics.R
SB=~/rstudio/rstudio.sbatch

mkdir -p "$OUT_HOST" "$LOGS_HOST"

# cancer | input basename (under $MB/<Cancer>/) | mem | time-hours
JOBS=(
  "Colorectal|Colorectal/Epithelial_Colorectal_merged_full_scanvi_final.RDS|512G|24"
  "Breast|Breast/Epithelial_Breast_subsample450k_merged_full_scanvi_final.RDS|512G|24"
  "Lung|Lung/Epithelial_Lung_merged_full_scanvi_final.RDS|256G|12"
  "Pancreas|Pancreas/Epithelial_Pancreas_merged_full_scanvi_final.RDS|256G|12"
  "Ovarian|Ovarian/Epithelial_Ovarian_merged_full_scanvi_final.RDS|256G|12"
  "Kidney|Kidney/Epithelial_Kidney_merged_full_scanvi_final.RDS|256G|12"
  "Skin|Skin/Epithelial_Skin_merged_full_scanvi_final.RDS|128G|6"
  "Liver_Biliary|Liver_Biliary/Epithelial_Liver_Biliary_merged_full_scanvi_final.RDS|128G|6"
  "Prostate|Prostate/Epithelial_Prostate_merged_full_scanvi_final.RDS|128G|6"
  # Sarcoma already run as pipeline test (jobs 20903675/76) -> Epithelial_Sarcoma_bench4.RDS + Sarcoma_validation/
)

echo "cancer,add_job,val_job" > "$LOGS_HOST/20_job_ids.csv"

for row in "${JOBS[@]}"; do
  IFS='|' read -r C REL MEM HRS <<< "$row"
  IN="$MB/$REL"
  BENCH="$OUT/Epithelial_${C}_bench4.RDS"
  VALDIR="$OUT/${C}_validation"

  JADD=$(sbatch --parsable \
    --job-name="epi_add_${C}" \
    --output="$LOGS_HOST/10_add_${C}.%j.out" \
    --error="$LOGS_HOST/10_add_${C}.%j.err" \
    --time="0-${HRS}:00:00" --mem="$MEM" --cpus-per-task=15 \
    "$SB" --script "$ADD" \
      --input_rds  "$IN" \
      --output_rds "$BENCH" \
      --batch_col  project_name \
      --label_col  scimilarity_celltype \
      --n_hvg 3000 --n_pcs 30 --min_batch 20)

  JVAL=$(sbatch --parsable \
    --job-name="epi_val_${C}" \
    --dependency=afterok:${JADD} \
    --output="$LOGS_HOST/30_val_${C}.%j.out" \
    --error="$LOGS_HOST/30_val_${C}.%j.err" \
    --time="0-04:00:00" --mem=128G --cpus-per-task=8 \
    "$SB" --script "$VAL" \
      --input_rds  "$BENCH" \
      --output_dir "$VALDIR" \
      --batch_col  project_name \
      --label_col  scimilarity_celltype \
      --n_cells 2000 --seed 42)

  echo "${C},${JADD},${JVAL}" | tee -a "$LOGS_HOST/20_job_ids.csv"
done

echo
echo "Submitted. Monitor: squeue -u hadasra"
echo "Job id map: $LOGS_HOST/20_job_ids.csv"
