#!/bin/bash
# Submit PCA+UMAP+clustering for every scANVI-corrected final object.
# Each object = one cancer type × one lineage, already batch-corrected with scANVI.
# Clustering uses the scANVI latent embedding (stored as reduction "scanvi").
#
# Output: /ofir_grp/users/hadasra/temp_objects/clustered_lineages/<lineage>_<cancer>_clustered.RDS
#         /mnt/.../clustered_lineages/umap_plots/<lineage>_<cancer>_UMAP_*.png
#         /mnt/.../clustered_lineages/<lineage>_<cancer>_cluster_summary.tsv

set -euo pipefail

SCRIPT=/ofir_grp/users/hadasra/create_sc_object_workflow/merge/cluster_lineage.R
SBATCH_R=~/rstudio/rstudio.sbatch
IN_BASE=/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query
OUT_BASE=/ofir_grp/users/hadasra/temp_objects/clustered_lineages
SUM_BASE=/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages
LOG_DIR="${SUM_BASE}/logs"

mkdir -p "${SUM_BASE}" "${LOG_DIR}"

# Memory tier by file size: <200M→32G, <1G→64G, <3G→128G, >=3G→256G
pick_mem() {
  local bytes=$1
  if   (( bytes < 209715200  )); then echo "32G"
  elif (( bytes < 1073741824 )); then echo "64G"
  elif (( bytes < 3221225472 )); then echo "128G"
  else                                echo "256G"
  fi
}

pick_time() {
  local bytes=$1
  if   (( bytes < 209715200  )); then echo "0-01:00:00"
  elif (( bytes < 1073741824 )); then echo "0-02:00:00"
  elif (( bytes < 3221225472 )); then echo "0-06:00:00"
  else                                echo "0-12:00:00"
  fi
}

echo "=============================================="
echo "Submitting clustering jobs for all scANVI final objects"
echo "=============================================="

submitted=0
skipped=0

while IFS= read -r rds_path; do
  fname=$(basename "$rds_path")
  # e.g. T-cells_Lung_merged_full_scanvi_final.RDS
  tag="${fname/_merged_full_scanvi_final.RDS/}"   # T-cells_Lung
  lineage=$(echo "$tag" | rev | cut -d_ -f2- | rev)  # T-cells  (all but last word)
  cancer=$(echo "$tag" | rev | cut -d_ -f1 | rev)    # Lung

  # Rebuild properly: lineage = everything before last underscore-separated cancer word
  # Use the directory name as cancer type (more reliable)
  dir=$(dirname "$rds_path")
  if [[ "$(basename $dir)" == "merged_by_query" ]]; then
    # Top-level (Breast objects stored at root of merged_by_query)
    cancer_tag="Breast"
  else
    cancer_tag="$(basename $dir)"
  fi
  lineage_tag="${fname/_${cancer_tag}_merged_full_scanvi_final.RDS/}"
  job_name="${lineage_tag}_${cancer_tag}"

  out_rds="${OUT_BASE}/${job_name}_clustered.RDS"
  sum_tsv="${OUT_BASE}/${job_name}_cluster_summary.tsv"   # must use /ofir_grp/ — container doesn't bind /mnt/new_groups/

  # Skip if already done
  out_rds_check="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages/${job_name}_clustered.RDS"
  if [[ -f "${out_rds_check}" ]]; then
    echo "  [DONE]       ${job_name} — already exists, skipping"
    skipped=$((skipped + 1))
    continue
  fi

  bytes=$(stat -c%s "$rds_path" 2>/dev/null || echo 0)
  mem=$(pick_mem "$bytes")
  wtime=$(pick_time "$bytes")

  JOB=$(sbatch --parsable \
    --job-name="cl_${job_name:0:20}" \
    --time="${wtime}" \
    --mem="${mem}" \
    --cpus-per-task=8 \
    --output="${LOG_DIR}/cl_${job_name}_%j.out" \
    --error="${LOG_DIR}/cl_${job_name}_%j.err" \
    "${SBATCH_R}" \
      --script  "${SCRIPT}" \
      --input   "/ofir_grp/users/hadasra/temp_objects/merged_by_query/${rds_path#${IN_BASE}/}" \
      --output  "${out_rds}" \
      --lineage "${job_name}" \
      --n_pcs   30 \
      --min_dist 0.3 \
      --summary_out "${sum_tsv}")

  echo "  [${mem} / ${wtime}] ${job_name} -> job ${JOB}"
  submitted=$((submitted + 1))

done < <(find "${IN_BASE}" -name "*_full_scanvi_final.RDS" | grep -v "/other_" | sort)

echo ""
echo "Submitted: ${submitted} jobs  (skipped: ${skipped} already done)"
echo "Logs:      ${LOG_DIR}/"
echo "Summaries: ${SUM_BASE}/*_cluster_summary.tsv"
echo "Outputs:   ${SUM_BASE}/"
echo ""
echo "When done, aggregate results:"
echo "  head -1 \$(ls ${SUM_BASE}/*_cluster_summary.tsv | head -1) > ${SUM_BASE}/all_summaries.tsv"
echo "  tail -n+2 -q ${SUM_BASE}/*_cluster_summary.tsv >> ${SUM_BASE}/all_summaries.tsv"
echo "  column -t ${SUM_BASE}/all_summaries.tsv"
