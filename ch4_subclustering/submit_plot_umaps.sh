#!/bin/bash
# submit_plot_umaps.sh
#
# Re-generate UMAP plots for ALL clustered lineage objects.
# Creates two PNGs per object in clustered_lineages/umap_plots/:
#   {lineage}_UMAP_cluster.png      — louvain_res0.50 cluster IDs
#   {lineage}_UMAP_scimilarity.png  — scimilarity_cluster_label (majority-vote)
#
# Usage:
#   bash submit_plot_umaps.sh [--dependency <jobid1,jobid2,...>]

set -euo pipefail

IN_BASE_MNT=/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages
IN_BASE_R=/ofir_grp/users/hadasra/temp_objects/clustered_lineages
SCRIPT_R=/ofir_grp/users/hadasra/create_sc_object_workflow/merge/plot_umaps.R
SBATCH_R=~/rstudio/rstudio.sbatch
LOG_DIR="${IN_BASE_MNT}/umap_logs"

DRY_RUN=false
DEPENDENCY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)    DRY_RUN=true; shift ;;
    --dependency) DEPENDENCY="$2"; shift 2 ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done

mkdir -p "${LOG_DIR}"

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*" >&2; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

dep_flag=""
[[ -n "${DEPENDENCY}" ]] && dep_flag="--dependency=afterok:${DEPENDENCY}"

# Memory: just needs to load RDS + render; ~3-4x file size
pick_mem() {
  local bytes=$1
  if   (( bytes < 524288000   )); then echo "32G"
  elif (( bytes < 2147483648  )); then echo "64G"
  elif (( bytes < 5368709120  )); then echo "128G"
  else                                echo "256G"
  fi
}

echo "=============================================="
echo " UMAP re-plots — all clustered lineages"
echo " Dry run   : ${DRY_RUN}"
echo " Dependency: ${DEPENDENCY:-none}"
echo "=============================================="
echo ""

submitted=0

while IFS= read -r rds_path; do
  fname=$(basename "$rds_path")
  lineage="${fname/_clustered.RDS/}"

  bytes=$(stat -c%s "$rds_path" 2>/dev/null || echo 0)
  mem=$(pick_mem "$bytes")

  JOB=$(maybe_sbatch \
    ${dep_flag} \
    --job-name="umap_${lineage:0:22}" \
    --time="0-02:00:00" \
    --mem="${mem}" \
    --cpus-per-task=2 \
    --output="${LOG_DIR}/umap_${lineage}_%j.out" \
    --error="${LOG_DIR}/umap_${lineage}_%j.err" \
    "${SBATCH_R}" \
      --script    "${SCRIPT_R}" \
      --input_dir "${IN_BASE_R}" \
      --lineage   "${lineage}")

  echo "  [${mem}] ${lineage}  job=${JOB}"
  submitted=$((submitted + 1))

done < <(find "${IN_BASE_MNT}" -maxdepth 1 -name "*_clustered.RDS" | sort)

echo ""
echo "Submitted: ${submitted} jobs"
echo "Logs:      ${LOG_DIR}/"
echo "Output:    ${IN_BASE_MNT}/umap_plots/"
echo ""
echo "Monitor:   squeue -u hadasra | grep umap_"
