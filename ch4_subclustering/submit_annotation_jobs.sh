#!/bin/bash
# submit_annotation_jobs.sh
#
# Submit one SLURM job per lineage for SCimilarity annotation.
# Memory is set per-object based on file size.
#
# Usage:
#   bash submit_annotation_jobs.sh [--input_dir <path>] [--output_dir <path>] [--dry-run]
#
# Defaults:
#   input_dir  : /mnt/.../temp_objects/clustered_lineages
#   output_dir : /mnt/.../temp_objects/annotated_lineages

set -euo pipefail

IN_BASE_PY=/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages
IN_BASE_R=/ofir_grp/users/hadasra/temp_objects/clustered_lineages

SCRIPT_R=/ofir_grp/users/hadasra/create_sc_object_workflow/merge/annotate_lineages.R
SBATCH_R=~/rstudio/rstudio.sbatch
LOG_DIR="${IN_BASE_PY}/annotation_logs"

DRY_RUN=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input_dir)  IN_BASE_PY="$2"; IN_BASE_R="$2"; shift 2 ;;
    --dry-run)    DRY_RUN=true; shift ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done

mkdir -p "${LOG_DIR}"

pick_mem() {
  local bytes=$1
  if   (( bytes < 524288000   )); then echo "32G"
  elif (( bytes < 2147483648  )); then echo "64G"
  elif (( bytes < 5368709120  )); then echo "128G"
  else                                echo "256G"
  fi
}

pick_time() {
  local bytes=$1
  if   (( bytes < 524288000   )); then echo "0-01:30:00"
  elif (( bytes < 2147483648  )); then echo "0-03:00:00"
  elif (( bytes < 5368709120  )); then echo "0-06:00:00"
  else                                echo "0-12:00:00"
  fi
}

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*"; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

echo "================================================="
echo " SCimilarity annotation — all clustered lineages"
echo "================================================="
echo " Input : ${IN_BASE_PY}"
echo " Dry   : ${DRY_RUN}"
echo ""

submitted=0; skipped=0

while IFS= read -r rds_path; do
  fname=$(basename "$rds_path")
  lineage="${fname/_clustered.RDS/}"
  lineage="${lineage/.RDS/}"

  # Skip if summary TSV already exists (means annotation already ran)
  out_tsv="${IN_BASE_PY}/${lineage}_scimilarity_summary.tsv"
  if [[ -f "${out_tsv}" ]]; then
    echo "  [DONE]   ${lineage}"
    skipped=$((skipped + 1))
    continue
  fi

  bytes=$(stat -c%s "$rds_path" 2>/dev/null || echo 0)
  mem=$(pick_mem "$bytes")
  wtime=$(pick_time "$bytes")

  JOB=$(maybe_sbatch \
    --job-name="annot_${lineage:0:22}" \
    --time="${wtime}" \
    --mem="${mem}" \
    --cpus-per-task=8 \
    --output="${LOG_DIR}/annot_${lineage}_%j.out" \
    --error="${LOG_DIR}/annot_${lineage}_%j.err" \
    "${SBATCH_R}" \
      --script    "${SCRIPT_R}" \
      --input_dir "${IN_BASE_R}" \
      --lineage   "${lineage}" \
      --overwrite)

  echo "  [SUBMIT] ${lineage}  mem=${mem} time=${wtime}  job=${JOB}"
  submitted=$((submitted + 1))

done < <(find "${IN_BASE_PY}" -maxdepth 1 -name "*_clustered.RDS" | sort)

echo ""
echo "Submitted : ${submitted}"
echo "Skipped   : ${skipped} (already done)"
echo "Logs      : ${LOG_DIR}/"
echo ""
echo "Monitor: squeue -u hadasra | grep annot_"
echo "Results: ls ${IN_BASE_PY}/*_scimilarity_summary.tsv 2>/dev/null | wc -l"
