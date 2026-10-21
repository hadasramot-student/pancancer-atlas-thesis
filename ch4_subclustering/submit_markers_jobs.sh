#!/bin/bash
# submit_markers_jobs.sh
#
# Submit one FindAllMarkers SLURM job per *_clustered.RDS file.
# Memory is higher than annotation — FindAllMarkers is CPU/RAM intensive.
#
# Usage:
#   bash submit_markers_jobs.sh [--dry-run] [--dependency <jobid1,jobid2,...>]

set -euo pipefail

IN_BASE_PY=/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages
IN_BASE_R=/ofir_grp/users/hadasra/temp_objects/clustered_lineages
SCRIPT_R=/ofir_grp/users/hadasra/create_sc_object_workflow/merge/find_markers.R
SBATCH_R=~/rstudio/rstudio.sbatch
LOG_DIR="${IN_BASE_PY}/markers_logs"

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

# Memory: FindAllMarkers needs more RAM than annotation (holds expression matrix + test results)
# Roughly 3–4× the RDS file size
pick_mem() {
  local bytes=$1
  if   (( bytes < 524288000   )); then echo "64G"
  elif (( bytes < 2147483648  )); then echo "128G"
  elif (( bytes < 8589934592  )); then echo "256G"
  else                                echo "512G"
  fi
}

pick_time() {
  local bytes=$1
  if   (( bytes < 524288000   )); then echo "0-02:00:00"
  elif (( bytes < 2147483648  )); then echo "0-06:00:00"
  elif (( bytes < 8589934592  )); then echo "0-12:00:00"
  else                                echo "1-00:00:00"
  fi
}

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*"; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

dep_flag=""
[[ -n "${DEPENDENCY}" ]] && dep_flag="--dependency=afterok:${DEPENDENCY}"

echo "================================================="
echo " FindAllMarkers — all clustered lineages"
echo "================================================="
echo " Input     : ${IN_BASE_PY}"
echo " Dependency: ${DEPENDENCY:-none}"
echo " Dry run   : ${DRY_RUN}"
echo ""

submitted=0; skipped=0

while IFS= read -r rds_path; do
  fname=$(basename "$rds_path")
  lineage="${fname/_clustered.RDS/}"

  out_markers="${IN_BASE_PY}/${lineage}_markers.csv"
  if [[ -f "${out_markers}" ]]; then
    echo "  [DONE]   ${lineage}"
    skipped=$((skipped + 1))
    continue
  fi

  bytes=$(stat -c%s "$rds_path" 2>/dev/null || echo 0)
  mem=$(pick_mem "$bytes")
  wtime=$(pick_time "$bytes")

  JOB=$(maybe_sbatch \
    ${dep_flag} \
    --job-name="mk_${lineage:0:22}" \
    --time="${wtime}" \
    --mem="${mem}" \
    --cpus-per-task=8 \
    --output="${LOG_DIR}/mk_${lineage}_%j.out" \
    --error="${LOG_DIR}/mk_${lineage}_%j.err" \
    "${SBATCH_R}" \
      --script    "${SCRIPT_R}" \
      --input_dir "${IN_BASE_R}" \
      --lineage   "${lineage}" \
      --n_cores   8)

  echo "  [SUBMIT] ${lineage}  mem=${mem} time=${wtime}  job=${JOB}"
  submitted=$((submitted + 1))

done < <(find "${IN_BASE_PY}" -maxdepth 1 -name "*_clustered.RDS" | sort)

echo ""
echo "Submitted : ${submitted}"
echo "Skipped   : ${skipped} (already done)"
echo "Logs      : ${LOG_DIR}/"
echo ""
echo "Monitor  : squeue -u hadasra | grep mk_"
echo "Results  : ls ${IN_BASE_PY}/*_markers.csv 2>/dev/null | wc -l"
