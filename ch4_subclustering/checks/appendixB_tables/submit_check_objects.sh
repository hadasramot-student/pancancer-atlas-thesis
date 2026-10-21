#!/bin/bash
# submit_check_objects.sh — one SLURM job per object (check_object.R), so each finished
# object's output is saved on its own. Skips objects whose outputs already exist.
#
# Usage: bash submit_check_objects.sh [object ...]     (default: all 75 in objects.tsv)

set -euo pipefail

HOST_BASE=/groups/ofircohen-group/users/hadasra
CONT_BASE=/ofir_grp/users/hadasra
CHECK_DIR=thesis_checks/appendixB_2026-10-21
CL_DIR=temp_objects/clustered_lineages
OUT_SUB=$CHECK_DIR/object_checks
SBATCH_R=~/rstudio/rstudio.sbatch

mkdir -p "$HOST_BASE/$OUT_SUB" "$HOST_BASE/$CHECK_DIR/logs"

pick_mem() {
  local bytes=$1
  if   (( bytes < 524288000   )); then echo "32G"
  elif (( bytes < 2147483648  )); then echo "64G"
  elif (( bytes < 8589934592  )); then echo "128G"
  else                                echo "256G"
  fi
}

if [[ $# -gt 0 ]]; then wanted=" $* "; else wanted=""; fi

tail -n +2 "$HOST_BASE/$CHECK_DIR/scripts/objects.tsv" | while IFS=$'\t' read -r sheet obj res v88; do
  [[ -n "$wanted" && "$wanted" != *" $obj "* ]] && continue
  if [[ -f "$HOST_BASE/$OUT_SUB/${obj}_marker_check.tsv" ]]; then
    echo "skip $obj (done)"; continue
  fi
  bytes=$(stat -c %s "$HOST_BASE/$CL_DIR/${obj}_clustered.RDS")
  mem=$(pick_mem "$bytes")
  sbatch --job-name "abchk_${obj}" --mem "$mem" --cpus-per-task 2 --time 0-06:00:00 \
    --output "$HOST_BASE/$CHECK_DIR/logs/%x-%j.out" --error "$HOST_BASE/$CHECK_DIR/logs/%x-%j.err" \
    "$SBATCH_R" --script "$CONT_BASE/$CHECK_DIR/scripts/check_object.R" \
    "$obj" "$res" "$CONT_BASE/$CL_DIR" "$CONT_BASE/$OUT_SUB"
done
