#!/bin/bash
# submit_gsea_hallmarks.sh
#
# Submit one SLURM job per lineage to run fgsea Hallmarks enrichment.
#
# Usage:
#   bash submit_gsea_hallmarks.sh [--dry-run]

set -euo pipefail

IN_BASE_R=/ofir_grp/users/hadasra/temp_objects/clustered_lineages
OUT_DIR_R=/ofir_grp/users/hadasra/temp_objects/clustered_lineages/FindMarkers_results/enrichment
SCRIPT_R=/ofir_grp/users/hadasra/create_sc_object_workflow/merge/gsea_hallmarks.R
SBATCH_TEMPLATE=~/rstudio/rstudio.sbatch
LOG_DIR=/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages/FindMarkers_results/enrichment/logs

DRY_RUN=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

mkdir -p "${LOG_DIR}"

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*"; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

# Collect all lineage names from marker files
LINEAGES=()
for f in /mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages/*_markers.csv; do
  base=$(basename "$f" _markers.csv)
  # skip *_top_markers.csv files (basename ends with "_top" after stripping _markers.csv)
  [[ "$base" == *_top ]] && continue
  LINEAGES+=("$base")
done

echo "Found ${#LINEAGES[@]} lineages to process"
echo ""

for lineage in "${LINEAGES[@]}"; do
  # Skip if output already exists
  out_file="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages/FindMarkers_results/enrichment/${lineage}_hallmarks_gsea.csv"
  if [[ -f "$out_file" ]]; then
    echo "  [SKIP] ${lineage} (output exists)"
    continue
  fi

  JOB=$(maybe_sbatch \
    --job-name="gsea_${lineage:0:22}" \
    --time=0-01:00:00 \
    --mem=16G \
    --cpus-per-task=4 \
    --output="${LOG_DIR}/gsea_${lineage}_%j.out" \
    --error="${LOG_DIR}/gsea_${lineage}_%j.err" \
    "${SBATCH_TEMPLATE}" \
      --script    "${SCRIPT_R}" \
      --lineage   "${lineage}" \
      --input_dir "${IN_BASE_R}" \
      --output_dir "${OUT_DIR_R}")

  echo "  [SUBMIT] ${lineage}  job=${JOB}"
done

echo ""
echo "Done submitting. Monitor with: squeue -u hadasra"
echo "Logs: ${LOG_DIR}/"
