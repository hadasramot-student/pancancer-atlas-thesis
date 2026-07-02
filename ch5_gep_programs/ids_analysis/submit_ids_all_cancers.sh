#!/bin/bash
# submit_ids_all_cancers.sh
#
# Ch.5 (GEP programs) — submits one IDS job per cancer type on SLURM. Each
# job scores GeneNMF meta-programs on the cNMF cells, then computes the
# IDS-max matrix between all cNMF + GeneNMF programs for that cancer.
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects ./submit_ids_all_cancers.sh   # submit all
#   ./submit_ids_all_cancers.sh --dry-run    # print only
#   ./submit_ids_all_cancers.sh --pvals      # also permutation p-values

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IDS_SCRIPT="${SCRIPT_DIR}/run_ids.py"
ATLAS_DATA_DIR="${ATLAS_DATA_DIR:-data/temp_objects}"
LOG_DIR="${ATLAS_DATA_DIR}/ids_output/logs"
PYTHON_ENV="scanvi-env"

CANCERS=(Lung Colorectal Kidney Liver_Biliary Ovarian Pancreas Prostate Sarcoma Skin Brain)

DRY_RUN=false
PVALS_FLAG=""
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=true ;;
    --pvals)   PVALS_FLAG="--pvals" ;;
  esac
done

mkdir -p "${LOG_DIR}"

maybe_sbatch() {
  if $DRY_RUN; then echo "[DRY-RUN] sbatch $*" >&2; echo "99999"; else sbatch --parsable "$@"; fi
}

for CANCER in "${CANCERS[@]}"; do
  echo "--- ${CANCER} ---"
  J=$(maybe_sbatch \
    --job-name="ids_${CANCER}" \
    --time="0-01:00:00" --mem="96G" --cpus-per-task=4 \
    --output="${LOG_DIR}/ids_${CANCER}-%j.out" \
    --error="${LOG_DIR}/ids_${CANCER}-%j.err" \
    --wrap="conda run -n ${PYTHON_ENV} python ${IDS_SCRIPT} --cancer ${CANCER} ${PVALS_FLAG}")
  echo "  IDS: job ${J}"
done

echo ""
echo "Submitted IDS jobs for ${#CANCERS[@]} cancers. Logs -> ${LOG_DIR}"
