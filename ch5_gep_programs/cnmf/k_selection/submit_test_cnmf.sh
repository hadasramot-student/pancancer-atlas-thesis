#!/bin/bash
# submit_test_cnmf.sh
#
# Test-run cNMF for a single cancer type (default: Lung)
# - Downsampled to 15K cells
# - K = 10, 15, 20, 25, 30
# - 100 replicates (vs 200 in full run)
# - Expected runtime: ~1-2h total
#
# Usage:
#   ./submit_test_cnmf.sh              # Lung, default
#   ./submit_test_cnmf.sh Colorectal   # other cancer type
#   ./submit_test_cnmf.sh --dry-run    # print only

set -euo pipefail

CANCER="${1:-Lung}"
SCRIPT="/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/cnmf/test_run_cnmf.py"
LOG_DIR="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_test/logs"
PYTHON_ENV="scanvi-env"

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true && CANCER="Lung"

mkdir -p "${LOG_DIR}"

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*" >&2; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

echo "=============================================="
echo " cNMF TEST RUN"
echo " Cancer: ${CANCER}"
echo " Cells: 15K downsampled"
echo " K: 10, 15, 20, 25, 30"
echo " Replicates: 100"
echo " Dry run: ${DRY_RUN}"
echo "=============================================="

# ── Step 1: Prepare ────────────────────────────────────────────────────────
J1=$(maybe_sbatch \
  --job-name="cnmf_test_prepare_${CANCER}" \
  --time="0-02:00:00" --mem="64G" --cpus-per-task=4 \
  --output="${LOG_DIR}/test_prepare_${CANCER}-%j.out" \
  --error="${LOG_DIR}/test_prepare_${CANCER}-%j.err" \
  --wrap="conda run -n ${PYTHON_ENV} python ${SCRIPT} --step prepare --cancer ${CANCER}")
echo "Step 1 (prepare): job $J1"

# ── Step 2: Factorize (2 parallel workers) ─────────────────────────────────
FACTORIZE_JOBS=()
for worker in 0 1; do
  JF=$(maybe_sbatch \
    --job-name="cnmf_test_factorize_${CANCER}_w${worker}" \
    --dependency=afterok:${J1} \
    --time="0-03:00:00" --mem="32G" --cpus-per-task=4 \
    --output="${LOG_DIR}/test_factorize_${CANCER}_w${worker}-%j.out" \
    --error="${LOG_DIR}/test_factorize_${CANCER}_w${worker}-%j.err" \
    --wrap="conda run -n ${PYTHON_ENV} python ${SCRIPT} \
      --step factorize --cancer ${CANCER} \
      --worker ${worker} --total_workers 2")
  FACTORIZE_JOBS+=($JF)
  echo "Step 2 (factorize worker ${worker}): job $JF"
done

# ── Step 3: Combine ────────────────────────────────────────────────────────
FACTORIZE_DEP="afterok:${FACTORIZE_JOBS[0]}:${FACTORIZE_JOBS[1]}"
J3=$(maybe_sbatch \
  --job-name="cnmf_test_combine_${CANCER}" \
  --dependency=${FACTORIZE_DEP} \
  --time="0-01:00:00" --mem="32G" --cpus-per-task=2 \
  --output="${LOG_DIR}/test_combine_${CANCER}-%j.out" \
  --error="${LOG_DIR}/test_combine_${CANCER}-%j.err" \
  --wrap="conda run -n ${PYTHON_ENV} python ${SCRIPT} --step combine --cancer ${CANCER}")
echo "Step 3 (combine): job $J3"

echo ""
echo "=============================================="
echo " Test jobs submitted."
echo " Monitor:  squeue -u hadasra | grep cnmf_test"
echo " Logs:     ${LOG_DIR}/"
echo ""
echo " After combine finishes, check stability:"
echo "   ls /mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_test/${CANCER}/"
echo ""
echo " Then run consensus manually (pick best K):"
echo "   conda run -n scanvi-env python ${SCRIPT} \\"
echo "     --step consensus --cancer ${CANCER} --k <K>"
echo "=============================================="
