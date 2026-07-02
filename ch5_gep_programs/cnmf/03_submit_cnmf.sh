#!/bin/bash
# 03_submit_cnmf.sh
#
# Ch.5 (GEP programs) — cNMF step 3/4. Submits the full cNMF pipeline to SLURM:
#   Step 1 (R):      export raw counts -> counts.mtx (01_export_counts_for_cnmf.R)
#   Step 2 (Python): cNMF prepare
#   Step 3 (Python): cNMF factorize  (parallelized across 7 K values x N_WORKERS)
#   Step 4 (Python): cNMF combine
#
# Some steps run inside the R apptainer container (which may see this repo
# and the data dir at different mount points than this host shell); set
# ATLAS_REPO_DIR_CTR if that differs from ATLAS_REPO_DIR.
#
# Usage:
#   ATLAS_REPO_DIR=/path/to/repo ATLAS_DATA_DIR=/path/to/temp_objects \
#     ./03_submit_cnmf.sh           # submit all steps
#   ./03_submit_cnmf.sh --dry-run   # print only

set -euo pipefail

ATLAS_REPO_DIR="${ATLAS_REPO_DIR:-.}"
ATLAS_REPO_DIR_CTR="${ATLAS_REPO_DIR_CTR:-$ATLAS_REPO_DIR}"
ATLAS_DATA_DIR="${ATLAS_DATA_DIR:-data/temp_objects}"
SCRIPT_DIR="${ATLAS_REPO_DIR}/cnmf"
LOG_DIR="${ATLAS_DATA_DIR}/cnmf_output/logs"
SBATCH_R=~/rstudio/rstudio.sbatch
SBATCH_PY=~/rstudio/python.sbatch
N_WORKERS=4   # parallel factorize workers per K

mkdir -p "${LOG_DIR}"

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*" >&2; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

echo "=============================================="
echo " Pan-cancer cNMF pipeline"
echo " Dry run: ${DRY_RUN}"
echo "=============================================="

# ── Step 1: Export raw counts from RDS objects ─────────────────────────────
J1=$(maybe_sbatch \
  --job-name="cnmf_01_export" \
  --time="1-00:00:00" --mem="128G" --cpus-per-task=10 \
  --output="${LOG_DIR}/cnmf_01_export-%j.out" \
  --error="${LOG_DIR}/cnmf_01_export-%j.err" \
  "${SBATCH_R}" \
    --script "${ATLAS_REPO_DIR_CTR}/cnmf/01_export_counts_for_cnmf.R")
echo "Step 1 (export): job $J1"

# ── Step 2: cNMF prepare ──────────────────────────────────────────────────
J2=$(maybe_sbatch \
  --job-name="cnmf_02_prepare" \
  --dependency=afterok:${J1} \
  --time="0-04:00:00" --mem="128G" --cpus-per-task=8 \
  --output="${LOG_DIR}/cnmf_02_prepare-%j.out" \
  --error="${LOG_DIR}/cnmf_02_prepare-%j.err" \
  "${SBATCH_PY}" \
    --script "${SCRIPT_DIR}/02_run_cnmf.py" \
    --step prepare)
echo "Step 2 (prepare): job $J2"

# ── Step 3: cNMF factorize (N_WORKERS parallel jobs) ─────────────────────
FACTORIZE_JOBS=()
for worker in $(seq 0 $((N_WORKERS - 1))); do
  JF=$(maybe_sbatch \
    --job-name="cnmf_03_factorize_w${worker}" \
    --dependency=afterok:${J2} \
    --time="2-00:00:00" --mem="64G" --cpus-per-task=8 \
    --output="${LOG_DIR}/cnmf_03_factorize_w${worker}-%j.out" \
    --error="${LOG_DIR}/cnmf_03_factorize_w${worker}-%j.err" \
    "${SBATCH_PY}" \
      --script "${SCRIPT_DIR}/02_run_cnmf.py" \
      --step factorize \
      --worker ${worker} \
      --total_workers ${N_WORKERS})
  FACTORIZE_JOBS+=($JF)
  echo "Step 3 (factorize worker ${worker}): job $JF"
done

# ── Step 4: cNMF combine (after all factorize workers done) ───────────────
FACTORIZE_DEP=$(IFS=:; echo "afterok:${FACTORIZE_JOBS[*]}" | sed 's/ /:/g')
J4=$(maybe_sbatch \
  --job-name="cnmf_04_combine" \
  --dependency=${FACTORIZE_DEP} \
  --time="0-02:00:00" --mem="64G" --cpus-per-task=4 \
  --output="${LOG_DIR}/cnmf_04_combine-%j.out" \
  --error="${LOG_DIR}/cnmf_04_combine-%j.err" \
  "${SBATCH_PY}" \
    --script "${SCRIPT_DIR}/02_run_cnmf.py" \
    --step combine)
echo "Step 4 (combine): job $J4"

echo ""
echo "=============================================="
echo " Jobs submitted. After combine finishes:"
echo " Inspect stability plots, then run consensus:"
echo ""
echo "   conda run -n scanvi-env python ${SCRIPT_DIR}/02_run_cnmf.py \\"
echo "     --step consensus --k <optimal_K> --density_threshold 0.1"
echo ""
echo " Monitor: squeue -u $USER | grep cnmf"
echo " Logs: ${LOG_DIR}/"
echo "=============================================="
