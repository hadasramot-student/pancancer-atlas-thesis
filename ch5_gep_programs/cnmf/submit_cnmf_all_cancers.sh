#!/bin/bash
# submit_cnmf_all_cancers.sh
#
# Submits cNMF pipeline for each cancer type in parallel on SLURM.
# Each cancer type gets its own independent chain: prepare → factorize → combine → consensus
#
# Fixed K across all cancer types for cross-cancer comparability.
#
# Usage:
#   ./submit_cnmf_all_cancers.sh              # submit all
#   ./submit_cnmf_all_cancers.sh --dry-run    # print only
#   ./submit_cnmf_all_cancers.sh --k 20       # override K (default: 20)

set -euo pipefail

SCRIPT_DIR="/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/cnmf"
CNMF_SCRIPT="${SCRIPT_DIR}/run_cnmf_per_cancer.py"
LOG_DIR="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output/logs"
PYTHON_ENV="scanvi-env"

# ── Parameters ────────────────────────────────────────────────────────────
FIXED_K=20          # fixed K for all cancer types
N_ITER=200          # replicates per K (full run)
N_CELLS=0           # 0 = use all cells (no downsampling)
N_WORKERS=4         # parallel factorize workers per cancer type

# ── Cancer types + h5ad paths ─────────────────────────────────────────────
declare -A H5AD_PATHS
H5AD_PATHS["Colorectal"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Colorectal/Epithelial_Colorectal_merged_full.h5ad"
H5AD_PATHS["Kidney"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Kidney/Epithelial_Kidney_merged_full.h5ad"
H5AD_PATHS["Liver_Biliary"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Liver_Biliary/Epithelial_Liver_Biliary_merged_full.h5ad"
H5AD_PATHS["Lung"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Lung/Epithelial_Lung_merged_full.h5ad"
H5AD_PATHS["Ovarian"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Ovarian/Epithelial_Ovarian_merged_full.h5ad"
H5AD_PATHS["Pancreas"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Pancreas/Epithelial_Pancreas_merged_full.h5ad"
H5AD_PATHS["Prostate"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Prostate/Epithelial_Prostate_merged_full.h5ad"
H5AD_PATHS["Sarcoma"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Sarcoma/Epithelial_Sarcoma_merged_full.h5ad"
H5AD_PATHS["Skin"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Skin/Epithelial_Skin_merged_full.h5ad"
H5AD_PATHS["Brain"]="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Brain/Epithelial_Brain_merged_full.h5ad"

mkdir -p "${LOG_DIR}"

DRY_RUN=false
for arg in "$@"; do
  [[ "$arg" == "--dry-run" ]] && DRY_RUN=true
  [[ "$arg" == "--k" ]] && shift && FIXED_K=$1
done

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*" >&2; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

echo "============================================================"
echo " Pan-cancer cNMF — all cancer types in parallel"
echo " Fixed K: ${FIXED_K}"
echo " Replicates: ${N_ITER}"
echo " Factorize workers per cancer: ${N_WORKERS}"
echo " Dry run: ${DRY_RUN}"
echo "============================================================"

# ── Submit one chain per cancer type ──────────────────────────────────────
for CANCER in "${!H5AD_PATHS[@]}"; do
  H5AD="${H5AD_PATHS[$CANCER]}"

  # Skip if h5ad doesn't exist
  if [[ ! -f "$H5AD" ]] && [[ "$DRY_RUN" == "false" ]]; then
    echo "WARNING: ${CANCER} — h5ad not found, skipping: ${H5AD}"
    continue
  fi

  echo ""
  echo "--- ${CANCER} ---"

  # Step 1: Prepare
  J_PREP=$(maybe_sbatch \
    --job-name="cnmf_prep_${CANCER}" \
    --time="0-04:00:00" --mem="128G" --cpus-per-task=8 \
    --output="${LOG_DIR}/cnmf_prep_${CANCER}-%j.out" \
    --error="${LOG_DIR}/cnmf_prep_${CANCER}-%j.err" \
    --wrap="conda run -n ${PYTHON_ENV} python ${CNMF_SCRIPT} \
      --step prepare --cancer ${CANCER} \
      --n_cells ${N_CELLS} --n_iter ${N_ITER} --k_values ${FIXED_K}")
  echo "  Prepare: job ${J_PREP}"

  # Steps 2a-2d: Factorize (N_WORKERS parallel)
  FACTORIZE_JOBS=()
  for worker in $(seq 0 $((N_WORKERS - 1))); do
    JF=$(maybe_sbatch \
      --job-name="cnmf_fact_${CANCER}_w${worker}" \
      --dependency=afterok:${J_PREP} \
      --time="2-00:00:00" --mem="64G" --cpus-per-task=8 \
      --output="${LOG_DIR}/cnmf_fact_${CANCER}_w${worker}-%j.out" \
      --error="${LOG_DIR}/cnmf_fact_${CANCER}_w${worker}-%j.err" \
      --wrap="conda run -n ${PYTHON_ENV} python ${CNMF_SCRIPT} \
        --step factorize --cancer ${CANCER} \
        --worker ${worker} --total_workers ${N_WORKERS}")
    FACTORIZE_JOBS+=($JF)
    echo "  Factorize worker ${worker}: job ${JF}"
  done

  # Step 3: Combine
  FACT_DEP=$(printf "afterok:%s" "${FACTORIZE_JOBS[0]}"; printf ":%s" "${FACTORIZE_JOBS[@]:1}")
  J_COMB=$(maybe_sbatch \
    --job-name="cnmf_comb_${CANCER}" \
    --dependency=${FACT_DEP} \
    --time="0-02:00:00" --mem="64G" --cpus-per-task=4 \
    --output="${LOG_DIR}/cnmf_comb_${CANCER}-%j.out" \
    --error="${LOG_DIR}/cnmf_comb_${CANCER}-%j.err" \
    --wrap="conda run -n ${PYTHON_ENV} python ${CNMF_SCRIPT} \
      --step combine --cancer ${CANCER}")
  echo "  Combine: job ${J_COMB}"

  # Step 4: Consensus at fixed K
  J_CONS=$(maybe_sbatch \
    --job-name="cnmf_cons_${CANCER}" \
    --dependency=afterok:${J_COMB} \
    --time="0-01:00:00" --mem="64G" --cpus-per-task=4 \
    --output="${LOG_DIR}/cnmf_cons_${CANCER}-%j.out" \
    --error="${LOG_DIR}/cnmf_cons_${CANCER}-%j.err" \
    --wrap="conda run -n ${PYTHON_ENV} python ${CNMF_SCRIPT} \
      --step consensus --cancer ${CANCER} --k ${FIXED_K}")
  echo "  Consensus K=${FIXED_K}: job ${J_CONS}"

done

echo ""
echo "============================================================"
echo " All chains submitted."
echo " Monitor: squeue -u hadasra | grep cnmf"
echo " Logs: ${LOG_DIR}/"
echo "============================================================"
