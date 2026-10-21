#!/bin/bash
# run_scanvi_all_cancertypes.sh
#
# Submit the 3-step scANVI pipeline (steps 5→6→7) for every multi-batch
# merged object across all cancer types.
#
# Prerequisites:
#   1. Run check_batch_counts.R first to produce batch_counts.csv
#   2. All merged RDS files must exist under:
#      /mnt/.../temp_objects/merged_by_query/{CancerType}/{CellType}_{CancerType}_merged.RDS
#
# Usage:
#   ./run_scanvi_all_cancertypes.sh [--batch-csv <path>] [--dry-run]
#
# Options:
#   --batch-csv <path>   Path to batch_counts.csv  [default: auto-detect]
#   --dry-run            Print commands without submitting
#
# Paths inside the R container use /ofir_grp/ (bind: /groups/ofircohen-group → /ofir_grp).
# Paths for the Python host use /mnt/new_groups/ofircohen-group/.

set -euo pipefail

# ---------------------------------------------------------------------------
# Configurable paths
# ---------------------------------------------------------------------------
BASE_R="/ofir_grp/users/hadasra/temp_objects/merged_by_query"
BASE_PY="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query"
BASE_CHECK="/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query"

SCRIPT_DIR_R="/ofir_grp/users/hadasra/create_sc_object_workflow/merge"
SCRIPT_DIR_PY="/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/merge"
SEURAT_TO_H5AD_R="/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/stoic-leavitt-0e73f8/seurat_to_h5ad.R"

SBATCH_R=~/rstudio/rstudio.sbatch
SBATCH_PY=~/rstudio/python.sbatch

LOG_DIR="${BASE_CHECK}/pipeline_logs/scanvi_all_cancertypes"
SKIP_LOG="${LOG_DIR}/skipped_objects.txt"
JOB_LOG="${LOG_DIR}/submitted_jobs.txt"

BATCH_CSV_DEFAULT="${BASE_CHECK}/batch_counts.csv"

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
DRY_RUN=false
BATCH_CSV="${BATCH_CSV_DEFAULT}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)   DRY_RUN=true;  shift ;;
    --batch-csv) BATCH_CSV="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ ! -f "${BATCH_CSV}" ]]; then
  echo "ERROR: batch_counts.csv not found at: ${BATCH_CSV}"
  echo "Run check_batch_counts.R first:"
  echo "  sbatch ~/rstudio/rstudio.sbatch \\"
  echo "    --script /ofir_grp/users/hadasra/create_sc_object_workflow/merge/check_batch_counts.R"
  exit 1
fi

mkdir -p "${LOG_DIR}"
> "${SKIP_LOG}"
> "${JOB_LOG}"

echo "========================================================"
echo " scANVI pipeline — all cancer types"
echo "========================================================"
echo " Batch CSV : ${BATCH_CSV}"
echo " Dry run   : ${DRY_RUN}"
echo " Log dir   : ${LOG_DIR}"
echo "========================================================"
echo ""

N_SUBMITTED=0
N_SKIPPED_SINGLE=0
N_SKIPPED_MISSING=0
N_SKIPPED_DONE=0
N_ERRORS=0

# ---------------------------------------------------------------------------
# Submit helper: echoes the command; runs it unless --dry-run
# ---------------------------------------------------------------------------
maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "    [DRY-RUN] sbatch $*"
    echo "99999"   # fake job ID
  else
    sbatch --parsable "$@"
  fi
}

# ---------------------------------------------------------------------------
# Main loop — read CSV (skip header)
# ---------------------------------------------------------------------------
while IFS=',' read -r cell_type cancer_type rds_path n_batches status; do
  # strip quotes that write.csv adds around fields
  cell_type="${cell_type//\"/}"
  cancer_type="${cancer_type//\"/}"
  rds_path="${rds_path//\"/}"
  n_batches="${n_batches//\"/}"
  status="${status//\"/}"

  OBJ="${cell_type}_${cancer_type}"
  SUBDIR="${cancer_type}"

  # ---- Skip non-multi-batch rows ----
  case "${status}" in
    missing)
      echo "[SKIP-MISSING] ${OBJ}"
      echo "MISSING|${OBJ}|${rds_path}" >> "${SKIP_LOG}"
      N_SKIPPED_MISSING=$(( N_SKIPPED_MISSING + 1 ))
      continue
      ;;
    single_batch)
      echo "[SKIP-SINGLE]  ${OBJ} (1 batch)"
      echo "SINGLE_BATCH|${OBJ}|${rds_path}" >> "${SKIP_LOG}"
      N_SKIPPED_SINGLE=$(( N_SKIPPED_SINGLE + 1 ))
      continue
      ;;
    error*)
      echo "[ERROR]        ${OBJ} — ${status}"
      echo "ERROR|${OBJ}|${status}" >> "${SKIP_LOG}"
      N_ERRORS=$(( N_ERRORS + 1 ))
      continue
      ;;
    multi_batch)
      : # fall through
      ;;
    *)
      echo "[UNKNOWN]      ${OBJ} — status: ${status}"
      continue
      ;;
  esac

  # ---- Derive paths ----
  INPUT_R="${BASE_R}/${SUBDIR}/${OBJ}_merged.RDS"
  H5AD_R="${BASE_R}/${SUBDIR}/${OBJ}_merged_full.h5ad"
  H5AD_PY="${BASE_PY}/${SUBDIR}/${OBJ}_merged_full.h5ad"
  PREFIX_PY="${BASE_PY}/${SUBDIR}/${OBJ}_merged_full_scanvi"
  LATENT_R="${BASE_R}/${SUBDIR}/${OBJ}_merged_full_scanvi_scanvi_latent.csv"
  LABELS_R="${BASE_R}/${SUBDIR}/${OBJ}_merged_full_scanvi_scanvi_labels.csv"
  FINAL_R="${BASE_R}/${SUBDIR}/${OBJ}_merged_full_scanvi_final.RDS"
  FINAL_CHECK="${BASE_CHECK}/${SUBDIR}/${OBJ}_merged_full_scanvi_final.RDS"

  # ---- Skip already-completed objects ----
  if [[ -f "${FINAL_CHECK}" ]]; then
    echo "[DONE]         ${OBJ} — final RDS already exists, skipping"
    echo "ALREADY_DONE|${OBJ}|${FINAL_CHECK}" >> "${SKIP_LOG}"
    N_SKIPPED_DONE=$(( N_SKIPPED_DONE + 1 ))
    continue
  fi

  echo "[SUBMIT]       ${OBJ} (${n_batches} batches)"

  # ---- Step 5: Seurat → h5ad (R container) ----
  JOB5=$(maybe_sbatch \
    --job-name="s5_h5ad_${OBJ}" \
    --time=0-04:00:00 \
    --mem=128G \
    --cpus-per-task=8 \
    --output="${LOG_DIR}/s5_h5ad_${OBJ}-%j.out" \
    --error="${LOG_DIR}/s5_h5ad_${OBJ}-%j.err" \
    "${SBATCH_R}" \
      --script "${SEURAT_TO_H5AD_R}" \
      --input_rds   "${INPUT_R}" \
      --output_h5ad "${H5AD_R}")
  echo "  Step 5 (h5ad):   job ${JOB5}"

  # ---- Step 6: scANVI training (Python, GPU) ----
  JOB6=$(maybe_sbatch \
    --job-name="s6_scanvi_${OBJ}" \
    --dependency=afterok:"${JOB5}" \
    --time=1-00:00:00 \
    --gres=gpu:rtx_3090:1 \
    --output="${LOG_DIR}/s6_scanvi_${OBJ}-%j.out" \
    --error="${LOG_DIR}/s6_scanvi_${OBJ}-%j.err" \
    "${SBATCH_PY}" \
      --script "${SCRIPT_DIR_PY}/export_scanvi_to_csv.py" \
      --h5ad                  "${H5AD_PY}" \
      --output_prefix         "${PREFIX_PY}" \
      --batch_key             project_name \
      --labels_key            scimilarity_celltype \
      --unlabeled_category    Unknown \
      --layer_counts          counts \
      --n_latent              30 \
      --max_epochs_scvi       200 \
      --max_epochs_scanvi     100 \
      --label_key_out         scanvi_pred)
  echo "  Step 6 (scANVI): job ${JOB6}"

  # ---- Step 7: Import scANVI back into Seurat (R container) ----
  JOB7=$(maybe_sbatch \
    --job-name="s7_seurat_${OBJ}" \
    --dependency=afterok:"${JOB6}" \
    --time=0-04:00:00 \
    --mem=128G \
    --cpus-per-task=8 \
    --output="${LOG_DIR}/s7_seurat_${OBJ}-%j.out" \
    --error="${LOG_DIR}/s7_seurat_${OBJ}-%j.err" \
    "${SBATCH_R}" \
      --script "${SCRIPT_DIR_R}/scanvi_into_seurat.R" \
      --seurat      "${INPUT_R}" \
      --latent_csv  "${LATENT_R}" \
      --labels_csv  "${LABELS_R}" \
      --output      "${FINAL_R}" \
      --label_col   scanvi_pred \
      --meta_col    scanvi_pred \
      --reduction   scanvi \
      --overwrite)
  echo "  Step 7 (Seurat): job ${JOB7}"

  echo "  Jobs: ${JOB5} -> ${JOB6} -> ${JOB7}" >> "${JOB_LOG}"
  echo "${OBJ}|step5=${JOB5}|step6=${JOB6}|step7=${JOB7}" >> "${JOB_LOG}"
  echo ""

  N_SUBMITTED=$(( N_SUBMITTED + 1 ))

done < <(tail -n +2 "${BATCH_CSV}")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "========================================================"
echo " Submission complete"
echo "========================================================"
echo "  Pipelines submitted : ${N_SUBMITTED}"
echo "  Already done        : ${N_SKIPPED_DONE}"
echo "  Skipped (1 batch)   : ${N_SKIPPED_SINGLE}"
echo "  Skipped (missing)   : ${N_SKIPPED_MISSING}"
echo "  Errors in CSV       : ${N_ERRORS}"
echo ""
echo "  Skip log : ${SKIP_LOG}"
echo "  Job log  : ${JOB_LOG}"
echo ""
echo "  Monitor  : squeue -u hadasra | grep -E 's[567]_'"
echo "========================================================"
