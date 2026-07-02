#!/bin/bash
# submit_genenmf_all_cancers.sh
#
# Ch.5 (GEP programs) — submit one GeneNMF job per cancer type (epithelial
# compartment), in parallel, on SLURM. Each job runs
# genenmf/run_genenmf_per_cancer.R inside the Rstudio apptainer image (via
# ~/rstudio/rstudio.sbatch) and is fully independent.
#
# Per cancer it produces meta-programs, their top-50 genes, and a Jaccard
# comparison against that cancer's cNMF K=20 consensus programs.
#
# Usage:
#   ATLAS_GROUP_ROOT_HOST=/mnt/group ATLAS_GROUP_ROOT_CTR=/ctr_bind \
#     bash submit_genenmf_all_cancers.sh                 # submit all
#   bash submit_genenmf_all_cancers.sh --dry-run       # print sbatch lines only
#   bash submit_genenmf_all_cancers.sh --install       # submit install job first, chain the rest after it
#   bash submit_genenmf_all_cancers.sh --cancers "Lung Skin"   # subset
#
# TWO PATH DOMAINS — the apptainer image binds ATLAS_GROUP_ROOT_HOST -> ATLAS_GROUP_ROOT_CTR:
#   * HOST_* paths are used by SLURM and this shell:
#       sbatch --output/--error, mkdir, stat, file-existence tests.
#   * CTR_*  paths are used INSIDE the container:
#       the --script path and every path the R script opens (--input, --output_dir...).
# Using the wrong domain is silent failure (e.g. logs never written, stat -> 0 bytes).

set -euo pipefail

# ── Path roots ────────────────────────────────────────────────────────────────
HOST_ROOT="${ATLAS_GROUP_ROOT_HOST:?Set ATLAS_GROUP_ROOT_HOST, e.g. /mnt/new_groups/ofircohen-group}"
CTR_ROOT="${ATLAS_GROUP_ROOT_CTR:?Set ATLAS_GROUP_ROOT_CTR, e.g. /ofir_grp (the apptainer bind target)}"

# Scripts: derive the container path from this script's own (host) location, so it
# works whether run from the worktree or a canonical checkout.
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CTR_SCRIPT_DIR="${SELF_DIR/#${HOST_ROOT}/${CTR_ROOT}}"
R_SCRIPT="${CTR_SCRIPT_DIR}/run_genenmf_per_cancer.R"
INSTALL_R="${CTR_SCRIPT_DIR}/install_genenmf.R"
SBATCH_R="${HOME}/rstudio/rstudio.sbatch"        # apptainer Rstudio wrapper (host path)

# Data (host view for checks, container view for R args)
ATLAS_DATA_DIR="${ATLAS_DATA_DIR:-${HOST_ROOT}/users/$USER/temp_objects}"
CTR_DATA_DIR="${ATLAS_DATA_DIR/#${HOST_ROOT}/${CTR_ROOT}}"
HOST_BASE="${ATLAS_DATA_DIR}/merged_by_query"
CTR_BASE="${CTR_DATA_DIR}/merged_by_query"
HOST_OUT="${ATLAS_DATA_DIR}/genenmf_output"
CTR_OUT="${CTR_DATA_DIR}/genenmf_output"
LOG_DIR="${HOST_OUT}/logs"                        # sbatch writes here -> host path

# ── Cancer types ──────────────────────────────────────────────────────────────
# The 9 epithelial cancers + Brain (Brain loads other_Brain and is filtered to
# cell_type == "Malignant" at runtime). Breast is excluded until its cNMF job
# finishes — add "Breast" here once cnmf_output/Breast exists.
CANCERS=(Lung Ovarian Kidney Liver_Biliary Pancreas Prostate Sarcoma Skin Colorectal Brain)

# ── Analysis parameters (passed through to the R script) ──────────────────────
N_MP=20             # meta-programs per cancer (matches cNMF K=20)
K_RANGE="4:9"       # per-sample NMF k range
TOP_N=50            # top genes per program for Jaccard
JACCARD_THR=0.1     # "recovered by both" threshold
SAMPLE_COL=""       # empty = auto (column with most samples >= min_cells; favors orig.ident)

# ── Arg parsing ───────────────────────────────────────────────────────────────
DRY_RUN=false
DO_INSTALL=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    --install) DO_INSTALL=true; shift ;;
    --cancers) shift; read -r -a CANCERS <<< "$1"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

[[ "${DRY_RUN}" == "false" ]] && mkdir -p "${LOG_DIR}"

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] sbatch $*" >&2; echo "99999"
  else
    sbatch --parsable "$@"
  fi
}

# Memory / time tier by object size (scANVI-final RDS).
pick_mem() {
  local b=$1
  if   (( b < 1073741824 )); then echo "128G"   # < 1 GB
  elif (( b < 3221225472 )); then echo "256G"   # < 3 GB
  else                            echo "384G"
  fi
}
pick_time() {
  local b=$1
  if   (( b < 1073741824 )); then echo "0-12:00:00"
  elif (( b < 3221225472 )); then echo "1-00:00:00"
  else                            echo "2-00:00:00"
  fi
}

echo "============================================================"
echo " GeneNMF — per-cancer meta-programs + cNMF Jaccard comparison"
echo " Cancers : ${CANCERS[*]}"
echo " n_mp=${N_MP}  k=${K_RANGE}  top_n=${TOP_N}  jaccard>${JACCARD_THR}"
echo " Dry run : ${DRY_RUN}   Install first: ${DO_INSTALL}"
echo "============================================================"

# ── Optional install job; everything else waits on it ─────────────────────────
INSTALL_DEP=""
if [[ "${DO_INSTALL}" == "true" ]]; then
  J_INSTALL=$(maybe_sbatch \
    --job-name="genenmf_install" \
    --time="0-01:00:00" --mem="16G" --cpus-per-task=4 \
    --output="${LOG_DIR}/genenmf_install-%j.out" \
    --error="${LOG_DIR}/genenmf_install-%j.err" \
    "${SBATCH_R}" --script "${INSTALL_R}")
  echo "Install job: ${J_INSTALL}"
  INSTALL_DEP="--dependency=afterok:${J_INSTALL}"
fi

# ── One job per cancer ────────────────────────────────────────────────────────
submitted=0
for CANCER in "${CANCERS[@]}"; do
  # Brain: load other_Brain and filter to Malignant at runtime; others use Epithelial.
  if [[ "${CANCER}" == "Brain" ]]; then
    REL="Brain/other_Brain_merged_full_scanvi_final.RDS"
    FILTER_CT="Malignant"
  else
    REL="${CANCER}/Epithelial_${CANCER}_merged_full_scanvi_final.RDS"
    FILTER_CT=""
  fi
  HOST_RDS="${HOST_BASE}/${REL}"   # for shell checks / stat
  CTR_RDS="${CTR_BASE}/${REL}"     # for the R --input (inside container)

  if [[ ! -f "$HOST_RDS" ]] && [[ "$DRY_RUN" == "false" ]]; then
    echo "  [SKIP] ${CANCER} — object not found: ${HOST_RDS}"
    continue
  fi

  out_check="${HOST_OUT}/${CANCER}/${CANCER}_highconfidence_programs.csv"
  if [[ -f "$out_check" ]]; then
    echo "  [DONE] ${CANCER} — output exists, skipping"
    continue
  fi

  bytes=$(stat -c%s "$HOST_RDS" 2>/dev/null || echo 0)
  mem=$(pick_mem "$bytes")
  wtime=$(pick_time "$bytes")

  ARGS=(--script "${R_SCRIPT}"
        --cancer "${CANCER}"
        --input "${CTR_RDS}"
        --base_dir "${CTR_BASE}"
        --output_dir "${CTR_OUT}"
        --n_mp "${N_MP}"
        --k_range "${K_RANGE}"
        --top_n "${TOP_N}"
        --jaccard_thr "${JACCARD_THR}")
  [[ -n "${FILTER_CT}" ]]  && ARGS+=(--filter_celltype "${FILTER_CT}")
  [[ -n "${SAMPLE_COL}" ]] && ARGS+=(--sample_col "${SAMPLE_COL}")

  JOB=$(maybe_sbatch \
    --job-name="gnmf_${CANCER}" \
    ${INSTALL_DEP} \
    --time="${wtime}" --mem="${mem}" --cpus-per-task=15 \
    --output="${LOG_DIR}/gnmf_${CANCER}-%j.out" \
    --error="${LOG_DIR}/gnmf_${CANCER}-%j.err" \
    "${SBATCH_R}" "${ARGS[@]}")

  echo "  [SUBMIT ${mem}/${wtime}] ${CANCER} -> job ${JOB}"
  submitted=$((submitted + 1))
done

echo ""
echo "============================================================"
echo " Submitted ${submitted} GeneNMF job(s)."
echo " Monitor : squeue -u $USER | grep gnmf"
echo " Logs    : ${LOG_DIR}/"
echo " Outputs : ${HOST_OUT}/<Cancer>/"
echo "============================================================"
