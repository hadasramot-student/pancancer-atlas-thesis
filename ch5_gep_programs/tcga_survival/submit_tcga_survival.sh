#!/bin/bash
# submit_tcga_survival.sh
#
# Submits run_tcga_survival.R inside the Rstudio apptainer image via
# ~/rstudio/rstudio.sbatch. See submit_genenmf_all_cancers.sh for the
# HOST/CTR path-translation rationale (apptainer binds
# ATLAS_GROUP_ROOT_HOST -> ATLAS_GROUP_ROOT_CTR).
#
# Usage:
#   ATLAS_GROUP_ROOT_HOST=/mnt/new_groups/ofircohen-group ATLAS_GROUP_ROOT_CTR=/ofir_grp \
#     bash submit_tcga_survival.sh

set -euo pipefail

HOST_ROOT="${ATLAS_GROUP_ROOT_HOST:?Set ATLAS_GROUP_ROOT_HOST, e.g. /mnt/new_groups/ofircohen-group}"
CTR_ROOT="${ATLAS_GROUP_ROOT_CTR:?Set ATLAS_GROUP_ROOT_CTR, e.g. /ofir_grp (the apptainer bind target)}"

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
CTR_SCRIPT_DIR="${SELF_DIR/#${HOST_ROOT}/${CTR_ROOT}}"
R_SCRIPT="${CTR_SCRIPT_DIR}/run_tcga_survival.R"
SBATCH_R="${HOME}/rstudio/rstudio.sbatch"

ATLAS_DATA_DIR="${ATLAS_DATA_DIR:-${HOST_ROOT}/users/$USER/temp_objects}"
CTR_DATA_DIR="${ATLAS_DATA_DIR/#${HOST_ROOT}/${CTR_ROOT}}"
LOG_DIR="${ATLAS_DATA_DIR}/tcga_clinical/logs"
mkdir -p "${LOG_DIR}"

sbatch --parsable \
  --job-name tcga_survival \
  --output "${LOG_DIR}/%x-%J.out" \
  --error "${LOG_DIR}/%x-%J.err" \
  --cpus-per-task=4 \
  --mem=32G \
  --time 0-02:00:00 \
  --export=ALL,ATLAS_DATA_DIR="${CTR_DATA_DIR}" \
  "${SBATCH_R}" --script "${R_SCRIPT}"
