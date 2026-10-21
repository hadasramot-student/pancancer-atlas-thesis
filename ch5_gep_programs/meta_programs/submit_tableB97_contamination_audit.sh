#!/bin/bash
# submit_tableB97_contamination_audit.sh
#
# Prompt 9c -- submit the Table B97 contamination audit as a single SLURM job.
# Pattern copied from ch5_gep_programs/tcga_survival/submit_tcga_survival.sh and
# ids_analysis/submit_ids_all_cancers.sh (conda env `scanvi-env`, no container;
# all input/output paths are host paths under /groups/ofircohen-group).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
AUDIT_PY="${SCRIPT_DIR}/tableB97_contamination_audit.py"
PYTHON_ENV="scanvi-env"

LOG_DIR="/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output/logs"
mkdir -p "${LOG_DIR}"

sbatch --parsable \
  --job-name="tableB97_audit" \
  --partition=cpu \
  --time="0-00:30:00" \
  --mem="16G" \
  --cpus-per-task=2 \
  --output="${LOG_DIR}/tableB97_audit-%j.out" \
  --error="${LOG_DIR}/tableB97_audit-%j.err" \
  --wrap="set -euo pipefail; source /storage/modules/packages/anaconda/etc/profile.d/conda.sh; conda run --no-capture-output -n ${PYTHON_ENV} python ${AUDIT_PY}"
