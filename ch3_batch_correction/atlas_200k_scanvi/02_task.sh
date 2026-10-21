#!/bin/bash
# 02_task.sh <classlist_file>  -- run one 02_assemble_counts.R task.
# Picks the project on line $SLURM_ARRAY_TASK_ID of the given class list.
set -euo pipefail
LIST="$1"
WT=/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/figure-1-2-scanvi-embedding-4e0a23/fig12
SIF=/groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif

p=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[[ -n "$p" ]] || { echo "no project at line $SLURM_ARRAY_TASK_ID of $LIST"; exit 1; }
echo "[$(date)] array task $SLURM_ARRAY_TASK_ID -> project $p"

apptainer exec --bind /groups/ofircohen-group:/ofir_grp "$SIF" \
  Rscript "$WT/02_assemble_counts.R" \
    --project "$p" \
    --root /ofir_grp/users/hadasra/temp_objects/figure1_2_scanvi \
    --to   /ofir_grp/users/hadasra/temp_objects
