#!/bin/bash
#SBATCH --partition=main
#SBATCH --time=02:00:00
#SBATCH --job-name=ch3_dup_overlap
#SBATCH --output=%x-%J.out
#SBATCH --error=%x-%J.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G

echo "start: $(date)  node: $(hostname)"
/groups/ofircohen-group/rstudio_images/apptainer-rstudio.sh v1.1 \
  --script /groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/chapter3-benchmark-duplicates-3ba85e/ch3_dup_overlap_check.R
echo "end: $(date)"
