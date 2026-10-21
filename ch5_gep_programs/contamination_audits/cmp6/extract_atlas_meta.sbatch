#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=cmp6_atlas_meta
#SBATCH --output=%x-%J.out
#SBATCH --error=%x-%J.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G

echo "start $(date)"
/groups/ofircohen-group/rstudio_images/apptainer-rstudio.sh v1.1 \
  --script /groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/cmp6-myeloid-contamination-83948d/cmp6_check/extract_atlas_meta.R
echo "end $(date)"
