#!/bin/bash
#SBATCH --job-name=3ca_inspect
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work
#SBATCH --output=01_inspect.out
#SBATCH --error=01_inspect.err
#SBATCH --time=00:30:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2

apptainer exec --bind /groups/ofircohen-group:/ofir_grp \
  /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif \
  Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work/01_inspect.R
