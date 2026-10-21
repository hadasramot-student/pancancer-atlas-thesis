#!/bin/bash
#SBATCH --job-name=3ca_validate
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work
#SBATCH --output=02_validate.out
#SBATCH --error=02_validate.err
#SBATCH --time=00:40:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2

apptainer exec --bind /groups/ofircohen-group:/ofir_grp \
  /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif \
  Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work/02_validate.R
