#!/bin/bash
#SBATCH --job-name=3ca_exclusions
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/atlas-metadata-consensus-validation-e53a85/work
#SBATCH --output=01_exclusions.out
#SBATCH --error=01_exclusions.err
#SBATCH --time=00:40:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2

apptainer exec --bind /groups/ofircohen-group:/ofir_grp \
  /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif \
  Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/atlas-metadata-consensus-validation-e53a85/work/01_exclusions.R
