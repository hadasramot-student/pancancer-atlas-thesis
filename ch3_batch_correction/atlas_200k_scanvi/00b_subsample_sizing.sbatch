#!/bin/bash
#SBATCH --job-name=atlas_subsample_sizing
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/atlas-subsample-sizing-ab0533/work
#SBATCH --output=01_subsample_sizing.out
#SBATCH --error=01_subsample_sizing.err
#SBATCH --time=00:30:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2

apptainer exec --bind /groups/ofircohen-group:/ofir_grp \
  /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif \
  Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/atlas-subsample-sizing-ab0533/work/01_subsample_sizing.R
