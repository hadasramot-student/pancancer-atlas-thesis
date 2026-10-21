#!/bin/bash
#SBATCH --job-name=fig12_01_subsample
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi/logs
#SBATCH --output=01_subsample.%j.out
#SBATCH --error=01_subsample.%j.err
#SBATCH --partition=main
#SBATCH --time=01:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=4

apptainer exec --bind /groups/ofircohen-group:/ofir_grp \
  /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif \
  Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/figure-1-2-scanvi-embedding-4e0a23/fig12/01_subsample.R
