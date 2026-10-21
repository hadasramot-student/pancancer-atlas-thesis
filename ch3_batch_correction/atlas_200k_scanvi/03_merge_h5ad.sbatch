#!/bin/bash
#SBATCH --job-name=fig12_03_merge
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi/logs
#SBATCH --output=03_merge.%j.out
#SBATCH --error=03_merge.%j.err
#SBATCH --partition=main
#SBATCH --time=04:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4

apptainer exec --bind /groups/ofircohen-group:/ofir_grp \
  /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif \
  Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/figure-1-2-scanvi-embedding-4e0a23/fig12/03_merge_h5ad.R
