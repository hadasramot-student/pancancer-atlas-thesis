#!/bin/bash
#SBATCH --job-name=3ca_derived
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work
#SBATCH --output=03_derived.out
#SBATCH --error=03_derived.err
#SBATCH --time=00:05:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
apptainer exec --bind /groups/ofircohen-group:/ofir_grp /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work/03_derived.R
