#!/bin/bash
#SBATCH --job-name=proj_audit
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/plots_and_outputs/project_count_audit/dump_projects-%j.log
apptainer exec --bind /groups/ofircohen-group:/ofir_grp /groups/ofircohen-group/rstudio_images/rstudio_v1.1.sif \
  Rscript /ofir_grp/users/hadasra/create_sc_object_workflow/plots_and_outputs/project_count_audit/dump_projects.R
