#!/bin/bash
#SBATCH --partition=main
#SBATCH --time=2:00:00
#SBATCH --job-name=kirc_consol_pt
#SBATCH --output=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/kirc-survival-consolidated-pt-da69f5/consolidated_pt/logs/%x-%J.out
#SBATCH --error=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/kirc-survival-consolidated-pt-da69f5/consolidated_pt/logs/%x-%J.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

echo "start: $(date)  node: $(hostname)"
/groups/ofircohen-group/rstudio_images/apptainer-rstudio.sh v1.1 --script \
  /ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/kirc-survival-consolidated-pt-da69f5/consolidated_pt/01_consolidated_pt_survival.R
echo "end: $(date)"
