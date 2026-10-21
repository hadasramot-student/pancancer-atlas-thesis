#!/bin/bash
#SBATCH --job-name=fig12_05_figs
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi/logs
#SBATCH --output=05_figs.%j.out
#SBATCH --error=05_figs.%j.err
#SBATCH --partition=main
#SBATCH --time=00:40:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

source /storage/modules/packages/anaconda/etc/profile.d/conda.sh
conda activate scanvi-env
python /groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/figure-1-2-scanvi-embedding-4e0a23/fig12/05_figures.py \
  --root /groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi
