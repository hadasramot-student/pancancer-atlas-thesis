#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=cmp10_check
#SBATCH --output=%x-%J.out
#SBATCH --error=%x-%J.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G

# activate the python env used for the cMP6 audit (needs anndata, scipy, pandas):
# source ~/miniconda3/etc/profile.d/conda.sh && conda activate <your_env>

echo "start $(date)"
/home/hadasra/.conda/envs/scanvi-env/bin/python /groups/ofircohen-group/users/hadasra/create_sc_object_workflow/cmp10_check/analyze.py
echo "end $(date)"
