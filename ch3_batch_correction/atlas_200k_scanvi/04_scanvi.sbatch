#!/bin/bash
#SBATCH --job-name=fig12_04_scanvi
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi/logs
#SBATCH --output=04_scanvi.%j.out
#SBATCH --error=04_scanvi.%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:rtx_4090:1
#SBATCH --exclude=cs-3090-06
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -x
echo "host=$(hostname)  date=$(date)  jobid=$SLURM_JOB_ID"
nvidia-smi || { echo "NO GPU VISIBLE - aborting"; exit 1; }

source /storage/modules/packages/anaconda/etc/profile.d/conda.sh
conda activate scanvi-env
echo "python: $(which python)   env: $CONDA_DEFAULT_ENV"
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else '-')"

ROOT=/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi
python /groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/figure-1-2-scanvi-embedding-4e0a23/fig12/04_scanvi.py \
  --h5ad       "$ROOT/fig12_200k.h5ad" \
  --out_prefix "$ROOT/fig12_200k" \
  --batch_key project_name \
  --labels_key Consensus_Cell_Type \
  --n_hvg 2000 --n_latent 30 --n_layers 2 \
  --scvi_epochs 400 --scanvi_epochs 200 \
  --n_neighbors 30 --seed 42
echo "exit=$? date=$(date)"
