#!/bin/bash
# scanvi_example.slurm
#
# Ch.3 (Batch correction) - example SLURM job running scANVI batch
# correction + latent export for one lineage x cancer-type object.
# Adjust ATLAS_REPO_DIR / H5AD / OUTPUT_PREFIX for your object before
# submitting (`sbatch scanvi_example.slurm`).

#SBATCH --job-name=scanvi_hspc
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=scanvi_hspc.%j.out
#SBATCH --error=scanvi_hspc.%j.err

# Load the cluster's Anaconda module
module load anaconda

# Enable 'conda activate' in this non-interactive shell
eval "$(conda shell.bash hook)"

# Activate your environment
conda activate scanvi-env

ATLAS_REPO_DIR="${ATLAS_REPO_DIR:-.}"
ATLAS_DATA_DIR="${ATLAS_DATA_DIR:-data/temp_objects}"

# Run scANVI latent export for one example object (HSPC lineage, Breast)
python "${ATLAS_REPO_DIR}/ch3_batch_correction/export_scanvi_to_csv.py" \
  --h5ad "${ATLAS_DATA_DIR}/merged_by_query/Hematopoietic_Stem_and_Progenitor_cells_Breast_merged_batchco_scANVI.h5ad" \
  --output_prefix "${ATLAS_DATA_DIR}/merged_by_query/Hematopoietic_Stem_and_Progenitor_cells_Breast" \
  --layer_counts counts \
  --n_latent 30 \
  --max_epochs_scvi 200 \
  --max_epochs_scanvi 100 \
  --unlabeled_category Unknown \
  --label_key_out scanvi_pred
