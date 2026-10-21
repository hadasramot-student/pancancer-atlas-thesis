#!/bin/bash
# Submit one SLURM job per cancer_type x cell_type combination — all 11 cancer types from scratch.
# Uses apptainer container via rstudio.sbatch + run_rscript_vanilla.sh wrapper
# to avoid .Rprofile devtools/fs version conflict.

set -euo pipefail

BASE_OUT_REMOTE="/ofir_grp/users/hadasra/temp_objects/merged_by_query"
META_RDS="/ofir_grp/users/hadasra/temp_objects/all_cells.RDS"
INPUT_DIR="/ofir_grp/users/hadasra/temp_objects"
SCRIPT="/ofir_grp/users/hadasra/create_sc_object_workflow/merge/merge_by_cancertype_celltype.R"
WRAPPER="/ofir_grp/users/hadasra/create_sc_object_workflow/merge/run_rscript_vanilla.sh"
SBATCH_R="${HOME}/rstudio/rstudio.sbatch"
LOG_DIR="/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/merge"

CELL_TYPES=(
  "B-cell"
  "Endothelial"
  "Fibroblasts"
  "Hematopoietic_Stem_and_Progenitor_cells"
  "Myeloid"
  "NK-cells"
  "T-cells"
  "Epithelial"
  "other"
)

CANCER_TYPES=(
  "Brain"
  "Colorectal"
  "Head_and_Neck"
  "Kidney"
  "Liver_Biliary"
  "Lung"
  "Ovarian"
  "Pancreas"
  "Prostate"
  "Sarcoma"
  "Skin"
)

submitted=0

for cancer in "${CANCER_TYPES[@]}"; do
  # Create output dir on login node so R doesn't need to
  mkdir -p "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/${cancer}"
  out_dir_remote="${BASE_OUT_REMOTE}/${cancer}"

  for celltype in "${CELL_TYPES[@]}"; do
    safe_cell=$(echo "${celltype}" | sed 's/[^A-Za-z0-9._-]/_/g')
    safe_cancer=$(echo "${cancer}" | sed 's/[^A-Za-z0-9._-]/_/g')
    job_name="merge_${safe_cell}_${safe_cancer}"
    job_name="${job_name:0:64}"

    jobid=$(sbatch --parsable \
      --job-name="${job_name}" \
      --output="${LOG_DIR}/${job_name}-%j.out" \
      --error="${LOG_DIR}/${job_name}-%j.err" \
      --time=1-12:00:00 \
      --mem=200G \
      --cpus-per-task=10 \
      "${SBATCH_R}" \
        --script "${WRAPPER}" \
        "${SCRIPT}" \
        --meta_rds "${META_RDS}" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${out_dir_remote}" \
        --cancer_type "${cancer}" \
        --cell_type "${celltype}")

    echo "SUBMITTED job ${jobid}: ${celltype} x ${cancer}"
    ((submitted++)) || true
  done
done

echo ""
echo "Done. Submitted: ${submitted} jobs total."
echo "Monitor with: squeue -u hadasra | grep merge"
