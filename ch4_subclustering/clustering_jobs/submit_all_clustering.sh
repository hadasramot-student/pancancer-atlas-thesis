#!/bin/bash
# submit_all_clustering.sh
#
# Ch.4 (Subclustering) — generates and submits one PCA+UMAP+clustering SLURM
# job per lineage (8 lineages: T-cells, Myeloid, Epithelial, B-cells,
# Fibroblasts, Endothelial, NK-cells, HSPC), each running run_clustering.R
# (30 PCs, 5 Louvain/Leiden resolutions) inside the R apptainer container.
#
# Inputs:  ATLAS_DATA_DIR/merged_by_cell_type/<lineage RDS, per RDS_MAP below>
# Outputs: ATLAS_DATA_DIR/clustered_lineages/<lineage>_clustered.RDS + summary
#
# Usage (run from this directory):
#   ATLAS_DATA_DIR=/path/to/temp_objects APPTAINER_R_IMAGE=/path/to/rstudio.sif \
#   APPTAINER_BIND_ARGS="--bind /cluster/path:/container/path" \
#     ./submit_all_clustering.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/run_clustering.R"
SIF="${APPTAINER_R_IMAGE:?Set APPTAINER_R_IMAGE to the rstudio apptainer image path}"
ATLAS_DATA_DIR="${ATLAS_DATA_DIR:-data/temp_objects}"
OBJ_BASE="${ATLAS_DATA_DIR}/merged_by_cell_type"
OUT_BASE="${ATLAS_DATA_DIR}/clustered_lineages"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "$OUT_BASE" "$LOG_DIR"

# ── Lineage definitions: name | RDS path | mem | time ─────────────────────
# mem/time scaled to object size
declare -A RDS_MAP=(
  [T-cells]="${OBJ_BASE}/Tcells_seperated/Tcell_all_cancers_merged.RDS"
  [Myeloid]="${OBJ_BASE}/Myeloid_merged.RDS"
  [Epithelial]="${OBJ_BASE}/Epithelial_merged.RDS"
  [B-cells]="${OBJ_BASE}/B-cell_merged.RDS"
  [Fibroblasts]="${OBJ_BASE}/Fibroblasts_merged.RDS"
  [Endothelial]="${OBJ_BASE}/Endothelial_merged.RDS"
  [NK-cells]="${OBJ_BASE}/NK-cells_merged.RDS"
  [HSPC]="${OBJ_BASE}/Hematopoietic_Stem_and_Progenitor_cells_merged.RDS"
)

declare -A MEM_MAP=(
  [T-cells]="128G"
  [Myeloid]="128G"
  [Epithelial]="300G"
  [B-cells]="64G"
  [Fibroblasts]="64G"
  [Endothelial]="64G"
  [NK-cells]="32G"
  [HSPC]="32G"
)

declare -A TIME_MAP=(
  [T-cells]="06:00:00"
  [Myeloid]="06:00:00"
  [Epithelial]="12:00:00"
  [B-cells]="03:00:00"
  [Fibroblasts]="03:00:00"
  [Endothelial]="03:00:00"
  [NK-cells]="02:00:00"
  [HSPC]="02:00:00"
)

LINEAGES=(T-cells Myeloid Epithelial B-cells Fibroblasts Endothelial NK-cells HSPC)

echo "Submitting clustering jobs for ${#LINEAGES[@]} lineages..."
echo ""

for LINEAGE in "${LINEAGES[@]}"; do
  RDS="${RDS_MAP[$LINEAGE]}"
  MEM="${MEM_MAP[$LINEAGE]}"
  TIME="${TIME_MAP[$LINEAGE]}"

  if [[ ! -f "$RDS" ]]; then
    echo "  SKIP $LINEAGE: RDS not found at $RDS"
    continue
  fi

  SBATCH_SCRIPT="${SCRIPT_DIR}/job_${LINEAGE}.sbatch"

  cat > "$SBATCH_SCRIPT" << SBATCH
#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=cluster_${LINEAGE}
#SBATCH --output=${LOG_DIR}/cluster_${LINEAGE}.%j.out
#SBATCH --error=${LOG_DIR}/cluster_${LINEAGE}.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}

echo "=== cluster_${LINEAGE} ==="
echo "Node: \$(hostname)"
echo "Start: \$(date)"
echo ""

apptainer exec \\
  ${APPTAINER_BIND_ARGS:-} \\
  ${SIF} \\
  Rscript "${R_SCRIPT}" \\
    --rds "${RDS}" \\
    --lineage "${LINEAGE}" \\
    --outdir "${OUT_BASE}" \\
    --n_pcs 30 \\
    --min_dist 0.3

echo ""
echo "End: \$(date)"
SBATCH

  JOB_ID=$(sbatch "$SBATCH_SCRIPT" | awk '{print $NF}')
  echo "  Submitted $LINEAGE → job $JOB_ID  (mem=$MEM, time=$TIME)"
  echo "    RDS:  $RDS"
  echo "    Log:  ${LOG_DIR}/cluster_${LINEAGE}.${JOB_ID}.out"
  echo ""
done

echo "All jobs submitted. Monitor with:"
echo "  squeue -u $USER"
echo "  tail -f ${LOG_DIR}/cluster_<lineage>.<jobid>.out"
echo ""
echo "Output RDS files will be written to: $OUT_BASE"
echo "Summary files: ${OUT_BASE}/<lineage>_cluster_summary.txt"
