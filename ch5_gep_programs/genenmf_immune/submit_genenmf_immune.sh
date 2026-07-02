#!/bin/bash
#============================================================
# submit_genenmf_immune.sh
#
# Ch.5 (GEP programs) — run GeneNMF on the SAME merged immune h5ads used by
# Spectra, one R job per lineage, for the Spectra-vs-GeneNMF immune
# cross-comparison (same cells).
#
# Requires <lineage>_immune.h5ad to already exist (produced by the Spectra
# pipeline's merge step). If the Spectra gene_factors CSV exists, a Jaccard
# MP-vs-Spectra comparison is also written.
#
# The R step runs inside an apptainer container that may see the host
# filesystem at a different mount point than this shell; set
# ATLAS_GROUP_ROOT_HOST / ATLAS_GROUP_ROOT_CTR / ATLAS_REPO_DIR_CTR if so.
#
# Usage:
#   ATLAS_GROUP_ROOT_HOST=/mnt/group ATLAS_REPO_DIR=/path/to/repo \
#     bash submit_genenmf_immune.sh [--dry-run] [--lineages "B-cell NK-cells"]
#============================================================
set -euo pipefail

G_MNT="${ATLAS_GROUP_ROOT_HOST:?Set ATLAS_GROUP_ROOT_HOST, e.g. /mnt/new_groups/ofircohen-group/users/\$USER}"
G_R="${ATLAS_GROUP_ROOT_CTR:-$G_MNT}"
ATLAS_REPO_DIR="${ATLAS_REPO_DIR:-.}"
ATLAS_REPO_DIR_CTR="${ATLAS_REPO_DIR_CTR:-$ATLAS_REPO_DIR}"
SPECTRA_OUT_MNT=${G_MNT}/temp_objects/clustered_lineages/spectra_out
SPECTRA_OUT_R=${G_R}/temp_objects/clustered_lineages/spectra_out
GN_DIR_R=${ATLAS_REPO_DIR_CTR}/ch5_gep_programs/genenmf_immune
OUT_R=${G_R}/temp_objects/genenmf_immune_output
LOGS=${ATLAS_REPO_DIR}/ch5_gep_programs/genenmf_immune/logs
SBATCH_R=~/rstudio/rstudio.sbatch

mkdir -p "${LOGS}"

LINEAGES="T-cells Myeloid B-cell NK-cells Hematopoietic_Stem_and_Progenitor_cells"
DRY=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)  DRY=true; shift ;;
    --lineages) LINEAGES="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done
sb() { if $DRY; then echo "[DRY] sbatch $*" >&2; echo "$RANDOM"; else sbatch --parsable "$@"; fi; }

for LIN in ${LINEAGES}; do
  H5AD_MNT=${SPECTRA_OUT_MNT}/${LIN}_immune.h5ad
  if [[ ! -f "${H5AD_MNT}" ]]; then
    echo "  [SKIP] ${LIN}: merged h5ad not found yet (${H5AD_MNT})"; continue
  fi
  SF_R=${SPECTRA_OUT_R}/${LIN}_spectra_gene_factors.csv   # may not exist yet -> Jaccard skipped
  J=$(sb --job-name=gn_${LIN:0:16} --partition=main --time=0-12:00:00 \
      --mem=160G --cpus-per-task=12 \
      --output=${LOGS}/genenmf_${LIN}_%j.out --error=${LOGS}/genenmf_${LIN}_%j.err \
      ${SBATCH_R} --script ${GN_DIR_R}/run_genenmf_per_lineage.R \
        --lineage ${LIN} \
        --h5ad ${SPECTRA_OUT_R}/${LIN}_immune.h5ad \
        --spectra_factors ${SF_R} \
        --output_dir ${OUT_R})
  echo "  [SUBMIT] GeneNMF ${LIN} job=${J}"
done
echo "Monitor: squeue -u $USER | grep gn_   |  outputs -> ${OUT_R}/"
