#!/bin/bash
#============================================================
# submit_spectra.sh
#
# Ch.5 (GEP programs) — end-to-end Spectra pipeline for the immune
# lineages, one dependency chain per lineage:
#
#   A. RDS -> h5ad   (R, per cancer file; huge lineages subsampled)
#   B. merge h5ad    (Python CPU, per lineage; + stratified subsample)
#   C. export Hallmark JSON   (R, once)
#   D. build priors  (Python CPU)
#   E. run Spectra   (Python GPU)
#
# The R steps run inside an apptainer container that may see the host
# filesystem at a different mount point than this shell (e.g. host
# /mnt/new_groups/... bound to /ofir_grp inside the container) — set
# ATLAS_GROUP_ROOT_HOST / ATLAS_GROUP_ROOT_CTR if that applies to you.
#
# Usage:
#   ATLAS_GROUP_ROOT_HOST=/mnt/group ATLAS_GROUP_ROOT_CTR=/ctr_bind \
#     bash submit_spectra.sh [--dry-run] [--lineages "T-cells Myeloid"]
#============================================================
set -euo pipefail

# ---- paths (R sees ATLAS_GROUP_ROOT_CTR, Python/shell see ATLAS_GROUP_ROOT_HOST; same files) ----
G_MNT="${ATLAS_GROUP_ROOT_HOST:?Set ATLAS_GROUP_ROOT_HOST, e.g. /mnt/new_groups/ofircohen-group/users/\$USER}"
G_R="${ATLAS_GROUP_ROOT_CTR:-$G_MNT}"
RDS_R=${G_R}/temp_objects/clustered_lineages
RDS_MNT=${G_MNT}/temp_objects/clustered_lineages
H5AD_MNT=${RDS_MNT}/h5ad
H5AD_R=${RDS_R}/h5ad
SP_MNT="${ATLAS_REPO_DIR:-.}/ch5_gep_programs/spectra"
SP_R="${ATLAS_REPO_DIR_CTR:-${ATLAS_REPO_DIR:-.}}/ch5_gep_programs/spectra"
OUT=${RDS_MNT}/spectra_out
PRIORS=${SP_MNT}/priors
LOGS=${SP_MNT}/logs

SBATCH_R=~/rstudio/rstudio.sbatch
SBATCH_PY=~/rstudio/python.sbatch
CELLTYPE_KEY=scimilarity_celltype

mkdir -p "${H5AD_MNT}" "${OUT}" "${PRIORS}" "${LOGS}"

# All immune lineages + their TME compartment + per-file/merge subsample caps.
# (small lineages -> full; T-cells / Myeloid -> subsampled)
declare -A COMPARTMENT=( [T-cells]=T [Myeloid]=Myeloid [B-cell]=B [NK-cells]=NK
                         [Hematopoietic_Stem_and_Progenitor_cells]=HSPC )
declare -A PERFILE_MAX=( [T-cells]=60000 [Myeloid]=60000 [B-cell]=0 [NK-cells]=0
                         [Hematopoietic_Stem_and_Progenitor_cells]=0 )
declare -A MERGE_MAX=(   [T-cells]=200000 [Myeloid]=200000 [B-cell]=0 [NK-cells]=0
                         [Hematopoietic_Stem_and_Progenitor_cells]=0 )
declare -A CONV_MEM=(    [T-cells]=384G [Myeloid]=384G [B-cell]=64G [NK-cells]=96G
                         [Hematopoietic_Stem_and_Progenitor_cells]=64G )

LINEAGES_DEFAULT="T-cells Myeloid B-cell NK-cells Hematopoietic_Stem_and_Progenitor_cells"
LINEAGES="${LINEAGES_DEFAULT}"
DRY=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)  DRY=true; shift ;;
    --lineages) LINEAGES="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

sb() { if $DRY; then echo "[DRY] sbatch $*" >&2; echo "$RANDOM"; else sbatch --parsable "$@"; fi; }

# ---- C. Hallmark export (once, shared) ----
HALLMARK_JSON_MNT=${PRIORS}/hallmark.json
HALLMARK_JSON_R=${SP_R}/priors/hallmark.json
if [[ -f "${HALLMARK_JSON_MNT}" ]]; then
  echo "[C] Hallmark JSON exists, skipping export"; JOB_C=""
else
  JOB_C=$(sb --job-name=sp_hallmark --partition=main --time=0-00:30:00 --mem=8G --cpus-per-task=2 \
    --output=${LOGS}/hallmark_%j.out --error=${LOGS}/hallmark_%j.err \
    ${SBATCH_R} --script ${SP_R}/export_hallmark_gmt.R --out ${HALLMARK_JSON_R})
  echo "[C] Hallmark export job=${JOB_C}"
fi
dep_C=""; [[ -n "${JOB_C}" ]] && dep_C="afterok:${JOB_C}"

for LIN in ${LINEAGES}; do
  echo "================ ${LIN} ================"
  comp=${COMPARTMENT[$LIN]}; pf=${PERFILE_MAX[$LIN]}; mm=${MERGE_MAX[$LIN]}; cmem=${CONV_MEM[$LIN]}

  # ---- A. convert each per-cancer RDS for this lineage ----
  conv_jobs=()
  shopt -s nullglob
  for rds in ${RDS_MNT}/${LIN}_*_clustered.RDS; do
    base=$(basename "$rds" .RDS)               # e.g. T-cells_Lung_clustered
    out_h5ad_R=${H5AD_R}/${base}.h5ad
    out_h5ad_MNT=${H5AD_MNT}/${base}.h5ad
    [[ -f "${out_h5ad_MNT}" ]] && { echo "  [A skip] ${base}"; continue; }
    J=$(sb --job-name=sp_h5ad_${base:0:18} --partition=main --time=0-08:00:00 \
        --mem=${cmem} --cpus-per-task=8 \
        --output=${LOGS}/conv_${base}_%j.out --error=${LOGS}/conv_${base}_%j.err \
        ${SBATCH_R} --script ${SP_R}/seurat_to_h5ad_spectra.R \
          --input_rds ${RDS_R}/${base}.RDS --output_h5ad ${out_h5ad_R} \
          --celltype_key ${CELLTYPE_KEY} --max_cells ${pf})
    conv_jobs+=("$J"); echo "  [A] ${base} job=${J}"
  done
  shopt -u nullglob

  dep_A=""
  if [[ ${#conv_jobs[@]} -gt 0 ]]; then
    dep_A="afterok:$(IFS=:; echo "${conv_jobs[*]}")"
  fi

  # ---- B. merge per-lineage h5ad (+ subsample) ----
  MERGED=${OUT}/${LIN}_immune.h5ad
  depB_parts=(); [[ -n "${dep_A}" ]] && depB_parts+=("${dep_A}")
  depB=""; [[ ${#depB_parts[@]} -gt 0 ]] && depB="--dependency=$(IFS=,; echo "${depB_parts[*]}")"
  JOB_B=$(sb --job-name=sp_merge_${LIN:0:14} --partition=main --time=0-06:00:00 \
      --mem=256G --cpus-per-task=8 ${depB} \
      --output=${LOGS}/merge_${LIN}_%j.out --error=${LOGS}/merge_${LIN}_%j.err \
      ${SBATCH_PY} --script ${SP_MNT}/merge_lineage_h5ad.py \
        --lineage ${LIN} --h5ad-dir ${H5AD_MNT} --out ${MERGED} \
        --cell-type-key ${CELLTYPE_KEY} --max-cells ${mm} --min-cells-per-type 50)
  echo "  [B] merge job=${JOB_B}"

  # ---- D. build priors ----
  GSD=${PRIORS}/${LIN}_gene_set_dictionary.json
  depD_parts=("afterok:${JOB_B}"); [[ -n "${dep_C}" ]] && depD_parts+=("${dep_C}")
  JOB_D=$(sb --job-name=sp_prior_${LIN:0:14} --partition=main --time=0-00:30:00 \
      --mem=16G --cpus-per-task=2 --dependency=$(IFS=,; echo "${depD_parts[*]}") \
      --output=${LOGS}/prior_${LIN}_%j.out --error=${LOGS}/prior_${LIN}_%j.err \
      ${SBATCH_PY} --script ${SP_MNT}/build_spectra_priors.py \
        --h5ad ${MERGED} --cell-type-key ${CELLTYPE_KEY} --compartment ${comp} \
        --hallmark-json ${HALLMARK_JSON_MNT} --tme-json ${SP_MNT}/tme_gene_sets.json \
        --out-json ${GSD})
  echo "  [D] priors job=${JOB_D}"

  # ---- E. run Spectra on GPU ----
  JOB_E=$(sb --job-name=sp_run_${LIN:0:16} --partition=main --time=1-12:00:00 \
      --gres=gpu:rtx_6000:1 --mem=128G --cpus-per-task=8 \
      --dependency=afterok:${JOB_D} \
      --output=${LOGS}/spectra_${LIN}_%j.out --error=${LOGS}/spectra_${LIN}_%j.err \
      ${SBATCH_PY} --script ${SP_MNT}/run_spectra.py \
        --h5ad ${MERGED} --gene-set-json ${GSD} --cell-type-key ${CELLTYPE_KEY} \
        --out-prefix ${OUT}/${LIN}_spectra --epochs 150 --n-hvg 4000)
  echo "  [E] spectra GPU job=${JOB_E}"
done

echo ""
echo "Submitted. Monitor: squeue -u $USER | grep sp_"
echo "Outputs -> ${OUT}/   |  logs -> ${LOGS}/"
