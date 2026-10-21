#!/bin/bash
# run_epithelial_breast_subsample_pipeline.sh
#
# Chapter 4 gap-fix for Epithelial_Breast: subsample to 450,000 cells, then run
# the *unchanged* Ch.4 pipeline used for every other lineage x cancer-type object.
#
# Chain (SLURM --dependency=afterok):
#   1. subsample   subsample_epithelial_breast.R      -> ..._subsample450k_..._scanvi_final.RDS
#   2. cluster     merge/cluster_lineage.R            -> clustered_lineages/Epithelial_Breast_clustered.RDS
#                                                        + Epithelial_Breast_cluster_summary.tsv
#   3. annotate    merge/annotate_lineages.R          -> Epithelial_Breast_scimilarity_summary.tsv
#   4. markers     merge/find_markers.R               -> Epithelial_Breast_markers.csv / _top_markers.csv
#   5. gsea        merge/gsea_hallmarks.R             -> FindMarkers_results/enrichment/Epithelial_Breast_hallmarks_gsea.csv
#
# Usage:  ./run_epithelial_breast_subsample_pipeline.sh [--dry-run]

set -euo pipefail

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# ── Paths (container view: /groups/ofircohen-group is bound to /ofir_grp) ────
REPO_R=/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/epithelial-breast-subsample-cluster-08231a
MERGE_R=/ofir_grp/users/hadasra/create_sc_object_workflow/merge
SBATCH_R=~/rstudio/rstudio.sbatch

BREAST_DIR_R=/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast
IN_FULL="${BREAST_DIR_R}/Epithelial_Breast_merged_full_scanvi_final.RDS"
SUB_RDS="${BREAST_DIR_R}/Epithelial_Breast_subsample450k_merged_full_scanvi_final.RDS"

CLUST_DIR_R=/ofir_grp/users/hadasra/temp_objects/clustered_lineages
CLUST_RDS="${CLUST_DIR_R}/Epithelial_Breast_clustered.RDS"
SUM_TSV="${CLUST_DIR_R}/Epithelial_Breast_cluster_summary.tsv"
GSEA_OUT_R="${CLUST_DIR_R}/FindMarkers_results/enrichment"

# SLURM --output/--error are resolved by slurmd in the HOST namespace, so they
# must be real host paths (/groups/...), NOT the container-only /ofir_grp bind.
LOG_HOST=/groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages/epithelial_breast_subsample_logs
LOG_R="${LOG_HOST}"
mkdir -p "${LOG_HOST}"

LINEAGE=Epithelial_Breast
N=450000
SEED=42

maybe_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    { echo "[DRY-RUN] sbatch $*"; } >&2
    echo "$(( RANDOM + 100000 ))"
  else
    sbatch --parsable "$@"
  fi
}

echo "======================================================================"
echo " Epithelial_Breast — subsample (${N}, seed ${SEED}) + Ch.4 pipeline"
echo " Dry run: ${DRY_RUN}"
echo " Logs   : ${LOG_HOST}/"
echo "======================================================================"

# ── 1. Subsample ───────────────────────────────────────────────────────────
J1=$(maybe_sbatch \
  --job-name="ebsub_1_subsample" \
  --partition=main --time="0-08:00:00" --mem=300G --cpus-per-task=4 \
  --output="${LOG_R}/1_subsample_%j.out" --error="${LOG_R}/1_subsample_%j.err" \
  "${SBATCH_R}" --script "${REPO_R}/subsample_epithelial_breast.R" \
    --input "${IN_FULL}" \
    --output "${SUB_RDS}" \
    --report_dir "${CLUST_DIR_R}" \
    --n "${N}" --seed "${SEED}" --strat_col project_name)
echo "  [1] subsample   -> job ${J1}"

# ── 2. Cluster (PCA on HVGs + UMAP + Louvain @ 0.25/0.30/0.50/0.75/1.00/2.00) ─
J2=$(maybe_sbatch \
  --job-name="ebsub_2_cluster" \
  --dependency=afterok:${J1} \
  --partition=main --time="1-00:00:00" --mem=384G --cpus-per-task=8 \
  --output="${LOG_R}/2_cluster_%j.out" --error="${LOG_R}/2_cluster_%j.err" \
  "${SBATCH_R}" --script "${MERGE_R}/cluster_lineage.R" \
    --input "${SUB_RDS}" \
    --output "${CLUST_RDS}" \
    --lineage "${LINEAGE}" \
    --n_pcs 30 --min_dist 0.3 \
    --summary_out "${SUM_TSV}")
echo "  [2] cluster     -> job ${J2}"

# ── 3. Scimilarity annotation (majority vote per cluster, reuse existing labels) ─
J3=$(maybe_sbatch \
  --job-name="ebsub_3_annotate" \
  --dependency=afterok:${J2} \
  --partition=main --time="0-12:00:00" --mem=384G --cpus-per-task=8 \
  --output="${LOG_R}/3_annotate_%j.out" --error="${LOG_R}/3_annotate_%j.err" \
  "${SBATCH_R}" --script "${MERGE_R}/annotate_lineages.R" \
    --input_dir "${CLUST_DIR_R}" \
    --lineage "${LINEAGE}" \
    --overwrite)
echo "  [3] annotate    -> job ${J3}"

# ── 4. FindAllMarkers (Wilcoxon @ louvain_res0.50) ─────────────────────────
J4=$(maybe_sbatch \
  --job-name="ebsub_4_markers" \
  --dependency=afterok:${J3} \
  --partition=main --time="1-12:00:00" --mem=512G --cpus-per-task=8 \
  --output="${LOG_R}/4_markers_%j.out" --error="${LOG_R}/4_markers_%j.err" \
  "${SBATCH_R}" --script "${MERGE_R}/find_markers.R" \
    --input_dir "${CLUST_DIR_R}" \
    --lineage "${LINEAGE}" \
    --n_cores 8)
echo "  [4] markers     -> job ${J4}"

# ── 5. Hallmark GSEA (fgsea, MSigDB H) per cluster ─────────────────────────
J5=$(maybe_sbatch \
  --job-name="ebsub_5_gsea" \
  --dependency=afterok:${J4} \
  --partition=main --time="0-03:00:00" --mem=48G --cpus-per-task=4 \
  --output="${LOG_R}/5_gsea_%j.out" --error="${LOG_R}/5_gsea_%j.err" \
  "${SBATCH_R}" --script "${MERGE_R}/gsea_hallmarks.R" \
    --lineage "${LINEAGE}" \
    --input_dir "${CLUST_DIR_R}" \
    --output_dir "${GSEA_OUT_R}")
echo "  [5] gsea        -> job ${J5}"

echo ""
echo "Submitted chain: ${J1} -> ${J2} -> ${J3} -> ${J4} -> ${J5}"
echo "Monitor : squeue -u hadasra | grep ebsub_"
echo "Logs    : ${LOG_HOST}/"
echo ""
echo "Deliverables when complete:"
echo "  ${CLUST_DIR_R}/Epithelial_Breast_cluster_summary.tsv"
echo "  ${CLUST_DIR_R}/Epithelial_Breast_scimilarity_summary.tsv"
echo "  ${CLUST_DIR_R}/Epithelial_Breast_markers.csv / _top_markers.csv"
echo "  ${GSEA_OUT_R}/Epithelial_Breast_hallmarks_gsea.csv"
echo "  ${CLUST_DIR_R}/Epithelial_Breast_subsample_composition.csv   (legend / methods)"
echo "  ${CLUST_DIR_R}/Epithelial_Breast_SUBSAMPLE_NOTE.txt"
