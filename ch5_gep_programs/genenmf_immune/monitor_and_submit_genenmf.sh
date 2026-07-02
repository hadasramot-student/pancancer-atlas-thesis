#!/bin/bash
#============================================================
# monitor_and_submit_genenmf.sh
#
# Ch.5 (GEP programs) — idempotent watcher for the per-lineage immune
# GeneNMF run (the GeneNMF side of the Spectra-vs-GeneNMF immune
# cross-method comparison).
#
# For each immune lineage, GeneNMF must consume the SAME merged h5ad that
# Spectra uses (spectra_out/<LIN>_immune.h5ad), which the Spectra pipeline's
# convert+merge steps produce. This script checks, per lineage, whether that
# merged h5ad is ready and — if GeneNMF hasn't already been run or queued —
# submits the per-lineage GeneNMF job via submit_genenmf_immune.sh.
#
# Safe to run repeatedly (e.g. from cron / a session loop): it skips lineages
# that are DONE (output rds exists), IN PROGRESS (gn_ job in squeue), or
# NOT READY (merged h5ad missing or still being written).
#
# Prints one status line per lineage and a final ALL_DONE / PENDING marker so a
# caller can stop looping once everything is submitted.
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects ATLAS_REPO_DIR=/path/to/repo \
#     bash monitor_and_submit_genenmf.sh
#============================================================
set -euo pipefail

ATLAS_DATA_DIR="${ATLAS_DATA_DIR:-data/temp_objects}"
ATLAS_REPO_DIR="${ATLAS_REPO_DIR:-.}"
SPECTRA_OUT=${ATLAS_DATA_DIR}/clustered_lineages/spectra_out
OUT=${ATLAS_DATA_DIR}/genenmf_immune_output
SUBMIT=${ATLAS_REPO_DIR}/ch5_gep_programs/genenmf_immune/submit_genenmf_immune.sh

# B-cell already done; these are the four still pending.
LINEAGES="T-cells Myeloid NK-cells Hematopoietic_Stem_and_Progenitor_cells"

# A merged h5ad is "stable" only if its mtime is older than this many seconds,
# to avoid reading a file the merge job is still writing.
STABLE_AGE=180

pending=0
for LIN in ${LINEAGES}; do
  h5ad=${SPECTRA_OUT}/${LIN}_immune.h5ad
  done_rds=${OUT}/${LIN}/${LIN}_genenmf_metaprograms.rds

  if [[ -f "${done_rds}" ]]; then
    echo "  [DONE]        ${LIN} (metaprograms.rds present)"
    continue
  fi

  # gn_ job names are truncated to 16 chars in submit_genenmf_immune.sh.
  jobname="gn_${LIN:0:16}"
  if squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -qx "${jobname}"; then
    echo "  [IN PROGRESS] ${LIN} (job ${jobname} in queue)"
    pending=1
    continue
  fi

  if [[ ! -f "${h5ad}" ]]; then
    echo "  [NOT READY]   ${LIN} (no merged h5ad yet)"
    pending=1
    continue
  fi

  # File exists but may still be flushing — require a stable mtime.
  age=$(( $(date +%s) - $(stat -c %Y "${h5ad}") ))
  if (( age < STABLE_AGE )); then
    echo "  [SETTLING]    ${LIN} (h5ad age ${age}s < ${STABLE_AGE}s; will pick up next pass)"
    pending=1
    continue
  fi

  echo "  [SUBMIT]      ${LIN} (merged h5ad ready, ${age}s old) -> launching GeneNMF"
  bash "${SUBMIT}" --lineages "${LIN}"
  pending=1   # treat as pending until its output rds lands
done

if (( pending == 0 )); then
  echo "ALL_DONE: all pending immune lineages have GeneNMF metaprograms.rds outputs."
else
  echo "PENDING: at least one lineage still waiting/running/submitting."
fi
