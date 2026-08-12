#!/usr/bin/env python
# run_cnmf_hallmark.py
#
# Ch.5 (GEP programs) — per-program MSigDB Hallmark hypergeometric enrichment
# for cNMF programs (epithelial/tumor cell compartment).
#
# Reconstructed methodology (reverse-engineered from the existing 200-row
# ids_output/figures/cnmf_top_hallmark.csv covering Brain, Colorectal, Kidney,
# Liver_Biliary, Lung, Ovarian, Pancreas, Prostate, Sarcoma, Skin — the
# original generating script could not be located on the cluster):
#   - top-N genes per program = top 50 genes by z-score from
#     cnmf_<cancer>.gene_spectra_score.k_20.dt_0_1.txt
#   - background = all genes scored by cNMF for that cancer (28,231 genes,
#     confirmed identical across all 11 cancer types incl. Breast)
#   - gene sets = MSigDB Hallmark (H) collection via msigdbr (50 sets)
#   - one-sided hypergeometric test: P(X >= overlap) = hypergeom.sf(overlap-1, N, K, n)
#   - top_hallmark = gene set with the lowest p-value; if that p >= 0.05 the
#     row is labeled "Unannotated" but retains the actual best-hit pval/overlap
#     (confirmed from the existing file: max named pval = 0.0497, min
#     Unannotated pval = 0.0539 — clean p<0.05 cutoff)
#
# CONFIRMED (via two independent consistency checks against the existing 200
# rows) that top_hallmark IS a genuine best/minimum-p-value match across all
# 50 Hallmark sets, not an arbitrary pick:
#   1. Whenever a set with higher raw overlap is achievable but NOT reported,
#      it always has a much larger K (gene-set size) than the reported set —
#      exactly what minimum-p-value (not maximum-overlap) selection predicts.
#   2. There is a hard floor: overlap < 3 is always reported as overlap=0,
#      pval=1.0, top_hallmark="Unannotated" (verified on all 65 overlap=0
#      rows: pval is EXACTLY 1.0 in every case, not just close to it).
#
# The `high_confidence` column is NOT derived from the hallmark test itself —
# validated exactly (0 mismatches / 200 existing rows) against
# genenmf_output/<cancer>/<cancer>_highconfidence_programs.csv: a cNMF program
# is high_confidence iff it is the `best_cnmf_match` for >=1 GeneNMF
# meta-program with `recovered_by_both == TRUE` (i.e. cross-method,
# cNMF-vs-GeneNMF Jaccard-validated, independent of Hallmark significance).
#
# CAVEAT: spot-checks of the Hallmark test itself against the existing
# 200-row file reproduce large/robust hits exactly (e.g. Kidney cNMF_6
# overlap=30) but differ on some small (2-4 gene) borderline overlaps — most
# likely because the original run used a different msigdbr/MSigDB Hallmark
# snapshot than the one currently installed (msigdbr 26.1.0, the only version
# present on this cluster). This script does NOT recompute the existing 200
# rows; it only adds Breast.
#
# Usage:
#   python run_cnmf_hallmark.py --cancer Breast

import argparse
import json
import os
import pandas as pd
from scipy.stats import hypergeom

CNMF_BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
GENENMF_BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output"
HALLMARK_JSON = "/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/merge/spectra/priors/hallmark.json"
K_PROGRAMS = 20
DT = "dt_0_1"
N_TOP_GENES = 50
PVAL_CUTOFF = 0.05
MIN_OVERLAP = 3  # overlaps below this are floored to overlap=0, pval=1.0 (confirmed convention)


def load_high_confidence_set(cancer):
    path = os.path.join(GENENMF_BASE, cancer, f"{cancer}_highconfidence_programs.csv")
    if not os.path.exists(path):
        # No GeneNMF run exists for this cancer type, so cross-method
        # validation is not possible. Return None (not an empty set) so the
        # caller can write NA rather than a false "False".
        print(f"  WARNING: no GeneNMF cross-validation file for {cancer} "
              f"({path}) -- high_confidence will be written as NA, not False")
        return None
    hc = pd.read_csv(path)
    return set(hc.loc[hc["recovered_by_both"] == True, "best_cnmf_match"])


def load_hallmark():
    with open(HALLMARK_JSON) as f:
        h = json.load(f)
    return {k.replace("HALLMARK_", ""): set(v) for k, v in h.items()}


def load_gene_spectra_score(cancer):
    path = os.path.join(CNMF_BASE, cancer, f"cnmf_{cancer}",
                         f"cnmf_{cancer}.gene_spectra_score.k_{K_PROGRAMS}.{DT}.txt")
    return pd.read_csv(path, sep="\t", index_col=0)


def run_cancer(cancer, hallmark):
    df = load_gene_spectra_score(cancer)
    n_background = df.shape[1]
    hc_set = load_high_confidence_set(cancer)
    rows = []
    for prog_idx in df.index:
        top_genes = set(df.loc[prog_idx].sort_values(ascending=False).head(N_TOP_GENES).index)
        best_name, best_pval, best_overlap = "Unannotated", 1.0, 0
        for hname, hset in hallmark.items():
            overlap = len(top_genes & hset)
            if overlap == 0:
                continue
            pval = hypergeom.sf(overlap - 1, n_background, len(hset), N_TOP_GENES)
            if pval < best_pval:
                best_name, best_pval, best_overlap = hname, pval, overlap
        if best_overlap < MIN_OVERLAP:
            best_name, best_pval, best_overlap = "Unannotated", 1.0, 0
        top_hallmark = best_name if best_pval < PVAL_CUTOFF else "Unannotated"
        program = f"cNMF_{prog_idx}"
        rows.append({
            "cancer": cancer,
            "program": program,
            "top_hallmark": top_hallmark,
            "pval": best_pval,
            "overlap": best_overlap,
            "high_confidence": (program in hc_set) if hc_set is not None else pd.NA,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    hallmark = load_hallmark()
    out = run_cancer(args.cancer, hallmark)
    out_path = args.out or f"{args.cancer}_cnmf_top_hallmark.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} rows -> {out_path}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
