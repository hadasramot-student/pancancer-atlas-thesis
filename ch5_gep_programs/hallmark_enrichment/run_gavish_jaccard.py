#!/usr/bin/env python
# run_gavish_jaccard.py
#
# Ch.5 (GEP programs) — Jaccard-similarity matching of cNMF programs
# (epithelial/tumor compartment) against the Gavish et al. 2023 Nature
# pan-cancer malignant meta-programs (MP_list[["Cancer"]], 41 MPs, top-50
# genes each), downloaded from the 3CA Zenodo archive
# (doi.org/10.5281/zenodo.7688626, ITH_hallmarks/MPs_distribution/MP_list.RDS).
#
# No prior script/output for this analysis could be located anywhere on the
# cluster for the other 10 cancer types, so this is a fresh implementation
# run uniformly across all 11 cancer types (not just Breast) for a
# comparable, apples-to-apples result.
#
# Methodology: top-50 genes per cNMF program (by z-score, same source as the
# Hallmark test) vs. top-50 genes per Gavish "Cancer" MP; Jaccard = |A∩B|/|A∪B|;
# threshold > 0.05 (as specified). Best (max-Jaccard) MP per program is kept.
#
# Usage:
#   python run_gavish_jaccard.py --cancer Breast
#   python run_gavish_jaccard.py --cancer all

import argparse
import os
import pandas as pd

CNMF_BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
GAVISH_CSV = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/gavish_2023/gavish_MP_gene_lists_long.csv"
K_PROGRAMS = 20
DT = "dt_0_1"
N_TOP_GENES = 50
JACCARD_CUTOFF = 0.05

CANCERS = ["Brain", "Colorectal", "Kidney", "Liver_Biliary", "Lung", "Ovarian",
           "Pancreas", "Prostate", "Sarcoma", "Skin", "Breast"]

# Mapping from Gavish "Cancer" MP name -> broad conserved category (the 5
# candidates named in this analysis). MPs not listed here don't map to any of
# the 5 (they're still matched/reported, just not bucketed into a category).
CATEGORY_MAP = {
    "MP6 Hypoxia": "hypoxia",
    "MP17 Interferon/MHC-II (I)": "MHC-II/IFN-gamma",
    "MP18 Interferon/MHC-II (II)": "MHC-II/IFN-gamma",
    "MP12 EMT-I": "EMT",
    "MP13 EMT-II": "EMT",
    "MP14 EMT-III": "EMT",
    "MP15 EMT IV": "EMT",
    # NOTE: Gavish's 41 malignant MPs have no MP literally named
    # "mitochondrial stress" or "NF-kB inflammatory signaling" -- those two
    # categories cannot be derived from this reference list. See README notes.
}


def load_gavish_cancer_mps():
    df = pd.read_csv(GAVISH_CSV)
    df = df[df["cell_type"] == "Cancer"]
    return {mp: set(sub["gene"]) for mp, sub in df.groupby("meta_program")}


def load_gene_spectra_score(cancer):
    path = os.path.join(CNMF_BASE, cancer, f"cnmf_{cancer}",
                         f"cnmf_{cancer}.gene_spectra_score.k_{K_PROGRAMS}.{DT}.txt")
    return pd.read_csv(path, sep="\t", index_col=0)


def jaccard(a, b):
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def run_cancer(cancer, gavish_mps):
    df = load_gene_spectra_score(cancer)
    rows = []
    for prog_idx in df.index:
        top_genes = set(df.loc[prog_idx].sort_values(ascending=False).head(N_TOP_GENES).index)
        best_mp, best_j = None, 0.0
        for mp, mp_genes in gavish_mps.items():
            j = jaccard(top_genes, mp_genes)
            if j > best_j:
                best_mp, best_j = mp, j
        matched = best_j > JACCARD_CUTOFF
        rows.append({
            "cancer": cancer,
            "program": f"cNMF_{prog_idx}",
            "best_gavish_mp": best_mp if matched else "Unmatched",
            "jaccard": best_j,
            "matched": matched,
            "conserved_category": CATEGORY_MAP.get(best_mp, "cancer-type-specific") if matched else "cancer-type-specific",
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer", default="all")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    gavish_mps = load_gavish_cancer_mps()
    cancers = CANCERS if args.cancer == "all" else [args.cancer]
    dfs = [run_cancer(c, gavish_mps) for c in cancers]
    out = pd.concat(dfs, ignore_index=True)
    out_path = args.out or ("all_cancers_gavish_jaccard.csv" if args.cancer == "all"
                             else f"{args.cancer}_gavish_jaccard.csv")
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} rows -> {out_path}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
