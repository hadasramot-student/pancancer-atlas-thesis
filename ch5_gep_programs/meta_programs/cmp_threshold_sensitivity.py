#!/usr/bin/env python
"""
Full-range Jaccard-threshold sensitivity for the cMP (consensus meta-program)
clustering of all 240 cNMF programs (12 cancer types x K=20), extending the
existing two-point check (0.05, 0.1) in meta_program_catalogue.py.

Reuses the exact same loading + greedy_cluster logic as
meta_program_catalogue.py, but skips the per-MP Hallmark enrichment and
Gavish-match annotation steps (not needed for counts) so all 7 thresholds run
quickly.
"""
import os
import numpy as np
import pandas as pd

CNMF_BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
K_PROGRAMS = 20
DT = "dt_0_1"
N_TOP_GENES = 50
MIN_CANCER_TYPES = 2
THRESHOLDS = [0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.5]

CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney",
           "Liver_Biliary", "Lung", "Ovarian", "Pancreas", "Prostate",
           "Sarcoma", "Skin"]

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def load_gene_spectra_score(cancer):
    path = os.path.join(CNMF_BASE, cancer, f"cnmf_{cancer}",
                         f"cnmf_{cancer}.gene_spectra_score.k_{K_PROGRAMS}.{DT}.txt")
    return pd.read_csv(path, sep="\t", index_col=0)


def load_all_programs():
    programs = {}
    for cancer in CANCERS:
        df = load_gene_spectra_score(cancer)
        for prog_idx in df.index:
            ranked = df.loc[prog_idx].sort_values(ascending=False)
            top_genes = set(ranked.head(N_TOP_GENES).index)
            key = f"{cancer}_cNMF_{prog_idx}"
            programs[key] = {"cancer": cancer, "top_set": top_genes}
    return programs


def jaccard(a, b):
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def build_jaccard_matrix(programs):
    keys = list(programs.keys())
    n = len(keys)
    J = np.zeros((n, n))
    sets = [programs[k]["top_set"] for k in keys]
    for i in range(n):
        for j in range(i + 1, n):
            j_val = jaccard(sets[i], sets[j])
            J[i, j] = j_val
            J[j, i] = j_val
    return keys, J


def greedy_cluster(keys, J, threshold):
    n = len(keys)
    pool = set(range(n))
    groups = []
    while pool:
        pool_list = list(pool)
        best_i, best_count, best_partners = None, -1, None
        for i in pool_list:
            partners = [j for j in pool_list if j != i and J[i, j] > threshold]
            if len(partners) > best_count:
                best_i, best_count, best_partners = i, len(partners), partners
        if best_count <= 0:
            break
        group = [best_i] + best_partners
        groups.append(group)
        pool -= set(group)
    singletons = sorted(pool)
    return groups, singletons


def max_group_size_hist(valid_mps):
    sizes = [len(g) for g, _ in valid_mps]
    return (min(sizes) if sizes else None,
            float(np.median(sizes)) if sizes else None,
            max(sizes) if sizes else None)


def main():
    print("Loading 240 cNMF programs (top-50 genes each)...")
    programs = load_all_programs()
    print(f"  {len(programs)} programs loaded")

    print("Computing pairwise Jaccard matrix (240x240)...")
    key_list, J = build_jaccard_matrix(programs)

    rows = []
    for thr in THRESHOLDS:
        groups, singleton_idx = greedy_cluster(key_list, J, thr)

        valid_mps = []
        n_single_cancer_clusters = 0
        for g in groups:
            cancers_in_group = sorted(set(programs[key_list[i]]["cancer"] for i in g))
            if len(cancers_in_group) >= MIN_CANCER_TYPES:
                valid_mps.append((g, cancers_in_group))
            else:
                n_single_cancer_clusters += 1

        n_pan_cancer = sum(1 for _, c in valid_mps if len(c) >= 7)
        minsz, medsz, maxsz = max_group_size_hist(valid_mps)

        row = {
            "threshold": thr,
            "n_cMPs": len(valid_mps),
            "n_pan_cancer_cMPs (>=7 cancers)": n_pan_cancer,
            "n_single_cancer_type_clusters": n_single_cancer_clusters,
            "n_singletons": len(singleton_idx),
            "n_raw_groups_before_recurrence_filter": len(groups),
            "cMP_size_min": minsz,
            "cMP_size_median": medsz,
            "cMP_size_max": maxsz,
        }
        rows.append(row)
        print(f"thr={thr}: cMPs={len(valid_mps)}, pan-cancer={n_pan_cancer}, "
              f"single-cancer-only={n_single_cancer_clusters}, singletons={len(singleton_idx)}")

    out = pd.DataFrame(rows)
    out_csv = os.path.join(OUTDIR, "cmp_threshold_sensitivity.csv")
    out.to_csv(out_csv, index=False)
    print(f"\nWrote: {out_csv}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
