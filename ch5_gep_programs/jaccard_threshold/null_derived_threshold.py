#!/usr/bin/env python
"""
Null/permutation-derived Jaccard threshold for cNMF-GeneNMF HC designation
(replaces the fixed 0.1 convention with a significance-based cutoff).

Two independent null distributions of "chance" top-50-gene Jaccard scores:

  A. Gene-label permutation null: draw two independent random 50-gene lists
     from the empirical gene pool (frequency-weighted by how often each gene
     appears across ALL observed cNMF + GeneNMF top-50 lists, not uniform
     over the full expressed background), compute Jaccard. Repeat N_PERM
     times. Captures "how much overlap do two random top-50 lists get just
     from shared usual-suspect genes (ribosomal/MT/cell-cycle) being common
     in the pool".

  B. Cross-cancer-type shuffled-pairing null: real (non-randomized) top-50
     lists, but cNMF program from cancer A paired with GeneNMF MP from a
     DIFFERENT cancer type B != A. These pairs have no a priori reason to
     represent the same biology. Captures the full real-data null including
     any structure permutation A destroys.

Convention (matching the existing pipeline / prior sensitivity sweep):
HC count = number of GeneNMF meta-programs (across the 11 cancer types with
a GeneNMF run) whose best Jaccard match to a same-cancer cNMF program clears
the threshold. This is the same metric that produced 44 at Jaccard>0.1 in
the prior sweep.
"""
import os
import numpy as np
import pandas as pd

BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output"
CANCERS = ["Kidney", "Brain", "Colorectal", "Pancreas", "Breast",
           "Liver_Biliary", "Prostate", "Lung", "Ovarian", "Sarcoma", "Skin"]
N_PERM = 20000
TOP_N = 50
RNG = np.random.default_rng(0)

OUTDIR = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output/jaccard_threshold_sensitivity"


def jaccard(a, b):
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


# ---------------------------------------------------------------- load data
def load_lists():
    cnmf = {}       # cancer -> {prog_id: set(genes)}
    genenmf = {}    # cancer -> {mp_id: set(genes)}
    all_gene_counts = {}  # gene -> count of (list) appearances across cNMF U GeneNMF

    for cancer in CANCERS:
        cdf = pd.read_csv(os.path.join(BASE, cancer, f"{cancer}_cnmf_top50_genes.csv"))
        cnmf[cancer] = {}
        for prog, sub in cdf.groupby("cnmf_program"):
            genes = set(sub.sort_values("rank")["gene"].head(TOP_N))
            cnmf[cancer][prog] = genes
            for g in genes:
                all_gene_counts[g] = all_gene_counts.get(g, 0) + 1

        mdf = pd.read_csv(os.path.join(BASE, cancer, f"{cancer}_genenmf_MP_top50_genes.csv"))
        genenmf[cancer] = {}
        for mp, sub in mdf.groupby("meta_program"):
            genes = set(sub.sort_values("rank")["gene"].head(TOP_N))
            genenmf[cancer][mp] = genes
            for g in genes:
                all_gene_counts[g] = all_gene_counts.get(g, 0) + 1

    return cnmf, genenmf, all_gene_counts


# ------------------------------------------------- Method A: permutation null
def weighted_sample_without_replacement(genes_arr, weights, k, rng):
    """Efraimidis-Spirakis weighted reservoir sampling w/o replacement."""
    u = rng.random(len(genes_arr))
    keys = u ** (1.0 / weights)
    top_idx = np.argpartition(-keys, k)[:k]
    return genes_arr[top_idx]


def permutation_null(all_gene_counts, n_perm, rng):
    genes_arr = np.array(list(all_gene_counts.keys()))
    weights = np.array([all_gene_counts[g] for g in genes_arr], dtype=float)
    jvals = np.empty(n_perm)
    for i in range(n_perm):
        listA = set(weighted_sample_without_replacement(genes_arr, weights, TOP_N, rng))
        listB = set(weighted_sample_without_replacement(genes_arr, weights, TOP_N, rng))
        jvals[i] = jaccard(listA, listB)
    return jvals


# --------------------------------------------- Method B: cross-cancer null
def cross_cancer_null(cnmf, genenmf):
    jvals = []
    for cancer_c in CANCERS:
        for prog, cgenes in cnmf[cancer_c].items():
            for cancer_m in CANCERS:
                if cancer_m == cancer_c:
                    continue  # exclude same-cancer-type pairs (the real test pairs)
                for mp, mgenes in genenmf[cancer_m].items():
                    jvals.append(jaccard(cgenes, mgenes))
    return np.array(jvals)


# --------------------------------------------------- real observed matrix
def real_best_jaccard(cnmf, genenmf):
    """Per (cancer, MP): best Jaccard among same-cancer cNMF programs.
    Matches the existing <Cancer>_highconfidence_programs.csv convention."""
    rows = []
    for cancer in CANCERS:
        for mp, mgenes in genenmf[cancer].items():
            best_j = max(jaccard(mgenes, cgenes) for cgenes in cnmf[cancer].values())
            rows.append({"cancer": cancer, "meta_program": mp, "best_jaccard": best_j})
    return pd.DataFrame(rows)


def main():
    print("Loading top-50 gene lists (11 cancers)...")
    cnmf, genenmf, all_gene_counts = load_lists()
    n_cnmf = sum(len(v) for v in cnmf.values())
    n_mp = sum(len(v) for v in genenmf.values())
    n_unique_genes = len(all_gene_counts)
    print(f"  {n_cnmf} cNMF programs, {n_mp} GeneNMF MPs, {n_unique_genes} unique genes in pool "
          f"(pool size = {sum(all_gene_counts.values())} gene-list memberships)")

    print("Computing real observed best-Jaccard per GeneNMF MP (sanity check vs prior pipeline)...")
    real_df = real_best_jaccard(cnmf, genenmf)
    real_df.to_csv(os.path.join(OUTDIR, "real_best_jaccard_recomputed.csv"), index=False)
    check_44 = int((real_df["best_jaccard"] > 0.1).sum())
    print(f"  recomputed HC count at Jaccard>0.1: {check_44} (prior pipeline reported 44)")

    print(f"\nMethod A: gene-label permutation null ({N_PERM} draws)...")
    null_A = permutation_null(all_gene_counts, N_PERM, RNG)
    print(f"  mean={null_A.mean():.4f}, median={np.median(null_A):.4f}, max={null_A.max():.4f}")

    print("\nMethod B: cross-cancer-type shuffled-pairing null (real gene lists, mismatched cancer types)...")
    null_B = cross_cancer_null(cnmf, genenmf)
    print(f"  n_pairs={len(null_B)}, mean={null_B.mean():.4f}, median={np.median(null_B):.4f}, max={null_B.max():.4f}")

    pct_levels = [95, 99, 99.9]
    summary_rows = []
    for name, arr in [("A_permutation", null_A), ("B_cross_cancer", null_B)]:
        for p in pct_levels:
            thr = float(np.percentile(arr, p))
            n_hc = int((real_df["best_jaccard"] > thr).sum())
            per_cancer = real_df[real_df["best_jaccard"] > thr].groupby("cancer").size().to_dict()
            summary_rows.append({
                "null_method": name,
                "percentile": p,
                "threshold": round(thr, 4),
                "n_HC": n_hc,
                "per_cancer_breakdown": "; ".join(f"{c}:{per_cancer.get(c,0)}" for c in CANCERS),
            })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(OUTDIR, "null_derived_threshold_summary.csv"), index=False)

    null_A_stats = pd.Series(null_A).describe(percentiles=[.5, .9, .95, .99, .999])
    null_B_stats = pd.Series(null_B).describe(percentiles=[.5, .9, .95, .99, .999])
    null_A_stats.to_csv(os.path.join(OUTDIR, "null_A_permutation_stats.csv"))
    null_B_stats.to_csv(os.path.join(OUTDIR, "null_B_crosscancer_stats.csv"))

    print("\n=== Null A (permutation) percentiles ===")
    print(null_A_stats)
    print("\n=== Null B (cross-cancer) percentiles ===")
    print(null_B_stats)
    print("\n=== Threshold x HC-count summary ===")
    print(summary.to_string(index=False))

    # where does 0.1 fall in each null?
    for name, arr in [("A_permutation", null_A), ("B_cross_cancer", null_B)]:
        pct_of_01 = float((arr < 0.1).mean() * 100)
        print(f"\n0.1 sits at the {pct_of_01:.2f}th percentile of null {name} "
              f"(i.e. {100-pct_of_01:.2f}% of null draws exceed 0.1)")

    print(f"\nWrote: {OUTDIR}/null_derived_threshold_summary.csv")
    print(f"Wrote: {OUTDIR}/null_A_permutation_stats.csv")
    print(f"Wrote: {OUTDIR}/null_B_crosscancer_stats.csv")
    print(f"Wrote: {OUTDIR}/real_best_jaccard_recomputed.csv")

    # save raw null arrays for plotting
    np.save(os.path.join(OUTDIR, "null_A_permutation_values.npy"), null_A)
    np.save(os.path.join(OUTDIR, "null_B_crosscancer_values.npy"), null_B)


if __name__ == "__main__":
    main()
