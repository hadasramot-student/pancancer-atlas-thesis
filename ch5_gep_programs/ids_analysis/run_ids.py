#!/usr/bin/env python
# run_ids.py
#
# Ch.5 (GEP programs) — IDS (InterDependence Score) implementation, the
# custom cross-method validation metric for GEP programs, used alongside
# Jaccard overlap for the epithelial cNMF-vs-GeneNMF comparison.
#
# ---------------------------------------------------------------------------
# ENV / RUNTIME NOTES (critical — do not lose):
#   * The `ids` package (github.com/aradha/interdependence_scores) MUST be
#     installed NON-EDITABLE (`pip install .`). An editable install points into
#     /tmp, which gets wiped, breaking `import ids`.
#   * This analysis directory must NOT be named "ids" — a dir named `ids` on the
#     worktree root shadows the installed `ids` package on sys.path (it becomes
#     an empty PEP-420 namespace package). Hence `ids_analysis/`.
#   * The interactive/login node has only ~23 GB RAM. IDS itself is cheap, but
#     loading + scoring the large h5ads is memory-heavy (Colorectal peaked at
#     ~41 GB). Run via SLURM (see submit_ids_all_cancers.sh, --mem=96G), not on
#     the interactive node, for the big cancers.
# ---------------------------------------------------------------------------
#
# Compute InterDependence Scores (IDS, Radhakrishnan et al. PNAS 2025,
# doi:10.1073/pnas.2509860122) between transcriptional programs identified by
# cNMF and GeneNMF, per cancer type.
#
# IDS measures linear + nonlinear statistical dependence between two variables
# given paired samples. Here:
#   - samples   = single cells (within one cancer type)
#   - variables = program activities per cell
#       * cNMF programs   -> consensus usage matrix (cells x K)
#       * GeneNMF MPs     -> sc.tl.score_genes of each MP's top genes on the
#                            same cells (log-normalized expression)
#
# Both methods are scored on the *same* cells, so the resulting square IDS
# matrix contains cNMF-cNMF, cNMF-GeneNMF and GeneNMF-GeneNMF dependence in one
# object -> directly usable for the cross-method comparison.
#
# Per the paper, each variable is min-max scaled to [0, 8] before IDS, and
# IDS-max (num_terms=6, bandwidth_term=1/2) is used.
#
# Usage:
#   python run_ids.py --cancer Lung          # one cancer
#   python run_ids.py --cancer all           # all cancers
#   python run_ids.py --cancer all --pvals   # also permutation p-values

import argparse
import os
import math
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from ids.numpy_dependence import compute_IDS_numpy

_ATLAS_DATA_DIR = os.environ.get("ATLAS_DATA_DIR", "data/temp_objects")
CNMF_BASE    = os.path.join(_ATLAS_DATA_DIR, "cnmf_output")
GENENMF_BASE = os.path.join(_ATLAS_DATA_DIR, "genenmf_output")
OUT_BASE     = os.path.join(_ATLAS_DATA_DIR, "ids_output")

CANCERS = ["Lung", "Colorectal", "Kidney", "Liver_Biliary", "Ovarian",
           "Pancreas", "Prostate", "Sarcoma", "Skin", "Brain"]

K = 20                 # cNMF components
DT = "dt_0_1"          # cNMF density threshold tag
N_TOP_MP_GENES = 50    # top genes per GeneNMF meta-program to score
SCORE_SEED = 0
SCALE_MAX = 8.0        # paper: min-max scale each variable to [0, 8]


def load_cnmf_usages(cancer):
    """cNMF consensus usage matrix -> DataFrame (cells x cNMF_1..cNMF_K)."""
    path = os.path.join(CNMF_BASE, cancer, f"cnmf_{cancer}",
                        f"cnmf_{cancer}.usages.k_{K}.{DT}.consensus.txt")
    u = pd.read_csv(path, sep="\t", index_col=0)
    u.columns = [f"cNMF_{c}" for c in u.columns]
    return u


def load_mp_genesets(cancer):
    """GeneNMF meta-program top genes -> dict {MP: [genes]}."""
    path = os.path.join(GENENMF_BASE, cancer, f"{cancer}_genenmf_MP_top50_genes.csv")
    df = pd.read_csv(path)
    df = df[df["rank"] <= N_TOP_MP_GENES]
    return {mp: sub.sort_values("rank")["gene"].tolist()
            for mp, sub in df.groupby("meta_program")}


def score_genenmf(cancer, cell_order):
    """Per-cell activity of each GeneNMF MP via signature scoring.

    Loads the cNMF h5ad (raw counts, same cells), log-normalizes, and runs
    sc.tl.score_genes per MP. Returns DataFrame (cells x MP...) aligned to
    cell_order.
    """
    h5ad = os.path.join(CNMF_BASE, cancer, f"cnmf_{cancer}.h5ad")
    adata = ad.read_h5ad(h5ad)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    genesets = load_mp_genesets(cancer)
    var_set = set(adata.var_names)
    scores = pd.DataFrame(index=adata.obs_names)
    for mp, genes in sorted(genesets.items(), key=lambda kv: kv[0]):
        present = [g for g in genes if g in var_set]
        if len(present) < 5:
            print(f"    [warn] {mp}: only {len(present)} genes present, skipping")
            continue
        sc.tl.score_genes(adata, present, score_name=f"_score_{mp}",
                          ctrl_size=50, random_state=SCORE_SEED, use_raw=False)
        scores[mp] = adata.obs[f"_score_{mp}"].values
    return scores.reindex(cell_order)


def minmax_scale(M, smax=SCALE_MAX):
    """Min-max scale each column to [0, smax]."""
    M = M.astype(np.float64)
    lo = M.min(axis=0, keepdims=True)
    hi = M.max(axis=0, keepdims=True)
    rng = hi - lo
    rng[rng == 0] = 1.0
    return (M - lo) / rng * smax


def long_format(ids_mat, labels, methods):
    """Square IDS matrix -> long table of unique program pairs (i<j)."""
    rows = []
    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            rows.append({
                "program_a": labels[i], "program_b": labels[j],
                "method_a": methods[i], "method_b": methods[j],
                "pair_type": "_".join(sorted([methods[i], methods[j]])),
                "ids": ids_mat[i, j],
            })
    return pd.DataFrame(rows).sort_values("ids", ascending=False)


def run_cancer(cancer, pvals=False):
    print(f"\n{'='*60}\n{cancer}\n{'='*60}")
    out_dir = os.path.join(OUT_BASE, cancer)
    os.makedirs(out_dir, exist_ok=True)

    print("  Loading cNMF usages...")
    cnmf = load_cnmf_usages(cancer)
    print(f"    cNMF: {cnmf.shape[0]} cells x {cnmf.shape[1]} programs")

    print("  Scoring GeneNMF meta-programs on the same cells...")
    gnmf = score_genenmf(cancer, cnmf.index)
    print(f"    GeneNMF: {gnmf.shape[1]} meta-programs scored")

    # Align to common cells (should already match; intersection is a safety net)
    common = cnmf.index.intersection(gnmf.index)
    common = common[gnmf.loc[common].notna().all(axis=1)]
    cnmf, gnmf = cnmf.loc[common], gnmf.loc[common]
    print(f"  Common cells used: {len(common)}")

    combined = pd.concat([cnmf, gnmf], axis=1)
    labels = list(combined.columns)
    methods = ["cNMF"] * cnmf.shape[1] + ["GeneNMF"] * gnmf.shape[1]

    M = minmax_scale(combined.values)
    print(f"  Computing IDS-max on {M.shape[0]} cells x {M.shape[1]} programs...")
    if pvals:
        ids_mat, p_mat = compute_IDS_numpy(M, p_norm="max", p_val=True, num_tests=100)
        pd.DataFrame(p_mat, index=labels, columns=labels).to_csv(
            os.path.join(out_dir, f"{cancer}_ids_pvals_matrix.csv"))
    else:
        ids_mat = compute_IDS_numpy(M, p_norm="max")

    # Save square matrix
    pd.DataFrame(ids_mat, index=labels, columns=labels).to_csv(
        os.path.join(out_dir, f"{cancer}_ids_matrix.csv"))

    # Save long-format pairs (with method annotation)
    lf = long_format(ids_mat, labels, methods)
    if pvals:
        pmap = {(labels[i], labels[j]): p_mat[i, j]
                for i in range(len(labels)) for j in range(len(labels))}
        lf["pval"] = [pmap[(a, b)] for a, b in zip(lf["program_a"], lf["program_b"])]
    lf.insert(0, "cancer", cancer)
    lf.to_csv(os.path.join(out_dir, f"{cancer}_ids_pairs.csv"), index=False)

    # Quick summary by pair type
    print("  Top cross-method (cNMF-GeneNMF) pairs:")
    cross = lf[lf["pair_type"] == "GeneNMF_cNMF"].head(5)
    for _, r in cross.iterrows():
        print(f"    {r['program_a']:>10} ~ {r['program_b']:<8} IDS={r['ids']:.3f}")
    print(f"  Saved -> {out_dir}")
    return lf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer", default="all")
    ap.add_argument("--pvals", action="store_true",
                    help="also compute permutation p-values (slower)")
    args = ap.parse_args()

    cancers = CANCERS if args.cancer == "all" else [args.cancer]
    all_lf = []
    for c in cancers:
        all_lf.append(run_cancer(c, pvals=args.pvals))

    if len(all_lf) > 1:
        os.makedirs(OUT_BASE, exist_ok=True)
        combined = pd.concat(all_lf, ignore_index=True)
        combined.to_csv(os.path.join(OUT_BASE, "all_cancers_ids_pairs.csv"), index=False)
        print(f"\nCombined pairs -> {os.path.join(OUT_BASE, 'all_cancers_ids_pairs.csv')}")


if __name__ == "__main__":
    main()
