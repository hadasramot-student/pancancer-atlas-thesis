#!/usr/bin/env python
"""
Gavish-style meta-program catalogue for Ch.5 cNMF programs (epithelial/tumor
compartment; K=20 per cancer type x 12 cancer types = 240 programs).

Adaptation of Gavish et al. 2023 (Nature) MP-building algorithm:
  1. pairwise Jaccard (top-50 gene overlap) across ALL 240 cNMF programs
     (all cancer-type pairs, incl. same-cancer-type pairs)
  2. iterative/greedy grouping: repeatedly take the program with the most
     Jaccard-similar partners (> threshold), group it + partners into one
     meta-program, remove from pool, repeat until no groupable programs remain
  3. recurrence filter: a group only counts as a genuine meta-program if it
     spans >= MIN_CANCER_TYPES distinct cancer types (Gavish's per-tumor
     recurrence filter, adapted to "per cancer type" since we have one
     NMF solution per cancer type, not per tumor)
  4. consensus signature per MP = genes ranked by frequency across member
     top-50 lists (ties broken by mean spectra-score rank)
  5. per-MP Hallmark hypergeometric enrichment on the consensus top-50 genes,
     with BH-FDR correction across the 50 Hallmark tests for that MP
  6. per-MP match to Gavish's own 41 "Cancer" meta-programs via Jaccard > 0.05
     between consensus top-50 and each Gavish MP's top-50 gene list

Run at two thresholds (0.1 primary, 0.05 for sensitivity) so threshold choice
can be compared.
"""

import json
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

CNMF_BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
GAVISH_CSV = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/gavish_2023/gavish_MP_gene_lists_long.csv"
HALLMARK_JSON = "/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/merge/spectra/priors/hallmark.json"

K_PROGRAMS = 20
DT = "dt_0_1"
N_TOP_GENES = 50
GAVISH_CUTOFF = 0.05
MIN_HALLMARK_OVERLAP = 3
HALLMARK_PVAL_CUTOFF = 0.05
MIN_CANCER_TYPES = 2          # recurrence filter (Gavish: min independent tumors -> here, min cancer types)
PANCANCER_FRACTION = 0.5      # "majority of the 12 cancer types" -> >=6 -> we'll use >=7 (strict majority)

CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney",
           "Liver_Biliary", "Lung", "Ovarian", "Pancreas", "Prostate",
           "Sarcoma", "Skin"]

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def load_gene_spectra_score(cancer):
    path = os.path.join(CNMF_BASE, cancer, f"cnmf_{cancer}",
                         f"cnmf_{cancer}.gene_spectra_score.k_{K_PROGRAMS}.{DT}.txt")
    return pd.read_csv(path, sep="\t", index_col=0)


def load_all_programs():
    """Returns dict: prog_key -> {'cancer', 'idx', 'top_genes' (list, ranked desc), 'top_set' (set)}"""
    programs = {}
    n_background = None
    for cancer in CANCERS:
        df = load_gene_spectra_score(cancer)
        if n_background is None:
            n_background = df.shape[1]
        else:
            assert df.shape[1] == n_background, f"{cancer} background gene count mismatch"
        for prog_idx in df.index:
            ranked = df.loc[prog_idx].sort_values(ascending=False)
            top_genes = list(ranked.head(N_TOP_GENES).index)
            key = f"{cancer}_cNMF_{prog_idx}"
            programs[key] = {
                "cancer": cancer,
                "idx": prog_idx,
                "top_genes": top_genes,
                "top_set": set(top_genes),
            }
    return programs, n_background


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
    """Gavish-style greedy iterative clustering.
    Returns: list of groups (each a list of indices into `keys`), plus list of
    singleton indices (no partner above threshold at the time they were
    considered)."""
    n = len(keys)
    pool = set(range(n))
    groups = []
    while pool:
        pool_list = list(pool)
        # count partners above threshold, restricted to current pool
        best_i, best_count, best_partners = None, -1, None
        for i in pool_list:
            partners = [j for j in pool_list if j != i and J[i, j] > threshold]
            if len(partners) > best_count:
                best_i, best_count, best_partners = i, len(partners), partners
        if best_count <= 0:
            break  # no more groupable programs; remainder are singletons
        group = [best_i] + best_partners
        groups.append(group)
        pool -= set(group)
    singletons = sorted(pool)
    return groups, singletons


def consensus_genes(programs, keys, group_idx, key_list):
    """Rank genes by (frequency across member top-50 lists desc, mean rank asc)."""
    freq = {}
    rank_sum = {}
    for idx in group_idx:
        key = key_list[idx]
        genes = programs[key]["top_genes"]
        for rank, g in enumerate(genes):
            freq[g] = freq.get(g, 0) + 1
            rank_sum[g] = rank_sum.get(g, 0) + rank
    genes_sorted = sorted(freq.keys(), key=lambda g: (-freq[g], rank_sum[g] / freq[g]))
    return genes_sorted, freq


def load_hallmark():
    with open(HALLMARK_JSON) as f:
        h = json.load(f)
    return {k.replace("HALLMARK_", ""): set(v) for k, v in h.items()}


def bh_correct(pvals):
    """Benjamini-Hochberg FDR correction. pvals: list -> returns list of qvals in same order."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.array(pvals)[order]
    q = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity (cumulative min from the largest rank down)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out.tolist()


def hallmark_enrichment(top_genes_set, hallmark, n_background):
    """Runs hypergeometric test vs all 50 hallmark sets, BH-corrects across
    them, returns best hit by min raw p-value."""
    names, pvals, overlaps = [], [], []
    for hname, hset in hallmark.items():
        overlap = len(top_genes_set & hset)
        if overlap < MIN_HALLMARK_OVERLAP:
            pval = 1.0
        else:
            pval = hypergeom.sf(overlap - 1, n_background, len(hset), len(top_genes_set))
        names.append(hname)
        pvals.append(pval)
        overlaps.append(overlap)
    qvals = bh_correct(pvals)
    best_i = int(np.argmin(pvals))
    best_name = names[best_i] if pvals[best_i] < HALLMARK_PVAL_CUTOFF else "Unannotated"
    return {
        "hallmark": best_name,
        "pval": pvals[best_i],
        "qval": qvals[best_i],
        "overlap": overlaps[best_i],
    }


def load_gavish_cancer_mps():
    df = pd.read_csv(GAVISH_CSV)
    df = df[df["cell_type"] == "Cancer"]
    return {mp: set(sub["gene"]) for mp, sub in df.groupby("meta_program")}


def gavish_match(top_genes_set, gavish_mps):
    best_mp, best_j = None, 0.0
    for mp, mp_genes in gavish_mps.items():
        j = jaccard(top_genes_set, mp_genes)
        if j > best_j:
            best_mp, best_j = mp, j
    matched = best_j > GAVISH_CUTOFF
    return (best_mp if matched else "Unmatched"), best_j, matched


def build_catalogue(programs, key_list, J, threshold, hallmark, n_background, gavish_mps, mp_prefix):
    groups, singleton_idx = greedy_cluster(key_list, J, threshold)

    valid_mps = []
    failed_recurrence = []  # groups that clustered but span < MIN_CANCER_TYPES

    for g in groups:
        cancers_in_group = sorted(set(programs[key_list[i]]["cancer"] for i in g))
        if len(cancers_in_group) >= MIN_CANCER_TYPES:
            valid_mps.append((g, cancers_in_group))
        else:
            failed_recurrence.append((g, cancers_in_group))

    rows = []
    for mp_num, (g, cancers_in_group) in enumerate(
            sorted(valid_mps, key=lambda x: -len(x[0])), start=1):
        mp_id = f"{mp_prefix}{mp_num}"
        members = [key_list[i] for i in g]
        genes_ranked, freq = consensus_genes(programs, key_list, g, key_list)
        consensus_top50 = set(genes_ranked[:N_TOP_GENES])
        hm = hallmark_enrichment(consensus_top50, hallmark, n_background)
        gm_mp, gm_j, gm_matched = gavish_match(consensus_top50, gavish_mps)
        n_cancers = len(cancers_in_group)
        classification = "pan-cancer" if n_cancers >= 7 else "cancer-type-specific"
        rows.append({
            "MP_ID": mp_id,
            "n_member_programs": len(members),
            "member_programs": ";".join(members),
            "n_distinct_cancer_types": n_cancers,
            "cancer_types_spanned": ";".join(cancers_in_group),
            "classification": classification,
            "consensus_top_genes": ";".join(genes_ranked[:15]),
            "best_hallmark": hm["hallmark"],
            "hallmark_pval": hm["pval"],
            "hallmark_qval_BH": hm["qval"],
            "hallmark_overlap": hm["overlap"],
            "gavish_match": gm_mp,
            "gavish_jaccard": round(gm_j, 4),
        })

    df = pd.DataFrame(rows)

    # non-recurring clusters (grouped by Jaccard but only 1 cancer type) + true singletons
    non_recurring_rows = []
    for g, cancers_in_group in failed_recurrence:
        members = [key_list[i] for i in g]
        non_recurring_rows.append({
            "programs": ";".join(members),
            "n_programs": len(members),
            "cancer_type": cancers_in_group[0],
            "reason": "clustered but within a single cancer type (fails recurrence filter)",
        })
    for i in singleton_idx:
        non_recurring_rows.append({
            "programs": key_list[i],
            "n_programs": 1,
            "cancer_type": programs[key_list[i]]["cancer"],
            "reason": "singleton (no Jaccard partner above threshold)",
        })
    non_recurring_df = pd.DataFrame(non_recurring_rows)

    return df, non_recurring_df, len(groups), len(singleton_idx)


def main():
    print("Loading 240 cNMF programs (top-50 genes each)...")
    programs, n_background = load_all_programs()
    key_list = list(programs.keys())
    print(f"  {len(key_list)} programs loaded; background = {n_background} genes")

    print("Computing pairwise Jaccard matrix (240x240)...")
    key_list, J = build_jaccard_matrix(programs)

    hallmark = load_hallmark()
    gavish_mps = load_gavish_cancer_mps()

    for threshold, prefix, tag in [(0.1, "cMP", "0.1"), (0.05, "cMP", "0.05")]:
        print(f"\n=== Threshold Jaccard > {threshold} ===")
        df, non_recurring_df, n_raw_groups, n_true_singletons = build_catalogue(
            programs, key_list, J, threshold, hallmark, n_background, gavish_mps, prefix)
        out_csv = os.path.join(OUTDIR, f"meta_program_catalogue_j{tag.replace('.', '')}.csv")
        df.to_csv(out_csv, index=False)
        out_nonrec = os.path.join(OUTDIR, f"non_recurring_programs_j{tag.replace('.', '')}.csv")
        non_recurring_df.to_csv(out_nonrec, index=False)
        print(f"  raw Jaccard-groups (before recurrence filter): {n_raw_groups}")
        print(f"  valid meta-programs (>= {MIN_CANCER_TYPES} cancer types): {len(df)}")
        print(f"  single-cancer-type clusters (fail recurrence filter): "
              f"{(non_recurring_df['reason'].str.contains('single cancer')).sum() if len(non_recurring_df) else 0}")
        print(f"  true singleton programs: {n_true_singletons}")
        if len(df):
            print(f"  MP size distribution: min={df.n_member_programs.min()}, "
                  f"median={df.n_member_programs.median()}, max={df.n_member_programs.max()}")
            print(f"  pan-cancer MPs (>=7 cancer types): {(df.classification=='pan-cancer').sum()}")
        print(f"  wrote: {out_csv}")
        print(f"  wrote: {out_nonrec}")


if __name__ == "__main__":
    main()
