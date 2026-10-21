#!/usr/bin/env python
# build_hc85_ids_summary.py
#
# Recomputes the median InterDependence Score (IDS) for the EXPANDED
# high-confidence (HC) cNMF-GeneNMF program pair set, epithelial compartment.
#
# Background: the thesis originally validated HC pairs at an arbitrary
# Jaccard > 0.1 threshold (44 pairs, median IDS 0.10-0.21 across cancers).
# Section 5.2.7 replaced that with a statistically-derived threshold,
# Jaccard > 0.0538 (99th percentile of a null distribution of chance-level
# cross-cancer-type Jaccard), which expands the HC set to 85 pairs. IDS for
# those 85 pairs had not yet been recomputed -- this script does that using
# IDS values already computed by run_ids.py (IDS-max, Radhakrishnan et al.
# PNAS 2025; per-cell program activity; unchanged methodology) for every
# cNMF-GeneNMF program pair per cancer type, so no new IDS computation is
# needed beyond adding Breast (which run_ids.py's default CANCERS list had
# omitted).
#
# HC pair identity = each GeneNMF meta-program's single best cNMF match, as
# already determined by the R GeneNMF pipeline (best_cnmf_match /
# best_jaccard columns in <cancer>_highconfidence_programs.csv), filtered to
# best_jaccard > 0.0538. This is the same pairing logic used for the
# original 44-pair HC designation (run_genenmf_per_cancer.R), just at the
# new threshold -- not a re-derivation from raw gene lists, since the
# per-MP best-match-by-Jaccard-to-cNMF computation already lives there and
# top-50 gene lists / best_jaccard values are unchanged.
#
# Inputs:
#   <ATLAS>/genenmf_output/<Cancer>/<Cancer>_highconfidence_programs.csv
#     (meta_program, best_cnmf_match, best_jaccard)
#   <ATLAS>/ids_output/<Cancer>/<Cancer>_ids_pairs.csv
#     (long-format IDS for every cNMF-GeneNMF / cNMF-cNMF / GeneNMF-GeneNMF
#     pair in that cancer; produced by run_ids.py)
#
# Output:
#   ids_analysis/hc85_ids_pairs.csv       -- one row per HC pair
#   ids_analysis/hc85_ids_summary.txt     -- overall + per-cancer stats,
#                                             old-44 vs new-85 comparison

import os
import numpy as np
import pandas as pd

ATLAS = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects"
GENENMF_BASE = os.path.join(ATLAS, "genenmf_output")
IDS_BASE = os.path.join(ATLAS, "ids_output")
OUT_DIR = "/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/ids_analysis"

CANCERS = ["Breast", "Kidney", "Lung", "Colorectal", "Skin", "Brain",
           "Liver_Biliary", "Ovarian", "Pancreas", "Prostate", "Sarcoma"]

JACCARD_NEW = 0.0538   # statistically-derived threshold (Sec. 5.2.7), strict >
JACCARD_OLD = 0.1      # original arbitrary threshold, strict > (matches
                        # run_genenmf_per_cancer.R / make_fig44 JACCARD_THR_EPI)


def load_hc_pairs(cancer, thr):
    path = os.path.join(GENENMF_BASE, cancer, f"{cancer}_highconfidence_programs.csv")
    df = pd.read_csv(path)
    df = df[df["best_jaccard"] > thr].copy()
    df["cancer"] = cancer
    return df[["cancer", "meta_program", "best_cnmf_match", "best_jaccard"]]


def load_ids_lookup(cancer):
    path = os.path.join(IDS_BASE, cancer, f"{cancer}_ids_pairs.csv")
    df = pd.read_csv(path)
    cross = df[df["pair_type"] == "GeneNMF_cNMF"]
    lookup = {}
    for _, r in cross.iterrows():
        a, b = r["program_a"], r["program_b"]
        lookup[frozenset((a, b))] = r["ids"]
    return lookup


def main():
    hc_rows = []
    missing = []
    for cancer in CANCERS:
        hc = load_hc_pairs(cancer, JACCARD_NEW)
        if hc.empty:
            continue
        ids_lookup = load_ids_lookup(cancer)
        for _, r in hc.iterrows():
            mp, cnmf_prog, jacc = r["meta_program"], r["best_cnmf_match"], r["best_jaccard"]
            key = frozenset((mp, cnmf_prog))
            ids_val = ids_lookup.get(key)
            if ids_val is None:
                missing.append((cancer, cnmf_prog, mp))
                continue
            hc_rows.append({
                "cancer_type": cancer,
                "cNMF_program_id": cnmf_prog,
                "GeneNMF_program_id": mp,
                "jaccard": jacc,
                "IDS": ids_val,
            })

    out = pd.DataFrame(hc_rows).sort_values(["cancer_type", "cNMF_program_id"])
    out.to_csv(os.path.join(OUT_DIR, "hc85_ids_pairs.csv"), index=False)

    if missing:
        print(f"WARNING: {len(missing)} HC pairs had no matching IDS entry:")
        for m in missing:
            print(f"  {m}")

    n = len(out)
    overall_median = out["IDS"].median()
    overall_min = out["IDS"].min()
    overall_max = out["IDS"].max()

    per_cancer = out.groupby("cancer_type")["IDS"].agg(["median", "min", "max", "count"])
    per_cancer = per_cancer.sort_values("median", ascending=False)

    # old-44 set for comparison / "new 41" flag
    old_rows = []
    for cancer in CANCERS:
        hc_old = load_hc_pairs(cancer, JACCARD_OLD)
        if hc_old.empty:
            continue
        ids_lookup = load_ids_lookup(cancer)
        for _, r in hc_old.iterrows():
            mp, cnmf_prog = r["meta_program"], r["best_cnmf_match"]
            key = frozenset((mp, cnmf_prog))
            ids_val = ids_lookup.get(key)
            if ids_val is not None:
                old_rows.append({"cancer_type": cancer, "cNMF_program_id": cnmf_prog,
                                  "GeneNMF_program_id": mp, "IDS": ids_val})
    old_df = pd.DataFrame(old_rows)

    out["pair_key"] = out["cancer_type"] + "|" + out["cNMF_program_id"] + "|" + out["GeneNMF_program_id"]
    old_df["pair_key"] = old_df["cancer_type"] + "|" + old_df["cNMF_program_id"] + "|" + old_df["GeneNMF_program_id"]
    new_only = out[~out["pair_key"].isin(old_df["pair_key"])]

    lines = []
    lines.append(f"HC pairs (Jaccard > {JACCARD_NEW}): n = {n}")
    lines.append(f"Overall median IDS = {overall_median:.4f}")
    lines.append(f"Overall min IDS    = {overall_min:.4f}")
    lines.append(f"Overall max IDS    = {overall_max:.4f}")
    lines.append("")
    lines.append("Median IDS by cancer type (n pairs):")
    for cancer, row in per_cancer.iterrows():
        lines.append(f"  {cancer:<15s} median={row['median']:.4f}  "
                      f"min={row['min']:.4f}  max={row['max']:.4f}  n={int(row['count'])}")
    lines.append("")
    lines.append(f"Old HC set (Jaccard > {JACCARD_OLD}): n = {len(old_df)}")
    lines.append(f"Old-set median IDS = {old_df['IDS'].median():.4f}")
    lines.append("")
    lines.append(f"Newly-added pairs (in 85-set, not in old {len(old_df)}-set): n = {len(new_only)}")
    if len(new_only):
        lines.append(f"  New-only median IDS = {new_only['IDS'].median():.4f}")
        lines.append(f"  New-only min/max    = {new_only['IDS'].min():.4f} / {new_only['IDS'].max():.4f}")
        # Mann-Whitney-style rank comparison without scipy dependency: just report medians + a simple check
        old_median = old_df["IDS"].median()
        new_median = new_only["IDS"].median()
        direction = "LOWER than" if new_median < old_median else ("HIGHER than" if new_median > old_median else "comparable to")
        lines.append(f"  New-only median IDS is {direction} the old-44 median "
                      f"({new_median:.4f} vs {old_median:.4f}).")

    summary_txt = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "hc85_ids_summary.txt"), "w") as f:
        f.write(summary_txt + "\n")

    print(summary_txt)
    print(f"\nWrote {os.path.join(OUT_DIR, 'hc85_ids_pairs.csv')}")
    print(f"Wrote {os.path.join(OUT_DIR, 'hc85_ids_summary.txt')}")


if __name__ == "__main__":
    main()
