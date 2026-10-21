#!/usr/bin/env python
# 01_hc_pairs_and_hallmark_bh.py
#
# Comment 12 (Figure 6.1, A.29, A.30, A.38) -- data step. No cNMF/GeneNMF/IDS
# re-run: reads saved outputs only.
#   1. HC pair tables at Jaccard > 0.0538 and > 0.1 (cNMF vs GeneNMF, top-50)
#      -> hc_pairs_0.0538.csv, hc_pairs_0.1.csv, hc_pairs_extra_0.0538_not_0.1.csv
#   2. Kidney redundancy check (85 -> 74 reconstruction)
#      -> kidney_redundancy_families.csv
#   3. Hallmark hypergeometric tests, all 50 sets per cNMF program, BH within
#      program -> hallmark_all_tests_bh.csv, hallmark_top_per_program_bh.csv
#
# Hallmark test reproduces run_cnmf_hallmark.py: top-50 genes by z-score from
# cnmf_<cancer>.gene_spectra_score.k_20.dt_0_1.txt, background = all scored
# genes (28,231), one-sided hypergeometric P(X >= overlap), Hallmark sets from
# merge/spectra/priors/hallmark.json.

import json
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

BASE = "/groups/ofircohen-group/users/hadasra"
CNMF = f"{BASE}/temp_objects/cnmf_output"
GNMF = f"{BASE}/temp_objects/genenmf_output"
HALLMARK_JSON = f"{BASE}/create_sc_object_workflow/merge/spectra/priors/hallmark.json"
OLD_TOP = f"{BASE}/temp_objects/ids_output/figures/cnmf_top_hallmark_new.csv"
OUT = f"{BASE}/thesis_figures/comment12_2026-10-21"

CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney", "Liver_Biliary",
           "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"]
THR_LOW, THR_HIGH = 0.0538, 0.1
N_TOP, MIN_OVERLAP, Q_CUT = 50, 3, 0.05
KIDNEY_MP_MERGE = 0.25  # MP-MP top-50 Jaccard used to collapse redundant Kidney MPs


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.minimum(r, 1.0)
    return q


def top50_spectra(cancer):
    df = pd.read_csv(f"{CNMF}/{cancer}/cnmf_{cancer}/cnmf_{cancer}.gene_spectra_score.k_20.dt_0_1.txt",
                     sep="\t", index_col=0)
    return df.shape[1], {f"cNMF_{i}": list(df.loc[i].sort_values(ascending=False).head(N_TOP).index)
                         for i in df.index}


# ---------------------------------------------------------------- 1. HC pairs
pairs = []
for c in CANCERS:
    f = f"{GNMF}/{c}/{c}_highconfidence_programs.csv"
    if not os.path.exists(f):
        continue
    hc = pd.read_csv(f)
    cg = pd.read_csv(f"{GNMF}/{c}/{c}_cnmf_top50_genes.csv")
    mg = pd.read_csv(f"{GNMF}/{c}/{c}_genenmf_MP_top50_genes.csv")
    cset = {k: set(v.gene) for k, v in cg.groupby("cnmf_program")}
    mset = {k: set(v.gene) for k, v in mg.groupby("meta_program")}
    for _, r in hc.iterrows():
        a, b = cset[r.best_cnmf_match], mset[r.meta_program]
        pairs.append(dict(cancer=c, cnmf_program=r.best_cnmf_match, genenmf_mp=r.meta_program,
                          jaccard=r.best_jaccard, n_shared_genes=len(a & b),
                          jaccard_recomputed=round(len(a & b) / len(a | b), 4),
                          n_mp_genes=r.n_mp_genes))
pairs = pd.DataFrame(pairs)
pairs["mp_num"] = pairs.genenmf_mp.str[2:].astype(int)
pairs = pairs.sort_values(["cancer", "mp_num"]).drop(columns="mp_num")
n_mismatch = (abs(pairs.jaccard - pairs.jaccard_recomputed) > 1e-3).sum()
cols = ["cancer", "cnmf_program", "genenmf_mp", "jaccard", "n_shared_genes"]
p_low = pairs[pairs.jaccard > THR_LOW]
p_high = pairs[pairs.jaccard > THR_HIGH]
p_low[cols].to_csv(f"{OUT}/hc_pairs_0.0538.csv", index=False)
p_high[cols].to_csv(f"{OUT}/hc_pairs_0.1.csv", index=False)
extra = p_low[p_low.jaccard <= THR_HIGH]
extra[cols].to_csv(f"{OUT}/hc_pairs_extra_0.0538_not_0.1.csv", index=False)

counts = pd.DataFrame({
    "pairs_gt_0.0538": p_low.groupby("cancer").size(),
    "pairs_gt_0.1": p_high.groupby("cancer").size(),
    "HC_cNMF_programs_gt_0.0538": p_low.groupby("cancer").cnmf_program.nunique(),
    "HC_cNMF_programs_gt_0.1": p_high.groupby("cancer").cnmf_program.nunique(),
}).reindex([c for c in CANCERS if c != "Head_and_Neck"]).fillna(0).astype(int)
counts.loc["Total"] = counts.sum()
counts.to_csv(f"{OUT}/hc_counts_by_cancer.csv")

# ------------------------------------------------------- 2. Kidney redundancy
mg = pd.read_csv(f"{GNMF}/Kidney/Kidney_genenmf_MP_top50_genes.csv")
mset = {k: set(v.gene) for k, v in mg.groupby("meta_program")}
met = pd.read_csv(f"{GNMF}/Kidney/Kidney_genenmf_MP_metrics.csv", index_col=0)
kid = p_low[p_low.cancer == "Kidney"].copy()
ms = list(kid.genenmf_mp)
parent = {m: m for m in ms}


def find(x):
    while parent[x] != x:
        x = parent[x]
    return x


for i, a in enumerate(ms):
    for b in ms[i + 1:]:
        if len(mset[a] & mset[b]) / len(mset[a] | mset[b]) >= KIDNEY_MP_MERGE:
            parent[find(a)] = find(b)
kid["family_root"] = [find(m) for m in ms]
kid["sampleCoverage"] = met.loc[kid.genenmf_mp, "sampleCoverage"].values
# representative = highest cross-method Jaccard, ties -> higher sample coverage
kid = kid.sort_values(["family_root", "jaccard", "sampleCoverage"], ascending=[True, False, False])
kid["kept"] = ~kid.duplicated("family_root")
fam_id = {r: f"F{i+1}" for i, r in enumerate(kid.family_root.unique())}
kid["family"] = kid.family_root.map(fam_id)
kid[["cancer", "genenmf_mp", "cnmf_program", "jaccard", "sampleCoverage", "family", "kept"]] \
    .to_csv(f"{OUT}/kidney_redundancy_families.csv", index=False)
n_kid_units = kid.kept.sum()
alt_units = p_high[p_high.cancer == "Kidney"].cnmf_program.nunique()

# ------------------------------------------------------- 3. Hallmark + BH
hallmark = {k.replace("HALLMARK_", ""): set(v) for k, v in json.load(open(HALLMARK_JSON)).items()}
hc_low = p_low.groupby("cancer").cnmf_program.apply(set).to_dict()
hc_high = p_high.groupby("cancer").cnmf_program.apply(set).to_dict()
rows = []
for c in CANCERS:
    N, tops = top50_spectra(c)
    for prog, genes in tops.items():
        g = set(genes)
        recs = []
        for h, hs in hallmark.items():
            ov = len(g & hs)
            p = hypergeom.sf(ov - 1, N, len(hs), N_TOP) if ov > 0 else 1.0
            recs.append(dict(cancer=c, program=prog, hallmark=h, overlap=ov, set_size=len(hs),
                             background=N, pval=p))
        q = bh([r["pval"] for r in recs])
        for r, qq in zip(recs, q):
            r["qval"] = qq
            r["enriched"] = (qq < Q_CUT) and (r["overlap"] >= MIN_OVERLAP)
            r["HC_0.0538"] = prog in hc_low.get(c, set())
            r["HC_0.1"] = prog in hc_high.get(c, set())
        rows += recs
allt = pd.DataFrame(rows)
allt.to_csv(f"{OUT}/hallmark_all_tests_bh.csv", index=False)

# top pathway per program = smallest q (ties -> larger overlap); Unannotated if not enriched
best = (allt.sort_values(["cancer", "program", "qval", "overlap"], ascending=[True, True, True, False])
        .groupby(["cancer", "program"], as_index=False).first())
best.loc[~best.enriched, "hallmark"] = "Unannotated"
best = best.rename(columns={"hallmark": "top_hallmark_bh"})

old = pd.read_csv(OLD_TOP).rename(columns={"top_hallmark": "top_hallmark_old", "pval": "pval_old",
                                           "overlap": "overlap_old", "high_confidence": "HC_old"})
best = best.merge(old, on=["cancer", "program"], how="left")
best.to_csv(f"{OUT}/hallmark_top_per_program_bh.csv", index=False)

# ------------------------------------------------------- summary printout
print("Jaccard recompute mismatches:", n_mismatch)
print(counts.to_string())
print(f"Extra pairs 0.0538-not-0.1: {len(extra)}")
print(f"Kidney units after MP-MP merge (J>={KIDNEY_MP_MERGE}): {n_kid_units} "
      f"-> total {len(p_low) - (len(kid) - n_kid_units)}")
print(f"Alt rule (Kidney >0.1, unique cNMF): {alt_units} -> total {len(p_low) - (len(kid) - alt_units)}")
print(kid[["genenmf_mp", "cnmf_program", "jaccard", "sampleCoverage", "family", "kept"]].to_string(index=False))
print("tests per program:", allt.groupby(["cancer", "program"]).size().unique())
print("background sizes:", allt.background.unique())
print("programs enriched (q<0.05, >=3):", (best.top_hallmark_bh != "Unannotated").sum(), "/", len(best))
print("old named (raw p<0.05):", (best.top_hallmark_old != "Unannotated").sum())
same = (best.top_hallmark_bh == best.top_hallmark_old)
print("top pathway identical old vs BH:", same.sum(), "/", len(best))
print(best.loc[~same, ["cancer", "program", "top_hallmark_old", "pval_old", "overlap_old",
                       "top_hallmark_bh", "pval", "qval", "overlap"]].to_string(index=False))
