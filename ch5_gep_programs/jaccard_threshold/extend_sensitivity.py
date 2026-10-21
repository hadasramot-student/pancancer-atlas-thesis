#!/usr/bin/env python
"""
Extends the cNMF-GeneNMF HC-designation Jaccard-threshold sensitivity to the
full range requested (adds 0.2 and 0.5 to the existing 0.05-0.15 run).

Re-reads each cancer's <Cancer>_highconfidence_programs.csv (best_jaccard per
GeneNMF meta-program vs its best-matching cNMF program; already computed by
genenmf/run_genenmf_per_cancer.R) and recounts how many meta-programs clear
each threshold -- same method as the original 0.05/0.075/0.1/0.125/0.15 run,
just re-thresholding the same best_jaccard values, no re-clustering needed.
"""
import pandas as pd
import os

BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output"
CANCERS = ["Kidney", "Brain", "Colorectal", "Pancreas", "Breast",
           "Liver_Biliary", "Prostate", "Lung", "Ovarian", "Sarcoma", "Skin"]
THRESHOLDS = [0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.5]

best_jaccard = {}
max_overall = 0.0
max_where = None
for cancer in CANCERS:
    path = os.path.join(BASE, cancer, f"{cancer}_highconfidence_programs.csv")
    df = pd.read_csv(path)
    best_jaccard[cancer] = df["best_jaccard"].values
    if len(df):
        i = df["best_jaccard"].idxmax()
        if df["best_jaccard"][i] > max_overall:
            max_overall = df["best_jaccard"][i]
            max_where = f"{cancer} {df['meta_program'][i]} vs {df['best_cnmf_match'][i]}"

rows = []
for thr in THRESHOLDS:
    row = {"threshold": thr}
    total = 0
    for cancer in CANCERS:
        n = int((best_jaccard[cancer] > thr).sum())
        row[cancer] = n
        total += n
    row["Total"] = total
    rows.append(row)

out = pd.DataFrame(rows)[["threshold"] + CANCERS + ["Total"]]
out_csv = os.path.join(BASE, "jaccard_threshold_sensitivity", "jaccard_threshold_sensitivity_full.csv")
out.to_csv(out_csv, index=False)
print(out.to_string(index=False))
print(f"\nHighest single best_jaccard observed anywhere: {max_overall:.4f} ({max_where})")
print(f"Wrote: {out_csv}")
