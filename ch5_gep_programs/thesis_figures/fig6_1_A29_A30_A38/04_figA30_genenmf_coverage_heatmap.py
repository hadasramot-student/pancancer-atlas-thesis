#!/usr/bin/env python
# 04_figA30_genenmf_coverage_heatmap.py
#
# Figure A.30 (appendix_a_fig40) -- fraction of samples in which each GeneNMF
# meta-program is active (sampleCoverage), per cancer type (11 types; Head and
# Neck has no GeneNMF run).
#   red border = HC meta-program (best cNMF match Jaccard > 0.0538)
#   white cell = meta-program index not produced for that cancer type
# Input: temp_objects/genenmf_output/<cancer>/<cancer>_genenmf_MP_metrics.csv
#        hc_pairs_0.0538.csv (script 01)
# Output: appendix_a_fig40.png (400 dpi) / .pdf, figA30_sample_coverage.csv

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = "/groups/ofircohen-group/users/hadasra"
OUT = f"{BASE}/thesis_figures/comment12_2026-10-21"
CANCERS = ["Brain", "Breast", "Colorectal", "Kidney", "Liver_Biliary",
           "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"]
LABEL = {"Liver_Biliary": "Liver/Biliary"}


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    cov = {}
    for c in CANCERS:
        m = pd.read_csv(f"{BASE}/temp_objects/genenmf_output/{c}/{c}_genenmf_MP_metrics.csv", index_col=0)
        cov[c] = m["sampleCoverage"]
    M = pd.DataFrame(cov)
    nmax = max(int(i[2:]) for i in M.index)
    mps = [f"MP{i}" for i in range(1, nmax + 1)]
    M = M.reindex(mps)
    M.to_csv(f"{OUT}/figA30_sample_coverage.csv")

    hc = pd.read_csv(f"{OUT}/hc_pairs_0.0538.csv")
    hcset = set(zip(hc.cancer, hc.genenmf_mp))

    fig = plt.figure(figsize=(5.4, 5.6))
    ax = fig.add_axes([0.14, 0.20, 0.68, 0.77])
    cax = fig.add_axes([0.86, 0.20, 0.03, 0.77])
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("white")
    im = ax.imshow(np.ma.masked_invalid(M.values.astype(float)), cmap=cmap, vmin=0, vmax=1, aspect="auto")
    n_hc = 0
    for j, c in enumerate(CANCERS):
        for i, p in enumerate(mps):
            if (c, p) in hcset:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#D00000", lw=1.6, zorder=3))
                n_hc += 1
    ax.set_xticks(range(len(CANCERS)))
    ax.set_xticklabels([LABEL.get(c, c) for c in CANCERS], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8.5)
    ax.set_yticks(range(len(mps)))
    ax.set_yticklabels(mps, fontsize=8)
    ax.set_ylabel("GeneNMF meta-program", fontsize=9)
    ax.tick_params(length=2, width=0.6)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Fraction of samples with program active", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.text(0.14, 0.045, "Red border = high-confidence program (cNMF and GeneNMF, Jaccard > 0.0538)",
             fontsize=8, ha="left", va="bottom")
    fig.text(0.14, 0.015, "White = meta-program not produced for that cancer type",
             fontsize=8, ha="left", va="bottom")
    fig.savefig(f"{OUT}/appendix_a_fig40.png", dpi=400)
    fig.savefig(f"{OUT}/appendix_a_fig40.pdf")
    print(f"HC borders: {n_hc}; MPs per cancer: {M.notna().sum().to_dict()}; "
          f"coverage range {np.nanmin(M.values):.3f}-{np.nanmax(M.values):.3f}")


if __name__ == "__main__":
    main()
