#!/usr/bin/env python
# 03_figA29_cnmf_usage_heatmap.py
#
# Figure A.29 (appendix_a_fig39) -- mean normalized usage of the 20 cNMF
# programs (rows) per cancer type (columns), 12 cancer types.
#   usage per cell normalized to sum 1 across the 20 programs, then averaged
#   over all cells of that cancer type.
#   red border = HC cNMF program (best match of >= 1 GeneNMF MP, Jaccard > 0.0538)
# Input: temp_objects/cnmf_output/<cancer>/cnmf_<cancer>/*.usages.k_20.dt_0_1.consensus.txt
#        hc_pairs_0.0538.csv (script 01)
# Output: appendix_a_fig39.png (400 dpi) / .pdf, figA29_mean_usage.csv

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = "/groups/ofircohen-group/users/hadasra"
OUT = f"{BASE}/thesis_figures/comment12_2026-10-21"
CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney", "Liver_Biliary",
           "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"]
LABEL = {"Head_and_Neck": "Head and Neck", "Liver_Biliary": "Liver/Biliary"}


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    mat = {}
    ncells = {}
    for c in CANCERS:
        u = pd.read_csv(f"{BASE}/temp_objects/cnmf_output/{c}/cnmf_{c}/cnmf_{c}.usages.k_20.dt_0_1.consensus.txt",
                        sep="\t", index_col=0)
        s = u.sum(axis=1)
        u = u[s > 0].div(s[s > 0], axis=0)
        mat[c] = u.mean(axis=0).values
        ncells[c] = len(u)
    progs = [f"cNMF_{i}" for i in range(1, 21)]
    M = pd.DataFrame(mat, index=progs)
    M.to_csv(f"{OUT}/figA29_mean_usage.csv")

    hc = pd.read_csv(f"{OUT}/hc_pairs_0.0538.csv")
    hcset = set(zip(hc.cancer, hc.cnmf_program))

    fig = plt.figure(figsize=(5.4, 5.6))
    ax = fig.add_axes([0.16, 0.20, 0.66, 0.77])
    cax = fig.add_axes([0.86, 0.20, 0.03, 0.77])
    vmax = np.ceil(M.values.max() * 20) / 20
    im = ax.imshow(M.values, cmap="YlGnBu", vmin=0, vmax=vmax, aspect="auto")
    n_hc = 0
    for j, c in enumerate(CANCERS):
        for i, p in enumerate(progs):
            if (c, p) in hcset:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#D00000", lw=1.6, zorder=3))
                n_hc += 1
    ax.set_xticks(range(len(CANCERS)))
    ax.set_xticklabels([LABEL.get(c, c) for c in CANCERS], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8.5)
    ax.set_yticks(range(20))
    ax.set_yticklabels([p.replace("_", " ") for p in progs], fontsize=8)
    ax.set_ylabel("cNMF program", fontsize=9)
    ax.tick_params(length=2, width=0.6)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Mean normalized usage", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.text(0.16, 0.025, "Red border = high-confidence program (cNMF and GeneNMF, Jaccard > 0.0538)",
             fontsize=8, ha="left", va="bottom", color="black")
    fig.savefig(f"{OUT}/appendix_a_fig39.png", dpi=400)
    fig.savefig(f"{OUT}/appendix_a_fig39.pdf")
    print(f"HC borders: {n_hc}; vmax {vmax}; max {M.values.max():.3f}; cells per cancer {ncells}")


if __name__ == "__main__":
    main()
