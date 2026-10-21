#!/usr/bin/env python
# 05_figA38_hc_program_associations.py
#
# Figure A.38 (appendix_a_fig48) -- adapted from
# create_sc_object_workflow/plots_and_outputs/make_fig44_hc_program_associations.py
# Changes vs the original: JACCARD_THR_EPI 0.1 -> 0.0538 (panel A); panel B
# unchanged (>= 0.1); no titles/suptitle; 5.4 in print width with panel B
# stacked under panel A; all text >= 8 pt; plain cancer labels.
#   A. epithelial: cNMF vs GeneNMF HC meta-programs (Jaccard > 0.0538), 11 cancer types
#      + Kidney zoom (19 programs in one column)
#   B. immune: GeneNMF vs Spectra HC meta-programs (Jaccard >= 0.1), 5 lineages
# y = sample coverage, color = Jaccard to best cross-method match, label = GeneNMF MP ID.
# Output: appendix_a_fig48.png (400 dpi) / .pdf

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

ATLAS = "/groups/ofircohen-group/users/hadasra/temp_objects"
EPI_DIR = os.path.join(ATLAS, "genenmf_output")
IMM_DIR = os.path.join(ATLAS, "genenmf_immune_output")
OUT = "/groups/ofircohen-group/users/hadasra/thesis_figures/comment12_2026-10-21"

CANCERS = ["Breast", "Kidney", "Lung", "Colorectal", "Skin", "Brain",
           "Liver_Biliary", "Ovarian", "Pancreas", "Prostate", "Sarcoma"]
CLABEL = {"Liver_Biliary": "Liver/Biliary"}
LINEAGES = ["B-cell", "T-cells", "NK-cells", "Myeloid",
            "Hematopoietic_Stem_and_Progenitor_cells"]
LINEAGE_LABEL = {"B-cell": "B cell", "T-cells": "T cell", "NK-cells": "NK",
                 "Myeloid": "Myeloid", "Hematopoietic_Stem_and_Progenitor_cells": "HSPC"}

JACCARD_THR_EPI = 0.0538  # strict >
JACCARD_THR_IMM = 0.1     # inclusive >=
FS = 8

rng = np.random.default_rng(0)


def load_epithelial():
    rows = []
    for c in CANCERS:
        met = pd.read_csv(os.path.join(EPI_DIR, c, f"{c}_genenmf_MP_metrics.csv"))
        met = met.rename(columns={met.columns[0]: "meta_program"})
        hc = pd.read_csv(os.path.join(EPI_DIR, c, f"{c}_highconfidence_programs.csv"))
        df = hc.merge(met[["meta_program", "sampleCoverage"]], on="meta_program", how="left")
        df = df[df["best_jaccard"] > JACCARD_THR_EPI].copy()
        df["group"] = CLABEL.get(c, c)
        df["mp_num"] = df["meta_program"].str.replace("MP", "", regex=False).astype(int)
        rows.append(df[["group", "mp_num", "sampleCoverage", "best_jaccard"]])
    return pd.concat(rows, ignore_index=True)


def load_immune():
    best = pd.read_csv(os.path.join(IMM_DIR, "Spectra_vs_GeneNMF_jaccard_best_match_all_lineages.csv"))
    rows = []
    for lin in LINEAGES:
        met = pd.read_csv(os.path.join(IMM_DIR, lin, f"{lin}_genenmf_MP_metrics.csv"))
        met = met.rename(columns={met.columns[0]: "meta_program"})
        df = best[best["lineage"] == lin].merge(met[["meta_program", "sampleCoverage"]],
                                               on="meta_program", how="left")
        df["best_jaccard"] = df["best_jaccard"].astype(float)
        df = df[df["best_jaccard"] >= JACCARD_THR_IMM].copy()
        df["group"] = LINEAGE_LABEL[lin]
        df["mp_num"] = df["meta_program"].str.replace("MP", "", regex=False).astype(int)
        rows.append(df[["group", "mp_num", "sampleCoverage", "best_jaccard"]])
    return pd.concat(rows, ignore_index=True)


def draw_panel(ax, df, groups, cmap, norm, skip_group=None, width=0.5):
    order = {g: i for i, g in enumerate(groups)}
    x = np.array([order[g] for g in df["group"]], float) + rng.uniform(-width / 2, width / 2, len(df))
    y = df["sampleCoverage"].values
    ax.scatter(x, y, c=df["best_jaccard"].values, cmap=cmap, norm=norm, s=22,
               edgecolor="black", linewidth=0.5, zorder=3)
    ax.axhline(0.5, color="black", linestyle="--", linewidth=0.7, zorder=1)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=40, ha="right", rotation_mode="anchor", fontsize=8.5)
    ax.set_ylabel("Sample coverage", fontsize=9)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.tick_params(labelsize=8, length=2, width=0.6)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, zorder=0)
    texts = [ax.text(xi, yi, str(n), fontsize=FS, ha="center", va="center", zorder=4)
             for xi, yi, n, g in zip(x, y, df["mp_num"], df["group"]) if g != skip_group]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
                expand=(1.15, 1.3), force_text=(0.3, 0.5), force_static=(0.2, 0.3))


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": FS,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    epi = load_epithelial()
    imm = load_immune()
    allj = pd.concat([epi["best_jaccard"], imm["best_jaccard"]])
    cmap = plt.get_cmap("viridis")
    norm = matplotlib.colors.Normalize(vmin=0.05, vmax=np.ceil(allj.max() * 10) / 10)
    groups_epi = [CLABEL.get(c, c) for c in CANCERS]

    W, H = 5.4, 7.6
    fig = plt.figure(figsize=(W, H))
    axA = fig.add_axes([0.55 / W, 4.30 / H, 3.45 / W, 3.05 / H])
    axK = fig.add_axes([4.40 / W, 4.30 / H, 0.92 / W, 3.05 / H])
    axB = fig.add_axes([0.55 / W, 0.75 / H, 2.6 / W, 2.75 / H])
    cax = fig.add_axes([3.75 / W, 0.75 / H, 0.13 / W, 2.75 / H])

    draw_panel(axA, epi, groups_epi, cmap, norm, skip_group="Kidney")
    kidx = groups_epi.index("Kidney")
    kid = epi[epi["group"] == "Kidney"]
    axA.add_patch(matplotlib.patches.Rectangle((kidx - 0.33, kid.sampleCoverage.min() - 0.03), 0.66,
                                               kid.sampleCoverage.max() - kid.sampleCoverage.min() + 0.06,
                                               fill=False, ec="black", lw=0.7, ls=(0, (2, 2)), zorder=5))

    xk = rng.uniform(-0.8, 0.8, len(kid))
    yk = kid["sampleCoverage"].values
    axK.scatter(xk, yk, c=kid["best_jaccard"].values, cmap=cmap, norm=norm, s=26,
                edgecolor="black", linewidth=0.5, zorder=3)
    kt = [axK.text(xi, yi, str(n), fontsize=FS, ha="center", va="center", zorder=4)
          for xi, yi, n in zip(xk, yk, kid["mp_num"])]
    axK.set_xlim(-1.1, 1.1)
    axK.set_ylim(yk.min() - 0.04, yk.max() + 0.04)
    adjust_text(kt, ax=axK, arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
                expand=(1.3, 1.6), force_text=(0.5, 0.8), force_static=(0.3, 0.4))
    axK.set_xticks([0])
    axK.set_xticklabels([f"Kidney (n = {len(kid)})"], fontsize=8.5)
    axK.tick_params(labelsize=8, length=2, width=0.6)
    axK.axhline(0.5, color="black", linestyle="--", linewidth=0.7, zorder=1)

    draw_panel(axB, imm, [LINEAGE_LABEL[l] for l in LINEAGES], cmap, norm)

    for ax, letter in [(axA, "A"), (axB, "B")]:
        ax.text(-0.13 if ax is axA else -0.18, 1.02, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="bottom", ha="left")

    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Jaccard similarity\n(best cross-method match)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    fig.text(0.55 / W, 0.12 / H,
             "A: cNMF and GeneNMF, Jaccard > 0.0538.  B: GeneNMF and Spectra, Jaccard ≥ 0.1.\n"
             "Numbers = GeneNMF meta-program ID; dashed line = 50% sample coverage.",
             fontsize=8, ha="left", va="bottom")

    fig.savefig(f"{OUT}/appendix_a_fig48.png", dpi=400)
    fig.savefig(f"{OUT}/appendix_a_fig48.pdf")
    print(f"Epithelial HC: {len(epi)} across {epi.group.nunique()} cancers "
          f"(Kidney {len(kid)}); per cancer {epi.group.value_counts().to_dict()}")
    print(f"Immune HC: {len(imm)} across {imm.group.nunique()} lineages; {imm.group.value_counts().to_dict()}")
    print(f"Jaccard range {allj.min():.3f}-{allj.max():.3f}")


if __name__ == "__main__":
    main()
