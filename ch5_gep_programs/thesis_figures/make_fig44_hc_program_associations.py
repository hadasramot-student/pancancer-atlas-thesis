#!/usr/bin/env python
# make_fig44_hc_program_associations.py
#
# Ch.5 Fig 44 -- "HC program associations across cancer types and immune
# lineages". Two panels:
#   A. Arm 1 (epithelial): cNMF vs GeneNMF HC programs (Jaccard > 0.1)
#      across 11 cancer types.
#   B. Arm 2 (immune): GeneNMF vs Spectra HC programs (Jaccard >= 0.1)
#      across 5 immune lineages.
# Each dot = one high-confidence meta-program. y = sample coverage (fraction
# of samples in which the program is active), color = Jaccard similarity to
# its best cross-method match, number label = GeneNMF meta-program ID.
# Dashed line at 50% sample coverage.
#
# Legibility redesign (label smear in the dense Kidney column, 16/19 HC
# programs in one x-position):
#   (a) label font significantly larger relative to dot size
#   (b) adjustText-based label repulsion + thin leader lines in both panels
#   (c) a zoomed inset callout for the Kidney column specifically, since
#       16 overlapping programs at ~0.5 cm of column width cannot be
#       resolved by repulsion alone at 6-inch print width
#
# Data is unchanged from the underlying pipeline outputs (no re-derivation
# of Jaccard/coverage values, no threshold changes) -- this script only
# assembles + plots them; dot positions and color encoding are exactly the
# values already computed by run_genenmf_per_cancer.R / run_genenmf_per_lineage.R.
#
# Inputs:
#   <ATLAS>/genenmf_output/<Cancer>/<Cancer>_genenmf_MP_metrics.csv        (sampleCoverage)
#   <ATLAS>/genenmf_output/<Cancer>/<Cancer>_highconfidence_programs.csv   (best_jaccard, recovered_by_both)
#   <ATLAS>/genenmf_immune_output/<Lineage>/<Lineage>_genenmf_MP_metrics.csv (sampleCoverage)
#   <ATLAS>/genenmf_immune_output/Spectra_vs_GeneNMF_jaccard_best_match_all_lineages.csv (best_jaccard)
#
# Output:
#   plots_and_outputs/Fig44_HC_program_associations.png

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle
from adjustText import adjust_text

ATLAS = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects"
EPI_DIR = os.path.join(ATLAS, "genenmf_output")
IMM_DIR = os.path.join(ATLAS, "genenmf_immune_output")
OUT_PNG = "/mnt/new_groups/ofircohen-group/users/hadasra/create_sc_object_workflow/plots_and_outputs/Fig44_HC_program_associations.png"

CANCERS = ["Breast", "Kidney", "Lung", "Colorectal", "Skin", "Brain",
           "Liver_Biliary", "Ovarian", "Pancreas", "Prostate", "Sarcoma"]
LINEAGES = ["B-cell", "T-cells", "NK-cells", "Myeloid",
            "Hematopoietic_Stem_and_Progenitor_cells"]
LINEAGE_LABEL = {"B-cell": "B-cell", "T-cells": "T-cell", "NK-cells": "NK",
                  "Myeloid": "Myeloid", "Hematopoietic_Stem_and_Progenitor_cells": "HSPC"}

JACCARD_THR_EPI = 0.1   # strict >   (matches run_genenmf_per_cancer.R --jaccard_thr default)
JACCARD_THR_IMM = 0.1   # inclusive >=

rng = np.random.default_rng(0)  # fixed seed -> reproducible jitter, not a data change


def load_epithelial():
    rows = []
    for c in CANCERS:
        met = pd.read_csv(os.path.join(EPI_DIR, c, f"{c}_genenmf_MP_metrics.csv"))
        met = met.rename(columns={met.columns[0]: "meta_program"})
        hc = pd.read_csv(os.path.join(EPI_DIR, c, f"{c}_highconfidence_programs.csv"))
        df = hc.merge(met[["meta_program", "sampleCoverage"]], on="meta_program", how="left")
        df = df[df["best_jaccard"] > JACCARD_THR_EPI].copy()
        df["group"] = c
        df["mp_num"] = df["meta_program"].str.replace("MP", "", regex=False).astype(int)
        rows.append(df[["group", "mp_num", "sampleCoverage", "best_jaccard"]])
    return pd.concat(rows, ignore_index=True)


def load_immune():
    best = pd.read_csv(os.path.join(IMM_DIR, "Spectra_vs_GeneNMF_jaccard_best_match_all_lineages.csv"))
    rows = []
    for lin in LINEAGES:
        met = pd.read_csv(os.path.join(IMM_DIR, lin, f"{lin}_genenmf_MP_metrics.csv"))
        met = met.rename(columns={met.columns[0]: "meta_program"})
        lin_best = best[best["lineage"] == lin]
        df = lin_best.merge(met[["meta_program", "sampleCoverage"]], on="meta_program", how="left")
        df["best_jaccard"] = df["best_jaccard"].astype(float)
        df = df[df["best_jaccard"] >= JACCARD_THR_IMM].copy()
        df["group"] = LINEAGE_LABEL[lin]
        df["mp_num"] = df["meta_program"].str.replace("MP", "", regex=False).astype(int)
        rows.append(df[["group", "mp_num", "sampleCoverage", "best_jaccard"]])
    return pd.concat(rows, ignore_index=True)


def jittered_x(df, groups, width=0.32):
    order = {g: i for i, g in enumerate(groups)}
    x = np.array([order[g] for g in df["group"]], dtype=float)
    x += rng.uniform(-width / 2, width / 2, size=len(df))
    return x


def draw_panel(ax, df, groups, cmap, norm, label_fontsize, dot_size, title,
                highlight_group=None, skip_label_for_highlight=False, title_fontsize=11):
    x = jittered_x(df, groups)
    y = df["sampleCoverage"].values
    c = df["best_jaccard"].values

    sc = ax.scatter(x, y, c=c, cmap=cmap, norm=norm, s=dot_size,
                     edgecolor="0.25", linewidth=0.5, zorder=3)

    ax.axhline(0.5, color="0.4", linestyle="--", linewidth=0.9, zorder=1)
    ax.text(len(groups) - 0.55, 0.5, "50%", fontsize=7.5, color="0.35",
            va="bottom", ha="right", style="italic")

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Sample coverage", fontsize=10)
    ax.set_ylim(-0.03, 1.05)
    ax.set_title(title, fontsize=title_fontsize, pad=8)
    ax.grid(axis="y", color="0.9", linewidth=0.6, zorder=0)

    texts = []
    for xi, yi, num, grp in zip(x, y, df["mp_num"], df["group"]):
        if skip_label_for_highlight and grp == highlight_group:
            continue  # labeled separately in the inset
        texts.append(ax.text(xi, yi, str(num), fontsize=label_fontsize,
                              fontweight="bold", ha="center", va="center", zorder=4))
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="0.35", lw=0.6),
                expand=(1.3, 1.6), force_text=(0.4, 0.6), force_static=(0.2, 0.3))
    return sc


def main():
    epi = load_epithelial()
    imm = load_immune()

    all_jaccard = pd.concat([epi["best_jaccard"], imm["best_jaccard"]])
    cmap = plt.get_cmap("viridis")
    norm = matplotlib.colors.Normalize(vmin=all_jaccard.min(), vmax=all_jaccard.max())

    fig = plt.figure(figsize=(17.0, 6.2), dpi=300)
    # Dedicated column for the Kidney zoom (a real subplot, not a floating
    # overlay) so it structurally cannot cover data in neighboring columns.
    gs = fig.add_gridspec(1, 3, width_ratios=[11, 3.2, 5], wspace=0.38)
    axA = fig.add_subplot(gs[0, 0])
    axK = fig.add_subplot(gs[0, 1])
    axB = fig.add_subplot(gs[0, 2])

    # ---- Panel A: epithelial, 11 cancer types (Kidney gets inset callout) --
    draw_panel(axA, epi, CANCERS, cmap, norm, label_fontsize=10.5, dot_size=70,
               title="A. Epithelial HC programs across cancer types\n"
                     "(cNMF ∩ GeneNMF, Jaccard > 0.1)",
               highlight_group="Kidney", skip_label_for_highlight=True)
    axA.text(0.01, 0.99, "A", transform=axA.transAxes, fontsize=14, fontweight="bold",
              va="top", ha="left")

    # Kidney zoom: a genuinely separate subplot (axK, its own gridspec cell),
    # not a floating overlay -- so it cannot cover data points in neighboring
    # columns (Colorectal/Brain), unlike an inset_axes positioned by
    # axes-fraction on top of axA. The main panel already drew the Kidney
    # dots (unlabeled); axK re-renders just that column, spread wider on x,
    # with its own labels -- same dot positions (coverage values), just
    # spread further apart for legibility.
    kidney = epi[epi["group"] == "Kidney"]
    kidney_x_idx = CANCERS.index("Kidney")
    yk = kidney["sampleCoverage"].values
    ck = kidney["best_jaccard"].values
    ymin, ymax = yk.min() - 0.04, yk.max() + 0.04

    xk = rng.uniform(-0.85, 0.85, size=len(kidney))
    axK.scatter(xk, yk, c=ck, cmap=cmap, norm=norm, s=110,
                edgecolor="0.25", linewidth=0.6, zorder=3)
    ktexts = [axK.text(xi, yi, str(num), fontsize=10, fontweight="bold",
                        ha="center", va="center", zorder=4)
              for xi, yi, num in zip(xk, yk, kidney["mp_num"])]
    adjust_text(ktexts, ax=axK,
                arrowprops=dict(arrowstyle="-", color="0.35", lw=0.6),
                expand=(1.5, 1.9), force_text=(0.6, 0.9), force_static=(0.3, 0.4))
    axK.set_xlim(-1.0, 1.0)
    axK.set_ylim(max(-0.02, ymin), min(1.05, ymax))
    axK.set_xticks([])
    axK.set_title(f"Kidney zoom\n(n={len(kidney)} HC programs)", fontsize=9.5, pad=6)
    axK.tick_params(axis="y", labelsize=8)
    for spine in axK.spines.values():
        spine.set_edgecolor("0.3")
        spine.set_linewidth(0.9)

    # dotted rectangle around the real Kidney column in axA + connector
    # lines to axK, drawn with figure-level ConnectionPatch so it works
    # across separate axes (no risk of drawing over other columns' data).
    rect_x0, rect_x1 = kidney_x_idx - 0.35, kidney_x_idx + 0.35
    rect_y0, rect_y1 = max(-0.02, ymin), min(1.05, ymax)
    axA.add_patch(Rectangle((rect_x0, rect_y0), rect_x1 - rect_x0, rect_y1 - rect_y0,
                             fill=False, ec="0.3", lw=0.9, ls=(0, (2, 2)), zorder=5))
    for (xa, ya), (xb, yb) in [((rect_x1, rect_y1), (-1.0, rect_y1)),
                                ((rect_x1, rect_y0), (-1.0, rect_y0))]:
        con = ConnectionPatch(xyA=(xa, ya), coordsA=axA.transData,
                               xyB=(xb, yb), coordsB=axK.transData,
                               color="0.3", lw=0.8, ls=(0, (2, 2)), zorder=1)
        fig.add_artist(con)

    # ---- Panel B: immune, 5 lineages -------------------------------------
    draw_panel(axB, imm, [LINEAGE_LABEL[l] for l in LINEAGES], cmap, norm,
               label_fontsize=10.5, dot_size=70,
               title="B. Immune HC programs\nacross lineages\n"
                     "(GeneNMF ∩ Spectra,\nJaccard ≥ 0.1)", title_fontsize=10)
    axB.text(0.02, 0.99, "B", transform=axB.transAxes, fontsize=14, fontweight="bold",
              va="top", ha="left")
    axB.set_ylabel("")

    # Finalize axes geometry BEFORE creating the colorbar. fig.colorbar(ax=[...])
    # shrinks the given axes in place to make room for itself, based on
    # whatever the figure margins are *at that moment*; calling
    # subplots_adjust afterward would reposition axA/axK/axB again without
    # moving the already-placed colorbar axes, desyncing the two and letting
    # axB's plot box grow back out under the colorbar (this previously cut
    # off the HSPC column's dots, not just the title).
    fig.subplots_adjust(top=0.74, bottom=0.16)

    # ---- shared colorbar ---------------------------------------------------
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[axA, axK, axB], orientation="vertical", fraction=0.02, pad=0.10)
    cbar.set_label("Jaccard similarity\n(best cross-method match)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.suptitle("HC program associations across cancer types and immune lineages",
                 fontsize=13, y=0.97)

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")
    print(f"Epithelial HC programs: {len(epi)} across {epi['group'].nunique()} cancer types "
          f"(Kidney: {len(kidney)})")
    print(f"Immune HC programs: {len(imm)} across {imm['group'].nunique()} lineages")


if __name__ == "__main__":
    main()
