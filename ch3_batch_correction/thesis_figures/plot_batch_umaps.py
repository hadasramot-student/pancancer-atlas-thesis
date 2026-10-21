#!/usr/bin/env python
"""Comment 18: print redraw of the Ch.3 batch-correction UMAPs (Breast benchmark).

Coordinates come from extract_batch_umap_data.R (same embeddings as
pancancer-atlas-thesis/ch3_batch_correction/plot_publication_figures.R; the
Uncorrected and scANVI UMAPs are recomputed there with the same settings:
30 neighbours, min_dist 0.3).

Kinds:
  overview  5 UMAPs (3 + 2) coloured by study/project, no legend
  subtypes  5 UMAPs (3 + 2) coloured by Scimilarity subtype (top 10 + Other),
            one legend in the empty 6th slot
  ch3_fig1  T cells, Uncorrected | scANVI by study, no legend (3.6 x 2.9 in)
  ch3_fig2  T cells, Uncorrected | scANVI by subtype, legend right (5.4 x 2.9 in)

Subtype colours follow the original script exactly: generic labels demoted to
"Other", top 10 by count (Other may take a slot), kept labels sorted
alphabetically (case-insensitive, as R sort() did), palette in that order,
Other = #BBBBBB. Labels are then shortened for display only (RENAME).

Usage:
  python plot_batch_umaps.py --lineage T-cells --kind subtypes --out <dir>/appendix_a_fig54
"""
import argparse
import os
import sys
import textwrap

import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from c18_style import plt, FIG_W, save  # noqa: E402

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
METHODS = ["Uncorrected", "MNN", "ComBat", "Harmony", "scANVI"]
TOP_N = 10
OTHER_COL = "#BBBBBB"

# discrete_colors() of plot_publication_figures.R
BASE_COLS = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
    "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB",
    "#E78AC3", "#A6D854", "#FFD92F", "#1B9E77", "#D95F02",
    "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
    "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#FDBF6F",
    "#CAB2D6", "#6A3D9A", "#B15928", "#FF7F00", "#FFFF99",
]
# study colours carry no legend, so drop the pale yellows for print
STUDY_COLS = [c for c in BASE_COLS if c not in ("#FFD92F", "#FFFF99")]

DEMOTE = {"Endothelial": ["endothelial cell"], "Fibroblasts": ["fibroblast"]}

APPROVED_TCELL = {
    "Activated CD4+ T cell": "#e41a1c", "Activated CD8+ T cell": "#377eb8",
    "CD4+ T cell": "#4daf4a", "CD8+ cytotoxic T cell": "#ff7f00",
    "CD8+ memory T cell": "#984ea3", "CD8+ T cell": "#a65628",
    "Gamma-delta T cell": "#f781bf", "Naive CD4+ T cell": "#66c2a5",
    "Regulatory T cell": "#fc8d62", "T cell (unspecified)": "#8da0cb",
    "Other": "#bbbbbb",
}

# display names (Cell Ontology -> short); anything not listed is only capitalised
RENAME = {
    # T cells
    "activated CD4-positive, alpha-beta T cell": "Activated CD4+ T cell",
    "activated CD8-positive, alpha-beta T cell": "Activated CD8+ T cell",
    "CD4-positive, alpha-beta T cell": "CD4+ T cell",
    "CD8-positive, alpha-beta cytotoxic T cell": "CD8+ cytotoxic T cell",
    "CD8-positive, alpha-beta memory T cell": "CD8+ memory T cell",
    "CD8-positive, alpha-beta T cell": "CD8+ T cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "Naive CD4+ T cell",
    "T cell": "T cell (unspecified)",
    # B cells
    "B cell": "B cell (unspecified)",
    # Endothelial
    "endothelial cell of artery": "Arterial endothelial cell",
    "endothelial cell of lymphatic vessel": "Lymphatic endothelial cell",
    "endothelial cell of vascular tree": "Vascular endothelial cell",
    # Fibroblasts
    "myofibroblast cell": "Myofibroblast",
    "stromal cell of ovary": "Ovarian stromal cell",
    # Myeloid
    "CD14-positive, CD16-positive monocyte": "CD14+ CD16+ monocyte",
    # NK cells
    "CD16-negative, CD56-bright natural killer cell, human": "CD56bright CD16− NK cell",
    "CD16-positive, CD56-dim natural killer cell, human": "CD56dim CD16+ NK cell",
    "natural killer cell": "NK cell",
    "mature NK T cell": "Mature NKT cell",
}


def display(lab):
    if lab in RENAME:
        return RENAME[lab]
    return lab[:1].upper() + lab[1:]


def ramp(cols, n):
    """R colorRampPalette(cols)(n) (linear RGB interpolation)."""
    if n <= len(cols):
        return cols[:n]
    rgb = np.array([mcolors.to_rgb(c) for c in cols])
    x = np.linspace(0, 1, len(cols))
    xi = np.linspace(0, 1, n)
    return [mcolors.to_hex([np.interp(xi[i], x, rgb[:, k]) for k in range(3)])
            for i in range(n)]


def subtype_colours(raw, lineage):
    v = raw.astype(str)
    dem = [d.lower() for d in DEMOTE.get(lineage, [])]
    v = v.where(~v.str.lower().isin(dem), "Other")
    # R: sort(table(v), decreasing = TRUE) is stable, so ties keep table()'s
    # alphabetical order
    vc = v.value_counts()
    keep = sorted(vc.index, key=lambda k: (-vc[k], k.lower(), k))[:TOP_N]
    v = v.where(v.isin(keep), "Other")
    shown = sorted([k for k in keep if k != "Other"], key=lambda s: (s.lower(), s))
    cmap = dict(zip(shown, BASE_COLS[:len(shown)]))
    cmap["Other"] = OTHER_COL
    return v, shown, cmap


ap = argparse.ArgumentParser()
ap.add_argument("--lineage", required=True)
ap.add_argument("--kind", choices=["overview", "subtypes", "ch3_fig1", "ch3_fig2"], required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

df = pd.read_csv(os.path.join(DATA_DIR, f"{a.lineage}_Breast_umapdata.csv.gz"))
n = len(df)
rng = np.random.default_rng(42)
order = rng.permutation(n)
df = df.iloc[order].reset_index(drop=True)

by_study = a.kind in ("overview", "ch3_fig1")
if by_study:
    studies = sorted(df["project_name"].astype(str).unique(), key=lambda s: (s.lower(), s))
    smap = dict(zip(studies, ramp(STUDY_COLS, len(studies))))
    colours = df["project_name"].astype(str).map(smap).values
    legend_items = None
else:
    v, shown, cmap = subtype_colours(df["scimilarity_celltype"], a.lineage)
    colours = v.map(cmap).values
    legend_items = [(display(k), cmap[k]) for k in shown] + [("Other", OTHER_COL)]
    if a.lineage == "T-cells":
        got = {lab: c.lower() for lab, c in legend_items}
        assert got == APPROVED_TCELL, f"T-cell colours differ from approved:\n{got}"
    print("subtypes shown: " + "; ".join(f"{display(k)} [{k}]" if display(k) != k else k
                                         for k in shown))

methods = METHODS if a.kind in ("overview", "subtypes") else ["Uncorrected", "scANVI"]

# point size (pt^2) / alpha by cell number: small objects get bigger dots
s_pt = float(np.clip(8000.0 / n, 0.3, 6.0))
alpha = 0.6 if n > 10000 else 0.8


def panel(ax, m):
    x, y = df[f"{m}_1"].values, df[f"{m}_2"].values
    ax.scatter(x, y, s=s_pt, c=colours, alpha=alpha, linewidths=0, marker="o",
               rasterized=True)
    px, py = 0.03 * np.ptp(x), 0.03 * np.ptp(y)
    ax.set_xlim(x.min() - px, x.max() + px)
    ax.set_ylim(y.min() - py, y.max() + py)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
        sp.set_color("#555555")
    ax.set_title(m, fontsize=9, pad=3)


def legend(ax_or_fig, items, **kw):
    handles = [Line2D([], [], linestyle="", marker="o", markersize=5,
                      markerfacecolor=c, markeredgewidth=0,
                      label=textwrap.fill(lab, 28)) for lab, c in items]
    return ax_or_fig.legend(handles=handles, title="Cell type", fontsize=8, title_fontsize=8,
                            frameon=False, handlelength=1.0, handletextpad=0.4, labelspacing=0.22,
                            borderaxespad=0, alignment="left", **kw)


if a.kind in ("overview", "subtypes"):
    left, right, top, bottom = 0.05, 0.05, 0.22, 0.05
    hgap, vgap = 0.10, 0.25
    pw = (FIG_W - left - right - 2 * hgap) / 3
    ph = pw * 0.92
    fig_h = top + 2 * ph + vgap + bottom
    fig = plt.figure(figsize=(FIG_W, fig_h))

    def add(col, row, x_off=0.0):
        x0 = left + x_off + col * (pw + hgap)
        y0 = fig_h - top - ph - row * (ph + vgap)
        return fig.add_axes([x0 / FIG_W, y0 / fig_h, pw / FIG_W, ph / fig_h])

    # overview: centre the two bottom panels; subtypes: legend takes slot 6
    x_off = (pw + hgap) / 2 if a.kind == "overview" else 0.0
    slots = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 1, x_off), (1, 1, x_off)]
    for m, (c, r, xo) in zip(methods, slots):
        panel(add(c, r, xo), m)
    if a.kind == "subtypes":
        lax = add(2, 1)
        lax.axis("off")
        # start level with the panel titles so 11 entries fit beside the panels
        legend(lax, legend_items, loc="upper left",
               bbox_to_anchor=(-0.02, 1 + (vgap - 0.02) / ph))
else:
    fw = 3.6 if a.kind == "ch3_fig1" else FIG_W
    fig_h = 2.9
    fig = plt.figure(figsize=(fw, fig_h))
    left, top, bottom, hgap = 0.05, 0.22, 0.05, 0.10
    pw = (3.6 - 2 * left - hgap) / 2
    ph = fig_h - top - bottom
    for i, m in enumerate(methods):
        x0 = left + i * (pw + hgap)
        ax = fig.add_axes([x0 / fw, bottom / fig_h, pw / fw, ph / fig_h])
        panel(ax, m)
    if a.kind == "ch3_fig2":
        lax = fig.add_axes([(3.6 + 0.05) / fw, bottom / fig_h,
                            (fw - 3.65) / fw, ph / fig_h])
        lax.axis("off")
        legend(lax, legend_items, loc="center left", bbox_to_anchor=(0.0, 0.5))

save(fig, a.out)
print(f"{a.lineage} {a.kind}: {n} cells, {df['project_name'].nunique()} projects, "
      f"point size {s_pt:.2f} pt^2, alpha {alpha}")
