#!/usr/bin/env python
"""Comment 18: print redraw of the per-cluster marker dot plots.

Rows = Louvain clusters (res 0.50), columns = marker genes; clusters are
split into vertically stacked panels with at most MAX_COLS gene columns
each, sharing one colour bar and one size legend at the bottom.

Two modes, matching what the current thesis figures show:
  --mode log2fc   (figs 23/26/29/32/35/38) values straight from the saved
                  <lineage>_top_markers.csv: colour = avg_log2FC, size = pct.1,
                  each cluster shows only its own top-N genes (a gene shared by
                  two clusters gets one column, as in the original matplotlib
                  figure, whose categorical x-axis merged repeated names).
  --mode dotplot  (figs 51/52) Seurat DotPlot values exported by
                  extract_dotplot_data.R: colour = avg.exp.scaled, size = pct.exp,
                  full cluster x gene grid within each panel.

Usage:
  python plot_marker_dotplot.py --mode log2fc --lineage Myeloid_Prostate --top_n 5 \
      --out <dir>/appendix_a_fig23
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from c18_style import plt, FIG_W, MAX_H, save  # noqa: E402

MARKER_DIR = "/groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages"
MAX_COLS = 45

ap = argparse.ArgumentParser()
ap.add_argument("--mode", choices=["log2fc", "dotplot"], required=True)
ap.add_argument("--lineage", required=True)
ap.add_argument("--top_n", type=int, default=5)
ap.add_argument("--data", help="DotPlot CSV (dotplot mode)")
ap.add_argument("--out", required=True)
a = ap.parse_args()

# ---------------------------------------------------------------- data
# panels: list of (cluster_ids, genes, df with columns cl, gene, colour, size)
panels = []
if a.mode == "log2fc":
    mk = pd.read_csv(os.path.join(MARKER_DIR, f"{a.lineage}_top_markers.csv"))
    mk["cluster"] = mk["cluster"].astype(int)
    top = (mk.sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
             .groupby("cluster").head(a.top_n))
    clusters = sorted(top["cluster"].unique())
    per = max(1, MAX_COLS // a.top_n)
    n_pan = math.ceil(len(clusters) / per)
    size = math.ceil(len(clusters) / n_pan)
    for p in range(n_pan):
        cl = clusters[p * size:(p + 1) * size]
        sub = top[top["cluster"].isin(cl)]
        genes = list(dict.fromkeys(sub["gene"]))
        d = pd.DataFrame({"cl": sub["cluster"], "gene": sub["gene"],
                          "colour": sub["avg_log2FC"], "size": sub["pct.1"] * 100})
        panels.append((cl, genes, d))
    cbar_label = "Average log$_2$ fold change"
    size_label = "Cells expressing (%)"
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "ylorrd_print", plt.get_cmap("YlOrRd")(np.linspace(0.25, 1.0, 256)))
    all_c = pd.concat([p[2]["colour"] for p in panels])
    vmin, vmax = 0.0, float(all_c.max())
    size_ticks = [10, 25, 50, 75]
else:
    dd = pd.read_csv(a.data)
    for p in sorted(dd["panel"].unique()):
        sub = dd[dd["panel"] == p].copy()
        genes = list(sub.sort_values("gene_order")["features.plot"].drop_duplicates())
        cl = sorted(sub["id"].astype(int).unique())
        d = pd.DataFrame({"cl": sub["id"].astype(int), "gene": sub["features.plot"],
                          "colour": sub["avg.exp.scaled"], "size": sub["pct.exp"]})
        panels.append((cl, genes, d))
    cbar_label = "Scaled average expression"
    size_label = "Cells expressing (%)"
    # Seurat DotPlot palette (lightgrey -> #B2182B) with a darker grey low end for print
    cmap = mcolors.LinearSegmentedColormap.from_list("grey_red", ["#A6A6A6", "#B2182B"])
    all_c = pd.concat([p[2]["colour"] for p in panels])
    vmin, vmax = float(all_c.min()), float(all_c.max())
    size_ticks = [20, 40, 60, 80]

n_cols_max = max(len(g) for _, g, _ in panels)
assert n_cols_max <= MAX_COLS, n_cols_max

# ---------------------------------------------------------------- geometry (inches)
left, right = 0.42, 0.08
col_w = (FIG_W - left - right) / max(n_cols_max, 20)   # cap spacing for small objects
row_h = min(0.16, col_w * 1.25)
# rotated 8-pt gene labels, sized per panel
lab_hs = [max(len(g) for g in gs) * 8 * 0.58 / 72 + 0.10 for _, gs, _ in panels]
gap = 0.10
legend_h = 0.62
n_rows = sum(len(cl) + 0.5 for cl, _, _ in panels)
fixed_h = 0.08 + sum(lab_hs) + gap * len(panels) + legend_h
if fixed_h + n_rows * row_h > MAX_H:          # shrink rows to fit one page
    row_h = (MAX_H - fixed_h) / n_rows
    assert row_h >= 0.115, f"rows too tight for 8-pt labels ({row_h:.3f} in)"
panel_hs = [(len(cl) + 0.5) * row_h for cl, _, _ in panels]
fig_h = fixed_h + sum(panel_hs)

fig = plt.figure(figsize=(FIG_W, fig_h))
dot_max = (col_w * 72 * 0.92) ** 2           # area (pt^2) of a 100 % dot fits one column
def area(pct):
    return np.clip(pct, 0, 100) / 100.0 * dot_max

y_top = fig_h - 0.08
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
for (cl, genes, d), ph, lab_h in zip(panels, panel_hs, lab_hs):
    y0 = y_top - ph
    ax = fig.add_axes([left / FIG_W, y0 / fig_h, col_w * len(genes) / FIG_W, ph / fig_h])
    gx = {g: i for i, g in enumerate(genes)}
    cy = {c: i for i, c in enumerate(cl)}
    ax.scatter(d["gene"].map(gx), d["cl"].map(cy), s=area(d["size"].values),
               c=d["colour"].values, cmap=cmap, norm=norm,
               edgecolors="#333333", linewidths=0.25, zorder=3)
    ax.set_xlim(-0.5, len(genes) - 0.5)
    ax.set_ylim(-0.75, len(cl) - 0.25)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=8)
    ax.set_yticks(range(len(cl)))
    ax.set_yticklabels([str(c) for c in cl], fontsize=8)
    ax.set_ylabel("Cluster", fontsize=9)
    ax.yaxis.set_label_coords(-0.30 / (col_w * len(genes)), 0.5)   # same x in every panel
    ax.grid(True, color="#d0d0d0", linewidth=0.4, zorder=0)
    ax.tick_params(length=2, width=0.5, pad=1.5)
    for s in ax.spines.values():
        s.set_linewidth(0.6)
    y_top = y0 - lab_h - gap

# ---------------------------------------------------------------- shared legends
cax = fig.add_axes([left / FIG_W, 0.42 / fig_h, 1.9 / FIG_W, 0.09 / fig_h])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_label(cbar_label, fontsize=8, labelpad=2)
cb.ax.tick_params(labelsize=8, length=2, width=0.5, pad=1)
cb.outline.set_linewidth(0.5)

handles = [Line2D([], [], linestyle="", marker="o", markersize=math.sqrt(area(t)),
                  markerfacecolor="#555555", markeredgecolor="#333333",
                  markeredgewidth=0.25, label=f"{t}") for t in size_ticks]
fig.legend(handles=handles, title=size_label, loc="lower left",
           bbox_to_anchor=((left + 2.6) / FIG_W, 0.18 / fig_h), ncol=len(size_ticks),
           frameon=False, fontsize=8, title_fontsize=8, handletextpad=0.3,
           columnspacing=1.0, borderaxespad=0)

save(fig, a.out)
print(f"{a.lineage}: {sum(len(c) for c, _, _ in panels)} clusters, {len(panels)} panels, "
      f"max {n_cols_max} gene columns, top {a.top_n}")
