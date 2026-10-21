"""ch2_fig2 (thesis Figure 4.4): 4-panel atlas scANVI UMAP, print version.

Redraw of pancancer-atlas-thesis/ch3_batch_correction/atlas_200k_scanvi/05_figures.py (Figure 1 there).
Input is the master table written by that script (same 200,000 cells, same scANVI UMAP
coordinates, same colours, same seed-42 plotting order). Changes: 5.4 in width, no suptitle,
plain legend titles/labels, no UMAP1/UMAP2 labels; "Other" grey darkened for print.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from printstyle import save, WIDTH, plain

ROOT = "/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi"
df = pd.read_csv(f"{ROOT}/figures/fig12_master_table.csv")
RNG = np.random.default_rng(42)

META_ORDER = ["Epithelial", "T-cells", "Myeloid", "B-cell", "Fibroblasts",
              "Endothelial", "NK-cells", "Hematopoietic_Stem_and_Progenitor_cells", "other"]
META_COL = dict(zip(META_ORDER, ["#3b7dd8", "#e4572e", "#f2b134", "#5bb867", "#9b5de5",
                                 "#00bbf9", "#ff85a1", "#8d99ae", "#9a9ca5"]))  # other: #c9cad1 -> #9a9ca5
REPO_ORDER = ["3CA", "ScBaseCount", "Breastscape"]
REPO_COL = {"3CA": "#e4572e", "ScBaseCount": "#3b7dd8", "Breastscape": "#5bb867"}
AGREE_ORDER = ["Full", "Partial", "Other/Unknown"]
AGREE_COL = {"Full": "#2a9d8f", "Partial": "#e9c46a", "Other/Unknown": "#9a9ca5"}  # was #c9cad1

cancers = sorted(df["cancer_type"].unique(), key=lambda c: -(df["cancer_type"] == c).sum())
base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
CAN_COL = {lv: base[i % len(base)] for i, lv in enumerate(cancers)}

xy = df[["scANVI_UMAP1", "scANVI_UMAP2"]].values
PANELS = [("Consensus_Cell_Type", META_ORDER, META_COL, "a", "Cell type"),
          ("cancer_type", cancers, CAN_COL, "b", "Cancer type"),
          ("source_repository", REPO_ORDER, REPO_COL, "c", "Source repository"),
          ("agreement_level", AGREE_ORDER, AGREE_COL, "d", "Agreement")]

H = 6.5
fig = plt.figure(figsize=(WIDTH, H))
S = 2.0                       # UMAP panel side (in); legend in 2 columns underneath
pos = [(0.35, 4.25), (3.05, 4.25), (0.35, 0.62), (3.05, 0.62)]   # lower-left corners (in)
for (col, order, cmap, letter, ltitle), (x0, y0) in zip(PANELS, pos):
    ax = fig.add_axes([x0 / WIDTH, y0 / H, S / WIDTH, S / H])
    idx = RNG.permutation(len(df))
    pxy, cats = xy[idx], df[col].values[idx]
    for lv in order:
        m = cats == lv
        if m.any():
            ax.scatter(pxy[m, 0], pxy[m, 1], s=0.2, c=[cmap[lv]], linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    cnt = df[col].value_counts().to_dict()
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=4, markerfacecolor=cmap[lv],
                      markeredgewidth=0, label=f"{plain(lv)} ({cnt[lv]:,})")
               for lv in order if cnt.get(lv, 0) > 0]
    ax.legend(handles=handles, title=ltitle, loc="upper left", bbox_to_anchor=(-0.19, 0.0),
              ncol=2 if len(handles) > 4 else 1, frameon=False, handletextpad=0.1, labelspacing=0.2, columnspacing=0.4,
              borderaxespad=0.1, alignment="left")
save(fig, "ch2_fig2")
