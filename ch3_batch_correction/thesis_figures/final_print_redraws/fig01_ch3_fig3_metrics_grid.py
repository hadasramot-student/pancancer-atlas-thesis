"""ch3_fig3 (thesis Figure 4.3): batch-correction metric grid, print version.

Redraw of pancancer-atlas-thesis/ch3_batch_correction/plot_fig17_metrics_grid_large_fonts.R
(same input CSV, same ranking and best-method rule), at 5.4 in width without title/subtitle.
Rows share one y scale per metric (as facet_grid(scales = "free_y") did).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from printstyle import save, WIDTH

CSV = "/groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/all_validation_results.csv"
df = pd.read_csv(CSV)

METHODS = ["mnn", "combat", "harmony", "scanvi"]
MLAB = {"mnn": "MNN", "combat": "ComBat", "harmony": "Harmony", "scanvi": "scANVI"}
MCOL = {"MNN": "#377EB8", "ComBat": "#4DAF4A", "Harmony": "#FF7F00", "scANVI": "#E41A1C"}
CT = ["B-cell", "Endothelial", "Fibroblasts", "Hematopoietic_Stem_and_Progenitor_cells",
      "Myeloid", "NK-cells", "T-cells"]
CTLAB = ["B cells", "Endothelial", "Fibroblasts", "HSPCs", "Myeloid", "NK cells", "T cells"]
METRICS = [("batch_lisi", "Batch LISI\n(higher =\nbetter)", "higher"),
           ("label_lisi", "Label LISI\n(lower =\nbetter)", "lower"),
           ("kbet_accept", "kBET\nacceptance\n(higher =\nbetter)", "higher"),
           ("graph_connectivity", "Graph\nconnectivity\n(higher =\nbetter)", "higher"),
           ("f1_label", "F1 label\nscore\n(higher =\nbetter)", "higher")]

fig, axs = plt.subplots(len(METRICS), len(CT), figsize=(WIDTH, 6.2), sharey="row")
fig.subplots_adjust(left=0.18, right=0.995, top=0.955, bottom=0.085, wspace=0.06, hspace=0.18)

for i, (m, mlab, direction) in enumerate(METRICS):
    ymax = df[m].max() * 1.08
    for j, ct in enumerate(CT):
        ax = axs[i, j]
        sub = df[df.cell_type == ct].set_index("method").loc[METHODS]
        vals = sub[m].values
        order = np.argsort(-vals if direction == "higher" else vals, kind="stable")
        best = vals.max() if direction == "higher" else vals.min()
        for r, k in enumerate(order):
            lab = MLAB[METHODS[k]]
            isb = vals[k] == best
            ax.bar(r, vals[k], width=0.75, color=MCOL[lab], edgecolor="black",
                   linewidth=1.6 if isb else 0.35, zorder=3 if isb else 2)
        ax.set_xticks([]); ax.set_xlim(-0.6, 3.6); ax.set_ylim(0, ymax)
        ax.grid(axis="y", color="#b0b0b0", linewidth=0.4, zorder=0)
        ax.tick_params(axis="y", length=2, pad=1)
        if i == 0:
            ax.set_title(CTLAB[j], fontsize=8, pad=3)
        if j == 0:
            ax.set_ylabel(mlab, rotation=0, ha="right", va="center", fontsize=8, labelpad=4)
    axs[i, 0].yaxis.set_major_locator(plt.MaxNLocator(3))

handles = [Patch(facecolor=MCOL[MLAB[m]], edgecolor="black", linewidth=0.35, label=MLAB[m]) for m in METHODS]
handles.append(Patch(facecolor="white", edgecolor="black", linewidth=1.6, label="Best method"))
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, title="Method",
           bbox_to_anchor=(0.58, 0.0), handlelength=1.4, columnspacing=1.0)
save(fig, "ch3_fig3")
