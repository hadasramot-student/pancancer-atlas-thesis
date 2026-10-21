"""Ch.2 Fig 1: consensus cell-type composition per cancer type, redrawn at print size.

Input: counts_long.csv (cancer_type x Consensus_Cell_Type counts) from count_composition.R,
computed on the full 8,107,665-cell deduplicated atlas. Only the plot styling differs from
../make_bar_chart.py.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})

label_map = {
    "Epithelial": "Epithelial/Cancer", "T-cells": "T cells", "Myeloid": "Myeloid",
    "B-cell": "B cells", "NK-cells": "NK cells", "Fibroblasts": "Fibroblasts/Stroma",
    "Endothelial": "Endothelial/Vascular",
    "Hematopoietic_Stem_and_Progenitor_cells": "HSC/Stem-like", "other": "Other/Unknown",
}
order = list(label_map.values())
colors = ["#e64b35", "#4dbbd5", "#00a087", "#3c5488", "#f39b7f",
          "#8491b4", "#91d1c2", "#dc0000", "#b09c85"]
cancer_map = {"Liver_Biliary": "Liver/Biliary", "Head_and_Neck": "Head and Neck",
              "Other_Models": "Other models", "unknown": "Unknown"}

expected = {"Epithelial/Cancer": 2443133, "T cells": 1510735, "Myeloid": 979202,
            "B cells": 371455, "Fibroblasts/Stroma": 264397, "Endothelial/Vascular": 257365,
            "NK cells": 90334, "HSC/Stem-like": 32492, "Other/Unknown": 2158552}
EXPECTED_TOTAL = 8107665

df = pd.read_csv(os.path.join(HERE, "counts_long.csv"))
unmapped = set(df["Consensus_Cell_Type"]) - set(label_map)
if unmapped or df["cancer_type"].isna().any():
    sys.exit(f"STOP: unmapped cell types {unmapped} or NA cancer_type")
df["cell_type"] = df["Consensus_Cell_Type"].map(label_map)

# --- totals check: stop on any mismatch ---
tot = df.groupby("cell_type")["n"].sum()
bad = {k: (int(tot.get(k, 0)), v) for k, v in expected.items() if int(tot.get(k, 0)) != v}
print("Cell-type totals (observed vs thesis):")
for k in order:
    print(f"  {k:22s} {int(tot.get(k, 0)):>10,} {expected[k]:>10,}")
print(f"  {'Total':22s} {int(tot.sum()):>10,} {EXPECTED_TOTAL:>10,}")
if bad or int(tot.sum()) != EXPECTED_TOTAL:
    sys.exit(f"STOP: totals differ from thesis: {bad}; total {int(tot.sum())}")

counts = (df.pivot_table(index="cancer_type", columns="cell_type", values="n", aggfunc="sum")
          .fillna(0).astype(int).reindex(columns=order, fill_value=0))
pct = counts.div(counts.sum(axis=1), axis=0) * 100
pct = pct.sort_values("Epithelial/Cancer", ascending=True)  # highest at top in barh
counts = counts.loc[pct.index]
labels = [cancer_map.get(c, c) for c in pct.index]

# --- composition table ---
tab = pd.DataFrame({"cancer_type": labels, "total_cells": counts.sum(axis=1).values})
for ct in order:
    tab[f"{ct}_n"] = counts[ct].values
    tab[f"{ct}_pct"] = pct[ct].round(2).values
tab.iloc[::-1].to_csv(os.path.join(HERE, "ch2_fig1_composition.csv"), index=False)

# --- figure at final print size ---
fig, ax = plt.subplots(figsize=(5.4, 4.2))
y = np.arange(len(pct))
left = np.zeros(len(pct))
for ct, col in zip(order, colors):
    ax.barh(y, pct[ct].values, left=left, color=col, height=0.78,
            edgecolor="white", linewidth=0.3)
    left += pct[ct].values

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.set_ylim(-0.6, len(pct) - 0.4)
ax.set_xlim(0, 100)
ax.set_xticks(range(0, 101, 20))
ax.tick_params(axis="x", labelsize=8.5, length=3, width=0.6)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("% of cells", fontsize=9.5)
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_linewidth(0.6)

handles = [mpatches.Patch(facecolor=c, label=l) for l, c in zip(order, colors)]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.42, -0.13), ncol=3,
          fontsize=8.5, title="Cell type", title_fontsize=9, frameon=False,
          handlelength=1.1, handleheight=0.9, columnspacing=1.0, handletextpad=0.4,
          labelspacing=0.35)

fig.subplots_adjust(left=0.2, right=0.98, top=0.98, bottom=0.30)
fig.savefig(os.path.join(HERE, "ch2_fig1.png"), dpi=400)
fig.savefig(os.path.join(HERE, "ch2_fig1.pdf"))

# --- spot check vs thesis text ---
lab_pct = pct.set_axis(labels)
print("\nSpot check (thesis text in brackets):")
for c, ref in [("Head and Neck", 85), ("Liver/Biliary", 50)]:
    print(f"  T cells in {c:14s} {lab_pct.loc[c, 'T cells']:5.1f}%  [~{ref}%]")
for c, ref in [("Prostate", 49), ("Colorectal", 40), ("Breast", 39), ("Pancreas", 35)]:
    print(f"  Epithelial/Cancer in {c:10s} {lab_pct.loc[c, 'Epithelial/Cancer']:5.1f}%  [~{ref}%]")
print("\nOrder top->bottom:", ", ".join(labels[::-1]))
