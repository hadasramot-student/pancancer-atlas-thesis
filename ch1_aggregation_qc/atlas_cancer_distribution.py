# atlas_cancer_distribution.py
#
# Ch.1 (Aggregation & QC) - atlas composition figure: cell counts per cancer
# type across the merged pan-cancer atlas, as a horizontal bar chart (panel
# A) + pie chart (panel B).
#
# Inputs: none (cell counts are the published per-cancer-type totals for
# this atlas; update the `data` list below if the atlas composition changes,
# or wire it up to Table3_metatype_composition.csv for a data-driven version).
# Outputs: <ATLAS_FIGURES_DIR>/Fig_cancer_distribution.{png,pdf}
#
# Usage:
#   ATLAS_FIGURES_DIR=/path/to/plots_and_outputs python atlas_cancer_distribution.py

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

data = [
    ("Breast", 2931922),
    ("Lung", 1058285),
    ("Colorectal", 982477),
    ("Unknown", 562518),
    ("Pancreas", 384677),
    ("Ovarian", 382551),
    ("Kidney", 371604),
    ("Skin", 369212),
    ("Other Models", 311451),
    ("Prostate", 210716),
    ("Liver/Biliary", 188218),
    ("Brain", 175383),
    ("Neuroendocrine", 129559),
    ("Sarcoma", 80682),
    ("Head & Neck", 71896),
]

df = pd.DataFrame(data, columns=["cancer_type", "n"])
total = df["n"].sum()
df["pct"] = df["n"] / total * 100
df = df.sort_values("n", ascending=True)

colors = [
    "#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
    "#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85",
    "#FFDC91","#80C680","#A6CEE3","#CAB2D6","#FDBF6F"
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7),
                                gridspec_kw={"width_ratios": [1.6, 1]})

# --- Panel A: horizontal bar chart ---
bars = ax1.barh(df["cancer_type"], df["n"] / 1e6, color=colors[:len(df)], height=0.7, edgecolor="white")
for bar, (_, row) in zip(bars, df.iterrows()):
    ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
             f"{row['n']/1e6:.2f}M ({row['pct']:.1f}%)",
             va="center", fontsize=9)
ax1.set_xlabel("Number of cells (millions)", fontsize=11)
ax1.set_title("A", fontsize=14, fontweight="bold", loc="left")
ax1.set_xlim(0, df["n"].max()/1e6 * 1.35)
ax1.spines[["top","right"]].set_visible(False)
ax1.tick_params(axis="y", labelsize=10)

# --- Panel B: pie chart ---
df_pie = df.sort_values("n", ascending=False)
wedge_colors = colors[:len(df_pie)]
wedges, texts, autotexts = ax2.pie(
    df_pie["n"],
    labels=None,
    colors=wedge_colors,
    autopct=lambda p: f"{p:.1f}%" if p > 2 else "",
    startangle=140,
    pctdistance=0.75,
    wedgeprops={"edgecolor": "white", "linewidth": 0.8}
)
for at in autotexts:
    at.set_fontsize(8)
ax2.set_title("B", fontsize=14, fontweight="bold", loc="left")
ax2.legend(wedges, df_pie["cancer_type"], loc="lower left",
           bbox_to_anchor=(-0.15, -0.15), fontsize=8.5, frameon=False, ncol=2)

fig.suptitle(f"Distribution of {total/1e6:.2f}M cells across cancer types",
             fontsize=13, fontweight="bold", y=1.01)

plt.tight_layout()
out_dir = os.environ.get("ATLAS_FIGURES_DIR", "outputs/figures")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "Fig_cancer_distribution.png")
plt.savefig(out, dpi=300, bbox_inches="tight")
plt.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("Saved:", out)
