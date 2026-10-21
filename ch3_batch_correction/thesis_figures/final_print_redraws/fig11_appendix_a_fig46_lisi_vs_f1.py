"""appendix_a_fig46 (thesis Figure A.36): batch mixing (batch LISI) vs biological conservation (F1 label).

Same input as Figure 4.3 (all_validation_results.csv; 4 methods x 7 cell types = 28 points),
same colours (method) and shapes (cell type). The Pearson r printed here goes to the caption.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr
from printstyle import save, WIDTH, plain

df = pd.read_csv("/groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/all_validation_results.csv")
MLAB = {"mnn": "MNN", "combat": "ComBat", "harmony": "Harmony", "scanvi": "scANVI"}
MCOL = {"mnn": "#1f77b4", "combat": "#2ca02c", "harmony": "#ff7f0e", "scanvi": "#d62728"}
CT = ["B-cell", "Endothelial", "Fibroblasts", "Hematopoietic_Stem_and_Progenitor_cells", "Myeloid", "NK-cells", "T-cells"]
MARK = dict(zip(CT, ["o", "s", "^", "D", "P", "X", "v"]))

r, p = pearsonr(df.batch_lisi, df.f1_label)
print(f"Pearson r = {r:.3f}, p = {p:.3f}, n = {len(df)}")

fig, ax = plt.subplots(figsize=(WIDTH, 3.9))
fig.subplots_adjust(left=0.12, right=0.74, top=0.98, bottom=0.12)
for _, row in df.iterrows():
    ax.scatter(row.batch_lisi, row.f1_label, s=45, marker=MARK[row.cell_type], c=MCOL[row.method],
               edgecolors="black", linewidths=0.5, zorder=3)
ax.set_xlabel("Batch LISI (higher = more batch mixing)")
ax.set_ylabel("F1 label score\n(higher = better biological conservation)")
ax.grid(color="#c8c8c8", linewidth=0.4); ax.set_axisbelow(True)

h_m = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=MCOL[m], markeredgecolor="black",
              markersize=6, label=MLAB[m]) for m in MLAB]
h_c = [Line2D([0], [0], marker=MARK[c], linestyle="", markerfacecolor="#808080", markeredgecolor="black",
              markersize=6, label=plain(c)) for c in CT]
fig.legend(handles=h_m, title="Method (colour)", loc="upper left", bbox_to_anchor=(0.755, 0.98),
           frameon=False, alignment="left")
fig.legend(handles=h_c, title="Cell type (shape)", loc="upper left", bbox_to_anchor=(0.755, 0.6),
           frameon=False, alignment="left")
save(fig, "appendix_a_fig46")
