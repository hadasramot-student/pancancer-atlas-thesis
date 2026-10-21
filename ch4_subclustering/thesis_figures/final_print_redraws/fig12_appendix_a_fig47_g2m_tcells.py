"""appendix_a_fig47 (thesis Figure A.37): best G2M-checkpoint-enriched T cell cluster per cancer type.

Data: clustered_lineages/FindMarkers_results/enrichment/T-cells_<cancer>_hallmarks_gsea.csv
(fgsea Hallmark results per cluster). Per cancer type the cluster with the highest NES for
HALLMARK_G2M_CHECKPOINT is shown; dot size = -log10(padj); colour = NES.
Colour scale starts at a mid-red instead of white so that low-NES dots (Brain, Sarcoma) print;
size legend sits outside the data area.
"""
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from printstyle import save, WIDTH, plain

D = "/groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages/FindMarkers_results/enrichment"
rows = []
for f in sorted(glob.glob(f"{D}/T-cells_*_hallmarks_gsea.csv")):
    d = pd.read_csv(f)
    g = d[d.pathway == "HALLMARK_G2M_CHECKPOINT"].sort_values("NES", ascending=False).iloc[0]
    rows.append(dict(cancer=f.split("T-cells_")[1].replace("_hallmarks_gsea.csv", ""), cluster=g.cluster,
                     NES=g.NES, padj=g.padj, n_clusters=d.cluster.nunique()))
res = pd.DataFrame(rows).sort_values("NES").reset_index(drop=True)
print(res.to_string())

def size(padj):
    return 12 + 6 * -np.log10(padj)

cmap = LinearSegmentedColormap.from_list("reds_print", plt.get_cmap("Reds")(np.linspace(0.35, 1.0, 256)))
fig, ax = plt.subplots(figsize=(WIDTH, 3.6))
fig.subplots_adjust(left=0.16, right=0.66, top=0.97, bottom=0.17)
sc = ax.scatter(res.NES, range(len(res)), s=size(res.padj), c=res.NES, cmap=cmap, vmin=1.8, vmax=2.7,
                edgecolors="black", linewidths=0.6, zorder=3)
ax.set_yticks(range(len(res))); ax.set_yticklabels([plain(c) for c in res.cancer])
ax.set_ylim(-0.7, len(res) - 0.3)
ax.set_xlim(-0.12, 2.75); ax.axvline(0, color="black", linewidth=0.6)
ax.set_xlabel("Normalized enrichment score (NES),\nHALLMARK_G2M_CHECKPOINT")
ax.grid(axis="x", color="#c8c8c8", linewidth=0.4); ax.set_axisbelow(True)

cax = fig.add_axes([0.69, 0.17, 0.025, 0.8])
cb = fig.colorbar(sc, cax=cax); cb.set_label("NES")
h = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor="#808080", markeredgecolor="black",
            markersize=size(p) ** 0.5, label=lab) for p, lab in [(1e-8, "10$^{-8}$"), (1e-20, "10$^{-20}$"), (1e-35, "10$^{-35}$")]]
fig.legend(handles=h, title="Adjusted p", loc="upper left", bbox_to_anchor=(0.80, 0.97), frameon=False,
           labelspacing=1.6, borderpad=0.6, alignment="left")
save(fig, "appendix_a_fig47")
