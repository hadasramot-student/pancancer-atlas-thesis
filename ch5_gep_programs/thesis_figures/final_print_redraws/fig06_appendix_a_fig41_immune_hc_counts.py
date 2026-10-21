"""appendix_a_fig41 (thesis Figure A.31): number of HC GeneNMF immune programs per lineage."""
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from printstyle import save, WIDTH
from immune_common import best, LINEAGES, LABEL, HC_THR

n_all = best.groupby("lineage").size().reindex(LINEAGES)
n_hc = best.groupby("lineage").hc.sum().reindex(LINEAGES).astype(int)
print(dict(zip(LINEAGES, zip(n_hc, n_all))), "total HC", n_hc.sum())

fig, ax = plt.subplots(figsize=(WIDTH, 2.7))
fig.subplots_adjust(left=0.11, right=0.99, top=0.86, bottom=0.17)
x = range(len(LINEAGES))
ax.bar(x, n_hc, width=0.7, color="#2171b5", edgecolor="black", linewidth=0.4,
       label=f"High-confidence (Jaccard > {HC_THR:g})")
ax.bar(x, n_all - n_hc, bottom=n_hc, width=0.7, color="#c6dbef", edgecolor="black", linewidth=0.4,
       label="Not recovered by Spectra")
for i, (h, a) in enumerate(zip(n_hc, n_all)):
    ax.text(i, h + 0.3, f"{h}/{a}", ha="center", va="bottom", fontsize=8, fontweight="bold")
ax.set_xticks(list(x)); ax.set_xticklabels([LABEL[l] for l in LINEAGES])
ax.set_xlabel("Immune lineage"); ax.set_ylabel("Number of GeneNMF\nmeta-programs")
ax.set_ylim(0, 20.5); ax.yaxis.set_major_locator(MultipleLocator(5))
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False)
save(fig, "appendix_a_fig41")
