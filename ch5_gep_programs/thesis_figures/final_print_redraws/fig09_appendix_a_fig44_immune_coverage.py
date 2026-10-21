"""appendix_a_fig44 (thesis Figure A.34): GeneNMF meta-program sample coverage, immune compartment.

Colour = sampleCoverage from <lineage>_genenmf_MP_metrics.csv (fraction of samples in which the
meta-program is active); star = HC in that lineage (best Spectra Jaccard > 0.1); red row label =
HC in at least one lineage (as in the original figure).
"""
import matplotlib.pyplot as plt
from printstyle import save, WIDTH
from immune_common import best, coverage, LINEAGES, LABEL, MPS, HC_RED

cov = coverage()
hc = {(l, m) for l, m in zip(best.lineage, best.meta_program) if (best.hc[(best.lineage == l) & (best.meta_program == m)].iloc[0])}
hc_any = {m for _, m in hc}
print("HC cells", len(hc), "| rows HC in >=1 lineage", len(hc_any))

fig, ax = plt.subplots(figsize=(WIDTH, 4.6))
fig.subplots_adjust(left=0.15, right=0.86, top=0.985, bottom=0.11)
im = ax.imshow(cov.values, cmap="Blues", vmin=0, vmax=1, aspect="auto", interpolation="nearest")
for j, lin in enumerate(LINEAGES):
    for i, mp in enumerate(MPS):
        if (lin, mp) in hc:
            ax.plot(j, i, marker="*", markersize=9, markerfacecolor="white", markeredgecolor="black",
                    markeredgewidth=0.6)
ax.set_xticks(range(len(LINEAGES))); ax.set_xticklabels([LABEL[l] for l in LINEAGES])
ax.set_yticks(range(len(MPS))); ax.set_yticklabels(MPS)
for t, mp in zip(ax.get_yticklabels(), MPS):
    if mp in hc_any:
        t.set_color(HC_RED); t.set_fontweight("bold")
ax.set_xlabel("Immune lineage"); ax.set_ylabel("GeneNMF meta-program")
cax = fig.add_axes([0.885, 0.11, 0.025, 0.875])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Sample coverage (fraction of samples)")
save(fig, "appendix_a_fig44")
