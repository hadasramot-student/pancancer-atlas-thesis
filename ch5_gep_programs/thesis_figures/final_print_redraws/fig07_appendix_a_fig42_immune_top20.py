"""appendix_a_fig42 (thesis Figure A.32): top 20 GeneNMF-Spectra program pairs (immune compartment)."""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from printstyle import save, WIDTH
from immune_common import best, LINEAGES, LABEL, LCOL, HC_THR

# 7 pairs tie at J = 0.1236 for the last 5 places; keep the 5 shown in the original figure
# (T cells MP15, B cells MP14, HSPCs MP9, T cells MP13, T cells MP2; T cells MP4 and MP8 not shown).
TIE_ORDER = [("T-cells", "MP15"), ("B-cell", "MP14"), ("Hematopoietic_Stem_and_Progenitor_cells", "MP9"),
             ("T-cells", "MP13"), ("T-cells", "MP2")]
srt = best.sort_values("best_jaccard", ascending=False, kind="stable")
above = srt[srt.best_jaccard > 0.1236]
ties = srt.set_index(["lineage", "meta_program"]).loc[TIE_ORDER].reset_index()
top = pd.concat([above, ties]).head(20).iloc[::-1]
print(top[["lineage", "meta_program", "best_spectra_program", "best_jaccard"]].iloc[::-1].to_string(index=False))

fig, ax = plt.subplots(figsize=(WIDTH, 4.6))
fig.subplots_adjust(left=0.19, right=0.98, top=0.99, bottom=0.1)
y = range(len(top))
ax.barh(y, top.best_jaccard, color=[LCOL[l] for l in top.lineage], height=0.75, edgecolor="black", linewidth=0.3)
ax.set_yticks(list(y))
ax.set_yticklabels([f"{LABEL[l]} {m}" for l, m in zip(top.lineage, top.meta_program)])
for i, (v, s) in enumerate(zip(top.best_jaccard, top.best_spectra_program)):
    ax.text(v + 0.004, i, f"{v:.3f} ({s.replace('Spectra_', 'Spectra ')})", va="center", fontsize=8)
ax.axvline(HC_THR, color="black", linestyle="--", linewidth=0.7)
ax.set_xlim(0, 0.47); ax.set_ylim(-0.6, len(top) - 0.4)
ax.set_xlabel("Best-match Jaccard similarity (GeneNMF meta-program vs. Spectra)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(handles=[Patch(facecolor=LCOL[l], edgecolor="black", linewidth=0.3, label=LABEL[l]) for l in LINEAGES],
          title="Lineage", loc="lower right", frameon=False)
save(fig, "appendix_a_fig42")
