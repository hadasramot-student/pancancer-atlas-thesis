"""ch5_fig2 (thesis Figure 6.2): high-confidence (HC) cNMF-GeneNMF pairs per cancer type.

Data: per-cancer best Jaccard of every GeneNMF meta-program against its best-matching
same-cancer cNMF program (real_best_jaccard_recomputed.csv, written by
ch5_gep_programs/jaccard_threshold/null_derived_threshold.py) and the null-B 99th-percentile
threshold (0.0538) from null_derived_threshold_summary.csv. Same numbers as
thesis_figures/comment12_2026-10-21/hc_counts_by_cancer.csv (85 HC in total).

What one bar counts: the GeneNMF meta-programs of that cancer type (12-20 per cancer,
210 in total); the dark part are those whose best cNMF match has Jaccard > 0.0538.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from printstyle import save, WIDTH, plain

D = "/groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output/jaccard_threshold_sensitivity"
real = pd.read_csv(f"{D}/real_best_jaccard_recomputed.csv")
summ = pd.read_csv(f"{D}/null_derived_threshold_summary.csv")
THR = float(summ.query("null_method == 'B_cross_cancer' and percentile == 99.0").threshold.iloc[0])

ORDER = ["Brain", "Colorectal", "Kidney", "Liver_Biliary", "Lung", "Ovarian",
         "Pancreas", "Prostate", "Sarcoma", "Skin", "Breast"]          # order of the thesis figure
g = real.groupby("cancer").best_jaccard
n_mp = g.size().reindex(ORDER)
n_hc = g.apply(lambda x: int((x > THR).sum())).reindex(ORDER)
print("threshold", THR, "| HC total", n_hc.sum(), "| MPs total", n_mp.sum())
print(pd.DataFrame({"MPs": n_mp, "HC": n_hc}).to_string())

fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
fig.subplots_adjust(left=0.105, right=0.99, top=0.83, bottom=0.2)
x = range(len(ORDER))
ax.bar(x, n_hc, width=0.75, color="#2e86ab", edgecolor="black", linewidth=0.4,
       label=f"HC (best cNMF match Jaccard > {THR:g})")
ax.bar(x, n_mp - n_hc, bottom=n_hc, width=0.75, color="#bdd5ea", edgecolor="black", linewidth=0.4,
       label="Not HC")
for i, v in enumerate(n_hc):
    if v > 0:
        ax.text(i, v / 2 if v >= 3 else v + 0.3, str(v), ha="center", va="center" if v >= 3 else "bottom",
                fontsize=8, fontweight="bold", color="white" if v >= 3 else "black")
ax.set_xticks(list(x))
ax.set_xticklabels([plain(c) for c in ORDER], rotation=35, ha="right", rotation_mode="anchor")
ax.set_ylabel("Number of GeneNMF\nmeta-programs")
ax.set_ylim(0, 20.5)
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.grid(axis="y", color="#b0b0b0", linewidth=0.4, zorder=0); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False)
save(fig, "ch5_fig2")
