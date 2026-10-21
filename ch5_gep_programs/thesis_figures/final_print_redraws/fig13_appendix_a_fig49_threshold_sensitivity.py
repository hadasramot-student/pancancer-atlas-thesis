"""appendix_a_fig49 (thesis Figure A.39): Jaccard-threshold sensitivity, 0.05-0.5, panels a and b.

a: total HC cNMF-GeneNMF programs per threshold (jaccard_threshold_sensitivity_full.csv,
   written by ch5_gep_programs/jaccard_threshold/extend_sensitivity.py).
b: number of cMPs and of pan-cancer cMPs (>= 7 cancer types) per threshold from clustering the
   240 cNMF programs (ch5_gep_programs/meta_programs/cmp_threshold_sensitivity.csv).
Same values, colours and markers as the original; dashed line = 0.1.
"""
import pandas as pd
import matplotlib.pyplot as plt
from printstyle import save, WIDTH

a = pd.read_csv("/groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output/jaccard_threshold_sensitivity/jaccard_threshold_sensitivity_full.csv")
b = pd.read_csv("/groups/ofircohen-group/users/hadasra/pancancer-atlas-thesis/ch5_gep_programs/meta_programs/cmp_threshold_sensitivity.csv")
print(a[["threshold", "Total"]].to_string(index=False)); print(b.iloc[:, :3].to_string(index=False))

fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.9))
fig.subplots_adjust(left=0.135, right=0.99, top=0.91, bottom=0.15, wspace=0.3)
ax = axs[0]
ax.plot(a.threshold, a.Total, "-o", color="#2b6cb0", markersize=4, linewidth=1.3)
for x, y in zip(a.threshold, a.Total):
    ax.annotate(str(y), (x, y), textcoords="offset points", xytext=(4, 4), fontsize=8)
ax.set_ylabel("Total HC programs\n(cNMF ∩ GeneNMF)")
ax.set_ylim(-5, 112)

ax = axs[1]
ax.plot(b.threshold, b.n_cMPs, "-o", color="#b5521b", markersize=4, linewidth=1.3, label="Total cMPs")
ax.plot(b.threshold, b.iloc[:, 2], "--s", color="#2e8b57", markersize=4, linewidth=1.3,
        label="Pan-cancer cMPs\n(≥ 7 cancer types)")
for k, (x, y) in enumerate(zip(b.threshold, b.n_cMPs)):
    # first label to the left so 23 and 24 (0.05, 0.075) do not overlap
    ax.annotate(str(y), (x, y), textcoords="offset points", xytext=(-4, -3) if k == 0 else (4, 4),
                ha="right" if k == 0 else "left", va="top" if k == 0 else "baseline", fontsize=8)
ax.set_ylabel("Count"); ax.set_ylim(-1.5, 28)
ax.legend(loc="upper right", frameon=False)

for ax, lab in zip(axs, "ab"):
    ax.axvline(0.1, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_xlim(0, 0.56); ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xlabel("Jaccard threshold")
    ax.grid(color="#c8c8c8", linewidth=0.4); ax.set_axisbelow(True)
    ax.set_title(lab, loc="left", fontsize=9, fontweight="bold")
save(fig, "appendix_a_fig49")
