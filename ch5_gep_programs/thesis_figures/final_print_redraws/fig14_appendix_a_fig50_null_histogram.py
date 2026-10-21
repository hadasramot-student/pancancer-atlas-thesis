"""appendix_a_fig50 (thesis Figure A.40): real cNMF-GeneNMF best-Jaccard distribution vs null-derived thresholds.

Data written by ch5_gep_programs/jaccard_threshold/null_derived_threshold.py:
null_A_permutation_values.npy (20,000 draws), null_B_crosscancer_values.npy (42,000 pairs),
real_best_jaccard_recomputed.csv (210 GeneNMF meta-programs), null_derived_threshold_summary.csv.
Density histograms on a log y axis; legend above the plot so it covers no bars.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from printstyle import save, WIDTH

D = "/groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output/jaccard_threshold_sensitivity"
nA = np.load(f"{D}/null_A_permutation_values.npy")
nB = np.load(f"{D}/null_B_crosscancer_values.npy")
real = pd.read_csv(f"{D}/real_best_jaccard_recomputed.csv").best_jaccard.values
summ = pd.read_csv(f"{D}/null_derived_threshold_summary.csv")
t99 = summ.query("null_method == 'B_cross_cancer' and percentile == 99.0").threshold.iloc[0]
t999 = summ.query("null_method == 'B_cross_cancer' and percentile == 99.9").threshold.iloc[0]
print(len(nA), len(nB), len(real), t99, t999)

bins = np.linspace(0, 0.55, 61)
fig, ax = plt.subplots(figsize=(WIDTH, 4.4))
fig.subplots_adjust(left=0.11, right=0.98, top=0.72, bottom=0.11)
ax.hist(nA, bins=bins, density=True, color="#8c9bb0", alpha=0.6, label=f"Null A: gene-frequency permutation (n = {len(nA):,})")
ax.hist(nB, bins=bins, density=True, color="#f0a050", alpha=0.6, label=f"Null B: cross-cancer-type real pairs (n = {len(nB):,})")
ax.hist(real, bins=bins, density=True, color="#4a90d9", alpha=0.6, label=f"Observed: same-cancer best match (n = {len(real)})")
ax.axvline(0.1, color="black", linewidth=1.6, label="Fixed threshold (0.1)")
ax.axvline(t99, color="#b5521b", linestyle="--", linewidth=1.3, label=f"Null B 99th percentile ({t99:.4f})")
ax.axvline(t999, color="#b5521b", linestyle=":", linewidth=1.5, label=f"Null B 99.9th percentile ({t999:.4f})")
ax.set_yscale("log")
ax.set_xlim(0, 0.55)
ax.set_xlabel("Jaccard similarity (top-50 genes)")
ax.set_ylabel("Density (log scale)")
ax.legend(loc="lower left", bbox_to_anchor=(-0.1, 1.02), ncol=1, frameon=False, borderaxespad=0)
save(fig, "appendix_a_fig50")
