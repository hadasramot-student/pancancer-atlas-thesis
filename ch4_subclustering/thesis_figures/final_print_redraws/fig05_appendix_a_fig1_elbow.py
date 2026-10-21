"""appendix_a_fig1 (thesis Figure A.2): elbow plot, T cells, Liver/Biliary.

Data: PCA standard deviations of the first 30 PCs, exported read-only from
T-cells_Liver_Biliary_clustered.RDS by fig05_extract_pca_stdev.R (same object and reduction as
ch4_subclustering/elbow_check.R). Elbow rule as in elbow_check.R: PC after the largest drop.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from printstyle import save, WIDTH, OUT

d = pd.read_csv(f"{OUT}/data/Tcells_LiverBiliary_pca_stdev.csv")
elbow = int(np.argmax(-np.diff(d.stdev.values))) + 2
print("PC1 sd %.3f, PC30 sd %.3f, elbow PC%d" % (d.stdev.iloc[0], d.stdev.iloc[-1], elbow))

fig, ax = plt.subplots(figsize=(WIDTH, 2.8))
fig.subplots_adjust(left=0.1, right=0.98, top=0.96, bottom=0.16)
ax.plot(d.PC, d.stdev, color="black", linewidth=0.8, zorder=2)
ax.scatter(d.PC, d.stdev, s=14, color="black", zorder=3)
ax.axvline(elbow, color="firebrick", linestyle="--", linewidth=0.9)
ax.text(elbow + 0.5, d.stdev.max() * 0.95, f"elbow ~PC{elbow}", color="firebrick", fontsize=8, ha="left")
ax.set_xlabel("Principal component"); ax.set_ylabel("Standard deviation")
ax.set_xticks([1, 5, 10, 15, 20, 25, 30])
ax.grid(color="#c8c8c8", linewidth=0.4); ax.set_axisbelow(True)
save(fig, "appendix_a_fig1")
