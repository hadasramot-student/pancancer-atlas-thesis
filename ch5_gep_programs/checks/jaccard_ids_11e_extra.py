import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

OUT = "/home/hadasra/.tmp/claude-343831894/-groups-ofircohen-group-users-hadasra-create-sc-object-workflow--claude-worktrees-jaccard-ids-hallmark-stats-acdc92/86cb1c60-51b8-4654-a3d6-7d7718d1bab3/scratchpad"
df = pd.read_csv(f"{OUT}/11e_all_cnmf_genenmf_pairs.csv")
d2 = df[~df.cancer.isin(["Sarcoma"])]          # drop the null cancer

for name, g in [("all 11 cancers", df), ("excl Sarcoma", d2)]:
    print(f"\n--- {name} (N={len(g)}) ---")
    z = (g.jaccard == 0).sum()
    print(f"  pairs with Jaccard exactly 0            : {z} ({z/len(g)*100:.0f}%)")
    print(f"    ...their IDS: median={g[g.jaccard==0].ids.median():.2f}  "
          f"p90={g[g.jaccard==0].ids.quantile(.9):.2f}  max={g[g.jaccard==0].ids.max():.2f}")
    hi = g[g.ids >= 0.5]
    print(f"  pairs with IDS >= 0.5                   : {len(hi)}")
    print(f"    ...of those, Jaccard < 0.05           : {(hi.jaccard<0.05).sum()} ({(hi.jaccard<0.05).mean()*100:.0f}%)")
    print(f"    ...of those, Jaccard = 0              : {(hi.jaccard==0).sum()} ({(hi.jaccard==0).mean()*100:.0f}%)")
    hij = g[g.jaccard > 0.1]
    print(f"  pairs with Jaccard > 0.1 (HC-ish)       : {len(hij)}  -> IDS median={hij.ids.median():.2f}, min={hij.ids.min():.2f}")

# binned figure
bins = [-1e-9, 1e-9, 0.02, 0.05, 0.1, 0.2, 1]
labs = ["J=0", "0-0.02", "0.02-0.05", "0.05-0.1", "0.1-0.2", ">0.2"]
d2 = d2.copy()
d2["band"] = pd.cut(d2.jaccard, bins, labels=labs)
fig, ax = plt.subplots(figsize=(8, 5))
data = [d2[d2.band == l].ids.values for l in labs]
bp = ax.boxplot(data, labels=[f"{l}\nn={len(x)}" for l, x in zip(labs, data)],
                showmeans=True, widths=.6)
ax.set_xlabel("pairwise Jaccard band (top-50 gene lists)")
ax.set_ylabel("pairwise IDS")
ax.set_title("IDS distribution by Jaccard band (cNMF-GeneNMF pairs, excl. Sarcoma)\n"
             f"Pearson r={pearsonr(d2.jaccard,d2.ids)[0]:.3f}  Spearman rho={spearmanr(d2.jaccard,d2.ids)[0]:.3f}")
ax.grid(alpha=.25, axis="y")
fig.tight_layout(); fig.savefig(f"{OUT}/11e_ids_by_jaccard_band.png", dpi=140)
print(f"\nwrote {OUT}/11e_ids_by_jaccard_band.png")
