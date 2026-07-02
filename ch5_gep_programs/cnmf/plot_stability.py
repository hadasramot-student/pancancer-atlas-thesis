#!/usr/bin/env python
# plot_stability.py
#
# Ch.5 (GEP programs) — generates the cNMF K-sweep stability + reconstruction
# error plot (used to pick the final K=20) from cNMF combine output.
# Usage: python plot_stability.py --cancer Lung --cnmf_dir <path>

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--cancer",     default="Lung")
parser.add_argument("--run_name",   default=None)
parser.add_argument("--cnmf_dir",   default=os.path.join(os.environ.get("ATLAS_DATA_DIR", "data/temp_objects"), "cnmf_test"))
parser.add_argument("--density_threshold", type=float, default=0.1)
args = parser.parse_args()

run_name  = args.run_name or f"test_{args.cancer}"
base_dir  = os.path.join(args.cnmf_dir, args.cancer, run_name)
tmp_dir   = os.path.join(base_dir, "cnmf_tmp")
out_plot  = os.path.join(base_dir, f"{run_name}.stability_plot.pdf")
out_txt   = os.path.join(base_dir, f"{run_name}.stability_scores.txt")

K_VALUES = sorted([
    int(f.split(".k_")[1].split(".")[0])
    for f in os.listdir(tmp_dir)
    if "local_density_cache" in f and "merged" in f
])
print(f"Found K values: {K_VALUES}")

stability   = {}
recon_error = {}

for k in K_VALUES:
    # Stability: fraction of NMF runs within density threshold
    density_f = os.path.join(tmp_dir, f"{run_name}.local_density_cache.k_{k}.merged.df.npz")
    d = np.load(density_f, allow_pickle=True)
    df = pd.DataFrame(**d)
    stability[k] = (df.values < args.density_threshold).mean()

    # Reconstruction error: from the merged spectra (NMF objective value)
    # cNMF stores the NMF error in the params file
    params_f = os.path.join(tmp_dir, f"{run_name}.nmf_params.df.npz")
    if os.path.exists(params_f):
        params = pd.DataFrame(**np.load(params_f, allow_pickle=True))
        if "err" in params.columns:
            recon_error[k] = params.loc[params["K"] == k, "err"].mean() if "K" in params.columns else np.nan

print("\nStability scores:")
for k in K_VALUES:
    print(f"  K={k:3d}  stability={stability[k]:.4f}")

# Save scores as text
with open(out_txt, "w") as f:
    f.write("K\tstability\n")
    for k in K_VALUES:
        f.write(f"{k}\t{stability[k]:.4f}\n")
print(f"\nSaved scores to: {out_txt}")

# ── Plot ───────────────────────────────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=(7, 4))

ks   = K_VALUES
stab = [stability[k] for k in ks]

ax1.plot(ks, stab, "o-", color="#2166ac", linewidth=2, markersize=8, label="Stability")
ax1.set_xlabel("Number of programs (K)", fontsize=13)
ax1.set_ylabel("Stability", color="#2166ac", fontsize=13)
ax1.tick_params(axis="y", labelcolor="#2166ac")
ax1.set_xticks(ks)
ax1.set_ylim(0, 1.05)

for k, s in zip(ks, stab):
    ax1.annotate(f"{s:.3f}", (k, s), textcoords="offset points",
                 xytext=(0, 8), ha="center", fontsize=9, color="#2166ac")

# Mark best K
best_k = max(stability, key=stability.get)
ax1.axvline(x=best_k, color="gray", linestyle="--", alpha=0.5, label=f"Best K={best_k}")

# Reconstruction error if available
if recon_error and any(not np.isnan(v) for v in recon_error.values()):
    ax2 = ax1.twinx()
    errs = [recon_error.get(k, np.nan) for k in ks]
    ax2.plot(ks, errs, "s--", color="#d6604d", linewidth=1.5, markersize=6, label="Recon. error")
    ax2.set_ylabel("Reconstruction error", color="#d6604d", fontsize=13)
    ax2.tick_params(axis="y", labelcolor="#d6604d")
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(*zip(*[(h, l) for h, l in
                      zip(*ax1.get_legend_handles_labels())] +
                      [(h, l) for h, l in zip(lines2, labels2)]),
               loc="lower left", fontsize=9)
else:
    ax1.legend(loc="lower left", fontsize=9)

plt.title(f"cNMF Stability — {args.cancer} (density threshold={args.density_threshold})",
          fontsize=13, pad=10)
plt.tight_layout()
plt.savefig(out_plot, dpi=150, bbox_inches="tight")
print(f"Saved plot to: {out_plot}")
