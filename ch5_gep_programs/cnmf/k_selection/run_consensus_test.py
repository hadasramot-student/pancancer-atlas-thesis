#!/usr/bin/env python
# run_consensus_test.py
#
# Runs consensus for all K values and prints stability scores,
# so you can pick the optimal K for the test run.

import os
import numpy as np
import pandas as pd
from cnmf import cNMF

CANCER     = "Lung"
OUTPUT_DIR = f"/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_test/{CANCER}"
RUN_NAME   = f"test_{CANCER}"
K_VALUES   = [10, 15, 20, 25, 30]
DENSITY_THRESHOLD = 0.1

cnmf_obj = cNMF(output_dir=OUTPUT_DIR, name=RUN_NAME)

print(f"Running consensus for K = {K_VALUES}")
print(f"density_threshold = {DENSITY_THRESHOLD}\n")

stability_scores = {}

for k in K_VALUES:
    print(f"{'='*40}")
    print(f"K = {k}")
    try:
        cnmf_obj.consensus(
            k=k,
            density_threshold=DENSITY_THRESHOLD,
            show_clustering=False,
        )
        # Read stability score from the stats file
        stats_file = os.path.join(OUTPUT_DIR, RUN_NAME,
                                  f"{RUN_NAME}.k_{k}.dt_{str(DENSITY_THRESHOLD).replace('.','_')}.consensus_stats.df.npz")
        if os.path.exists(stats_file):
            stats = pd.DataFrame(**np.load(stats_file, allow_pickle=True))
            print(stats)
            if "stability" in stats.columns:
                stability_scores[k] = stats["stability"].mean()
        print(f"K={k} done.\n")
    except Exception as e:
        print(f"K={k} failed: {e}\n")

print("\n" + "="*40)
print("STABILITY SUMMARY:")
for k, s in stability_scores.items():
    print(f"  K={k:3d}  stability={s:.4f}")

if stability_scores:
    best_k = max(stability_scores, key=stability_scores.get)
    print(f"\nBest K by stability: {best_k}")
