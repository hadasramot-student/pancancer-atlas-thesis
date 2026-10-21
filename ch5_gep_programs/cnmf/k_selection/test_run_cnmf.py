#!/usr/bin/env python
# test_run_cnmf.py
#
# Test run: single cancer type, downsampled, reduced K range
# Goal: verify pipeline end-to-end in ~1-2h before full submission
#
# Usage:
#   python test_run_cnmf.py --step prepare  --cancer Lung
#   python test_run_cnmf.py --step factorize --cancer Lung
#   python test_run_cnmf.py --step combine   --cancer Lung
#   python test_run_cnmf.py --step consensus --cancer Lung --k 20

import argparse
import os
import scipy.io
import scipy.sparse as sp
import numpy as np
import pandas as pd
import anndata as ad
from cnmf import cNMF

# ── Config ─────────────────────────────────────────────────────────────────
H5AD_PATHS = {
    "Lung"       : "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Lung/Epithelial_Lung_merged_full.h5ad",
    "Colorectal" : "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Colorectal/Epithelial_Colorectal_merged_full.h5ad",
    "Ovarian"    : "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Ovarian/Epithelial_Ovarian_merged_full.h5ad",
    "Pancreas"   : "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Pancreas/Epithelial_Pancreas_merged_full.h5ad",
}

OUTPUT_BASE = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/cnmf_test"
N_CELLS     = 15000   # downsample to this many cells
N_HVG       = 2000
N_ITER      = 100     # fewer replicates for test (still enough for stability)
SEED        = 42
K_VALUES    = list(range(10, 35, 5))  # [10, 15, 20, 25, 30]


def load_and_prepare(cancer, n_cells=N_CELLS):
    path = H5AD_PATHS[cancer]
    print(f"Loading {cancer} from: {path}")

    import anndata as ad
    adata = ad.read_h5ad(path)
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Use counts layer if available, else X
    if "counts" in adata.layers:
        mat = adata.layers["counts"]
    else:
        mat = adata.X

    # Back-transform if log-normalized
    if hasattr(mat, 'data'):
        sample = mat.data[:200]
    else:
        sample = np.array(mat[:10, :10]).flatten()
        sample = sample[sample > 0]

    if len(sample) > 0 and not np.all(sample == np.floor(sample)):
        print("Counts are log-normalized — back-transforming with expm1() + round()")
        mat = mat.copy()
        mat.data = np.round(np.expm1(mat.data)).astype(np.float32)
        mat.data[mat.data < 0] = 0
        print(f"  Sample after back-transform: {mat.data[:5]}")

    # Downsample
    if adata.shape[0] > n_cells:
        print(f"Downsampling {adata.shape[0]} → {n_cells} cells...")
        np.random.seed(SEED)
        idx = np.random.choice(adata.shape[0], n_cells, replace=False)
        idx = np.sort(idx)
        mat = mat[idx]
        obs = adata.obs.iloc[idx].copy()
        var = adata.var.copy()
    else:
        obs = adata.obs.copy()
        var = adata.var.copy()
        print(f"Using all {adata.shape[0]} cells (< {n_cells} threshold)")

    # Build clean AnnData with counts layer
    adata_sub = ad.AnnData(X=mat)
    adata_sub.obs_names = obs.index
    adata_sub.var_names = var.index
    adata_sub.layers["counts"] = mat.copy()

    # Copy useful metadata
    for col in ["cancer_type", "project_name", "orig.ident"]:
        if col in obs.columns:
            adata_sub.obs[col] = obs[col].values

    print(f"Final object: {adata_sub.shape[0]} cells x {adata_sub.shape[1]} genes")
    if "cancer_type" in adata_sub.obs.columns:
        print(f"Cancer type: {adata_sub.obs['cancer_type'].value_counts().to_dict()}")

    return adata_sub


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step",    required=True,
                        choices=["prepare", "factorize", "combine", "consensus"])
    parser.add_argument("--cancer",  default="Lung",
                        choices=list(H5AD_PATHS.keys()))
    parser.add_argument("--k",       type=int, default=None)
    parser.add_argument("--worker",  type=int, default=0)
    parser.add_argument("--total_workers", type=int, default=1)
    parser.add_argument("--density_threshold", type=float, default=0.1)
    args = parser.parse_args()

    run_name   = f"test_{args.cancer}"
    output_dir = os.path.join(OUTPUT_BASE, args.cancer)
    os.makedirs(output_dir, exist_ok=True)

    cnmf_obj = cNMF(output_dir=output_dir, name=run_name)

    if args.step == "prepare":
        print(f"\n{'='*50}\nTEST RUN: PREPARE\nCancer: {args.cancer}\nCells: {N_CELLS}\nK: {K_VALUES}\nReplicates: {N_ITER}\n{'='*50}")
        adata = load_and_prepare(args.cancer)

        h5ad_path = os.path.join(output_dir, f"{run_name}.h5ad")
        adata.write_h5ad(h5ad_path)
        print(f"Saved to: {h5ad_path}")

        cnmf_obj.prepare(
            counts_fn         = h5ad_path,
            components        = K_VALUES,
            n_iter            = N_ITER,
            seed              = SEED,
            num_highvar_genes = N_HVG,
        )
        print("Prepare complete.")

    elif args.step == "factorize":
        print(f"\n{'='*50}\nTEST RUN: FACTORIZE\nWorker: {args.worker}/{args.total_workers}\n{'='*50}")
        cnmf_obj.factorize(
            worker_i      = args.worker,
            total_workers = args.total_workers,
        )
        print("Factorize complete.")

    elif args.step == "combine":
        print(f"\n{'='*50}\nTEST RUN: COMBINE\n{'='*50}")
        cnmf_obj.combine()
        print("Combine complete.")
        print(f"\nNow inspect stability plot and choose K, then run:")
        print(f"  python test_run_cnmf.py --step consensus --cancer {args.cancer} --k <K>")

    elif args.step == "consensus":
        if args.k is None:
            raise ValueError("--k required for consensus")
        print(f"\n{'='*50}\nTEST RUN: CONSENSUS  K={args.k}\n{'='*50}")
        cnmf_obj.consensus(
            k                 = args.k,
            density_threshold = args.density_threshold,
            show_clustering   = False,
        )
        print(f"Consensus complete. Output: {output_dir}/{run_name}/")


if __name__ == "__main__":
    main()
