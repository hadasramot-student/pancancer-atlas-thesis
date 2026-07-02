#!/usr/bin/env python
# run_cnmf_per_cancer.py
#
# Ch.5 (GEP programs) — alternate cNMF entrypoint: runs cNMF for a single
# cancer type at a time (vs. 02_run_cnmf.py's atlas-wide run), used for both
# a downsampled K-sweep test run and the full-cells fixed-K production run.
#
# Usage — test:
#   python run_cnmf_per_cancer.py --step prepare --cancer Lung
#   python run_cnmf_per_cancer.py --step factorize --cancer Lung --worker 0 --total_workers 2
#   python run_cnmf_per_cancer.py --step combine --cancer Lung
#   python run_cnmf_per_cancer.py --step consensus --cancer Lung --k 20
#
# Usage — full run (all cells, fixed K, more replicates):
#   python run_cnmf_per_cancer.py --step prepare --cancer Lung \
#       --n_cells 0 --n_iter 200 --k_values "20"
#   python run_cnmf_per_cancer.py --step factorize --cancer Lung \
#       --worker 0 --total_workers 4
#   python run_cnmf_per_cancer.py --step combine --cancer Lung
#   python run_cnmf_per_cancer.py --step consensus --cancer Lung --k 20

import argparse
import os
import numpy as np
import pandas as pd
import anndata as ad
from cnmf import cNMF

ATLAS_DATA_DIR = os.environ.get("ATLAS_DATA_DIR", "data/temp_objects")
_MERGED_BASE = os.path.join(ATLAS_DATA_DIR, "merged_by_query")

# Cancer -> Epithelial/tumor-compartment h5ad (see ../filter_brain_malignant.R
# for how Brain's is derived; the rest are the Ch.4 lineage-split objects).
H5AD_PATHS = {
    cancer: os.path.join(_MERGED_BASE, cancer, f"Epithelial_{cancer}_merged_full.h5ad")
    for cancer in [
        "Lung", "Colorectal", "Kidney", "Liver_Biliary", "Ovarian", "Pancreas",
        "Prostate", "Sarcoma", "Skin", "Brain", "Breast",
    ]
}

N_HVG = 2000
SEED  = 42


def load_counts(cancer, n_cells):
    path = H5AD_PATHS[cancer]
    print(f"Loading {cancer}: {path}")
    adata = ad.read_h5ad(path)
    print(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    mat = adata.layers["counts"] if "counts" in adata.layers else adata.X

    # Back-transform log-normalized counts if needed
    sample = mat.data[:200] if hasattr(mat, "data") else np.array(mat[:5, :5]).flatten()
    sample = sample[sample > 0]
    if len(sample) > 0 and not np.all(sample == np.floor(sample)):
        print("  Back-transforming log-normalized counts with expm1() + round()")
        mat = mat.copy()
        mat.data = np.round(np.expm1(mat.data)).astype(np.float32)
        mat.data[mat.data < 0] = 0

    # Downsample
    if n_cells > 0 and adata.shape[0] > n_cells:
        print(f"  Downsampling {adata.shape[0]} → {n_cells} cells")
        np.random.seed(SEED)
        idx = np.sort(np.random.choice(adata.shape[0], n_cells, replace=False))
        mat = mat[idx]
        obs = adata.obs.iloc[idx].copy()
    else:
        obs = adata.obs.copy()
        print(f"  Using all {adata.shape[0]} cells")

    out = ad.AnnData(X=mat)
    out.obs_names = obs.index
    out.var_names = adata.var_names
    out.layers["counts"] = mat.copy()
    for col in obs.columns:
        try:
            out.obs[col] = obs[col].values
        except Exception:
            pass

    print(f"  Final: {out.shape[0]} cells x {out.shape[1]} genes")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step",     required=True,
                        choices=["prepare", "factorize", "combine", "consensus"])
    parser.add_argument("--cancer",   required=True, choices=list(H5AD_PATHS.keys()))
    parser.add_argument("--output_base", default=os.path.join(ATLAS_DATA_DIR, "cnmf_output"))
    parser.add_argument("--k",        type=int, default=None)
    parser.add_argument("--k_values", type=str, default="10,15,20,25,30",
                        help="Comma-separated K values for prepare step")
    parser.add_argument("--n_iter",   type=int, default=200)
    parser.add_argument("--n_cells",  type=int, default=0,
                        help="Downsample to N cells (0 = use all)")
    parser.add_argument("--worker",   type=int, default=0)
    parser.add_argument("--total_workers", type=int, default=1)
    parser.add_argument("--density_threshold", type=float, default=0.1)
    args = parser.parse_args()

    k_values   = [int(x) for x in args.k_values.split(",")]
    run_name   = f"cnmf_{args.cancer}"
    output_dir = os.path.join(args.output_base, args.cancer)
    os.makedirs(output_dir, exist_ok=True)

    cnmf_obj = cNMF(output_dir=output_dir, name=run_name)

    if args.step == "prepare":
        print(f"\n{'='*50}\nPREPARE  cancer={args.cancer}  K={k_values}  n_iter={args.n_iter}  n_cells={args.n_cells}\n{'='*50}")
        adata     = load_counts(args.cancer, args.n_cells)
        h5ad_path = os.path.join(output_dir, f"{run_name}.h5ad")
        adata.write_h5ad(h5ad_path)
        cnmf_obj.prepare(
            counts_fn         = h5ad_path,
            components        = k_values,
            n_iter            = args.n_iter,
            seed              = SEED,
            num_highvar_genes = N_HVG,
        )
        print("Prepare complete.")

    elif args.step == "factorize":
        print(f"\n{'='*50}\nFACTORIZE  cancer={args.cancer}  worker={args.worker}/{args.total_workers}\n{'='*50}")
        cnmf_obj.factorize(worker_i=args.worker, total_workers=args.total_workers)
        print("Factorize complete.")

    elif args.step == "combine":
        print(f"\n{'='*50}\nCOMBINE  cancer={args.cancer}\n{'='*50}")
        cnmf_obj.combine()
        print("Combine complete.")

    elif args.step == "consensus":
        if args.k is None:
            raise ValueError("--k required for consensus")
        print(f"\n{'='*50}\nCONSENSUS  cancer={args.cancer}  K={args.k}\n{'='*50}")
        cnmf_obj.consensus(
            k=args.k,
            density_threshold=args.density_threshold,
            show_clustering=False,
        )
        print(f"Consensus complete. Output: {output_dir}/{run_name}/")


if __name__ == "__main__":
    main()
