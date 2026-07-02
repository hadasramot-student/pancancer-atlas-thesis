#!/usr/bin/env python
# 02_run_cnmf.py
#
# Ch.5 (GEP programs) — cNMF step 2/4. K-sweep (K=10..40) + final consensus.
# Kotliar et al. eLife 2019
#
# Steps:
#   1. Load raw counts -> build AnnData
#   2. cNMF prepare  (HVG selection, TPM normalization, write input files)
#   3. cNMF factorize (NMF runs per K, parallelized)
#   4. cNMF combine  (aggregate NMF runs per K)
#   5. cNMF consensus (select K, compute GEPs and usages)
#
# Inputs:  ATLAS_DATA_DIR/cnmf_input/{counts.mtx,genes.tsv,barcodes.tsv} (from 01_export_counts_for_cnmf.R)
# Outputs: ATLAS_DATA_DIR/cnmf_output/<RUN_NAME>/ (GEP spectra + usage matrices per K)
#
# Usage:
#   python 02_run_cnmf.py --step prepare
#   python 02_run_cnmf.py --step factorize --k 20 --worker 0 --total_workers 4
#   python 02_run_cnmf.py --step combine
#   python 02_run_cnmf.py --step consensus --k 20 --density_threshold 0.1

import argparse
import os
import scipy.io
import scipy.sparse as sp
import numpy as np
import pandas as pd
import anndata as ad
from cnmf import cNMF

# ── Config ─────────────────────────────────────────────────────────────────
ATLAS_DATA_DIR = os.environ.get("ATLAS_DATA_DIR", "data/temp_objects")
MTX_DIR    = os.path.join(ATLAS_DATA_DIR, "cnmf_input")
OUTPUT_DIR = os.path.join(ATLAS_DATA_DIR, "cnmf_output")
RUN_NAME   = "pancancer_epithelial"
N_HVG      = 2000       # HVGs for cNMF (Kotliar recommends 2000)
N_ITER     = 200        # NMF iterations per run
SEED       = 42
K_VALUES   = list(range(10, 45, 5))  # [10, 15, 20, 25, 30, 35, 40]
N_RUNS     = 200        # replicates per K for stability assessment

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_adata():
    print("Loading count matrix...")
    mat = scipy.io.mmread(os.path.join(MTX_DIR, "counts.mtx")).tocsr()
    genes    = pd.read_csv(os.path.join(MTX_DIR, "genes.tsv"),    header=None)[0].tolist()
    barcodes = pd.read_csv(os.path.join(MTX_DIR, "barcodes.tsv"), header=None)[0].tolist()
    meta     = pd.read_csv(os.path.join(MTX_DIR, "cell_metadata.csv"), index_col=0)

    # counts.mtx from R writeMM is genes x cells — transpose to cells x genes
    mat = mat.T.tocsr()
    print(f"Matrix: {mat.shape[0]} cells x {mat.shape[1]} genes")

    adata = ad.AnnData(X=mat)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.layers["counts"] = mat.copy()

    # attach metadata
    common_cells = meta.index.intersection(adata.obs_names)
    print(f"Cells with metadata: {len(common_cells)} / {adata.n_obs}")
    for col in meta.columns:
        adata.obs[col] = meta.reindex(adata.obs_names)[col].values

    # validate counts are integers
    sample = mat.data[:200]
    if not np.all(sample == np.floor(sample)):
        print("WARNING: counts layer has normalized floats — back-transforming with expm1() + round()")
        print("  (standard approach when raw counts unavailable)")
        # reverse log1p normalization → approximate integer pseudo-counts
        mat_pseudo = mat.copy()
        mat_pseudo.data = np.round(np.expm1(mat_pseudo.data)).astype(np.float32)
        mat_pseudo.data[mat_pseudo.data < 0] = 0
        adata.X = mat_pseudo
        adata.layers["counts"] = mat_pseudo.copy()
        print(f"  Pseudo-counts sample: {mat_pseudo.data[:5]}")
        print(f"  Are integers now: {np.all(mat_pseudo.data == np.floor(mat_pseudo.data))}")

    print(f"Cancer type distribution:")
    if "cancer_type" in adata.obs.columns:
        print(adata.obs["cancer_type"].value_counts().to_string())

    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", required=True,
                        choices=["prepare", "factorize", "combine", "consensus", "all"])
    parser.add_argument("--k",              type=int, default=None)
    parser.add_argument("--worker",         type=int, default=0)
    parser.add_argument("--total_workers",  type=int, default=1)
    parser.add_argument("--density_threshold", type=float, default=0.1)
    args = parser.parse_args()

    cnmf_obj = cNMF(output_dir=OUTPUT_DIR, name=RUN_NAME)

    if args.step in ("prepare", "all"):
        print(f"\n{'='*50}\nStep: PREPARE\nK values: {K_VALUES}\nHVGs: {N_HVG}\n{'='*50}")
        adata = load_adata()

        # Save h5ad for cNMF input
        h5ad_path = os.path.join(OUTPUT_DIR, f"{RUN_NAME}.h5ad")
        adata.write_h5ad(h5ad_path)
        print(f"Saved AnnData to: {h5ad_path}")

        cnmf_obj.prepare(
            counts_fn    = h5ad_path,
            components   = K_VALUES,
            n_iter       = N_RUNS,
            seed         = SEED,
            num_highvar_genes = N_HVG,
        )
        print("Prepare complete.")

    if args.step in ("factorize", "all"):
        k_list = [args.k] if args.k else K_VALUES
        print(f"\n{'='*50}\nStep: FACTORIZE\nK: {k_list}\nWorker: {args.worker}/{args.total_workers}\n{'='*50}")
        cnmf_obj.factorize(
            worker_i      = args.worker,
            total_workers = args.total_workers,
        )
        print("Factorize complete.")

    if args.step in ("combine", "all"):
        print(f"\n{'='*50}\nStep: COMBINE\n{'='*50}")
        cnmf_obj.combine()
        print("Combine complete.")

    if args.step in ("consensus", "all"):
        k = args.k
        if k is None:
            raise ValueError("--k required for consensus step")
        print(f"\n{'='*50}\nStep: CONSENSUS  K={k}\ndensity_threshold={args.density_threshold}\n{'='*50}")
        cnmf_obj.consensus(
            k                  = k,
            density_threshold  = args.density_threshold,
            show_clustering    = False,
        )
        print(f"Consensus complete for K={k}.")
        print(f"Output: {OUTPUT_DIR}/{RUN_NAME}/")


if __name__ == "__main__":
    main()
