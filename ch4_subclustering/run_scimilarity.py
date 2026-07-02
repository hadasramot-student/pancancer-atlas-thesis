#!/usr/bin/env python
"""
run_scimilarity.py
------------------
Query a Seurat-exported h5ad against the SCimilarity reference model and
write per-cell predictions to a CSV file.

Expected h5ad layout (produced by seurat_to_h5ad.R):
  adata.X          : log-normalised counts  (cells × genes)
  adata.layers["counts"] : raw integer counts

Usage:
  python run_scimilarity.py \
      --h5ad   <path/to/object.h5ad> \
      --out_csv <path/to/predictions.csv> \
      --model_path <path/to/scimilarity/model_v1.1> \
      [--k 50] [--threshold 0.60] [--use_gpu]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad",        required=True,  help="Input h5ad (from seurat_to_h5ad.R)")
    p.add_argument("--out_csv",     required=True,  help="Output CSV: cell | scimilarity_label | confidence")
    p.add_argument("--model_path",  required=True,  help="Path to SCimilarity model directory (model_v1.1)")
    p.add_argument("--k",           type=int, default=50,   help="KNN neighbours [50]")
    p.add_argument("--threshold",   type=float, default=0.60, help="vs2nd confidence threshold [0.60]")
    p.add_argument("--use_gpu",     action="store_true",    help="Use GPU for embedding")
    p.add_argument("--layer_counts", default="counts",      help="Layer name for raw counts [counts]")
    args = p.parse_args()

    try:
        from scimilarity import CellAnnotation
        from scimilarity.utils import align_dataset
    except ImportError:
        sys.exit("ERROR: scimilarity is not installed in this environment.\n"
                 "Install it or activate the correct conda env before running this script.")

    print(f"[scimilarity] Loading h5ad: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"[scimilarity] {adata.n_obs} cells × {adata.n_vars} genes")

    # Ensure X contains log-normalised counts (adata.X from seurat_to_h5ad.R is already lognorm)
    # If raw counts are in a layer, they are ignored here — scimilarity needs lognorm in X
    if args.layer_counts in adata.layers:
        print(f"[scimilarity] Raw counts found in layer '{args.layer_counts}' — X kept as lognorm")

    print(f"[scimilarity] Loading model from: {args.model_path}")
    ca = CellAnnotation(model_path=args.model_path, use_gpu=args.use_gpu)

    print("[scimilarity] Aligning genes to model gene order ...")
    adata = align_dataset(adata, ca.gene_order)

    print("[scimilarity] Computing embeddings ...")
    adata.obsm["X_scimilarity"] = ca.get_embeddings(adata.X)

    print(f"[scimilarity] KNN prediction (k={args.k}) ...")
    pred, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(
        adata.obsm["X_scimilarity"],
        k=args.k,
        weighting=True,
    )

    df = pd.DataFrame({
        "cell":               adata.obs_names,
        "scimilarity_label":  pred.values,
        "confidence":         nn_stats["vs2nd_weighted"].values,
    })
    df["scimilarity_label_filtered"] = np.where(
        df["confidence"] >= args.threshold,
        df["scimilarity_label"],
        "Unknown"
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"[scimilarity] Saved {len(df)} predictions → {args.out_csv}")
    print("[scimilarity] Done.")


if __name__ == "__main__":
    main()
