#!/usr/bin/env python
# export_scanvi_to_csv.py
#
# Ch.3 (Batch correction) - runs scVI + scANVI on one h5ad (converted from
# Seurat via seurat_to_h5ad.R) and exports the latent embedding + predicted
# labels to CSV, so they can be re-imported into the Seurat object
# (scanvi_into_seurat.R) alongside the CCA/MNN/ComBat/Harmony reductions.
#
# Inputs:  --h5ad  Seurat-exported h5ad with a raw-counts layer
# Outputs: <output_prefix>_scanvi_latent.csv, <output_prefix>_scanvi_labels.csv
#
# Usage: see scanvi_example.slurm for a full invocation.

import os
import argparse
import numpy as np
import pandas as pd

import scanpy as sc
import scvi


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run SCANVI on an h5ad file (from Seurat) and export latent space "
            "+ labels to CSV for Seurat."
        )
    )
    parser.add_argument(
        "--h5ad",
        required=True,
        help="Input h5ad file converted from Seurat."
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix for output CSV files (latent and labels)."
    )
    parser.add_argument(
        "--batch_key",
        default="project_name",
        help="Column in adata.obs with batch info [default: project_name]."
    )
    parser.add_argument(
        "--labels_key",
        default="Meta.scimilarity",
        help="Column in adata.obs with initial labels for SCANVI [default: Meta.scimilarity]."
    )
    parser.add_argument(
        "--unlabeled_category",
        default="Unknown",
        help="Category in labels_key treated as unlabeled [default: Unknown]."
    )
    parser.add_argument(
        "--layer_counts",
        default=None,
        help="Name of layer containing raw counts (if None, use adata.X)."
    )
    parser.add_argument(
        "--n_latent",
        type=int,
        default=30,
        help="Latent dimensionality for SCVI/SCANVI [default: 30]."
    )
    parser.add_argument(
        "--max_epochs_scvi",
        type=int,
        default=200,
        help="Max epochs for SCVI pretraining [default: 200]."
    )
    parser.add_argument(
        "--max_epochs_scanvi",
        type=int,
        default=100,
        help="Max epochs for SCANVI [default: 100]."
    )
    parser.add_argument(
        "--label_key_out",
        default="scanvi_pred",
        help="Name of obs/CSV column for predicted labels [default: scanvi_pred]."
    )
    parser.add_argument(
        "--n_hvgs",
        type=int,
        default=4000,
        help="Number of highly variable genes to select before training [default: 4000]. Set to 0 to skip HVG selection."
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    print(f"Loading AnnData from: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    print(adata)
    print("obs columns:", list(adata.obs.columns))
    print("layers:", list(adata.layers.keys()) if hasattr(adata, "layers") else None)

    if args.layer_counts is not None and args.layer_counts not in adata.layers:
        raise ValueError(f"layer_counts '{args.layer_counts}' not found in adata.layers")

    if args.batch_key not in adata.obs:
        raise ValueError(f"batch_key '{args.batch_key}' not found in adata.obs")

    if args.labels_key not in adata.obs:
        raise ValueError(f"labels_key '{args.labels_key}' not found in adata.obs")

    # ------------------------------------------------------------------
    # 1b. HVG selection (seurat_v3 works on raw counts)
    # ------------------------------------------------------------------
    if args.n_hvgs > 0:
        n_hvgs = min(args.n_hvgs, adata.n_vars)
        print(f"Selecting top {n_hvgs} highly variable genes (seurat_v3, from {adata.n_vars} total)...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvgs,
            flavor="seurat_v3",
            layer=args.layer_counts,
            subset=True
        )
        print(f"Genes after HVG selection: {adata.n_vars}")

    # ------------------------------------------------------------------
    # 2. Setup and train SCVI
    # ------------------------------------------------------------------
    print("Setting up SCVI anndata...")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=args.layer_counts,          # None => use adata.X
        batch_key=args.batch_key,
        labels_key=args.labels_key
    )

    print("Initializing SCVI model...")
    scvi_model = scvi.model.SCVI(adata, n_latent=args.n_latent)

    print(f"Training SCVI for up to {args.max_epochs_scvi} epochs (early stopping enabled)...")
    scvi_model.train(max_epochs=args.max_epochs_scvi, early_stopping=True)

    # ------------------------------------------------------------------
    # 3. Setup and train SCANVI
    # ------------------------------------------------------------------
    print("Initializing SCANVI model from SCVI...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=args.unlabeled_category
    )

    print(f"Training SCANVI for up to {args.max_epochs_scanvi} epochs (early stopping enabled)...")
    scanvi_model.train(max_epochs=args.max_epochs_scanvi, early_stopping=True)

    # ------------------------------------------------------------------
    # 4. Get latent representation
    # ------------------------------------------------------------------
    print("Computing SCANVI latent representation...")
    latent = scanvi_model.get_latent_representation(adata)
    print("Latent shape:", latent.shape)

    latent_df = pd.DataFrame(
        latent,
        index=adata.obs_names,
        columns=[f"SCANVI_{i+1}" for i in range(latent.shape[1])]
    )
    latent_df.insert(0, "cell", latent_df.index)

    latent_csv = f"{args.output_prefix}_scanvi_latent.csv"
    os.makedirs(os.path.dirname(latent_csv), exist_ok=True)
    latent_df.to_csv(latent_csv, index=False)
    print(f"Wrote latent CSV to: {latent_csv}")

    # ------------------------------------------------------------------
    # 5. Predict labels
    # ------------------------------------------------------------------
    print("Predicting labels with SCANVI...")
    pred_labels = scanvi_model.predict(adata)
    adata.obs[args.label_key_out] = pred_labels

    labels_df = pd.DataFrame({
        "cell": adata.obs_names,
        args.label_key_out: adata.obs[args.label_key_out].astype(str).values
    })

    labels_csv = f"{args.output_prefix}_scanvi_labels.csv"
    os.makedirs(os.path.dirname(labels_csv), exist_ok=True)
    labels_df.to_csv(labels_csv, index=False)
    print(f"Wrote labels CSV to: {labels_csv}")

    print("Done.")


if __name__ == "__main__":
    main()
