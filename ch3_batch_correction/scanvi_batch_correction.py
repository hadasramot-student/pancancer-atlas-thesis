#!/usr/bin/env python3
# scanvi_batch_correction.py
#
# Ch.3 (Batch correction) - scVI/scANVI batch correction. Trains scVI on raw
# counts (unsupervised), then scANVI (semi-supervised, seeded by an existing
# cell-type label) to produce a batch-corrected latent embedding.
#
# Inputs:  --input  h5ad with raw counts and a batch/label obs column
# Outputs: --output h5ad with X_scVI / X_scANVI in obsm and scanvi predictions
#
# Usage (run on a GPU node, conda env with scvi-tools):
#   python scanvi_batch_correction.py --input obj.h5ad --output obj_scanvi.h5ad

import os
import scanpy as sc
import scvi
from optparse import OptionParser


def log(msg):
    import datetime
    print(f"[{datetime.datetime.now()}] {msg}", flush=True)


def main():
    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-i", "--input",
                      dest="input",
                      help="Input h5ad file path (required)")
    parser.add_option("-o", "--output",
                      dest="output",
                      help="Output h5ad file path (required)")
    parser.add_option("-b", "--batch-col",
                      dest="batch_col",
                      default="project_name",
                      help="Column in adata.obs to use as batch key (default: project_name)")
    parser.add_option("-l", "--labels-col",
                      dest="labels_col",
                      default=None,
                      help="Column in adata.obs to use as cell-type labels for scANVI semi-supervised learning. "
                           "If not provided, all cells are treated as 'Unknown'.")
    parser.add_option("--n-hvg",
                      dest="n_hvg",
                      type="int",
                      default=3000,
                      help="Number of highly variable genes to use (default: 3000)")
    parser.add_option("--n-latent",
                      dest="n_latent",
                      type="int",
                      default=20,
                      help="Dimensionality of the latent space (default: 20)")
    parser.add_option("--n-layers",
                      dest="n_layers",
                      type="int",
                      default=2,
                      help="Number of hidden layers in the encoder/decoder (default: 2)")
    parser.add_option("--scvi-epochs",
                      dest="scvi_epochs",
                      type="int",
                      default=200,
                      help="Max training epochs for scVI (default: 200)")
    parser.add_option("--scanvi-epochs",
                      dest="scanvi_epochs",
                      type="int",
                      default=20,
                      help="Max training epochs for scANVI (default: 20)")
    parser.add_option("--n-neighbors",
                      dest="n_neighbors",
                      type="int",
                      default=30,
                      help="Number of neighbors for UMAP/clustering (default: 30)")
    parser.add_option("--resolution",
                      dest="resolution",
                      type="float",
                      default=0.8,
                      help="Leiden clustering resolution (default: 0.8)")
    parser.add_option("--use-gpu",
                      dest="use_gpu",
                      action="store_true",
                      default=False,
                      help="Use GPU for training (default: CPU)")

    (options, args) = parser.parse_args()

    if not options.input:
        parser.error("--input is required")
    if not options.output:
        parser.error("--output is required")

    if not options.use_gpu:
        scvi.settings.dl_num_workers = 0
        
    # --------------------------------------------------
    # Configure accelerator for scvi-tools
    # --------------------------------------------------
    if options.use_gpu:
        accelerator = "gpu"
        devices = "auto"
    else:
        accelerator = "cpu"
        devices = "auto"
        scvi.settings.dl_num_workers = 0
    # --------------------------------------------------
    # Load data
    # --------------------------------------------------
    log(f"Loading: {options.input}")
    adata = sc.read_h5ad(options.input)
    log(f"Loaded object: {adata.n_obs} cells x {adata.n_vars} genes")

    if options.batch_col not in adata.obs.columns:
        raise ValueError(
            f"Batch column '{options.batch_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # --------------------------------------------------
    # Ensure counts layer exists
    # --------------------------------------------------
    if "counts" in adata.layers:
        log("Using 'counts' layer for HVG selection and scVI/scANVI")
        count_layer = "counts"
    else:
        log("No 'counts' layer found — assuming adata.X contains raw counts")
        count_layer = None

    # --------------------------------------------------
    # HVG selection on raw counts
    # --------------------------------------------------
    log(f"Selecting {options.n_hvg} highly variable genes")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=options.n_hvg,
        subset=False,
        flavor="seurat_v3",
        batch_key=options.batch_col,
        layer=count_layer,
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    log(f"HVG subset: {adata_hvg.n_obs} cells x {adata_hvg.n_vars} genes")

    # --------------------------------------------------
    # Set up labels column for scANVI
    # --------------------------------------------------
    if options.labels_col and options.labels_col in adata_hvg.obs.columns:
        log(f"Using '{options.labels_col}' as cell-type labels for scANVI")
        adata_hvg.obs[options.labels_col] = adata_hvg.obs[options.labels_col].astype(str)
        adata_hvg.obs[options.labels_col] = adata_hvg.obs[options.labels_col].fillna("Unknown")
        labels_key = options.labels_col
    else:
        log("No valid labels column provided — assigning all cells as 'Unknown'")
        adata_hvg.obs["_scanvi_labels"] = "Unknown"
        labels_key = "_scanvi_labels"

    # --------------------------------------------------
    # scVI setup + training
    # --------------------------------------------------
    log("Setting up scVI model")
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer=count_layer,
        batch_key=options.batch_col,
    )

    scvi_model = scvi.model.SCVI(
        adata_hvg,
        n_latent=options.n_latent,
        n_layers=options.n_layers,
    )

    log(f"Training scVI (max {options.scvi_epochs} epochs)")
    scvi_model.train(
        max_epochs=options.scvi_epochs,
        early_stopping=True,
        accelerator=accelerator,
        devices=devices,
    )
    log("scVI training complete")

    # --------------------------------------------------
    # scANVI setup + training
    # --------------------------------------------------
    log("Setting up scANVI model from trained scVI")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",
        labels_key=labels_key,
        adata=adata_hvg,
    )

    log(f"Training scANVI (max {options.scanvi_epochs} epochs)")
    scanvi_model.train(
        max_epochs=options.scanvi_epochs,
        accelerator=accelerator,
        devices=devices,
    )
    log("scANVI training complete")

    # --------------------------------------------------
    # Extract latent embedding + predictions
    # --------------------------------------------------
    log("Extracting scANVI latent representation")
    adata_hvg.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
    adata_hvg.obs["scanvi_predicted_labels"] = scanvi_model.predict(adata_hvg)

    # --------------------------------------------------
    # Downstream analysis on scANVI embedding
    # --------------------------------------------------
    log(f"Computing neighbors (k={options.n_neighbors}) on scANVI embedding")
    sc.pp.neighbors(adata_hvg, use_rep="X_scANVI", n_neighbors=options.n_neighbors)

    log("Running UMAP")
    sc.tl.umap(adata_hvg)
    adata_hvg.obsm["X_umap_scANVI"] = adata_hvg.obsm["X_umap"].copy()

    log(f"Clustering with Leiden (resolution={options.resolution})")
    sc.tl.leiden(adata_hvg, resolution=options.resolution, key_added="scanvi_clusters")

    # --------------------------------------------------
    # Transfer results back to full AnnData
    # --------------------------------------------------
    log("Transferring results to full AnnData")
    adata.obsm["X_scANVI"] = adata_hvg.obsm["X_scANVI"]
    adata.obsm["X_umap_scANVI"] = adata_hvg.obsm["X_umap_scANVI"]
    adata.obs["scanvi_clusters"] = adata_hvg.obs["scanvi_clusters"]
    adata.obs["scanvi_predicted_labels"] = adata_hvg.obs["scanvi_predicted_labels"]
    adata.var["highly_variable_scANVI"] = adata.var_names.isin(adata_hvg.var_names)

    if labels_key == "_scanvi_labels":
        adata_hvg.obs.drop(columns=["_scanvi_labels"], inplace=True)

    # --------------------------------------------------
    # Saves
    # --------------------------------------------------
    log(f"Saving to: {options.output}")
    os.makedirs(os.path.dirname(os.path.abspath(options.output)), exist_ok=True)
    adata.write_h5ad(options.output)
    log("Done")


if __name__ == "__main__":
    main()