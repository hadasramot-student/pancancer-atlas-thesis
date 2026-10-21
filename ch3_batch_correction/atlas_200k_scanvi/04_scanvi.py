#!/usr/bin/env python3
"""
04_scanvi.py -- pan-atlas 200k-cell integration for Figure 1 / Figure 2.

One GPU job produces BOTH coordinate systems on the SAME cells / SAME HVGs:
  * Figure 2 "before" : raw PCA (30 PCs) -> UMAP           [X_pca,   X_umap_pca]
  * Figure 1 / Fig 2 "after" : scVI -> scANVI latent (30d) -> UMAP  [X_scANVI, X_umap_scanvi]

Procedure mirrors the Chapter 3 benchmarked / selected setup (scANVI), run once
across all 200k assembled cells instead of per lineage x cancer-type.

Input  h5ad (from 03_merge_and_h5ad.R):
  .X            = LogNormalize(scale.factor 1e4) expression  (Seurat NormalizeData)
  .layers[counts] = raw counts
  .obs          : project_name, Consensus_Cell_Type, cancer_type, source_repo,
                  agreement, cell_id
Outputs (--out_prefix):
  <p>_coords.csv        cell_id + PCA_UMAP1/2 + scANVI_UMAP1/2 + PC1..30 + scANVI_1..30
  <p>_scanvi_pred.csv   cell_id + scanvi_predicted_label
  <p>_integrated.h5ad   full object with obsm + hvg flag
  <p>_run_params.json   exact parameters + realised epochs + early-stop info
"""
import os, json, argparse, datetime
import numpy as np
import pandas as pd
import scanpy as sc
import scvi


def log(m): print(f"[{datetime.datetime.now():%H:%M:%S}] {m}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--batch_key", default="project_name")
    ap.add_argument("--labels_key", default="Consensus_Cell_Type")
    ap.add_argument("--n_hvg", type=int, default=2000)
    ap.add_argument("--n_latent", type=int, default=30)
    ap.add_argument("--n_layers", type=int, default=2)
    ap.add_argument("--scvi_epochs", type=int, default=400)
    ap.add_argument("--scanvi_epochs", type=int, default=200)
    ap.add_argument("--n_neighbors", type=int, default=30)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    scvi.settings.seed = args.seed
    sc.settings.verbosity = 1
    params = vars(args).copy()
    params["started"] = str(datetime.datetime.now())

    # ---------------------------------------------------------------- load ----
    log(f"loading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    log(f"{adata.n_obs} cells x {adata.n_vars} genes")
    assert "counts" in adata.layers, "expected raw 'counts' layer"
    for k in (args.batch_key, args.labels_key, "cell_id"):
        assert k in adata.obs, f"missing obs column {k}"
    adata.obs[args.batch_key] = adata.obs[args.batch_key].astype("category")
    adata.obs[args.labels_key] = adata.obs[args.labels_key].astype(str)

    params["n_cells"] = int(adata.n_obs)
    params["n_genes_in"] = int(adata.n_vars)
    params["n_batches"] = int(adata.obs[args.batch_key].nunique())
    params["label_levels"] = sorted(adata.obs[args.labels_key].unique().tolist())

    # ------------------------------------------------- HVG (shared, seurat_v3) ----
    log(f"HVG: top {args.n_hvg} (seurat_v3 on counts, batch={args.batch_key})")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=args.n_hvg, flavor="seurat_v3",
        batch_key=args.batch_key, layer="counts", subset=False,
    )
    hvg = adata[:, adata.var["highly_variable"]].copy()
    log(f"HVG subset: {hvg.n_vars} genes")

    # =========================================================== raw PCA ====
    # Figure 2 "before": standard Seurat-style PCA on scaled LogNormalized HVGs.
    log("raw PCA branch: scale -> PCA(30) -> neighbors -> UMAP")
    pca_ad = hvg.copy()                       # .X is LogNormalize(1e4) from Seurat
    sc.pp.scale(pca_ad, max_value=10)
    sc.tl.pca(pca_ad, n_comps=args.n_latent, svd_solver="arpack", random_state=args.seed)
    sc.pp.neighbors(pca_ad, use_rep="X_pca", n_neighbors=args.n_neighbors, random_state=args.seed)
    sc.tl.umap(pca_ad, random_state=args.seed)
    X_pca = pca_ad.obsm["X_pca"].copy()
    U_pca = pca_ad.obsm["X_umap"].copy()
    del pca_ad

    # =========================================================== scVI/scANVI ====
    log("scVI setup")
    scvi.model.SCVI.setup_anndata(hvg, layer="counts", batch_key=args.batch_key)
    m_scvi = scvi.model.SCVI(hvg, n_latent=args.n_latent, n_layers=args.n_layers)
    log(f"train scVI (<= {args.scvi_epochs} epochs, early stopping)")
    m_scvi.train(max_epochs=args.scvi_epochs, early_stopping=True,
                 accelerator="gpu", devices="auto")

    def n_epochs(model):
        for k in ("elbo_train", "elbo_validation", "train_loss_epoch",
                  "reconstruction_loss_train"):
            if k in model.history and len(model.history[k]):
                return int(len(model.history[k]))
        return -1
    params["scvi_epochs_run"] = n_epochs(m_scvi)

    log("scANVI from scVI")
    m_scan = scvi.model.SCANVI.from_scvi_model(
        m_scvi, adata=hvg, labels_key=args.labels_key, unlabeled_category="Unknown")
    log(f"train scANVI (<= {args.scanvi_epochs} epochs, early stopping)")
    m_scan.train(max_epochs=args.scanvi_epochs, early_stopping=True,
                 accelerator="gpu", devices="auto")
    params["scanvi_epochs_run"] = n_epochs(m_scan)

    log("scANVI latent + predictions")
    X_scanvi = m_scan.get_latent_representation()
    pred = m_scan.predict(hvg)

    log("scANVI neighbors -> UMAP")
    hvg.obsm["X_scANVI"] = X_scanvi
    sc.pp.neighbors(hvg, use_rep="X_scANVI", n_neighbors=args.n_neighbors, random_state=args.seed)
    sc.tl.umap(hvg, random_state=args.seed)
    U_scanvi = hvg.obsm["X_umap"].copy()

    # ------------------------------------------------------------- write out ----
    uid = adata.obs_names.astype(str).values          # row_uid = project_name|cell_id
    coords = pd.DataFrame({"row_uid": uid,
                           "PCA_UMAP1": U_pca[:, 0], "PCA_UMAP2": U_pca[:, 1],
                           "scANVI_UMAP1": U_scanvi[:, 0], "scANVI_UMAP2": U_scanvi[:, 1]})
    for i in range(args.n_latent):
        coords[f"PC{i+1}"] = X_pca[:, i]
    for i in range(X_scanvi.shape[1]):
        coords[f"scANVI_{i+1}"] = X_scanvi[:, i]
    coords.to_csv(f"{args.out_prefix}_coords.csv", index=False)
    log(f"wrote {args.out_prefix}_coords.csv")

    pd.DataFrame({"row_uid": uid, "scanvi_predicted_label": pred}).to_csv(
        f"{args.out_prefix}_scanvi_pred.csv", index=False)

    adata.obsm["X_pca"] = X_pca
    adata.obsm["X_umap_pca"] = U_pca
    adata.obsm["X_scANVI"] = X_scanvi
    adata.obsm["X_umap_scanvi"] = U_scanvi
    adata.obs["scanvi_predicted_label"] = pred
    adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)
    adata.write_h5ad(f"{args.out_prefix}_integrated.h5ad")
    log(f"wrote {args.out_prefix}_integrated.h5ad")

    params["finished"] = str(datetime.datetime.now())
    params["scvi_version"] = scvi.__version__
    params["scanpy_version"] = sc.__version__
    with open(f"{args.out_prefix}_run_params.json", "w") as fh:
        json.dump(params, fh, indent=2)
    log("DONE 04_scanvi")


if __name__ == "__main__":
    main()
