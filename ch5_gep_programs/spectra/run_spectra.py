#!/usr/bin/env python
"""
run_spectra.py

Ch.5 (GEP programs) — Spectra pipeline step E. Fit Spectra (prior-guided
topic model) on ONE immune lineage AnnData, on GPU.

Inputs:
  --h5ad          merged lineage h5ad (X = log-normalised, layers['counts'])
  --gene-set-json nested gene_set_dictionary from build_spectra_priors.py
  --cell-type-key obs column with immune cell types (matches dict keys)

Outputs (written under --out-prefix):
  *_cell_scores.csv     cells  x factors   (program loadings per cell  -> Task 3 vs cNMF)
  *_gene_factors.csv    factors x genes     (program gene weights)
  *_markers.csv         top genes per factor
  *_factor_meta.csv     factor -> cell type + matched prior set
  *_spectra.h5ad        adata with SPECTRA results in .obsm/.uns
  *_model.pickle        trained SPECTRA model (best effort)

Run via python.sbatch (activates scanvi-env) with --gres=gpu:...
"""
import argparse
import json
import os
# Reduce CUDA fragmentation before torch initialises the CUDA context.
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import pickle
import datetime
import numpy as np
import pandas as pd
import scanpy as sc
import torch

from Spectra import Spectra_gpu
from Spectra import Spectra_util


def log(msg):
    print(f"[{datetime.datetime.now()}] {msg}", flush=True)


def to_numpy(x):
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--gene-set-json", required=True)
    ap.add_argument("--cell-type-key", default="scimilarity_celltype")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-hvg", type=int, default=4000,
                    help="HVGs to union with prior genes for the Spectra vocab")
    ap.add_argument("--epochs", type=int, default=150,
                    help="Spectra_gpu training epochs (its default is ~100; "
                         "minibatched optimizer converges in this range — NOT the "
                         "CPU Spectra's 10000)")
    ap.add_argument("--lr", type=float, default=0.001)
    ap.add_argument("--lam", type=float, default=0.008)
    ap.add_argument("--delta", type=float, default=0.001)
    ap.add_argument("--rho", type=float, default=0.05)
    ap.add_argument("--no-cell-types", action="store_true",
                    help="Ignore cell types (single global factorisation)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    log(f"torch {torch.__version__} | CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        log(f"GPU: {torch.cuda.get_device_name(0)}")
    else:
        log("WARNING: no GPU visible -> Spectra_gpu will fall back to CPU (slow).")

    use_cell_types = not args.no_cell_types
    ck = args.cell_type_key

    # ---------------- load ----------------
    log(f"Loading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    log(f"{adata.n_obs} cells x {adata.n_vars} genes")
    if use_cell_types and ck not in adata.obs.columns:
        raise KeyError(f"cell_type_key '{ck}' not in obs")
    if use_cell_types:
        adata.obs[ck] = adata.obs[ck].astype(str)

    with open(args.gene_set_json) as fh:
        gsd = json.load(fh)

    # ---------------- vocab = HVG ∪ prior genes ----------------
    log(f"Selecting {args.n_hvg} HVGs (seurat_v3 on counts)")
    layer = "counts" if "counts" in adata.layers else None
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg,
                                flavor="seurat_v3", layer=layer, subset=False)
    hvg = set(adata.var_names[adata.var["highly_variable"]])

    prior_genes = set()
    for ct, sets in gsd.items():
        for _, genes in sets.items():
            prior_genes.update(genes)
    prior_genes &= set(adata.var_names)

    vocab = sorted(hvg | prior_genes)
    log(f"Vocab: {len(vocab)} genes ({len(hvg)} HVG ∪ {len(prior_genes)} prior genes)")
    adata = adata[:, vocab].copy()
    # Spectra_gpu has a bug indexing sparse X with a pandas Series when
    # use_highly_variable=True. Since we already subset adata to exactly the
    # model vocab, pass use_highly_variable=False and feed a dense X.
    import scipy.sparse as _sp
    if _sp.issparse(adata.X):
        adata.X = np.asarray(adata.X.todense(), dtype=np.float32)
    adata.var["spectra_vocab"] = True   # whole subsetted gene space is the vocab

    # ---------------- validate / clean gene set dict ----------------
    log("Validating gene set dictionary against data")
    gsd = Spectra_util.check_gene_set_dictionary(
        adata, gsd, obs_key=ck, global_key="global",
        return_dict=True, min_len=3, use_cell_types=use_cell_types)

    # report sizes
    n_global = len(gsd.get("global", {}))
    n_perct = {k: len(v) for k, v in gsd.items() if k != "global"}
    log(f"Dictionary after check: global={n_global} sets; "
        f"{len(n_perct)} cell types; total sets="
        f"{n_global + sum(n_perct.values())}")

    # ---------------- fit Spectra ----------------
    log(f"Fitting Spectra_gpu (epochs={args.epochs}, lam={args.lam}, "
        f"rho={args.rho}, use_cell_types={use_cell_types})")
    model = Spectra_gpu.est_spectra(
        adata=adata,
        gene_set_dictionary=gsd,
        L=None,                       # match #factors to #gene sets (recommended)
        use_highly_variable=False,    # adata already subset to the vocab above
        cell_type_key=ck if use_cell_types else None,
        use_cell_types=use_cell_types,
        use_weights=True,
        lam=args.lam, delta=args.delta, rho=args.rho,
        n_top_vals=50,
        num_epochs=args.epochs,
        lr_schedule=[args.lr] * args.epochs,
        verbose=True,
    )
    log("Spectra fit complete")

    # ---------------- collect + save ----------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)

    cell_scores = to_numpy(adata.obsm["SPECTRA_cell_scores"])
    factors = to_numpy(adata.uns["SPECTRA_factors"])          # factors x vocab
    vocab_genes = adata.var_names[adata.var["spectra_vocab"]].tolist() \
        if "spectra_vocab" in adata.var.columns else vocab
    factor_ids = [f"Spectra_{i+1}" for i in range(cell_scores.shape[1])]

    pd.DataFrame(cell_scores, index=adata.obs_names, columns=factor_ids) \
        .to_csv(f"{args.out_prefix}_cell_scores.csv")
    log(f"cell_scores: {cell_scores.shape} -> {args.out_prefix}_cell_scores.csv")

    pd.DataFrame(factors, index=factor_ids[:factors.shape[0]], columns=vocab_genes) \
        .to_csv(f"{args.out_prefix}_gene_factors.csv")
    log(f"gene_factors: {factors.shape} -> {args.out_prefix}_gene_factors.csv")

    # markers
    try:
        markers = adata.uns.get("SPECTRA_markers")
        md = pd.DataFrame(markers)
        md.index = [f"Spectra_{i+1}" for i in range(md.shape[0])]
        md.to_csv(f"{args.out_prefix}_markers.csv")
        log(f"markers -> {args.out_prefix}_markers.csv")
    except Exception as e:
        log(f"markers save skipped: {e}")

    # factor -> cell type mapping
    try:
        if use_cell_types:
            fct = Spectra_gpu.get_factor_celltypes(adata, ck) \
                if hasattr(Spectra_gpu, "get_factor_celltypes") else None
            if fct is not None:
                pd.Series(fct, name="cell_type").to_csv(
                    f"{args.out_prefix}_factor_meta.csv")
                log(f"factor_meta -> {args.out_prefix}_factor_meta.csv")
    except Exception as e:
        log(f"factor_meta skipped: {e}")

    # adata with results (coerce uns tensors to numpy first)
    try:
        adata.obsm["SPECTRA_cell_scores"] = cell_scores
        adata.uns["SPECTRA_factors"] = factors
        adata.write_h5ad(f"{args.out_prefix}_spectra.h5ad")
        log(f"adata -> {args.out_prefix}_spectra.h5ad")
    except Exception as e:
        log(f"h5ad save skipped: {e}")

    # model pickle
    try:
        with open(f"{args.out_prefix}_model.pickle", "wb") as fh:
            pickle.dump(model, fh)
        log(f"model -> {args.out_prefix}_model.pickle")
    except Exception as e:
        log(f"model pickle skipped: {e}")

    log("ALL DONE")


if __name__ == "__main__":
    main()
