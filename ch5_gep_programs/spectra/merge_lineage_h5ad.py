#!/usr/bin/env python
"""
merge_lineage_h5ad.py

Ch.5 (GEP programs) — Spectra pipeline step B. Concatenate the per-cancer
h5ad files for one immune lineage into a single AnnData, with optional
stratified subsampling (used for the huge lineages T-cells / Myeloid).
Output is ready for Spectra.

Per-file h5ads are produced by ../spectra/seurat_to_h5ad_spectra.R
(X = log-normalised, layers['counts'] = raw counts, full obs).

Subsampling (when --max-cells is set) is stratified by (cancer_type x cell_type_key)
so rare immune populations and small cancers are not lost.

Usage:
  python merge_lineage_h5ad.py \
    --lineage T-cells \
    --h5ad-dir /path/to/temp_objects/clustered_lineages/h5ad \
    --out /path/to/temp_objects/clustered_lineages/spectra_out/T-cells_immune.h5ad \
    --cell-type-key scimilarity_celltype \
    --max-cells 200000        # omit / 0 -> keep all cells
"""
import argparse
import glob
import os
import re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc


def log(msg):
    import datetime
    print(f"[{datetime.datetime.now()}] {msg}", flush=True)


def stratified_subsample(adata, group_cols, max_cells, seed=0):
    """Proportional stratified subsample down to ~max_cells across group_cols."""
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    obs = adata.obs
    grp = obs.groupby(group_cols, observed=True).indices  # dict -> positional idx
    frac = max_cells / adata.n_obs
    keep = []
    for _, idx in grp.items():
        n = len(idx)
        # at least 1 (if group non-empty), proportional otherwise
        k = max(1, int(round(n * frac))) if n > 0 else 0
        k = min(k, n)
        keep.append(rng.choice(idx, size=k, replace=False))
    keep = np.sort(np.concatenate(keep))
    log(f"Stratified subsample: {adata.n_obs} -> {len(keep)} cells "
        f"(target {max_cells}, {len(grp)} strata)")
    return adata[keep].copy()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lineage", required=True,
                   help="Lineage prefix, e.g. T-cells, Myeloid, B-cell, NK-cells, "
                        "Hematopoietic_Stem_and_Progenitor_cells")
    p.add_argument("--h5ad-dir", required=True,
                   help="Directory holding per-cancer <lineage>_<cancer>_clustered.h5ad")
    p.add_argument("--out", required=True, help="Output merged h5ad path")
    p.add_argument("--cell-type-key", required=True,
                   help="obs column to use as Spectra cell_type_key")
    p.add_argument("--max-cells", type=int, default=0,
                   help="Stratified subsample target; 0 = keep all")
    p.add_argument("--min-cells-per-type", type=int, default=50,
                   help="Drop cell types with fewer than this many cells (after merge)")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    # Per-cancer files for this lineage. Anchor the prefix so 'B-cell' does not
    # accidentally match other names; cancer token sits between prefix and suffix.
    pattern = os.path.join(args.h5ad_dir, f"{args.lineage}_*_clustered.h5ad")
    files = sorted(glob.glob(pattern))
    # Guard: exact prefix match (avoid e.g. NK-cells vs ...).
    rgx = re.compile(rf"^{re.escape(args.lineage)}_.+_clustered\.h5ad$")
    files = [f for f in files if rgx.match(os.path.basename(f))]
    if not files:
        raise FileNotFoundError(f"No h5ad files matched {pattern}")
    log(f"{args.lineage}: {len(files)} per-cancer h5ad files")

    parts = []
    for f in files:
        cancer = re.sub(rf"^{re.escape(args.lineage)}_", "",
                        os.path.basename(f)).replace("_clustered.h5ad", "")
        a = sc.read_h5ad(f)
        if args.cell_type_key not in a.obs.columns:
            raise KeyError(f"{f}: missing cell_type_key '{args.cell_type_key}'. "
                           f"Available: {list(a.obs.columns)[:40]}")
        a.obs["cancer_type"] = cancer
        a.obs["source_file"] = os.path.basename(f)
        # make barcodes unique across cancers
        a.obs_names = [f"{cancer}_{bc}" for bc in a.obs_names]
        parts.append(a)
        log(f"  loaded {cancer}: {a.n_obs} cells x {a.n_vars} genes")

    log("Concatenating (inner join on genes shared by all cancers)")
    adata = ad.concat(parts, join="inner", index_unique=None, merge="same")
    del parts
    log(f"Merged: {adata.n_obs} cells x {adata.n_vars} genes")

    # Clean the cell-type label column
    ck = args.cell_type_key
    adata.obs[ck] = adata.obs[ck].astype(str).fillna("Unknown").replace(
        {"nan": "Unknown", "NA": "Unknown", "": "Unknown"})

    # Drop tiny / unknown cell types that would give Spectra unstable per-type factors
    vc = adata.obs[ck].value_counts()
    keep_types = vc[vc >= args.min_cells_per_type].index.tolist()
    dropped = vc[vc < args.min_cells_per_type]
    if len(dropped):
        log(f"Dropping {len(dropped)} cell types < {args.min_cells_per_type} cells: "
            f"{dropped.to_dict()}")
    adata = adata[adata.obs[ck].isin(keep_types)].copy()
    log(f"After cell-type filter: {adata.n_obs} cells, "
        f"{adata.obs[ck].nunique()} cell types")

    # Optional stratified subsample
    if args.max_cells and args.max_cells > 0:
        adata = stratified_subsample(
            adata, group_cols=["cancer_type", ck],
            max_cells=args.max_cells, seed=args.seed)

    log(f"Final: {adata.n_obs} cells x {adata.n_vars} genes")
    log(f"cell_type_key '{ck}' distribution:\n{adata.obs[ck].value_counts()}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    adata.write_h5ad(args.out)
    log(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
