#!/usr/bin/env python
import scipy.io, scipy.sparse as sp
import anndata as ad
import pandas as pd
import numpy as np

MTX_DIR = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Brain/brain_mtx"
OUT     = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Brain/Epithelial_Brain_merged_full.h5ad"

print("Loading mtx...")
mat      = scipy.io.mmread(f"{MTX_DIR}/matrix.mtx").tocsr().T  # transpose: cells x genes
genes    = pd.read_csv(f"{MTX_DIR}/genes.tsv",   header=None)[0].tolist()
barcodes = pd.read_csv(f"{MTX_DIR}/barcodes.tsv",header=None)[0].tolist()
meta     = pd.read_csv(f"{MTX_DIR}/metadata.csv",index_col=0, low_memory=False)

print(f"Matrix: {mat.shape[0]} cells x {mat.shape[1]} genes")

# Back-transform if log-normalized
sample = mat.data[:200]
if len(sample) > 0 and not np.all(sample == np.floor(sample)):
    print("Back-transforming log-normalized counts...")
    mat = mat.copy()
    mat.data = np.round(np.expm1(mat.data)).astype(np.float32)
    mat.data[mat.data < 0] = 0

adata = ad.AnnData(X=mat)
adata.obs_names = barcodes
adata.var_names = genes
adata.layers["counts"] = mat.copy()

# Attach all metadata, fixing reserved columns and mixed types
reserved = {"_index"}
for col in meta.columns:
    if col in reserved:
        continue
    vals = meta.reindex(adata.obs_names)[col]
    if vals.dtype == object:
        vals = vals.fillna("").astype(str)
    adata.obs[col] = vals.values

print(f"obs columns: {list(adata.obs.columns)}")
print(f"cell_type: {adata.obs['cell_type'].value_counts().to_dict()}")

print(f"Writing to: {OUT}")
adata.write_h5ad(OUT)
print("Done.")
