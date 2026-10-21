#!/usr/bin/env python
# stack_breast_chunks.py
#
# Stacks per-chunk .mtx files into a single AnnData h5ad with raw counts.
# Uses scipy sparse (int32-safe via chunked vstack) and anndata.

import os
import scipy.io
import scipy.sparse as sp
import anndata as ad
import numpy as np
import pandas as pd

MTX_DIR  = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Breast/epithelial_mtx_chunks"
OUT_H5AD = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Breast/Epithelial_Breast_merged_full.h5ad"

# Load genes
genes = pd.read_csv(os.path.join(MTX_DIR, "genes.tsv"), header=None, names=["gene"])
print(f"Genes: {len(genes)}")

# Find chunk dirs
chunk_dirs = sorted([
    os.path.join(MTX_DIR, d)
    for d in os.listdir(MTX_DIR)
    if os.path.isdir(os.path.join(MTX_DIR, d)) and d.startswith("chunk_")
])
print(f"Chunks: {len(chunk_dirs)}")

# Stack chunks one by one — use lil_matrix accumulation to avoid cbind overflow
all_barcodes = []
matrices = []

for i, cdir in enumerate(chunk_dirs):
    print(f"  Loading chunk {i+1}/{len(chunk_dirs)}: {os.path.basename(cdir)}")
    mat = scipy.io.mmread(os.path.join(cdir, "matrix.mtx")).tocsc()
    barcodes = pd.read_csv(os.path.join(cdir, "barcodes.tsv"), header=None)[0].tolist()
    print(f"    {mat.shape[1]} cells x {mat.shape[0]} genes")
    # matrix is genes x cells (from writeMM in R), transpose to cells x genes
    matrices.append(mat.T.tocsr())
    all_barcodes.extend(barcodes)

print(f"Stacking {len(matrices)} matrices (cells x genes)...")
# vstack works on CSR (cells x genes) — no 32-bit overflow since each chunk is small
combined = sp.vstack(matrices, format="csr")
print(f"Combined: {combined.shape[0]} cells x {combined.shape[1]} genes")

# Build AnnData
print("Building AnnData...")
adata = ad.AnnData(X=combined)
adata.obs_names = all_barcodes
adata.var_names = genes["gene"].tolist()

# Store as raw counts layer (required for cNMF and scANVI)
adata.layers["counts"] = combined.copy()

# Add cancer_type metadata
adata.obs["cancer_type"] = "Breast"

# Extract project_name from barcode prefix (format: ProjectName_barcode)
adata.obs["project_name"] = [b.rsplit("_", 1)[0] if "_" in b else "unknown" for b in all_barcodes]

print(f"Writing to: {OUT_H5AD}")
adata.write_h5ad(OUT_H5AD)
print("Done!")
print(f"Final object: {adata.shape[0]} cells x {adata.shape[1]} genes")
