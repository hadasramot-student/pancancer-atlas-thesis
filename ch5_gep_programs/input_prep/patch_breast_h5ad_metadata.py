#!/usr/bin/env python
# patch_breast_h5ad_metadata.py
#
# Joins metadata from breast_epithelial_metadata.csv onto
# Epithelial_Breast_merged_full.h5ad and rewrites the file.

import anndata as ad
import pandas as pd
import numpy as np

H5AD = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Breast/Epithelial_Breast_merged_full.h5ad"
META = "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query/Breast/breast_epithelial_metadata.csv"

print(f"Loading metadata CSV...")
meta = pd.read_csv(META, index_col="cell_barcode", low_memory=False)
print(f"  {len(meta)} cells, columns: {list(meta.columns)}")

# Convert all object columns to string, replacing NaN with empty string
for col in meta.select_dtypes(include="object").columns:
    meta[col] = meta[col].fillna("").astype(str)

print(f"Loading h5ad...")
adata = ad.read_h5ad(H5AD)
print(f"  {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  obs columns before: {list(adata.obs.columns)}")

# Join metadata by barcode
n_before = adata.n_obs
common = meta.index.intersection(adata.obs_names)
print(f"\nBarcodes matched: {len(common)} / {adata.n_obs}")

# Drop reserved h5ad column names
reserved = {"_index"}
for col in meta.columns:
    if col in reserved:
        print(f"  Skipping reserved column: {col}")
        continue
    adata.obs[col] = meta.reindex(adata.obs_names)[col].values

print(f"obs columns after: {list(adata.obs.columns)}")

# Check scimilarity_celltype coverage
if "scimilarity_celltype" in adata.obs.columns:
    missing = adata.obs["scimilarity_celltype"].isna().sum()
    print(f"\nscimilarity_celltype: {adata.n_obs - missing} cells labeled, {missing} missing")
    print("\nTop 15 values:")
    print(adata.obs["scimilarity_celltype"].value_counts().head(15).to_string())

print(f"\nWriting back to: {H5AD}")
adata.write_h5ad(H5AD)
print("Done.")
