# Spectra on immune lineages (Aim 3.2, GEP method 3)

Prior-guided topic modeling (Spectra) on the immune lineages of the pan-cancer
atlas, run **per immune lineage** (T-cells, Myeloid, B-cell, NK-cells, HSPC),
each merged across cancer types. Outputs program loadings per cell to compare
against cNMF / GeneNMF in Task 3.

## Pipeline (one SLURM dependency chain per lineage)

| Step | Script | Where | Notes |
|------|--------|-------|-------|
| A. RDS → h5ad | `seurat_to_h5ad_spectra.R` | R container, per cancer file | huge lineages (T/Myeloid) subsampled to 60k cells/file |
| B. merge + subsample | `merge_lineage_h5ad.py` | scanvi-env CPU | stratified by cancer×cell-type; T/Myeloid capped at 200k |
| C. Hallmark export | `export_hallmark_gmt.R` | R container (once) | MSigDB collection H via `msigdbr` → `priors/hallmark.json` |
| D. build priors | `build_spectra_priors.py` | scanvi-env CPU | nested dict: global = Hallmark + TME-global; per-celltype = compartment TME sets |
| E. run Spectra | `run_spectra.py` | scanvi-env **GPU** | `Spectra_gpu.est_spectra`, RTX 3090 |

Launch everything:
```bash
cd ch5_gep_programs/spectra   # from the repo root
ATLAS_DATA_DIR=/path/to/temp_objects bash submit_spectra.sh                 # all 5 immune lineages
bash submit_spectra.sh --dry-run       # print job graph only
bash submit_spectra.sh --lineages "B-cell NK-cells"   # subset
```
Monitor: `squeue -u $USER | grep sp_` · logs in `logs/` · outputs in
`$ATLAS_DATA_DIR/clustered_lineages/spectra_out/`.

## Key configuration decisions
- **cell_type_key = `scimilarity_celltype`** — only column with clean,
  fine-grained *within-lineage* immune subtypes (controlled Cell-Ontology
  vocabulary); same key used for scANVI labels. Lineage-name columns
  (`Consensus_Cell_Type`, `maj2_cell`) are constant per object → useless here.
- **Features are gene symbols** (verified) → match Hallmark/TME priors directly.
- **Subsampling**: small lineages (B-cell, NK, HSPC) run in full; T-cells and
  Myeloid subsampled (per-file 60k, merged 200k), stratified by cancer×cell-type.
- **Priors**: `tme_gene_sets.json` is literature-curated — **review/edit it**.
  Compartment block (T/Myeloid/B/NK/HSPC) is attached to every cell type in the
  lineage; `global` (Hallmark + global TME) applies to all cells.

## Outputs per lineage (`<lineage>_spectra_*`)
- `_cell_scores.csv` — cells × programs (the loadings for the cNMF comparison)
- `_gene_factors.csv` — programs × genes (gene weights)
- `_markers.csv` — top genes per program
- `_spectra.h5ad` — AnnData with `obsm['SPECTRA_cell_scores']`, `uns['SPECTRA_factors']`
- `_model.pickle` — trained model

## Gotchas (already handled)
- `Spectra_gpu.est_spectra` crashes indexing a sparse `X` with a pandas Series
  when `use_highly_variable=True`. `run_spectra.py` pre-subsets to the vocab,
  densifies X, and passes `use_highly_variable=False`.
- Use `temp_objects/seurat_to_h5ad.R` semantics, **not** the repo copy at
  `create_sc_object_workflow/seurat_to_h5ad.R` (that one has a debug `quit()`).
- `apptainer` only exists on compute nodes → all R steps go through `sbatch`.
