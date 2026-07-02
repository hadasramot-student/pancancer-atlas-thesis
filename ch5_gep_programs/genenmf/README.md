# GeneNMF meta-programs (per cancer type)

Cross-validation of cNMF programs using **GeneNMF** (Carmona lab 2023,
[carmonalab/GeneNMF](https://github.com/carmonalab/GeneNMF)). GeneNMF runs NMF
*per sample* and clusters the resulting programs across samples into
**meta-programs (MPs)** — cross-sample stable transcriptional programs within
each cancer type. Each cancer's MPs are then compared to its cNMF K=20
consensus programs by Jaccard overlap of the top-50 genes.

## Files

| File | Purpose |
|------|---------|
| `install_genenmf.R`            | One-off `install_github("carmonalab/GeneNMF")` into the user library. |
| `run_genenmf_per_cancer.R`     | Full pipeline for **one** cancer type (run NMF → meta-programs → top-50 genes → Jaccard vs cNMF). |
| `submit_genenmf_all_cancers.sh`| Submits one independent SLURM job per cancer, in parallel. |

## Run

All jobs run inside the apptainer Rstudio image via `~/rstudio/rstudio.sbatch`.
If your apptainer bind mount maps the host filesystem to a different path
inside the container (as on the original cluster this was developed on, where
`/mnt/new_groups/ofircohen-group` on the host is bound to `/ofir_grp` inside
the container), set `ATLAS_GROUP_ROOT_HOST` / `ATLAS_GROUP_ROOT_CTR` (see
`submit_genenmf_all_cancers.sh`) so paths resolve correctly on both sides.

```bash
# from the HPC login node, with scripts placed under this repo's genenmf/

# 1. Install GeneNMF once, then submit all cancers chained after it:
bash submit_genenmf_all_cancers.sh --install

# (or, if GeneNMF is already installed, just submit:)
bash submit_genenmf_all_cancers.sh

# inspect what would be submitted:
bash submit_genenmf_all_cancers.sh --dry-run

# subset:
bash submit_genenmf_all_cancers.sh --cancers "Lung Skin"
```

## Cancer types

The 9 epithelial cancers (Lung, Ovarian, Kidney, Liver_Biliary, Pancreas,
Prostate, Sarcoma, Skin, Colorectal) **plus Brain**. Brain's tumour object is
`other_Brain` filtered to `cell_type == "Malignant"` (49,688 cells), which was
already saved as `Epithelial_Brain_merged_full_scanvi_final.RDS`, so it loads
identically to the others. **Breast is excluded** until its cNMF job finishes —
add `Breast` to the `CANCERS=(...)` array in the submit script once
`cnmf_output/Breast/` exists.

## Per-sample split

GeneNMF needs a per-sample grouping. The R script auto-selects
`sample` → `orig.ident` → `project_name` (first present), or pass
`--sample_col`. Samples with `< --min_cells` (default 50) cells are dropped;
at least `--min_samples` (default 5) must remain.

## Outputs — `$ATLAS_DATA_DIR/genenmf_output/<Cancer>/`

| File | Contents |
|------|----------|
| `<Cancer>_genenmf_programs.rds`       | Raw per-sample NMF programs (`multiNMF`). |
| `<Cancer>_genenmf_metaprograms.rds`   | Full `getMetaPrograms` object. |
| `<Cancer>_genenmf_MP_top50_genes.csv` | Long table: cancer, meta_program, rank, gene. |
| `<Cancer>_genenmf_MP_metrics.csv`     | Per-MP QC (sample coverage, silhouette, …). |
| `<Cancer>_cnmf_top50_genes.csv`       | Top-50 genes per cNMF K=20 program (from `gene_spectra_score`). |
| `<Cancer>_jaccard_MP_vs_cNMF.csv`     | Jaccard matrix, MPs × cNMF programs. |
| `<Cancer>_highconfidence_programs.csv`| Per MP: best cNMF match, best Jaccard, `recovered_by_both` (Jaccard > 0.1). |

## Key parameters (defaults)

`--n_mp 20` (match cNMF K), `--k_range 4:9` (per-sample NMF), `--nfeatures 2000`,
`--top_n 50`, `--jaccard_thr 0.1`, `--min_cells 50`. Edit the matching variables
at the top of `submit_genenmf_all_cancers.sh` to change them for all cancers.
