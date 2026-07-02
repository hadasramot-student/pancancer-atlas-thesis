# Pan-Cancer Single-Cell Atlas — Analysis Pipeline

Analysis code for a pan-cancer single-cell RNA-seq atlas spanning 13 cancer
types: aggregation & QC, multi-method cell-type annotation with consensus
voting, batch correction & evaluation, lineage subclustering, and NMF-based
transcriptional program (GEP) discovery. This code accompanies the MSc
thesis *[thesis title]* — see **Code Availability** in the thesis for the
link back to this repository.

Mixed R (Seurat, SingleR, SCType, ComBat/sva, batchelor/MNN, harmony,
GeneNMF, fgsea, msigdbr) and Python (scANVI/scvi-tools, Scimilarity, cNMF,
Spectra) pipeline, developed for a SLURM cluster with apptainer containers
(some steps as job arrays; Spectra requires a GPU).

## Pipeline overview

```
Ch.1 Aggregation & QC          Ch.2 Annotation              Ch.3 Batch correction        Ch.4 Subclustering            Ch.5 GEP programs
─────────────────────          ───────────────              ──────────────────────       ──────────────────            ─────────────────
ScBaseCount ┐                  SingleR                      ComBat                       split into 8 lineages         tumor/malignant cell
3CA         ├─▶ QC filter  ──▶ SCType            ──▶        MNN               ──▶         (T/Myeloid/Epithelial/       selection (per cancer)
Breastscape ┘   + schema       Scimilarity                  Harmony                       B/Fibroblast/Endothelial/         │
                                original-metadata            scANVI                        NK/HSPC)                   ┌────┼────┐
                                cross-validation (3CA)             │                             │                    cNMF GeneNMF Spectra
                                       │                     LISI/kBET/F1/                  per-lineage PCA→          (epithelial) (immune)
                                consensus vote ──▶            graph connectivity            UMAP→Louvain (5 res)            │
                                Consensus_Cell_Type                                               │                   Hallmark priors /
                                                                                             Scimilarity cluster            IDS / Jaccard
                                                                                             annotation + markers      cross-method comparison
                                                                                             + Hallmark (fgsea)
```

Each stage's scripts live in the correspondingly-named directory below.
`common/functions.R` is the shared R helper library for Ch.1/Ch.2/Ch.4
(QC/normalization, SingleR/SCType/Scimilarity annotation, consensus voting,
lineage-splitting merge). `config/` holds the shared parameters (QC
thresholds, cancer types, lineages, K values, clustering resolutions).

## Repository layout

```
ch1_aggregation_qc/   Data retrieval (ScBaseCount/3CA/Breastscape), QC, metadata
                       standardization, atlas composition figures/tables
ch2_annotation/        SingleR, SCType, Scimilarity, 3CA cross-validation,
                       consensus voting, method-agreement figures
ch3_batch_correction/ ComBat, MNN, Harmony, scANVI, LISI/kBET/F1/graph
                       connectivity evaluation, method-comparison figures
ch4_subclustering/     Lineage splitting, per-lineage PCA/UMAP/clustering
                       (SLURM job array), Scimilarity cluster annotation,
                       marker identification, Hallmark enrichment
ch5_gep_programs/       Tumor cell selection, cNMF, GeneNMF (epithelial +
    cnmf/                immune), Spectra (+ TME priors, GPU), IDS,
    genenmf/              cross-method comparison (Jaccard + IDS)
    genenmf_immune/
    ids_analysis/
    spectra/
common/                 Shared R helper library (functions.R)
config/                 Shared parameters (QC thresholds, cancer types,
                        lineages, K values, resolutions) — see config/README.md
```

## Reproduction

1. **Environment.** R: see `r_packages.txt` (run `dump_session_info.R`
   inside your R environment to capture exact versions / generate a
   `renv.lock`). Python: `conda env create -f environment.yml` for the main
   environment (scANVI/cNMF/GeneNMF/IDS/Spectra); Scimilarity needs its own
   isolated environment, see `requirements-scimilarity.txt`.
2. **Environment variables.** Every script takes paths via CLI args or
   environment variables rather than hardcoding a filesystem layout — set
   `ATLAS_DATA_DIR` (intermediate/output data root) and the other variables
   documented in `config/README.md` before running anything.
3. **Source data.** ScBaseCount, 3CA (Tirosh lab Curated Cancer Cell Atlas),
   and Breastscape raw/pre-processed objects are not redistributed here —
   obtain them from their respective sources and point the Ch.1 ingestion
   scripts at them.
4. **Run order.** Follow the chapter numbering: Ch.1 → Ch.2 (per-object,
   via `ch1_aggregation_qc/00_main_dispatch.R`) → merge into the atlas → Ch.3
   (per lineage x cancer-type object) → Ch.4 (per lineage, split from the
   Ch.2 consensus types) → Ch.5 (per cancer type for cNMF/GeneNMF, per
   immune lineage for Spectra/GeneNMF-immune). Each script's header
   documents its specific inputs/outputs and a runnable usage example; SLURM
   submission wrappers (`submit_*.sh`) exist for the array-style steps
   (Ch.4 clustering, Ch.5 cNMF/GeneNMF/Spectra/IDS).

This repository intentionally excludes large data products (raw counts,
`.RDS`/`.h5ad` objects, model checkpoints) — see `.gitignore`. Only scripts,
configs, small reference tables, and documentation are tracked.

## Link to thesis

[thesis title / link / DOI — fill in before publishing]
