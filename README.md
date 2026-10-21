# Pan-Cancer Single-Cell Atlas — Analysis Pipeline

Analysis code for a pan-cancer single-cell RNA-seq atlas spanning 13 cancer
types: aggregation & QC, multi-method cell-type annotation with consensus
voting, batch correction & evaluation, lineage subclustering, and NMF-based
transcriptional program (GEP) discovery. This code accompanies the MSc
thesis *Pan-Cancer Single-cell Atlas of Tumors and their Microenvironment:
Harmonization, Cell-Type Curation, and Transcriptional Programs
Characterization* (Hadas Ramot, Ben-Gurion University of the Negev, 2026) —
see **Appendix C (Code Availability)** in the thesis for the link back to
this repository.

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
                                                                                             + Hallmark (fgsea)              │
                                                                                                                        Hallmark GSEA + Gavish/
                                                                                                                        Tirosh meta-program match
                                                                                                                             │
                                                                                                                        TCGA survival
                                                                                                                        (BRCA/KIRC outcome)
```

Each stage's scripts live in the correspondingly-named directory below (stage labels
Ch.1–Ch.5 = thesis Chapters 2–6, see *Chapter mapping*).
`common/functions.R` is the shared R helper library for Ch.1/Ch.2/Ch.4
(QC/normalization, SingleR/SCType/Scimilarity annotation, consensus voting,
lineage-splitting merge). `config/` holds the shared parameters (QC
thresholds, cancer types, lineages, K values, clustering resolutions).

## Chapter mapping

The folder names and the `Ch.N` labels inside the script headers follow the
original working numbering, which is one lower than the chapter numbers in the
submitted thesis:

| Repository folder | Script headers | Thesis chapter |
|---|---|---|
| `ch1_aggregation_qc/` | Ch.1 | **Chapter 2** — Data aggregation and QC (3CA, Breastscape, ScBaseCount; nFeature_RNA ≥ 200; removal of 5 duplicate Breast studies; 8,107,665 cells, 198 projects) |
| `ch2_annotation/` | Ch.2 | **Chapter 3** — Cell-type annotation (SingleR, SCType, Scimilarity, consensus meta-types) |
| `ch3_batch_correction/` | Ch.3 | **Chapter 4** — Batch correction (ComBat, MNN, Harmony, scANVI; LISI, kBET, graph connectivity, F1; epithelial benchmark; 200,000-cell atlas scANVI embedding) |
| `ch4_subclustering/` | Ch.4 | **Chapter 5** — Subclustering (75 lineage × cancer-type objects, Louvain at several resolutions, Scimilarity cluster labels, FindMarkers, fgsea Hallmark) |
| `ch5_gep_programs/` | Ch.5 | **Chapter 6** — Gene expression programs (cNMF K = 20, GeneNMF, IDS, Spectra, Hallmark, cNMF–GeneNMF Jaccard + null threshold 0.0538, consensus meta-programs, Gavish et al. comparison, TCGA survival) |

Figure scripts written for the thesis revision follow the same rule (e.g.
`ch5_fig1` = thesis Figure 6.1).

## Repository layout

```
ch1_aggregation_qc/            (thesis Ch.2)
    00_main_dispatch.R, 01-03b_*   ingest + QC per source (ScBaseCount / 3CA / Breastscape)
    create_tsvs.R, rbind_tsvs*.R   per-object metadata -> atlas metadata table (all_cells)
    rebuild_all_cells.R            rebuild all_cells from the *_step4.RDS.tsv files
    fix_59219_gap.R                fix the 59,219-cell gap -> "gapfix" all_cells table
    duplicate_studies/             cross-repository duplicate-study detection (5 Breast studies)
    project_count_audit/           project counts on the deduplicated atlas (198 projects)
    repository_inventory/          raw cell counts per repository and cancer type
    atlas_*.py                     atlas composition figures
    thesis_figures/                Figures 2.1 and A.1 (cells per cancer type after / before QC)
    checks/                        one-off sanity checks kept for traceability
ch2_annotation/                (thesis Ch.3)
    04_object_annotation.R         SingleR + SCType (+ 3CA original labels)
    my_scimilarity.py              Scimilarity annotation
    05_consensus_annotation.R      consensus vote -> Consensus_Cell_Type
    agreement_plots_methods.R      method agreement figures
    table4_agreement.R             Table 4: per-meta-type inter-method agreement
    validation_3ca_author_labels/  consensus vs held-out 3CA author labels (88.8% / 70.6%)
    thesis_figures/fig3_1_composition/  consensus cell-type composition per cancer type
    checks/
ch3_batch_correction/          (thesis Ch.4)
    merge_by_cancertype_celltype.R + submit_all_cancertypes.sh
                                   build the lineage x cancer-type objects from all_cells
    create_dev_subset.R            benchmark subsample (batch x subtype stratified)
    batch1-4_*.R                   CCA / fastMNN / ComBat / Harmony
    seurat_to_h5ad.R, scanvi_batch_correction.py, export_scanvi_to_csv.py,
    scanvi_into_seurat.R           scANVI round trip
    batch_correction_validation_metrics.R, compare_batch_methods.R,
    batch_correction_methods_summary.R   LISI / kBET / graph connectivity / F1
    pipeline/                      SLURM chains that ran the steps above (benchmark + all objects)
    epithelial_benchmark/          epithelial batch-correction benchmark
    benchmark_sensitivity/         duplicate-study overlap, dedup re-run, seed control
    atlas_200k_scanvi/             200,000-cell atlas scANVI embedding + UMAP figures
    plot_*.R, umap_batch_comparison.R   figures
    thesis_figures/                Figures 4.1, 4.2, 4.5 and the study-of-origin panel
    checks/
ch4_subclustering/             (thesis Ch.5)
    merge_by_cell_type_from_all_cells.R  lineage split
    elbow_check.R                  PCA dimensionality check
    clustering_jobs/               cluster_lineage.R + submit_clustering_per_cancertype.sh
                                   (75 objects, Louvain at 0.25/0.30/0.50/0.75/1.0/2.0);
                                   run_clustering.R + submit_all_clustering.sh (per-lineage run)
    annotate_lineages.R, run_scimilarity.py + submit_annotation_jobs.sh   Scimilarity cluster labels
    find_markers.R + submit_markers_jobs.sh                               FindMarkers
    gsea_hallmarks.R + submit_gsea_hallmarks.sh                           fgsea Hallmark
    plot_umaps.R + submit_plot_umaps.sh, plots_for_clusters.R, umap_by_scimilarity.R, fig19/fig23_*
    epithelial_breast_subsample/   450,000-cell stratified subsample of Epithelial_Breast
    report_epithelial_brain.R      Epithelial_Brain report
    audits/                        Kidney c5/c9 ambient check, cell-count cross-check
    thesis_figures/umap_redraws/   print-size UMAP redraws
    checks/
ch5_gep_programs/              (thesis Ch.6)
    input_prep/, filter_brain_malignant.R   epithelial / tumour h5ad inputs (Breast chunks, Brain)
    cnmf/                          cNMF (run_cnmf_per_cancer.py + submit_cnmf_all_cancers.sh, K = 20)
        k_selection/               K = 10-30 stability sweep on Lung (choice of K = 20)
    genenmf/, genenmf_immune/      GeneNMF (epithelial per cancer type; immune per lineage)
    spectra/                       Spectra (GPU) with Hallmark + TME priors
    ids_analysis/                  IDS + cross-method comparison
    jaccard_threshold/             null distribution for the cNMF-GeneNMF Jaccard threshold (0.0538)
    hallmark_enrichment/           Hallmark enrichment + Gavish et al. (2023) MP comparison
    meta_programs/                 consensus meta-programs, threshold sensitivity, Table B97 audit
    contamination_audits/          cMP6 (myeloid), cMP8 (fibroblast), cMP10 (MHC-II) checks
    tcga_survival/                 TCGA-BRCA/KIRC survival (run_tcga_survival.R)
        kirc_mp11/                 KIRC MP11 overall survival on consolidated primary tumours
        purity_adjustment/         purity-adjusted Cox models (03_purity_survival_v2.R = final)
    thesis_figures/                Figures 6.1, 6.3, 6.4, 6.5, A.29, A.30, A.38
    checks/
common/functions.R             shared R helpers
config/                        shared parameters (see config/README.md)
```

`checks/` folders hold one-off verification scripts. They are not part of the
main pipeline, but are kept because they may have produced numbers quoted in
the thesis.

## Environments

- **R** (Seurat v5, SingleR, SCType, sva/ComBat, batchelor, harmony, GeneNMF,
  fgsea, msigdbr, survival, ...): `r_packages.txt`. Most R steps ran inside
  the lab's RStudio apptainer image (`APPTAINER_R_IMAGE`). Run
  `dump_session_info.R` inside your R environment to record exact versions.
- **Python** (scanpy, scvi-tools/scANVI, cNMF, Spectra, IDS):
  `conda env create -f environment.yml` (or `pip install -r requirements.txt`).
- **Scimilarity** needs its own isolated environment:
  `requirements-scimilarity.txt`.
- scANVI and Spectra need a GPU.

## Run order

Each script's header documents its inputs, outputs and usage. In order:

**Chapter 2 — aggregation & QC (`ch1_aggregation_qc/`)**
1. `01_ingest_scbasecount.R`, `02b_run_3ca_batch.R` (→ `02_ingest_3ca.R`),
   `03b_run_breastscape_batch.R` (→ `03_ingest_breastscape.R`), each calling
   `00_main_dispatch.R` per object (ingest, QC with nFeature_RNA ≥ 200, then
   the Chapter 3 annotation steps below).
2. `create_tsvs.R` → `rbind_tsvs.R` / `rbind_tsvs_breastscape.R` →
   `rebuild_all_cells.R` → `fix_59219_gap.R` (atlas metadata table).
3. `duplicate_studies/*` (identify the 5 duplicate Breast studies),
   `project_count_audit/` (`dump_projects.R` → `build_tables.py`),
   `repository_inventory/`, `atlas_cancer_distribution.py`,
   `atlas_composition_barchart.py`.

**Chapter 3 — annotation (`ch2_annotation/`)**, per object, called from `00_main_dispatch.R`
1. `04_object_annotation.R` → `my_scimilarity.py` → `05_consensus_annotation.R`.
2. `agreement_plots_methods.R`, `table4_agreement.R`.
3. `validation_3ca_author_labels/01_inspect.R` → `02_validate.R` → `03_derived.R` → `04_exclusions.R`.
4. `thesis_figures/fig3_1_composition/count_composition.R` → `make_ch2_fig1.py`.

**Chapter 4 — batch correction (`ch3_batch_correction/`)**
1. `merge_by_cancertype_celltype.R` (via `submit_all_cancertypes.sh`) → lineage × cancer-type objects.
2. Benchmark: `pipeline/run_batch_correction_pipeline.sh` (or `run_batch_correction_no_cca.sh`) runs
   `create_dev_subset.R` → `batch1_cca_integration.R` → `batch2_fastmnn.R` → `batch3_combat.R` →
   `batch4_harmony.R` → `seurat_to_h5ad.R` → `scanvi_batch_correction.py` / `export_scanvi_to_csv.py`
   → `scanvi_into_seurat.R` → `batch_correction_validation_metrics.R`; then
   `compare_batch_methods.R`, `batch_correction_methods_summary.R` and the `plot_*.R` scripts.
3. scANVI on all objects: `pipeline/run_scanvi_all_cancertypes.sh` / `run_scanvi_full.sh`
   (`check_batch_counts.R` first).
4. `epithelial_benchmark/00_probe_metadata.R` → `20_submit_all.sh` (→ `10_add_batch_methods.R`).
5. `benchmark_sensitivity/ch3_dup_overlap_check.R`, `ch3_dedup_rerun.R`, `ch3_seed_control.R`.
6. `atlas_200k_scanvi/00_probe.R` → `00b_subsample_sizing.R` → `01_subsample.R` → `02_submit.sh`
   (→ `02_task.sh` → `02_assemble_counts.R`) → `03_merge_h5ad.R` → `04_scanvi.py` → `05_figures.py`.

**Chapter 5 — subclustering (`ch4_subclustering/`)**
1. `elbow_check.R`; `clustering_jobs/submit_clustering_per_cancertype.sh` (→ `cluster_lineage.R`).
   Epithelial_Breast: `epithelial_breast_subsample/run_epithelial_breast_subsample_pipeline.sh`.
2. `submit_annotation_jobs.sh` (→ `annotate_lineages.R`, `run_scimilarity.py`) →
   `submit_markers_jobs.sh` (→ `find_markers.R`) → `submit_gsea_hallmarks.sh` (→ `gsea_hallmarks.R`).
3. `submit_plot_umaps.sh` (→ `plot_umaps.R`), `plots_for_clusters.R`, `umap_by_scimilarity.R`,
   `fig19_*`, `fig23_*`, `report_epithelial_brain.R`, `audits/*`,
   `thesis_figures/umap_redraws/submit_extract.sh` (→ `extract_umap_data.R`) → `plot_umaps_print.R`.

**Chapter 6 — gene expression programs (`ch5_gep_programs/`)**
1. Inputs: `input_prep/*`, `filter_brain_malignant.R`.
2. cNMF: `cnmf/k_selection/submit_test_cnmf.sh` (→ `test_run_cnmf.py`, `run_consensus_test.py`) for
   the K sweep, then `cnmf/submit_cnmf_all_cancers.sh` (→ `run_cnmf_per_cancer.py`, K = 20),
   `cnmf/plot_stability.py`. (`01_export_counts_for_cnmf.R` / `02_run_cnmf.py` /
   `03_submit_cnmf.sh` are the earlier atlas-wide entry point.)
3. GeneNMF: `genenmf/submit_genenmf_all_cancers.sh` (→ `run_genenmf_per_cancer.R`);
   immune: `genenmf_immune/submit_genenmf_immune.sh`.
4. Spectra: `spectra/submit_spectra.sh` (see `spectra/README.md`).
5. IDS: `ids_analysis/submit_ids_all_cancers.sh` (→ `run_ids.py`), `plot_ids_crossmethod.py`,
   `build_hc85_ids_summary.py`, `jaccard_ids_stats.py`.
6. Jaccard threshold: `jaccard_threshold/null_derived_threshold.py`, `extend_sensitivity.py`.
7. Hallmark + Gavish: `hallmark_enrichment/inspect_and_export_mp_list.R` → `run_cnmf_hallmark.py`
   → `run_gavish_jaccard.py` → `make_gavish_mp_dotplot.py`.
8. Meta-programs: `meta_programs/meta_program_catalogue.py`, `cmp_threshold_sensitivity.py`,
   `submit_tableB97_contamination_audit.sh` (→ `tableB97_contamination_audit.py`);
   `contamination_audits/*`.
9. TCGA: `tcga_survival/submit_tcga_survival.sh` (→ `run_tcga_survival.R`);
   `tcga_survival/kirc_mp11/01_consolidated_pt_survival.R`;
   `tcga_survival/purity_adjustment/01_inspect_kidney_epi.R` → `01b_reanalyze_labels.R` →
   `02_purity_survival.R` → `03_purity_survival_v2.R` → `04_figures.R`.
10. Figures: `thesis_figures/*` (see the index below).

## Figure / table → script index

| Thesis item | Script |
|---|---|
| Figure 2.1 — cells per cancer type after QC; Figure A.1 — before QC | `ch1_aggregation_qc/thesis_figures/make_fig2_1_and_A1_cell_bars.py` (²) |
| Figure 3.1 — consensus cell-type composition per cancer type | `ch2_annotation/thesis_figures/fig3_1_composition/make_ch2_fig1.py` |
| Table 4 — inter-method agreement per meta-type | `ch2_annotation/table4_agreement.R` |
| 3CA author-label validation (88.8% / 70.6%) | `ch2_annotation/validation_3ca_author_labels/` |
| 5 duplicate Breast studies, 8,107,665 cells | `ch1_aggregation_qc/duplicate_studies/` |
| 198 projects | `ch1_aggregation_qc/project_count_audit/` |
| Batch-correction benchmark metrics and figures | `ch3_batch_correction/batch_correction_validation_metrics.R`, `compare_batch_methods.R`, `plot_*.R`, `thesis_figures/plot_fig03_study_only.R` |
| Benchmark duplicate overlap / seed control | `ch3_batch_correction/benchmark_sensitivity/` |
| Epithelial benchmark | `ch3_batch_correction/epithelial_benchmark/` |
| Atlas-wide 200,000-cell scANVI UMAPs | `ch3_batch_correction/atlas_200k_scanvi/05_figures.py` |
| Figures 4.1 and 4.2 — Uncorrected vs scANVI UMAP panels (from Figures A.15 / A.16) | `ch3_batch_correction/thesis_figures/make_fig4_1_4_2_two_panel_umaps.py` (²) |
| Figure 4.5 — Uncorrected (PCA) vs scANVI, by source repository | `ch3_batch_correction/thesis_figures/make_fig4_5_pca_vs_scanvi.py` (²) |
| Subclustering UMAPs (print redraws) | `ch4_subclustering/thesis_figures/umap_redraws/plot_umaps_print.R` |
| Subclustering marker dot plots | `ch4_subclustering/fig19_tcells_liver_biliary_dotplot.R`, `fig23_epithelial_kidney_dotplot.R` |
| Figure 6.1 — cNMF Hallmark enrichment (BH-corrected) | `ch5_gep_programs/thesis_figures/fig6_1_A29_A30_A38/01_hc_pairs_and_hallmark_bh.py` → `02_fig6_1_hallmark_dotplot.py` |
| Figure A.29 — cNMF usage heatmap | `.../fig6_1_A29_A30_A38/03_figA29_cnmf_usage_heatmap.py` |
| Figure A.30 — GeneNMF sample coverage | `.../fig6_1_A29_A30_A38/04_figA30_genenmf_coverage_heatmap.py` |
| Figure A.38 — HC program associations | `.../fig6_1_A29_A30_A38/05_figA38_hc_program_associations.py` (adapted from `thesis_figures/make_fig44_hc_program_associations.py`) |
| HC pairs, Kidney redundancy families, numbers in the text | `.../fig6_1_A29_A30_A38/01_hc_pairs_and_hallmark_bh.py`, `06_report_numbers.py` |
| Jaccard null threshold 0.0538 | `ch5_gep_programs/jaccard_threshold/null_derived_threshold.py` |
| Choice of K = 20 | `ch5_gep_programs/cnmf/k_selection/`, `cnmf/plot_stability.py` |
| Consensus meta-programs, Table B97 contamination audit | `ch5_gep_programs/meta_programs/` |
| Figure 6.3 — median IDS heatmap | `ch5_gep_programs/thesis_figures/make_fig6_3_ids_heatmap.py` (²) |
| Figure 6.4 — KIRC MP11 survival; Part 3 Cox models | `ch5_gep_programs/thesis_figures/fig6_4_fig6_5/ch5_fig4_and_part3.R` |
| Figure 6.5 — Gavish et al. meta-program comparison | `ch5_gep_programs/thesis_figures/fig6_4_fig6_5/make_ch5_fig5_gavish_dotplot.py` / `make_ch5_fig5.py` |
| KIRC MP11 purity-adjusted Cox | `ch5_gep_programs/tcga_survival/purity_adjustment/03_purity_survival_v2.R` |

(²) Made outside the cluster during the thesis revision (October 2026; Python 3 +
matplotlib + Pillow). These redraw figures from numbers or images that already
existed and do no new analysis: Figures 2.1 / A.1 and 6.3 use values typed in
from the original figure labels (Breast post-QC corrected to 59 projects after
the project-count audit); Figures 4.1, 4.2 and 4.5 crop panels from the
original images (Figures A.15 / A.16 and the 30 September Figure 4.5 image)
and redraw titles and legends.

## External inputs (not included)

- **Single-cell sources:** ScBaseCount, 3CA (Tirosh lab Curated Cancer Cell
  Atlas) and Breastscape. Obtain them from their providers and point the
  Chapter 2 ingestion scripts at them.
- **TCGA:** GDC STAR gene-count files (`tpm_unstranded`) for TCGA-KIRC and
  TCGA-BRCA, and the TCGA Clinical Data Resource (TCGA-CDR; Liu et al., 2018,
  *Cell*) for survival endpoints. Patient-level clinical tables are not
  redistributed here.
- **Gavish et al. (2023)** meta-programs (*Nature* 618:598–606): the
  `MP_list.RDS` file from the Tirosh lab repository
  [tiroshlab/3ca](https://github.com/tiroshlab/3ca) (release v1.0.0,
  `ITH_hallmarks/MPs_distribution/MP_list.RDS`).
  `ch5_gep_programs/hallmark_enrichment/inspect_and_export_mp_list.R` turns
  it into a long gene list. The Tirosh lab code is not copied here.
- **MSigDB Hallmark** gene sets, via the R package `msigdbr`
  (`spectra/export_hallmark_gmt.R`).
- **Scimilarity** pretrained model (`SCIMILARITY_MODEL_PATH`).
- Lab-internal `Utils.R` / `Utils.MS.R`, the 3CA project index and the
  Ensembl lookup table (see `config/README.md`).

## Paths to adapt

The core pipeline scripts read paths from CLI arguments or environment
variables (`ATLAS_DATA_DIR`, ... — see `config/README.md`). The scripts below
were added as they ran on the BGU cluster and still contain **hard-coded
absolute paths** (`/ofir_grp/users/hadasra/...`, which is the apptainer bind
of `/groups/ofircohen-group/...` = `/mnt/new_groups/ofircohen-group/...`, and
`/home/hadasra/R/...` for the user R library). They were left unchanged so the
code matches what produced the thesis results. Edit the path variables at the
top of each script before running:

- `ch1_aggregation_qc/`: `fix_59219_gap.R`, `rebuild_all_cells.R`, `duplicate_studies/*`,
  `project_count_audit/dump_projects.*`, `repository_inventory/*.R`, `checks/*`
- `ch2_annotation/`: `table4_agreement.R`, `validation_3ca_author_labels/*`,
  `thesis_figures/fig3_1_composition/count_composition.*`, `checks/*`
- `ch3_batch_correction/`: `merge_by_cancertype_celltype.R`, `submit_all_cancertypes.sh`,
  `pipeline/*`, `atlas_200k_scanvi/*`, `epithelial_benchmark/*`, `benchmark_sensitivity/*`,
  `thesis_figures/*`, `checks/*`
- `ch4_subclustering/`: `fig19_*`, `fig23_*`, `elbow_check.R`, `submit_*.sh`,
  `clustering_jobs/submit_clustering_per_cancertype.sh`, `epithelial_breast_subsample/*`,
  `audits/*`, `thesis_figures/umap_redraws/submit_extract.sh`, `checks/*`
- `ch5_gep_programs/`: `input_prep/*`, `cnmf/submit_cnmf_all_cancers.sh`, `cnmf/k_selection/*`,
  `hallmark_enrichment/*`, `meta_programs/*`, `jaccard_threshold/*`, `ids_analysis/build_hc85_ids_summary.py`,
  `ids_analysis/jaccard_ids_stats.py`, `contamination_audits/*`, `tcga_survival/kirc_mp11/*`,
  `tcga_survival/purity_adjustment/*`, `thesis_figures/*`, `checks/*`

The revision figure scripts marked (²) in the index above
(`ch1_aggregation_qc/thesis_figures/make_fig2_1_and_A1_cell_bars.py`,
`ch3_batch_correction/thesis_figures/make_fig4_1_4_2_two_panel_umaps.py`,
`ch3_batch_correction/thesis_figures/make_fig4_5_pca_vs_scanvi.py`,
`ch5_gep_programs/thesis_figures/make_fig6_3_ids_heatmap.py`) use **relative**
input/output paths that point to the working copy used during the revision
(e.g. `figs/`, `c5/`, `c8/`, `c13/`, `v123/content/images/`,
`v117docx/word/media/image5.png`). Create those folders, or edit the paths,
before running. The source images that the 4.x scripts crop are not included.

This repository intentionally excludes data products (raw counts, `.RDS` /
`.h5ad` objects, figures, model checkpoints, patient-level clinical data) —
see `.gitignore`. Only scripts, configs, small reference/result tables and
documentation are tracked.

## Thesis and code version

*Pan-Cancer Single-cell Atlas of Tumors and their Microenvironment:
Harmonization, Cell-Type Curation, and Transcriptional Programs
Characterization* — Hadas Ramot, MSc thesis, Ben-Gurion University of the
Negev, 2026.

The code version used for the thesis is the tagged release
[`v1.0-thesis`](https://github.com/hadasramot-student/pancancer-atlas-thesis/releases/tag/v1.0-thesis).
