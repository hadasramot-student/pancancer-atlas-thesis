# Configuration

- **`params.yaml`** — QC thresholds, the 13 cancer types (pointer), the 8
  lineages (pointer), clustering resolutions, and per-method GEP parameters
  (cNMF K-sweep, GeneNMF, Spectra, IDS). Reference values; most scripts
  expose the same numbers as CLI flags with matching defaults rather than
  reading this file directly.
- **`cancer_types.txt`** — the 13 cancer types in the atlas, one per line.
- **`lineages.txt`** — the 8 lineage compartments used for Ch.4 subclustering
  and (for the immune 5 of them) Ch.5 Spectra/GeneNMF.
- **`SC.3CA.Tirosh.project_file_index.tsv`** (not included — see below) —
  the 3CA project index used by Ch.1/Ch.2 3CA scripts.
- **`lab_utils/`** (not included — see below) — lab-internal R generic
  utilities (`Utils.R`, `Utils.MS.R`) sourced by `common/functions.R` and
  `ch1_aggregation_qc/02_ingest_3ca.R`.

## Environment variables

Scripts read paths from environment variables (with repo-relative
fallbacks) instead of hardcoding a filesystem layout. Set the ones relevant
to the stage you're running:

| Variable | Used by | Meaning |
|---|---|---|
| `ATLAS_DATA_DIR` | most scripts | Root of all intermediate/output data (the old `temp_objects/`) |
| `ATLAS_REPO_DIR` | orchestrators, SLURM submit scripts | Path to this repo checkout |
| `R_USER_LIBRARY` | R scripts | Personal R library path (`.libPaths()` prepend) |
| `SCIMILARITY_MODEL_PATH` | `ch2_annotation/my_scimilarity.py`, `ch4_subclustering/run_scimilarity.py`, `ch4_subclustering/annotate_lineages.R` | Local Scimilarity model checkpoint directory |
| `SC3CA_PROJECT_INDEX` | Ch.1/Ch.2 3CA scripts | Path to `SC.3CA.Tirosh.project_file_index.tsv` |
| `DROPBOX_ROOT` / `LAB_UTILS_DIR` | `common/functions.R`, `02_ingest_3ca.R` | Directory containing the lab-internal `Utils.R`/`Utils.MS.R` |
| `RAW_SCBASECOUNT_DIR` / `RAW_BREASTSCAPE_DIR` | Ch.1 batch drivers | Raw per-source input directories |
| `APPTAINER_R_IMAGE` | Ch.4/Ch.5 SLURM submit scripts | Path to the Rstudio apptainer `.sif` image |
| `ATLAS_GROUP_ROOT_HOST` / `ATLAS_GROUP_ROOT_CTR` | GeneNMF/Spectra SLURM submit scripts | If your apptainer bind-mounts the data dir to a different path inside the container than outside, set the host and container roots so R (inside) and Python/shell (outside) resolve the same files |

Reproduction from scratch requires the original source data (ScBaseCount,
3CA, Breastscape) and compute (SLURM + apptainer + a GPU node for scANVI and
Spectra) — see the root README for the full walkthrough.
