"""Shared data for the Arm 2 (immune, Spectra vs GeneNMF) appendix figures A.31-A.34.

Inputs written by ch5_gep_programs/genenmf_immune/compare_spectra_vs_genenmf.py and
run_genenmf_per_lineage.R. HC = GeneNMF meta-program whose best Spectra match (top-50 genes)
has Jaccard > 0.1 (same rule as the original figures; >= 0.1 gives the same 29 programs).
"""
import re
import pandas as pd

D = "/groups/ofircohen-group/users/hadasra/temp_objects/genenmf_immune_output"
LINEAGES = ["T-cells", "Myeloid", "B-cell", "NK-cells", "Hematopoietic_Stem_and_Progenitor_cells"]
LABEL = {"T-cells": "T cells", "Myeloid": "Myeloid", "B-cell": "B cells", "NK-cells": "NK cells",
         "Hematopoietic_Stem_and_Progenitor_cells": "HSPCs"}
LCOL = {"T-cells": "#e41a1c", "Myeloid": "#ff7f00", "B-cell": "#984ea3", "NK-cells": "#4daf4a",
        "Hematopoietic_Stem_and_Progenitor_cells": "#377eb8"}
HC_THR = 0.1
HC_RED = "#c00000"

best = pd.read_csv(f"{D}/Spectra_vs_GeneNMF_jaccard_best_match_all_lineages.csv")
full = pd.read_csv(f"{D}/Spectra_vs_GeneNMF_jaccard_full_all_lineages.csv")
best["hc"] = best.best_jaccard > HC_THR


def mp_num(mp):
    return int(re.sub(r"\D", "", mp))


MPS = [f"MP{i}" for i in range(1, 20)]


def coverage():
    rows = []
    for lin in LINEAGES:
        m = pd.read_csv(f"{D}/{lin}/{lin}_genenmf_MP_metrics.csv", index_col=0)
        rows.append(m.sampleCoverage.rename(lin))
    return pd.concat(rows, axis=1).reindex(MPS)
