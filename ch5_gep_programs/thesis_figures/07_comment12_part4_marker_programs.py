#!/usr/bin/env python
# 07_comment12_part4_marker_programs.py
#
# Comment 12, part 4 (read-only, no figures). For five marker-defined programs,
# find in each cancer the cNMF K=20 program(s) whose top-50 genes (same
# definition as script 01) hold >= 2 markers, and look up the target Hallmark
# set(s) in hallmark_all_tests_bh.csv (enriched = q < 0.05 and overlap >= 3).
# Output: part4_marker_programs.csv (stdout is used for REPORT_part4.txt)

import pandas as pd

BASE = "/groups/ofircohen-group/users/hadasra"
CNMF = f"{BASE}/temp_objects/cnmf_output"
OUT = f"{BASE}/thesis_figures/comment12_2026-10-21"
CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney", "Liver_Biliary",
           "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"]
CATS = {
    "Hypoxia": (["VEGFA", "SLC2A1", "LDHA", "PGK1", "DDIT4", "BNIP3", "NDRG1"], ["HYPOXIA"]),
    "MHC-II/IFNG": (["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74"],
                    ["INTERFERON_GAMMA_RESPONSE", "ALLOGRAFT_REJECTION"]),
    "NF-kB/TNFA": (["CXCL8", "CXCL1", "CXCL2", "IL1B", "ICAM1", "SOD2", "NAMPT"],
                   ["TNFA_SIGNALING_VIA_NFKB"]),
    "EMT": (["VIM", "FN1", "COL1A1", "SPARC", "POSTN", "THBS2"], ["EPITHELIAL_MESENCHYMAL_TRANSITION"]),
    "G2M": (["TOP2A", "CCNB1", "CCNB2", "UBE2C", "BUB1", "TPX2", "NUSAP1", "BIRC5"],
            ["G2M_CHECKPOINT", "E2F_TARGETS"]),
}

bh = pd.read_csv(f"{OUT}/hallmark_all_tests_bh.csv")
bh["enr"] = (bh.qval < 0.05) & (bh.overlap >= 3)
rows = []
for c in CANCERS:
    df = pd.read_csv(f"{CNMF}/{c}/cnmf_{c}/cnmf_{c}.gene_spectra_score.k_20.dt_0_1.txt", sep="\t", index_col=0)
    top = {f"cNMF_{i}": set(df.loc[i].sort_values(ascending=False).head(50).index) for i in df.index}
    for cat, (markers, sets) in CATS.items():
        missing = [m for m in markers if m not in df.columns]
        for p, g in top.items():
            hit = [m for m in markers if m in g]
            if len(hit) < 2:
                continue
            t = bh[(bh.cancer == c) & (bh.program == p)]
            for s in sets:
                r = t[t.hallmark == s].iloc[0]
                rows.append(dict(category=cat, cancer=c, program=p, n_markers=len(hit), markers=",".join(hit),
                                 hallmark=s, overlap=r.overlap, qval=r.qval, enriched=bool(r.enr),
                                 all_enriched=";".join(t[t.enr].sort_values("qval").hallmark),
                                 markers_not_in_panel=",".join(missing)))
res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/part4_marker_programs.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 120, "display.max_rows", 500)
print(res.drop(columns=["all_enriched", "markers_not_in_panel"]).to_string(index=False))
print()
print(res.drop_duplicates(["category", "cancer", "program"])[["category", "cancer", "program", "all_enriched"]].to_string(index=False))
print(res[["category", "markers_not_in_panel"]].drop_duplicates().to_string(index=False))
