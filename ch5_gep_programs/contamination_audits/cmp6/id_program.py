#!/usr/bin/env python
"""Identify which cNMF program (K=20) in each affected cancer type is the
MHC-II / interferon-response program that feeds into consensus cMP6."""
import pandas as pd, numpy as np

BASE = "/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
CANCERS = ["Brain", "Head_and_Neck", "Liver_Biliary", "Ovarian", "Prostate"]

# cMP6 consensus gene list (Table B97) + close MHC-II / IFN relatives
CMP6 = ["LAPTM5", "CD74", "HLA-DRA", "HLA-DPB1", "HLA-DQA1", "HLA-DRB1", "PTPRC", "C1QB"]
MHC2 = ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1",
        "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB"]
MYELOID = ["PTPRC", "C1QB", "C1QA", "C1QC", "LYZ", "CD68", "AIF1", "TYROBP",
           "FCER1G", "LAPTM5", "SRGN", "CD14", "ITGAM", "CSF1R"]

def load_score(c):
    f = f"{BASE}/{c}/cnmf_{c}/cnmf_{c}.gene_spectra_score.k_20.dt_0_1.txt"
    return pd.read_csv(f, sep="\t", index_col=0)  # 20 programs x genes

for c in CANCERS:
    s = load_score(c)                      # programs (rows) x genes (cols)
    # per-gene rank of each program's loading (1 = highest-loading program for that gene)
    rank = s.rank(axis=0, ascending=False)
    present = [g for g in CMP6 if g in s.columns]
    missing = [g for g in CMP6 if g not in s.columns]
    sub_rank = rank[present]
    sub_z = s[present]
    # score each program: how many cMP6 genes it is the top program for, mean rank, mean z
    top1 = (sub_rank == 1).sum(axis=1)
    mean_rank = sub_rank.mean(axis=1)
    mean_z = sub_z.mean(axis=1)
    best = mean_z.idxmax()
    print("=" * 78)
    print(f"{c}: cMP6 genes present={present}  missing={missing}")
    tbl = pd.DataFrame({"n_cMP6_genes_top1": top1, "mean_rank_of_cMP6": mean_rank.round(2),
                        "mean_z_of_cMP6": mean_z.round(3)}).sort_values("mean_z_of_cMP6", ascending=False)
    print(tbl.head(6).to_string())
    print(f"  --> BEST MATCH: program {best}")
    # show that program's top 25 genes + its rank/z on each signature gene
    top_genes = s.loc[best].sort_values(ascending=False).head(25)
    print(f"  program {best} top-25 genes: {', '.join(top_genes.index)}")
    detail = pd.DataFrame({"z": s.loc[best, present].round(2),
                           "rank_among_20_programs": sub_rank.loc[best].astype(int)})
    print(detail.to_string())
    # also MHC-II / myeloid extended footprint for that program
    for name, gl in [("MHC-II", MHC2), ("myeloid", MYELOID)]:
        gg = [g for g in gl if g in s.columns]
        r = rank.loc[best, gg].astype(int)
        print(f"  {name}: median rank of {best} on {len(gg)} genes = {int(r.median())}; "
              f"genes where prog {best} is top-1: {[g for g in gg if r[g] == 1]}")
