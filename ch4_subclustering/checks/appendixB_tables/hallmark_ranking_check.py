#!/usr/bin/env python3
"""Which ranking of the saved fgsea results reproduces v88 Tables B10-B17 (top 3 Hallmark per cluster)?
Read-only. Usage: hallmark_ranking_check.py  -> prints per table, and writes hallmark_ranking_check.csv"""
import os, re
import pandas as pd
CHK = "/groups/ofircohen-group/users/hadasra/thesis_checks/appendixB_2026-10-21"
ENR = "/groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages/FindMarkers_results/enrichment"
T = {"B10": "T-cells_Liver_Biliary", "B11": "Myeloid_Prostate", "B12": "Epithelial_Kidney",
     "B13": "B-cell_Liver_Biliary", "B14": "Fibroblasts_Brain", "B15": "Endothelial_Kidney",
     "B17": "Hematopoietic_Stem_and_Progenitor_cells_Breast"}
key = lambda p: re.sub(r"[^a-z0-9]", "", str(p).lower().replace("hallmark_", ""))
RANK = {"absNES desc": lambda g: g.assign(k=-g.NES.abs()).sort_values("k"),
        "NES desc": lambda g: g.sort_values("NES", ascending=False),
        "padj asc": lambda g: g.sort_values(["padj", "pval"]),
        "pval asc": lambda g: g.sort_values("pval"),
        "padj asc, absNES desc": lambda g: g.assign(k=-g.NES.abs()).sort_values(["padj", "k"])}
rows = []
for tb, obj in T.items():
    t = pd.read_csv(f"{CHK}/thesis_v88_tables/{tb}.csv", dtype=str).iloc[:, :4]
    t.columns = ["Cluster", "Pathway", "NES", "padj"]
    g = pd.read_csv(f"{ENR}/{obj}_hallmarks_gsea.csv")
    g["cluster"] = g.cluster.astype(str); g["key"] = g.pathway.map(key)
    # do the thesis (pathway, NES) pairs exist in the saved fgsea file?
    m = t.assign(key=t.Pathway.map(key)).merge(g, left_on=["Cluster", "key"], right_on=["cluster", "key"], how="left")
    nes_ok = int((abs(m.NES_x.astype(float) - m.NES_y) < 0.0015).sum())
    sig = g[g.padj < 0.05]
    th = t.assign(key=t.Pathway.map(key)).groupby("Cluster").key.apply(list).to_dict()
    for name, f in RANK.items():
        re_ = {c: list(f(d).key.head(3)) for c, d in sig.groupby("cluster")}
        same_order = sum(th.get(c) == re_.get(c) for c in set(th) | set(re_))
        same_set = sum(set(th.get(c, [])) == set(re_.get(c, [])) for c in set(th) | set(re_))
        rows.append(dict(table=tb, object=obj, ranking=name, clusters=len(set(th) | set(re_)),
                         same_order=same_order, same_set=same_set,
                         thesis_rows=len(t), thesis_NES_found_in_saved_fgsea=nes_ok))
R = pd.DataFrame(rows)
R.to_csv(f"{CHK}/hallmark_ranking_check.csv", index=False)
print(R.to_string(index=False))

# Clusters still different under padj ranking: are the swapped pathways tied in padj?
out = []
for tb, obj in T.items():
    t = pd.read_csv(f"{CHK}/thesis_v88_tables/{tb}.csv", dtype=str).iloc[:, :4]
    t.columns = ["Cluster", "Pathway", "NES", "padj"]
    g = pd.read_csv(f"{ENR}/{obj}_hallmarks_gsea.csv")
    g["cluster"] = g.cluster.astype(str); g["key"] = g.pathway.map(key)
    sig = g[g.padj < 0.05].sort_values(["padj", "pval"])
    th = t.assign(key=t.Pathway.map(key)).groupby("Cluster").key.apply(list).to_dict()
    for c, d in sig.groupby("cluster"):
        re_ = list(d.key.head(3))
        if th.get(c) == re_:
            continue
        cut = d.padj.iloc[min(2, len(d) - 1)]                       # padj of the 3rd-ranked set
        tied = d[d.padj == cut].key.tolist()
        out.append(dict(table=tb, object=obj, cluster=c, thesis=", ".join(th.get(c, [])),
                        padj_ranked=", ".join(re_), padj_at_rank3=cut,
                        n_sets_tied_at_that_padj=len(tied),
                        explained_by_tie=set(th.get(c, [])) <= set(d[d.padj <= cut].key)))
O = pd.DataFrame(out)
O.to_csv(f"{CHK}/hallmark_padj_tie_check.csv", index=False)
print(O.to_string(index=False))
