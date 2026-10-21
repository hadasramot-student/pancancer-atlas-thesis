#!/usr/bin/env python3
"""Tables B30 / B38 are at res0.30 but markers and fgsea exist only at res0.50. Do the thesis
marker / pathway columns equal the res0.50 values of the SAME cluster ID? And how well does each
res0.30 cluster overlap the res0.50 cluster with that ID? Read-only; uses saved files only.
Writes res030_marker_source_check.csv."""
import re
import pandas as pd
CHK = "/groups/ofircohen-group/users/hadasra/thesis_checks/appendixB_2026-10-21"
CL = "/groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages"
key = lambda p: re.sub(r"[^a-z0-9]", "", str(p).lower())
rows = []
for tb, obj in [("B30", "Endothelial_Kidney"), ("B38", "Epithelial_Kidney")]:
    t = pd.read_csv(f"{CHK}/thesis_v88_tables/{tb}.csv", dtype=str)
    t.columns = ["Cluster", "label", "n", "conf", "markers", "hallmark"]
    mk = pd.read_csv(f"{CL}/{obj}_top_markers.csv").sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
    top3 = {str(c): ", ".join(g.gene.head(3)) for c, g in mk.groupby("cluster")}
    g = pd.read_csv(f"{CL}/FindMarkers_results/enrichment/{obj}_hallmarks_gsea.csv")
    g = g[g.padj < 0.05].assign(a=lambda d: -d.NES.abs()).sort_values(["cluster", "a"])
    toph = {str(c): d.pathway.iloc[0].replace("HALLMARK_", "") for c, d in g.groupby("cluster")}
    s = pd.read_csv(f"{CL}/{obj}_scimilarity_summary.tsv", sep="\t", dtype={"cluster": str})
    n50 = s[s.resolution == "louvain_res0.50"].set_index("cluster")
    for _, r in t.iterrows():
        c = r.Cluster.strip()
        rows.append(dict(table=tb, object=obj, cluster=c, n_res030=r.n,
                         n_res050_same_id=int(n50.n_cells.get(c, 0)),
                         label_res050_same_id=n50.majority_label.get(c, ""),
                         conf_res050_same_id=n50.confidence_pct.get(c, ""),
                         thesis_markers=r.markers, res050_markers_same_id=top3.get(c, "none"),
                         markers_equal_res050_same_id=r.markers.strip() == top3.get(c, "none"),
                         thesis_hallmark=r.hallmark, res050_hallmark_same_id=toph.get(c, "—"),
                         hallmark_equal=key(r.hallmark) == key(toph.get(c, "—"))))
R = pd.DataFrame(rows)
R.to_csv(f"{CHK}/res030_marker_source_check.csv", index=False)
print(R.groupby("table")[["markers_equal_res050_same_id", "hallmark_equal"]].agg(["sum", "count"]))
pd.set_option("display.width", 300, "display.max_colwidth", 60)
print(R[R.table == "B38"][["cluster", "n_res030", "n_res050_same_id", "label_res050_same_id", "conf_res050_same_id",
                          "markers_equal_res050_same_id", "hallmark_equal"]].to_string(index=False))
