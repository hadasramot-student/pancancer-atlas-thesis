#!/usr/bin/env python
# 06_report_numbers.py -- numbers quoted in REPORT.txt (old vs new Figure 6.1 etc.)
import pandas as pd

OUT = "/groups/ofircohen-group/users/hadasra/thesis_figures/comment12_2026-10-21"
a = pd.read_csv(f"{OUT}/hallmark_all_tests_bh.csv")
b = pd.read_csv(f"{OUT}/hallmark_top_per_program_bh.csv")

g = a.groupby(["cancer", "program"])
print("programs with >=1 test raw p<0.05 (any overlap):", g.apply(lambda d: (d.pval < 0.05).any()).sum())
print("programs with >=1 test raw p<0.05 and >=3 genes:", g.apply(lambda d: ((d.pval < 0.05) & (d.overlap >= 3)).any()).sum())
print("programs with >=1 test BH q<0.05 (any overlap):", g.apply(lambda d: (d.qval < 0.05).any()).sum())
print("programs with >=1 test BH q<0.05 and >=3 genes:", g.apply(lambda d: d.enriched.any()).sum())
print("programs with >=2 enriched pathways:", g.apply(lambda d: d.enriched.sum() >= 2).sum())

# old figure definition: old top_hallmark (raw p<0.05), HC_old (Jaccard > 0.1), filled cells only
old = b[(b.top_hallmark_old != "Unannotated") & (b.HC_old.astype(str) == "True")]
oc = old.groupby(["top_hallmark_old", "cancer"]).size()
print(f"OLD: HC programs w/ named pathway {len(old)}; HC pathway-cancer cells {len(oc)}; "
      f"pathway rows {old.top_hallmark_old.nunique()}; cancers {old.cancer.nunique()}; "
      f"cells with >=2 HC programs {(oc >= 2).sum()}")
print("OLD HC_old by cancer:", b[b.HC_old.astype(str) == "True"].groupby("cancer").size().to_dict())

e = b[b.top_hallmark_bh != "Unannotated"]
for thr in ["HC_0.1", "HC_0.0538"]:
    h = e[e[thr]]
    hc = h.groupby(["top_hallmark_bh", "cancer"]).size()
    print(f"NEW stats, {thr}: HC programs w/ enriched pathway {len(h)}; filled cells {len(hc)}; "
          f"rows with filled dot {h.top_hallmark_bh.nunique()}; >=2 HC in cell {(hc >= 2).sum()}")
print("NEW all cells:", e.groupby(["top_hallmark_bh", "cancer"]).ngroups,
      "rows:", e.top_hallmark_bh.nunique(), "enriched programs per cancer:",
      e.groupby("cancer").size().to_dict())
print("HC_0.0538 programs:", int(b["HC_0.0538"].sum()), "of which unannotated:",
      int((b["HC_0.0538"] & (b.top_hallmark_bh == "Unannotated")).sum()))
cells = pd.read_csv(f"{OUT}/fig6_1_cells.csv")
print("x>=2 cells:")
print(cells[cells.n_hc >= 2][["top_hallmark_bh", "cancer", "n_hc", "programs"]].to_string(index=False))
