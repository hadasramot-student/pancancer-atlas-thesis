#!/usr/bin/env python3
"""regenerate_and_compare.py — rebuild the Appendix B cluster tables from saved outputs and
compare them with the thesis tables. Nothing is re-clustered, re-annotated or re-tested.

Inputs (read-only)
  object_checks/<obj>_clusters_from_object.tsv   check_object.R: per-cluster n / majority label /
                                                 confidence recomputed from <obj>_clustered.RDS
  object_checks/<obj>_marker_check.tsv           check_object.R: saved markers vs object pct.1/pct.2
  clustered_lineages/<obj>_top_markers.csv       saved FindMarkers (louvain_res0.50 only)
  clustered_lineages/<obj>_scimilarity_summary.tsv   saved votes at every resolution (source search)
  clustered_lineages/FindMarkers_results/enrichment/<obj>_hallmarks_gsea.csv   saved fgsea (res0.50)
  thesis_v88_tables/B<nn>.csv                    tables extracted from thesis/The_Thesis_v88.docx
  Point 16.1 tables for Epithelial Breast / Colorectal (not in v88)

Usage
  regenerate_and_compare.py regen <object>     one object -> regenerated/B<nn>_<object>_regenerated.csv
  regenerate_and_compare.py summary            xlsx + comparisons + Part 3 + REPORT.txt
"""
import glob, os, re, sys
import pandas as pd

BASE = "/groups/ofircohen-group/users/hadasra"
CL = f"{BASE}/temp_objects/clustered_lineages"
ENR = f"{CL}/FindMarkers_results/enrichment"
CHK = f"{BASE}/thesis_checks/appendixB_2026-10-21"
OBJ_CHK = f"{CHK}/object_checks"
REGEN = f"{CHK}/regenerated"
THESIS = f"{CHK}/thesis_v88_tables"
P161 = (f"{BASE}/create_sc_object_workflow/.claude/worktrees/"
        "epithelial-colorectal-breast-tables-b36ee3/appendixB_tables")
COLS = ["Cluster", "Cell Type (Scimilarity)", "n Cells", "Conf. %",
        "Top Markers (log₂FC)", "Top Hallmark Pathway"]
NA_MK = "0.50 markers not comparable"

objects = pd.read_csv(f"{CHK}/scripts/objects.tsv", sep="\t", dtype=str).fillna("")


def pretty_pathway(p):
    return p.replace("HALLMARK_", "").replace("_", " ").title().replace("Il2 ", "IL2 ") \
        .replace("Il6 ", "IL6 ").replace("Jak ", "JAK ").replace("Stat3", "STAT3") \
        .replace("Stat5", "STAT5").replace("Tnfa ", "TNFA ").replace("Nfkb", "NFkB") \
        .replace("Mtorc1", "mTORC1").replace("Pi3K Akt Mtor", "PI3K AKT mTOR") \
        .replace("Dna ", "DNA ").replace("Uv ", "UV ").replace("E2F", "E2F").replace("G2M", "G2M") \
        .replace("Kras ", "KRAS ").replace("Tgf ", "TGF ").replace("Myc ", "MYC ").replace("P53", "p53") \
        .replace("Wnt ", "WNT ").replace("Il2", "IL2").replace("Il6", "IL6")


def norm_pathway(p):
    """Case/spacing-insensitive key, so 'Myc Targets V1' == 'MYC Targets V1'."""
    return re.sub(r"[^a-z0-9]", "", str(p).lower().replace("hallmark_", ""))


def load_markers(obj):
    f = f"{CL}/{obj}_top_markers.csv"
    mk = pd.read_csv(f)
    mk["cluster"] = mk["cluster"].astype(str)
    return mk.sort_values(["cluster", "avg_log2FC"], ascending=[True, False])


def load_gsea(obj):
    f = f"{ENR}/{obj}_hallmarks_gsea.csv"
    if not os.path.exists(f):
        return None
    g = pd.read_csv(f)
    g["cluster"] = g["cluster"].astype(str)
    g = g[g["padj"] < 0.05].copy()
    g["absNES"] = g["NES"].abs()
    return g.sort_values(["cluster", "absNES"], ascending=[True, False])


def regen(obj):
    row = objects[objects.object == obj].iloc[0]
    sheet, res = row.sheet, row.resolution
    votes = pd.read_csv(f"{OBJ_CHK}/{obj}_clusters_from_object.tsv", sep="\t", dtype={"cluster": str})
    v = votes[votes.resolution == f"louvain_res{res}"].copy()
    v["cl_int"] = v.cluster.astype(int)
    v = v.sort_values("cl_int")
    comparable = res == "0.50"          # markers and fgsea exist only at louvain_res0.50
    mk = load_markers(obj)
    gs = load_gsea(obj)
    out = []
    for _, r in v.iterrows():
        lab = r.majority_label + (" *" if r.confidence_pct < 50 else "")
        if comparable:
            top = mk[mk.cluster == r.cluster].head(3)
            markers = ", ".join(top.gene) if len(top) else "no markers saved"
            if gs is None:
                hall = "not available"
            else:
                h = gs[gs.cluster == r.cluster]
                hall = pretty_pathway(h.pathway.iloc[0]) if len(h) else "—"
        else:
            markers, hall = NA_MK, NA_MK
        out.append([r.cl_int, lab, int(r.n_cells), r.confidence_pct, markers, hall])
    df = pd.DataFrame(out, columns=COLS)
    os.makedirs(REGEN, exist_ok=True)
    df.to_csv(f"{REGEN}/{sheet}_{obj}_regenerated.csv", index=False)
    print(f"{sheet} {obj} res{res}: {len(df)} clusters, {df['n Cells'].sum():,} cells")


# ----------------------------------------------------------------- comparison ----
def read_thesis(sheet, obj, v88):
    if v88:
        f = f"{THESIS}/{v88}.csv"
        src = f"v88 {v88}"
    else:
        part = obj.split("_", 1)[1]
        f = f"{P161}/AppendixB_Epithelial_{part}.csv"
        src = "Point 16.1 table (not in v88)"
    if not os.path.exists(f):
        return None, src
    t = pd.read_csv(f, dtype=str)
    t.columns = COLS
    t["Cluster"] = t["Cluster"].astype(int)
    t["n_int"] = t["n Cells"].str.replace(",", "").astype(int)
    t["conf_num"] = t["Conf. %"].str.replace("%", "").astype(float)
    t["label_plain"] = t["Cell Type (Scimilarity)"].str.replace(r"\s*\*$", "", regex=True).str.strip()
    return t, src


SUMM = {}


def sci_summary(obj):
    if obj not in SUMM:
        s = pd.read_csv(f"{CL}/{obj}_scimilarity_summary.tsv", sep="\t", dtype={"cluster": str})
        s["cl_int"] = s.cluster.astype(int)
        SUMM[obj] = s
    return SUMM[obj]


def guess_source(t):
    """Best (object, resolution) whose (cluster, n_cells) pairs match the thesis rows."""
    pairs = set(zip(t.Cluster, t.n_int))
    best = []
    for obj in objects.object:
        s = sci_summary(obj)
        for res, g in s.groupby("resolution"):
            hit = len(pairs & set(zip(g.cl_int, g.n_cells)))
            if hit:
                best.append((hit, obj, res.replace("louvain_res", ""), len(g)))
    best.sort(reverse=True)
    return best[:2]


MKIDX = None


def marker_index():
    """top-3 marker triple -> list of (object, cluster), over all saved res0.50 marker files."""
    global MKIDX
    if MKIDX is None:
        MKIDX = {}
        for obj in objects.object:
            mk = load_markers(obj)
            for cl, g in mk.groupby("cluster"):
                MKIDX.setdefault(", ".join(g.gene.head(3)), []).append((obj, cl))
    return MKIDX


def compare(sheet, obj, res, v88):
    reg = pd.read_csv(f"{REGEN}/{sheet}_{obj}_regenerated.csv", dtype={"Top Markers (log₂FC)": str})
    t, src = read_thesis(sheet, obj, v88)
    rec = dict(sheet=sheet, object=obj, resolution=f"res{res}", thesis_source=src,
               n_clusters_thesis=None if t is None else len(t), n_clusters_regen=len(reg))
    if t is None:
        rec.update(status="NO THESIS TABLE", notes="")
        return rec, pd.DataFrame()
    m = t.merge(reg, on="Cluster", how="outer", suffixes=("_th", "_re"), indicator=True)
    diffs, minor = [], 0
    for _, r in m.iterrows():
        if r._merge != "both":
            diffs.append(dict(cluster=r.Cluster, field="row",
                              thesis="present" if r._merge == "left_only" else "missing",
                              regenerated="present" if r._merge == "right_only" else "missing"))
            continue
        lab_re = re.sub(r"\s*\*$", "", r["Cell Type (Scimilarity)_re"]).strip()
        if r.label_plain != lab_re:
            diffs.append(dict(cluster=r.Cluster, field="label", thesis=r.label_plain, regenerated=lab_re))
        if r.n_int != r["n Cells_re"]:
            diffs.append(dict(cluster=r.Cluster, field="n_cells", thesis=r.n_int, regenerated=r["n Cells_re"]))
        dc = abs(r.conf_num - r["Conf. %_re"])
        if dc > 1.0:
            diffs.append(dict(cluster=r.Cluster, field="conf", thesis=r.conf_num, regenerated=r["Conf. %_re"]))
        elif dc > 0:
            minor += 1
        if r["Top Markers (log₂FC)_re"] != NA_MK and \
                str(r["Top Markers (log₂FC)_th"]).strip() != str(r["Top Markers (log₂FC)_re"]).strip():
            diffs.append(dict(cluster=r.Cluster, field="markers", thesis=r["Top Markers (log₂FC)_th"],
                              regenerated=r["Top Markers (log₂FC)_re"]))
        hr = r["Top Hallmark Pathway_re"]
        if hr not in (NA_MK, "not available") and norm_pathway(r["Top Hallmark Pathway_th"]) != norm_pathway(hr):
            diffs.append(dict(cluster=r.Cluster, field="hallmark", thesis=r["Top Hallmark Pathway_th"],
                              regenerated=hr))
    d = pd.DataFrame(diffs)
    if len(d):
        d.insert(0, "object", obj)
        d.insert(0, "sheet", sheet)
    if not len(d):
        status = "MINOR" if minor else "IDENTICAL"
        notes = f"{minor} clusters differ in Conf. % by rounding only" if minor else ""
        if res != "0.50":
            notes = (notes + "; " if notes else "") + "markers/hallmark not compared (saved only at res0.50)"
    else:
        status = "DIFFERENT"
        notes = diagnose(t, reg, d, obj, res)
    rec.update(status=status, notes=notes)
    return rec, d


def diagnose(t, reg, d, obj, res):
    bits = []
    by = d.groupby("field").cluster.nunique().to_dict()
    bits.append("differing clusters by field: " + ", ".join(f"{k}={v}" for k, v in sorted(by.items())))
    if "row" in by:
        bits.append(f"cluster sets differ (thesis {len(t)} rows vs {len(reg)} clusters)")
    if "n_cells" in by or "label" in by:
        src = guess_source(t)
        if src:
            hit, so, sr, nc = src[0]
            frac = hit / len(t)
            if so == obj and sr == res:
                bits.append(f"cells best match this object at res{sr} ({hit}/{len(t)} rows)")
            elif so == obj:
                bits.append(f"n cells match THIS object at res{sr} ({hit}/{len(t)} rows) -> different resolution")
            elif frac >= 0.5:
                bits.append(f"n cells match {so} res{sr} ({hit}/{len(t)} rows) -> rows copied from another object")
            else:
                bits.append(f"best n-cell match only {hit}/{len(t)} rows ({so} res{sr}) -> scrambled/mixed rows")
        else:
            bits.append("no object/resolution matches the thesis n cells -> scrambled/unknown source")
    if "markers" in by:
        idx = marker_index()
        mk_rows = d[d.field == "markers"]
        srcs = {}
        for tri in mk_rows.thesis:
            for so, sc in idx.get(str(tri).strip(), []):
                srcs[so] = srcs.get(so, 0) + 1
        if srcs:
            top = sorted(srcs.items(), key=lambda x: -x[1])[:2]
            bits.append("thesis marker triples found in: " + ", ".join(f"{k} ({v})" for k, v in top))
        else:
            bits.append("thesis marker triples not found in any saved marker file")
    return "; ".join(bits)


# ----------------------------------------------------------------- Part 3 ----
PART3 = [  # (marker table, hallmark table, object, figure-table resolution)
    ("B02", "B10", "T-cells_Liver_Biliary"), ("B03", "B11", "Myeloid_Prostate"),
    ("B04", "B12", "Epithelial_Kidney"), ("B05", "B13", "B-cell_Liver_Biliary"),
    ("B06", "B14", "Fibroblasts_Brain"), ("B07", "B15", "Endothelial_Kidney"),
    ("B08", "B16", "NK-cells_Skin"), ("B09", "B17", "Hematopoietic_Stem_and_Progenitor_cells_Breast"),
]


def part3():
    recs, all_d = [], []
    for mt, ht, obj in PART3:
        mk = load_markers(obj)
        reg_mk = pd.concat([g.head(5) for _, g in mk.groupby("cluster")])
        reg_mk["cl_int"] = reg_mk.cluster.astype(int)
        reg_mk = reg_mk.sort_values(["cl_int", "avg_log2FC"], ascending=[True, False])
        out_mk = reg_mk[["cl_int", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]].rename(
            columns={"cl_int": "Cluster", "gene": "Gene", "avg_log2FC": "avg log₂FC", "p_val_adj": "Adj. p-value"})
        out_mk.to_csv(f"{REGEN}/{mt}_{obj}_top5_markers_regenerated.csv", index=False)
        # markers
        f = f"{THESIS}/{mt}.csv"
        t = pd.read_csv(f, dtype=str)
        t.columns = ["Cluster", "Gene", "lfc", "pct1", "pct2", "padj"]
        t["Cluster"] = t.Cluster.astype(int)
        th = t.groupby("Cluster").Gene.apply(list).to_dict()
        re_ = out_mk.groupby("Cluster").Gene.apply(list).to_dict()
        dif = [c for c in sorted(set(th) | set(re_)) if th.get(c) != re_.get(c)]
        lfc_dif = 0
        mm = t.merge(out_mk, on=["Cluster", "Gene"])
        lfc_dif = int((abs(mm.lfc.astype(float) - mm["avg log₂FC"]) > 0.0015).sum())
        status = "IDENTICAL" if not dif and not lfc_dif else ("MINOR" if not dif else "DIFFERENT")
        recs.append(dict(sheet=mt, object=obj, resolution="res0.50 (markers)", status=status,
                         n_clusters_thesis=len(th), n_clusters_regen=len(re_),
                         notes=(f"clusters with different top-5 genes: {dif}" if dif else "") +
                               (f"; {lfc_dif} log2FC values differ beyond rounding" if lfc_dif else "")))
        for c in dif:
            all_d.append(dict(sheet=mt, object=obj, cluster=c, field="top5 markers",
                              thesis=", ".join(th.get(c, [])), regenerated=", ".join(re_.get(c, []))))
        # hallmarks
        gs = load_gsea(obj)
        rows = []
        if gs is not None:
            for cl, g in gs.groupby("cluster"):
                for _, r in g.head(3).iterrows():
                    rows.append([int(cl), pretty_pathway(r.pathway), round(r.NES, 3), r.padj])
        reg_h = pd.DataFrame(rows, columns=["Cluster", "Hallmark Pathway", "NES", "Adj. p-value"]).sort_values(
            ["Cluster"], kind="stable")
        reg_h.to_csv(f"{REGEN}/{ht}_{obj}_top3_hallmarks_regenerated.csv", index=False)
        f = f"{THESIS}/{ht}.csv"
        if os.path.exists(f):
            t = pd.read_csv(f, dtype=str)
            t = t.iloc[:, :3]
            t.columns = ["Cluster", "Pathway", "NES"]
            t["Cluster"] = t.Cluster.astype(int)
            thh = t.groupby("Cluster").Pathway.apply(lambda x: [norm_pathway(p) for p in x]).to_dict()
        else:
            thh = {}          # v88 B16 has no table: "No pathways reached significance (padj < 0.05)"
        reh = reg_h.groupby("Cluster")["Hallmark Pathway"].apply(lambda x: [norm_pathway(p) for p in x]).to_dict()
        dif = [c for c in sorted(set(thh) | set(reh)) if thh.get(c) != reh.get(c)]
        status = "IDENTICAL" if not dif else "DIFFERENT"
        note = "thesis has no table (caption: no pathway padj < 0.05); " if not os.path.exists(f) else ""
        recs.append(dict(sheet=ht, object=obj, resolution="res0.50 (fgsea)", status=status,
                         n_clusters_thesis=len(thh), n_clusters_regen=len(reh),
                         notes=note + (f"clusters with different top-3 pathways: {dif}" if dif else
                                       "same pathways per cluster")))
        for c in dif:
            all_d.append(dict(sheet=ht, object=obj, cluster=c, field="top3 hallmarks",
                              thesis=", ".join(thh.get(c, [])), regenerated=", ".join(reh.get(c, []))))
    return pd.DataFrame(recs), pd.DataFrame(all_d)


def summary():
    recs, diffs = [], []
    for _, o in objects.iterrows():
        r, d = compare(o.sheet, o.object, o.resolution, o.v88_table)
        recs.append(r)
        diffs.append(d)
    # Representative-object duplicates in v88 (B91-B96) — same format, extra check
    for tb, obj, res in [("B91", "Myeloid_Prostate", "0.50"), ("B92", "B-cell_Liver_Biliary", "0.50"),
                         ("B93", "Fibroblasts_Brain", "0.50"), ("B94", "Endothelial_Kidney", "0.30"),
                         ("B95", "NK-cells_Skin", "0.50"),
                         ("B96", "Hematopoietic_Stem_and_Progenitor_cells_Breast", "0.50")]:
        sheet = objects[objects.object == obj].sheet.iloc[0]
        r, d = compare(sheet, obj, res, tb)
        r["sheet"] = f"{tb} (v88 duplicate of {sheet})"
        recs.append(r)
        if len(d):
            d["sheet"] = tb
        diffs.append(d)
    S = pd.DataFrame(recs)
    D = pd.concat([d for d in diffs if len(d)], ignore_index=True)
    S3, D3 = part3()
    S.to_csv(f"{CHK}/summary_B18_B92.csv", index=False)
    D.to_csv(f"{CHK}/differences_B18_B92.csv", index=False)
    S3.to_csv(f"{CHK}/summary_B02_B17.csv", index=False)
    D3.to_csv(f"{CHK}/differences_B02_B17.csv", index=False)
    with pd.ExcelWriter(f"{CHK}/Supplementary_Tables_AppendixB_regenerated.xlsx") as xw:
        for _, o in objects.iterrows():
            pd.read_csv(f"{REGEN}/{o.sheet}_{o.object}_regenerated.csv").to_excel(
                xw, sheet_name=o.sheet, index=False)
    write_report(S, D, S3, D3)
    print(S.status.value_counts().to_string())
    print(S3.status.value_counts().to_string())


def write_report(S, D, S3, D3):
    W = []
    W.append("APPENDIX B TABLES CHECK — 2026-10-21 (read-only; no clustering / Scimilarity / FindMarkers re-run)\n")
    W.append("Thesis tables: extracted from thesis/The_Thesis_v88.docx (no Supplementary_Tables_AppendixB.xlsx")
    W.append("and no script writing it exists on the cluster). v88 numbers the per-object tables B18-B90 (73")
    W.append("objects, no Epithelial Breast/Colorectal); the new numbering B18-B92 is used below and the v88")
    W.append("table actually compared is given in 'thesis_source'. Epithelial Breast/Colorectal (B39/B40) are")
    W.append("compared with the Point 16.1 tables (worktree epithelial-colorectal-breast-tables-b36ee3).\n")
    W.append("Regenerated values: n cells / majority label / confidence recomputed from each <object>_clustered.RDS")
    W.append("(check_object.R, scimilarity_celltype majority vote); top 3 markers by avg_log2FC from")
    W.append("<object>_top_markers.csv; top Hallmark = highest |NES| with padj < 0.05 in")
    W.append("FindMarkers_results/enrichment/<object>_hallmarks_gsea.csv. Markers and fgsea exist only at")
    W.append("louvain_res0.50, so for the two res0.30 tables (B30, B38) both columns read")
    W.append(f"'{NA_MK}' and are not compared.")
    W.append("Status: IDENTICAL; MINOR = only Conf. % differs by <= 1 point (rounding); DIFFERENT = any cluster")
    W.append("with a different label, n cells, confidence (> 1 point), markers, hallmark, or a missing/extra row.\n")
    W.append("COUNTS (B18-B92): " + ", ".join(f"{k} {v}" for k, v in
                                             S[~S.sheet.str.contains('duplicate')].status.value_counts().items()))
    W.append("COUNTS (B02-B17): " + ", ".join(f"{k} {v}" for k, v in S3.status.value_counts().items()) + "\n")
    pd.set_option("display.width", 400, "display.max_colwidth", 400, "display.max_rows", 500)
    W.append("SUMMARY B18-B92 (+ v88 B91-B96 duplicates)")
    W.append(S[["sheet", "object", "resolution", "status", "n_clusters_thesis", "n_clusters_regen",
                "thesis_source", "notes"]].to_string(index=False))
    W.append("\nSUMMARY B02-B17 (8 representative objects)")
    W.append(S3[["sheet", "object", "resolution", "status", "n_clusters_thesis", "n_clusters_regen",
                 "notes"]].to_string(index=False))
    W.append("\nFull per-cluster differences: differences_B18_B92.csv, differences_B02_B17.csv")
    W.append("Marker-file sanity: every <object>_top_markers.csv was checked against its object")
    W.append("(object_checks/<object>_marker_check.tsv; recomputed pct.1/pct.2 of the top 3 genes).")
    W.append(marker_sanity())
    W.append("\nFILES: regenerated/B<nn>_<object>_regenerated.csv (75), regenerated/B02-B17 top-5/top-3 CSVs,")
    W.append("Supplementary_Tables_AppendixB_regenerated.xlsx, summary_*.csv, differences_*.csv,")
    W.append("thesis_v88_tables/ (extracted thesis tables), object_checks/, scripts/.")
    open(f"{CHK}/REPORT.txt", "w").write("\n".join(W) + "\n")


def marker_sanity():
    bad = []
    for obj in objects.object:
        m = pd.read_csv(f"{OBJ_CHK}/{obj}_marker_check.tsv", sep="\t")
        dd = max((m["pct.1"] - m.pct1_obj).abs().max(), (m["pct.2"] - m.pct2_obj).abs().max())
        if not dd <= 0.002:
            bad.append(f"{obj} (max |diff| {dd:.3f})")
    return "  all 75 marker files match their objects (|diff| <= 0.002)" if not bad else \
        "  marker files NOT matching their object: " + "; ".join(bad)


if __name__ == "__main__":
    if sys.argv[1] == "regen":
        regen(sys.argv[2])
    else:
        summary()
