#!/usr/bin/env python
"""
cMP8 (EMT) fibroblast-contamination check -- same design as the cMP6 audit
(cmp6-myeloid-contamination-83948d/cmp6_check/analyze.py), with fibroblast /
stromal fraction in place of myeloid / immune fraction.
"""
import os
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from scipy.stats import spearmanr, pearsonr

BASE = "/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
CMP6_OUT = ("/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/"
            "worktrees/cmp6-myeloid-contamination-83948d/cmp6_check/out")
ATLAS = f"{CMP6_OUT}/atlas_meta_slim.csv.gz"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

PROG = {"Brain": "5", "Pancreas": "14", "Prostate": "11"}   # cMP8 member programs
CANCERS = list(PROG)

FIB = ["COL1A1", "DCN", "PDGFRB", "LUM", "PDGFRA"]
FIB_ID = ["DCN", "PDGFRB", "LUM", "PDGFRA"]          # fibroblast-identity (non-EMT) markers
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
EMT = ["COL1A2", "SPARC", "FN1", "VIM", "THBS2", "COL5A1", "COL6A3"]
GLIO = ["SOX2", "OLIG2", "GFAP"]
GENES = FIB + EPI + EMT + GLIO

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 80)
pd.set_option("display.max_rows", 500)


def sample_key_from_cols(sample, orig, proj):
    s = pd.Series(sample, dtype="object").where(
        lambda x: x.notna() & ~x.isin(["NA", "nan", "", "NaN"]))
    oi = pd.Series(orig, dtype="object").astype(str)
    pr = pd.Series(proj, dtype="object").astype(str)
    fb = np.where((oi != pr) & (~oi.isin(["SeuratProject", "nan", "NA", ""])),
                  pr + "::" + oi, pr)
    return s.fillna(pd.Series(fb, index=s.index)).astype(str).values


def corr(x, y, label, c):
    if len(x) < 4:
        return dict(cancer_type=c, comparison=label, n=len(x), spearman_r=np.nan,
                    spearman_p=np.nan, pearson_r=np.nan, pearson_p=np.nan)
    sr, sp_ = spearmanr(x, y); pr, pp = pearsonr(x, y)
    return dict(cancer_type=c, comparison=label, n=len(x), spearman_r=sr, spearman_p=sp_,
                pearson_r=pr, pearson_p=pp)


print("loading atlas meta ...", flush=True)
atlas = pd.read_csv(ATLAS, low_memory=False)
atlas = atlas[atlas.cancer_type.isin(CANCERS)].copy()
print("Consensus_Cell_Type levels in Brain/Pancreas/Prostate:")
print(atlas.Consensus_Cell_Type.value_counts(dropna=False).to_string())
atlas["skey"] = sample_key_from_cols(atlas["sample"], atlas["orig.ident"], atlas["project_name"])
id2skey = dict(zip(atlas.cell_id, atlas.skey))

levels = set(atlas.Consensus_Cell_Type.dropna().astype(str))
FIB_LVL = [l for l in levels if "fibro" in l.lower()]
STROMA_LVL = [l for l in levels if any(k in l.lower() for k in
              ["fibro", "endothel", "smooth", "pericyte", "mesench", "stroma"])]
print("fibroblast levels:", FIB_LVL, " stromal levels:", STROMA_LVL, flush=True)

comp = []
for (ct, sk), g in atlas.groupby(["cancer_type", "skey"]):
    n = len(g); vc = g.Consensus_Cell_Type.astype(str).value_counts()
    nf = sum(vc.get(k, 0) for k in FIB_LVL)
    comp.append({"cancer_type": ct, "skey": sk, "n_cells_total": n,
                 "n_fibroblast": int(nf), "fibroblast_frac": nf / n,
                 "stromal_frac": sum(vc.get(k, 0) for k in STROMA_LVL) / n,
                 "epithelial_frac": vc.get("Epithelial", 0) / n})
comp = pd.DataFrame(comp)

per_sample_all, corr_rows, dec_rows, conc_rows, byct_rows, prog_top = [], [], [], [], [], []

for c in CANCERS:
    pk = PROG[c]
    print("=" * 78, f"\n{c}  (cMP8 member cNMF program = {pk})", flush=True)
    spec = pd.read_csv(f"{BASE}/{c}/cnmf_{c}/cnmf_{c}.gene_spectra_score.k_20.dt_0_1.txt",
                       sep="\t", index_col=0)
    spec.index = spec.index.astype(str)
    top = list(spec.loc[pk].sort_values(ascending=False).head(25).index)
    print("  program top-25:", ", ".join(top))
    prog_top.append({"cancer_type": c, "program": pk, "top25": ";".join(top)})

    A = ad.read_h5ad(f"{BASE}/{c}/cnmf_{c}.h5ad")
    obs = A.obs.copy()
    U = pd.read_csv(f"{BASE}/{c}/cnmf_{c}/cnmf_{c}.usages.k_20.dt_0_1.consensus.txt",
                    sep="\t", index_col=0).loc[obs.index]
    U.columns = U.columns.astype(str)
    share = U.div(U.sum(1), axis=0)[pk].values
    obs["share"] = share

    cn = obs["cell_name"].astype("object")
    srx = obs["sample"].astype("object") if "sample" in obs else pd.Series(np.nan, index=obs.index)
    sk = cn.map(id2skey)
    fb = sample_key_from_cols(srx, obs["orig.ident"], obs["project_name"])
    sk = sk.where(sk.notna(), pd.Series(fb, index=obs.index))
    obs["skey"] = sk.astype(str).values

    if c == "Brain":
        epi_mask = obs["cell_type"].astype(str).eq("Malignant").values
        epi_def = "author cell_type == 'Malignant'"
    else:
        epi_mask = obs["Consensus_Cell_Type"].astype(str).eq("Epithelial").values
        epi_def = "Consensus_Cell_Type == 'Epithelial'"
    print(f"  cnmf cells={A.n_obs}  epithelial compartment ({epi_def}) = {int(epi_mask.sum())}", flush=True)

    gidx = [g for g in GENES if g in A.var_names]
    X = A[:, gidx].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    cnt = pd.DataFrame(X, index=obs.index, columns=gidx)
    print(f"  X integer-valued on marker genes: {bool(np.all(np.mod(X, 1) == 0))}; "
          f"max={X.max():.3g}; missing genes: {[g for g in GENES if g not in gidx]}", flush=True)

    # usage by meta-type inside the object
    lab = obs["Consensus_Cell_Type"].astype(str).values
    if c == "Brain":
        lab = np.concatenate([lab, "author:" + obs["cell_type"].astype(str).values])
        sh2 = np.concatenate([share, share])
    else:
        sh2 = share
    t = pd.DataFrame({"ct": lab, "share": sh2}).groupby("ct")["share"] \
        .agg(n="size", mean_share="mean", median_share="median").reset_index()
    t.insert(0, "cancer_type", c); byct_rows.append(t)

    # per-sample usage vs fibroblast fraction
    e = obs.loc[epi_mask, ["skey", "share"]]
    ps = e.groupby("skey")["share"].agg(mean_prog_share="mean", median_prog_share="median",
                                        n_epi_cells="size").reset_index()
    ps.insert(0, "cancer_type", c)
    ps = ps.merge(comp[comp.cancer_type == c].drop(columns="cancer_type"), on="skey", how="left")

    # top vs bottom usage decile, epithelial compartment
    ei = np.where(epi_mask)[0]
    order = np.argsort(share[ei]); nd = max(1, len(ei) // 10)
    lo, hi = ei[order[:nd]], ei[order[-nd:]]
    row = {"cancer_type": c, "program": pk, "n_epi": len(ei), "n_decile": nd,
           "hi_share_min": float(share[hi].min()), "hi_share_med": float(np.median(share[hi])),
           "lo_share_max": float(share[lo].max())}
    fid = [g for g in FIB_ID if g in gidx]; eid = [g for g in EPI if g in gidx]
    for grp, idx in [("hi", hi), ("lo", lo)]:
        s = cnt.iloc[idx] > 0
        for g in gidx:
            row[f"{grp}_{g}"] = float(s[g].mean())
        col1 = s["COL1A1"]; fib_any = s[fid].any(axis=1); epi_any = s[eid].any(axis=1)
        row[f"{grp}_anyFibID(DCN|PDGFRB|LUM|PDGFRA)"] = float(fib_any.mean())
        row[f"{grp}_anyEpi(EPCAM|KRT8|KRT18|KRT19|CDH1)"] = float(epi_any.mean())
        row[f"{grp}_n_COL1A1+"] = int(col1.sum())
        for g in fid + ["EPCAM", "KRT8", "KRT18"]:
            if g in s:
                row[f"{grp}_{g}|COL1A1+"] = float(s.loc[col1, g].mean()) if col1.any() else np.nan
        row[f"{grp}_anyFibID|COL1A1+"] = float(fib_any[col1].mean()) if col1.any() else np.nan
        row[f"{grp}_anyEpi|COL1A1+"] = float(epi_any[col1].mean()) if col1.any() else np.nan
        row[f"{grp}_COL1A1&DCN&PDGFRB&LUM"] = float((col1 & s["DCN"] & s["PDGFRB"] & s["LUM"]).mean())
        row[f"{grp}_COL1A1+&noFibID"] = float((col1 & ~fib_any).mean())
        row[f"{grp}_COL1A1+&noFibID&Epi+"] = float((col1 & ~fib_any & epi_any).mean())
        row[f"{grp}_COL1A1+&FibID+&Epi+"] = float((col1 & fib_any & epi_any).mean())
        row[f"{grp}_COL1A1+&FibID+&noEpi"] = float((col1 & fib_any & ~epi_any).mean())
    dec_rows.append(row)

    # per-sample top-decile / fibroblast-marker+ counts
    hb = obs.iloc[hi][["skey"]].copy()
    sh = cnt.iloc[hi] > 0
    hb["fibID"] = sh[fid].any(axis=1).values
    hb["col1_fibID"] = (sh["COL1A1"] & sh[fid].any(axis=1)).values
    hb["col1_noFibID"] = (sh["COL1A1"] & ~sh[fid].any(axis=1)).values
    agg = hb.groupby("skey").agg(n_hi_decile=("fibID", "size"), n_hi_fibID=("fibID", "sum"),
                                 n_hi_COL1A1_fibID=("col1_fibID", "sum"),
                                 n_hi_COL1A1_noFibID=("col1_noFibID", "sum")).reset_index()
    ps = ps.merge(agg, on="skey", how="left")
    for k in ["n_hi_decile", "n_hi_fibID", "n_hi_COL1A1_fibID", "n_hi_COL1A1_noFibID"]:
        ps[k] = ps[k].fillna(0).astype(int)
    ps["frac_epi_in_hi_decile"] = ps.n_hi_decile / ps.n_epi_cells
    per_sample_all.append(ps)

    d = ps.dropna(subset=["fibroblast_frac"]); d = d[d.n_epi_cells >= 10]
    for y in ["fibroblast_frac", "stromal_frac", "epithelial_frac"]:
        corr_rows.append(corr(d.mean_prog_share.values, d[y].values, f"{c}:usage~{y}", c))

    vc = hb["skey"].value_counts(); vf = hb.loc[hb.fibID, "skey"].value_counts()
    nfib = int(hb.fibID.sum())
    conc_rows.append({
        "cancer_type": c, "n_samples_total": len(ps),
        "n_samples_epi>=10": int((ps.n_epi_cells >= 10).sum()),
        "n_hi_decile": len(hb), "n_samples_with_hi_decile": int(len(vc)),
        "top_sample_hi": vc.index[0], "top_sample_share_hi": float(vc.iloc[0] / len(hb)),
        "top3_share_hi": float(vc.iloc[:3].sum() / len(hb)),
        "n_hi_fibID+": nfib, "n_samples_hosting_fibID+": int(len(vf)),
        "top_sample_fibID+": vf.index[0] if len(vf) else None,
        "top_sample_share_fibID+": float(vf.iloc[0] / nfib) if nfib else np.nan,
        "top3_share_fibID+": float(vf.iloc[:3].sum() / nfib) if nfib else np.nan,
        "n_samples(epi>=10)_with>=5pct_epi_in_hi_decile":
            int(((ps.frac_epi_in_hi_decile >= 0.05) & (ps.n_epi_cells >= 10)).sum()),
    })
    del A, X

per_sample = pd.concat(per_sample_all, ignore_index=True)
cod = pd.DataFrame(corr_rows)
pool = per_sample.dropna(subset=["fibroblast_frac"]); pool = pool[pool.n_epi_cells >= 10].copy()
extra = []
for y in ["fibroblast_frac", "stromal_frac"]:
    extra.append(corr(pool.mean_prog_share.values, pool[y].values, f"POOLED:usage~{y}", "POOLED"))
    ur = pool.groupby("cancer_type")["mean_prog_share"].rank()
    yr = pool.groupby("cancer_type")[y].rank()
    sr, sp_ = spearmanr(ur, yr)
    extra.append(dict(cancer_type="POOLED", comparison=f"POOLED(within-cancer-ranked):usage~{y}",
                      n=len(pool), spearman_r=sr, spearman_p=sp_))
cod = pd.concat([cod, pd.DataFrame(extra)], ignore_index=True)
dec = pd.DataFrame(dec_rows); conc = pd.DataFrame(conc_rows)
byct = pd.concat(byct_rows, ignore_index=True)

per_sample.to_csv(f"{OUT}/cmp8_task1_per_sample.csv", index=False)
cod.to_csv(f"{OUT}/cmp8_task1_correlations.csv", index=False)
dec.set_index("cancer_type").T.to_csv(f"{OUT}/cmp8_task2_decile_coexpression.csv")
conc.to_csv(f"{OUT}/cmp8_task3_concentration.csv", index=False)
byct.to_csv(f"{OUT}/cmp8_usage_by_metatype.csv", index=False)
pd.DataFrame(prog_top).to_csv(f"{OUT}/cmp8_member_program_top25.csv", index=False)

print("\n##### per-sample (epithelial usage vs fibroblast fraction) #####")
print(per_sample.sort_values(["cancer_type", "mean_prog_share"], ascending=[True, False]).to_string(index=False))
print("\n##### correlations #####"); print(cod.to_string(index=False))
print("\n##### decile co-expression (fraction count>0) #####"); print(dec.set_index("cancer_type").T.to_string())
print("\n##### concentration #####"); print(conc.set_index("cancer_type").T.to_string())
print("\n##### usage by meta-type #####"); print(byct.to_string(index=False))
print("\nwrote", OUT)
