#!/usr/bin/env python
"""
cMP10 immune/MHC-II contamination check -- same method as the cMP6 audit
(cmp6-myeloid-contamination-83948d/cmp6_check/analyze.py + supp.py), applied
to cMP10 = Breast_cNMF_11 ; Pancreas_cNMF_11 ; Skin_cNMF_5
(meta_program_catalogue_J_gt_0.10.csv).

  Task 1  per-sample cNMF program usage (epithelial cells) vs per-sample
          myeloid / immune / lymphoid / B / T fraction (Chapter-2 consensus),
          per cancer, pooled, pooled ranked-within-cancer.  n_epi>=10, n>=4.
  Task 2  top vs bottom usage decile among consensus-Epithelial cells:
          detection (count>0, >=2) of MHC-II/myeloid genes vs epithelial markers.
  Task 4  concentration of top-decile cells across samples.
  Supp    mean usage by consensus meta-type inside each cNMF object.

Differences from the cMP6 script (forced by object size, not method):
  * X is read in row chunks from a backed h5ad (Breast cnmf h5ad is 63 GB).
  * composition is recomputed from the same atlas_meta_slim.csv.gz for the
    three cMP10 cancers.
"""
import os
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from scipy.stats import spearmanr, pearsonr

BASE  = "/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
CMP6  = "/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/cmp6-myeloid-contamination-83948d/cmp6_check"
ATLAS = f"{CMP6}/out/atlas_meta_slim.csv.gz"
OUT   = "/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/cmp10_check/out"
os.makedirs(OUT, exist_ok=True)

PROG = {"Breast": "11", "Pancreas": "11", "Skin": "5"}
CANCERS = list(PROG)

# cMP6 audit genes (kept for direct comparability)
MARKERS  = ["PTPRC", "CD68", "LYZ"]
# cMP10-specific MHC-II / myeloid-leaning genes
MHC2     = ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "SPP1"]
EPI_CTRL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
# supp.py extra lineage markers + fibroblast (cMP10 also carries BGN/VCAN/IGFBP7)
XMARK    = ["C1QB", "AIF1", "CD3D", "CD3E", "CD79A", "MS4A1", "CD19",
            "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"]
ALLG = MARKERS + MHC2 + EPI_CTRL + XMARK

IMMUNE   = ["Myeloid", "T-cells", "B-cell", "NK-cells"]
LYMPHOID = ["T-cells", "B-cell", "NK-cells"]
YS = ["myeloid_frac", "immune_frac", "lymphoid_frac", "bcell_frac", "tcell_frac",
      "fibro_frac"]

def sample_key_from_cols(sample, orig, proj):
    """coalesce: SRX/ERX accession -> project::orig.ident -> project (as cMP6)"""
    s = pd.Series(sample, dtype="object").where(
        lambda x: x.notna() & ~x.isin(["NA", "nan", "", "NaN"]))
    oi = pd.Series(orig, dtype="object").astype(str)
    pr = pd.Series(proj, dtype="object").astype(str)
    fb = np.where((oi != pr) & (~oi.isin(["SeuratProject", "nan", "NA", ""])),
                  pr + "::" + oi, pr)
    return s.fillna(pd.Series(fb, index=s.index)).astype(str).values

def corr_block(x, y, label):
    if len(x) < 4:
        return dict(comparison=label, n=len(x), spearman_r=np.nan, spearman_p=np.nan,
                    pearson_r=np.nan, pearson_p=np.nan)
    sr, sp_ = spearmanr(x, y); pr, pp = pearsonr(x, y)
    return dict(comparison=label, n=len(x), spearman_r=sr, spearman_p=sp_,
                pearson_r=pr, pearson_p=pp)

def read_gene_counts(path, genes, chunk=200_000):
    """dense (n_cells x len(genes)) counts, reading X row-chunked from backed h5ad"""
    A = ad.read_h5ad(path, backed="r")
    have = [g for g in genes if g in A.var_names]
    miss = [g for g in genes if g not in A.var_names]
    if miss: print("    genes missing from var:", miss, flush=True)
    cols = np.array([A.var_names.get_loc(g) for g in have])
    out = np.zeros((A.n_obs, len(have)), dtype=np.float32)
    for i in range(0, A.n_obs, chunk):
        blk = A.X[i:i + chunk]
        blk = blk[:, cols]
        out[i:i + chunk] = blk.toarray() if sp.issparse(blk) else np.asarray(blk)
        print(f"    X rows {i:,}-{min(i+chunk, A.n_obs):,} / {A.n_obs:,}", flush=True)
    obs = A.obs.copy()
    A.file.close()
    return obs, pd.DataFrame(out, index=obs.index, columns=have)

# ----------------------------------------------------------------------------
print("loading atlas meta ...", flush=True)
atlas = pd.read_csv(ATLAS, low_memory=False)
atlas = atlas[atlas.cancer_type.isin(CANCERS)].copy()
atlas["skey"] = sample_key_from_cols(atlas["sample"], atlas["orig.ident"], atlas["project_name"])
id2skey = dict(zip(atlas.cell_id, atlas.skey))

comp_rows = []
for (ct, sk), g in atlas.groupby(["cancer_type", "skey"]):
    n = len(g); vc = g.Consensus_Cell_Type.value_counts()
    comp_rows.append({"cancer_type": ct, "skey": sk, "n_cells_total": n,
        "n_myeloid": int(vc.get("Myeloid", 0)),
        "myeloid_frac": vc.get("Myeloid", 0) / n,
        "immune_frac": sum(vc.get(k, 0) for k in IMMUNE) / n,
        "lymphoid_frac": sum(vc.get(k, 0) for k in LYMPHOID) / n,
        "bcell_frac": vc.get("B-cell", 0) / n,
        "tcell_frac": vc.get("T-cells", 0) / n,
        "fibro_frac": vc.get("Fibroblasts", 0) / n,
        "epithelial_frac": vc.get("Epithelial", 0) / n})
comp = pd.DataFrame(comp_rows)
comp.to_csv(f"{OUT}/per_sample_composition.csv", index=False)
print(f"  {len(comp)} atlas samples across {CANCERS}", flush=True)
del atlas

per_sample_all, cod_all, t2_rows, conc_rows, byct_all = [], [], [], [], []

for c in CANCERS:
    pk = PROG[c]
    print("=" * 78, f"\n{c}  (cMP10 member cNMF program = {pk})", flush=True)
    obs, cntdf = read_gene_counts(f"{BASE}/{c}/cnmf_{c}.h5ad", ALLG)
    U = pd.read_csv(f"{BASE}/{c}/cnmf_{c}/cnmf_{c}.usages.k_20.dt_0_1.consensus.txt",
                    sep="\t", index_col=0)
    U.columns = U.columns.astype(str)
    U = U.loc[obs.index]
    Unorm = U.div(U.sum(1), axis=0)
    prog_share = Unorm[pk].values
    prog_raw = U[pk].values

    cn = obs["cell_name"].astype("object") if "cell_name" in obs else pd.Series(obs.index, index=obs.index)
    srx = obs["sample"].astype("object") if "sample" in obs else pd.Series(np.nan, index=obs.index)
    sk = cn.map(id2skey)
    need = sk.isna()
    sk_fb = sample_key_from_cols(srx, obs["orig.ident"], obs["project_name"])
    sk = sk.where(~need, pd.Series(sk_fb, index=obs.index))
    obs["skey"] = sk.astype(str).values
    obs["prog_share"] = prog_share; obs["prog_raw"] = prog_raw
    print(f"  sample keys: {int((~need).sum()):,} via atlas cell_id join, "
          f"{int(need.sum()):,} via obs fallback", flush=True)

    # ---- Supp: usage by consensus meta-type -----------------------------
    ctv = obs["Consensus_Cell_Type"].astype(str)
    tab = (pd.DataFrame({"ct": ctv.values, "share": prog_share}).groupby("ct")
           .agg(n=("share", "size"), mean_share=("share", "mean"),
                median_share=("share", "median")).reset_index())
    tab.insert(0, "cancer_type", c); tab.insert(1, "program", pk)
    byct_all.append(tab)
    if "cell_type" in obs:
        t2b = (pd.DataFrame({"ct": "author:" + obs["cell_type"].astype(str).values,
                             "share": prog_share}).groupby("ct")
               .agg(n=("share", "size"), mean_share=("share", "mean"),
                    median_share=("share", "median")).reset_index())
        t2b = t2b[t2b.n >= 100]
        t2b.insert(0, "cancer_type", c); t2b.insert(1, "program", pk)
        byct_all.append(t2b)

    epi_mask = ctv.eq("Epithelial").values
    epi_def = "Consensus_Cell_Type == 'Epithelial'"
    n_epi = int(epi_mask.sum())
    print(f"  cnmf cells={len(obs):,}  epithelial ({epi_def}) = {n_epi:,}", flush=True)

    # ================= TASK 1 ===========================================
    e = obs.loc[epi_mask, ["skey", "prog_share", "prog_raw"]]
    ps = (e.groupby("skey").agg(mean_prog_share=("prog_share", "mean"),
                                mean_prog_raw=("prog_raw", "mean"),
                                n_epi_cells=("prog_share", "size")).reset_index())
    ps["cancer_type"] = c
    ps = ps.merge(comp[comp.cancer_type == c].drop(columns="cancer_type"), on="skey", how="left")
    per_sample_all.append(ps)
    d = ps.dropna(subset=["myeloid_frac", "mean_prog_share"])
    d = d[d.n_epi_cells >= 10]
    for y in YS:
        r = corr_block(d["mean_prog_share"].values, d[y].values, f"{c}:usage~{y}")
        r["cancer_type"] = c; cod_all.append(r)
    print(f"  Task1 n_samples usable = {len(d)}  (epi samples total {len(ps)}, "
          f"unmatched to atlas {int(ps.myeloid_frac.isna().sum())})", flush=True)

    # ================= TASK 2 ===========================================
    ei = np.where(epi_mask)[0]
    order = np.argsort(prog_share[ei], kind="stable")
    ndec = max(1, len(ei) // 10)
    lo, hi = ei[order[:ndec]], ei[order[-ndec:]]
    row = {"cancer_type": c, "program": pk, "n_epi": n_epi, "n_decile": ndec,
           "epi_def": epi_def,
           "high_share_min": float(prog_share[hi].min()),
           "high_share_med": float(np.median(prog_share[hi])),
           "low_share_max": float(prog_share[lo].max()),
           "low_share_med": float(np.median(prog_share[lo]))}
    for grp, idx in {"high": hi, "low": lo}.items():
        sub = cntdf.iloc[idx]
        for g in sub.columns:
            row[f"{grp}_{g}_pos"] = float((sub[g] > 0).mean())
            row[f"{grp}_{g}_ge2"] = float((sub[g] >= 2).mean())
            row[f"{grp}_{g}_mean"] = float(sub[g].mean())
        for setname, gl in [("MYEL3", MARKERS), ("MHC2", MHC2), ("EPI", EPI_CTRL)]:
            gl = [g for g in gl if g in sub]
            m = sub[gl].values > 0
            row[f"{grp}_{setname}_ANY_pos"] = float(m.any(1).mean())
            row[f"{grp}_{setname}_ALL_pos"] = float(m.all(1).mean())
        mh = (sub[[g for g in MHC2 if g in sub]].values > 0).any(1)
        ep = (sub[[g for g in EPI_CTRL if g in sub]].values > 0).any(1)
        row[f"{grp}_MHC2any_and_EPIany"] = float((mh & ep).mean())
    t2_rows.append(row)

    # ================= TASK 4 ===========================================
    hb = obs.iloc[hi][["skey"]].copy()
    sub = cntdf.iloc[hi]
    hb["myel3"] = (sub[[g for g in MARKERS if g in sub]].values > 0).any(1)
    hb["mhc2"]  = (sub[[g for g in MHC2 if g in sub]].values > 0).any(1)
    vc_all = hb["skey"].value_counts()
    vc_epi = obs.loc[epi_mask, "skey"].value_counts()   # baseline spread of all epi cells
    r = {"cancer_type": c, "n_high_decile": len(hb),
         "n_samples_with_epi_cells": int(len(vc_epi)),
         "n_samples_with_high_decile_cells": int(len(vc_all)),
         "top_sample_share_of_high_decile": float(vc_all.iloc[0] / len(hb)),
         "top3_sample_share_of_high_decile": float(vc_all.iloc[:3].sum() / len(hb)),
         "top_sample_high_decile": vc_all.index[0],
         "top_sample_share_of_all_epi": float(vc_epi.get(vc_all.index[0], 0) / n_epi),
         "top3_samples_share_of_all_epi": float(vc_epi.reindex(vc_all.index[:3]).fillna(0).sum() / n_epi),
         "largest_epi_sample_share": float(vc_epi.iloc[0] / n_epi),
         "n_samples_for_50pct_high": int(np.searchsorted(vc_all.cumsum().values, 0.5 * len(hb)) + 1),
         "n_samples_for_50pct_epi": int(np.searchsorted(vc_epi.cumsum().values, 0.5 * n_epi) + 1)}
    for lab, col in [("marker_pos", "myel3"), ("mhc2_pos", "mhc2")]:
        dr = hb[hb[col]]; vc = dr["skey"].value_counts()
        r[f"n_high_decile_{lab}"] = int(len(dr))
        r[f"n_samples_hosting_{lab}"] = int(len(vc))
        r[f"top_sample_share_of_{lab}"] = float(vc.iloc[0] / len(dr)) if len(dr) else np.nan
        r[f"top3_sample_share_of_{lab}"] = float(vc.iloc[:3].sum() / len(dr)) if len(dr) else np.nan
        r[f"top_sample_{lab}"] = vc.index[0] if len(dr) else None
    conc_rows.append(r)
    psb = pd.DataFrame({"n_high_decile": vc_all, "n_epi": vc_epi.reindex(vc_all.index)})
    psb["high_over_expected"] = psb.n_high_decile / (psb.n_epi * 0.1)
    psb.insert(0, "cancer_type", c)
    psb.to_csv(f"{OUT}/task4_high_decile_by_sample_{c}.csv")
    del cntdf, obs, U, Unorm

per_sample = pd.concat(per_sample_all, ignore_index=True)
per_sample.to_csv(f"{OUT}/task1_per_sample.csv", index=False)
cod = pd.DataFrame(cod_all); cod.to_csv(f"{OUT}/task1_correlations_by_cancer.csv", index=False)
t2 = pd.DataFrame(t2_rows); t2.to_csv(f"{OUT}/task2_coexpression.csv", index=False)
conc = pd.DataFrame(conc_rows); conc.to_csv(f"{OUT}/task4_concentration.csv", index=False)
bc = pd.concat(byct_all, ignore_index=True); bc.to_csv(f"{OUT}/supp_usage_by_metatype.csv", index=False)

# ---- POOLED correlations (as cMP6) ----------------------------------------
pool = per_sample.dropna(subset=["myeloid_frac", "mean_prog_share"]).copy()
pool = pool[pool.n_epi_cells >= 10]
pooled_rows = []
for y in YS:
    pooled_rows.append(corr_block(pool["mean_prog_share"].values, pool[y].values,
                                  f"POOLED:usage~{y}"))
for y in ["myeloid_frac", "immune_frac"]:
    p2 = pool.copy()
    p2["u_rank"] = p2.groupby("cancer_type")["mean_prog_share"].rank()
    p2["m_rank"] = p2.groupby("cancer_type")[y].rank()
    sr, sp_ = spearmanr(p2["u_rank"], p2["m_rank"])
    pooled_rows.append(dict(comparison=f"POOLED(within-cancer-ranked):usage~{y}",
                            n=len(p2), spearman_r=sr, spearman_p=sp_,
                            pearson_r=np.nan, pearson_p=np.nan))
pooled = pd.DataFrame(pooled_rows)
pooled.to_csv(f"{OUT}/task1_correlations_pooled.csv", index=False)

# ---- print summary --------------------------------------------------------
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 80)
pd.set_option("display.max_rows", 300)
print("\n\n########## TASK 1: correlations by cancer ##########")
print(cod[["cancer_type", "comparison", "n", "spearman_r", "spearman_p", "pearson_r", "pearson_p"]].to_string(index=False))
print("\n########## TASK 1: pooled ##########")
print(pooled.to_string(index=False))
print("\n########## TASK 2: detection (count>0), high / low decile ##########")
for g in ALLG + ["MYEL3_ANY", "MHC2_ANY", "MHC2_ALL", "EPI_ANY", "EPI_ALL"]:
    hc, lc = f"high_{g}_pos", f"low_{g}_pos"
    if hc in t2:
        print(f"  {g:10s} " + "  ".join(f"{c}: {h:6.1%} / {l:6.1%}"
              for c, h, l in zip(t2.cancer_type, t2[hc], t2[lc])))
print("  MHC2any&EPIany " + "  ".join(f"{c}: {h:6.1%} / {l:6.1%}" for c, h, l in
      zip(t2.cancer_type, t2["high_MHC2any_and_EPIany"], t2["low_MHC2any_and_EPIany"])))
print("\n########## TASK 2: count>=2, high / low decile ##########")
for g in MARKERS + MHC2 + EPI_CTRL:
    hc, lc = f"high_{g}_ge2", f"low_{g}_ge2"
    if hc in t2:
        print(f"  {g:10s} " + "  ".join(f"{c}: {h:6.1%} / {l:6.1%}"
              for c, h, l in zip(t2.cancer_type, t2[hc], t2[lc])))
print("\n########## TASK 4: concentration across samples ##########")
print(conc.T.to_string())
print("\n########## SUPP: usage by meta-type ##########")
print(bc.to_string(index=False))
print("\nwrote:", OUT)
