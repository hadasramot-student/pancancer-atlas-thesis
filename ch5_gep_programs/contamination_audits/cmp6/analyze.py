#!/usr/bin/env python
"""
cMP6 myeloid-contamination check.

Two cheap tests, using existing cNMF usages + Chapter-2 consensus meta-types:
  (1) per-sample cNMF MHC-II program usage  vs  per-sample myeloid fraction
  (2) macrophage/leukocyte marker co-expression in top- vs bottom-decile
      high-usage epithelial-compartment cells
"""
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from scipy.stats import spearmanr, pearsonr

BASE = "/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
OUT  = "/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/cmp6-myeloid-contamination-83948d/cmp6_check/out"
ATLAS = f"{OUT}/atlas_meta_slim.csv.gz"

# MHC-II/IFN cNMF program index per cancer type (K=20), identified from
# gene_spectra_score by top-loading of the cMP6 gene set (see id_program.py).
PROG = {"Brain": "10", "Head_and_Neck": "8", "Liver_Biliary": "7",
        "Ovarian": "4", "Prostate": "8"}
CANCERS = list(PROG)

MARKERS = ["PTPRC", "CD68", "LYZ"]
EPI_CTRL = ["EPCAM", "KRT8", "KRT18"]
IMMUNE = ["Myeloid", "T-cells", "B-cell", "NK-cells"]
LYMPHOID = ["T-cells", "B-cell", "NK-cells"]

# ----------------------------------------------------------------------------
def sample_key_from_cols(sample, orig, proj):
    """coalesce: SRX/ERX accession  ->  project::orig.ident  ->  project"""
    s = pd.Series(sample, dtype="object").where(
        lambda x: x.notna() & ~x.isin(["NA", "nan", "", "NaN"]))
    oi = pd.Series(orig, dtype="object").astype(str)
    pr = pd.Series(proj, dtype="object").astype(str)
    fb = np.where((oi != pr) & (~oi.isin(["SeuratProject", "nan", "NA", ""])),
                  pr + "::" + oi, pr)
    return s.fillna(pd.Series(fb, index=s.index)).astype(str).values

# ----------------------------------------------------------------------------
print("loading atlas meta ...")
atlas = pd.read_csv(ATLAS, low_memory=False)
atlas = atlas[atlas.cancer_type.isin(CANCERS)].copy()
atlas["skey"] = sample_key_from_cols(atlas["sample"], atlas["orig.ident"], atlas["project_name"])
# cell_id -> skey lookup (named studies join here)
id2skey = dict(zip(atlas.cell_id, atlas.skey))

# per-sample composition from the atlas (Chapter-2 consensus meta-types)
comp_rows = []
for (ct, sk), g in atlas.groupby(["cancer_type", "skey"]):
    n = len(g)
    vc = g.Consensus_Cell_Type.value_counts()
    row = {"cancer_type": ct, "skey": sk, "n_cells_total": n,
           "n_myeloid": int(vc.get("Myeloid", 0)),
           "myeloid_frac": vc.get("Myeloid", 0) / n,
           "immune_frac": sum(vc.get(k, 0) for k in IMMUNE) / n,
           "lymphoid_frac": sum(vc.get(k, 0) for k in LYMPHOID) / n,
           "bcell_frac": vc.get("B-cell", 0) / n,
           "tcell_frac": vc.get("T-cells", 0) / n,
           "epithelial_frac": vc.get("Epithelial", 0) / n}
    comp_rows.append(row)
comp = pd.DataFrame(comp_rows)
comp.to_csv(f"{OUT}/per_sample_composition.csv", index=False)
print(f"  {len(comp)} atlas samples across {CANCERS}")

# ----------------------------------------------------------------------------
per_sample_all = []
codetail_all = []
task2_rows = []
conc_rows = []

for c in CANCERS:
    pk = PROG[c]
    print("=" * 78, f"\n{c}  (MHC-II cNMF program = {pk})")
    A = ad.read_h5ad(f"{BASE}/{c}/cnmf_{c}.h5ad")
    obs = A.obs.copy()
    U = pd.read_csv(f"{BASE}/{c}/cnmf_{c}/cnmf_{c}.usages.k_20.dt_0_1.consensus.txt",
                    sep="\t", index_col=0)
    U = U.loc[obs.index]                      # align to obs order
    Unorm = U.div(U.sum(1), axis=0)           # per-cell program shares (sum=1)
    prog_share = Unorm[pk].values
    prog_raw   = U[pk].values

    # ---- sample key per cnmf cell -----------------------------------------
    cn = obs["cell_name"].astype("object")
    srx = obs["sample"].astype("object") if "sample" in obs else pd.Series(np.nan, index=obs.index)
    sk = cn.map(id2skey)                                   # named-study join by cell_id
    need = sk.isna()
    sk_fallback = sample_key_from_cols(srx, obs["orig.ident"], obs["project_name"])
    sk = sk.where(~need, pd.Series(sk_fallback, index=obs.index))
    obs["skey"] = sk.astype(str).values
    obs["prog_share"] = prog_share
    obs["prog_raw"] = prog_raw

    # ---- epithelial-compartment mask ------------------------------------
    if c == "Brain":
        epi_mask = obs["cell_type"].astype(str).eq("Malignant").values   # consensus has no epi for glioma
        epi_def = "author cell_type == 'Malignant' (no consensus-Epithelial cells exist for Brain)"
    else:
        epi_mask = obs["Consensus_Cell_Type"].astype(str).eq("Epithelial").values
        epi_def = "Consensus_Cell_Type == 'Epithelial'"
    n_epi = int(epi_mask.sum())
    print(f"  cnmf cells={A.n_obs}  epithelial-compartment ({epi_def}) = {n_epi}")

    # ---- marker counts from raw X ---------------------------------------
    genes = {g: A.var_names.get_loc(g) for g in MARKERS + EPI_CTRL if g in A.var_names}
    X = A.X.tocsc() if sp.issparse(A.X) else A.X
    cnt = {}
    for g, j in genes.items():
        col = X[:, j]
        cnt[g] = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
    cntdf = pd.DataFrame(cnt, index=obs.index)

    # ================= TASK 1 : per-sample usage vs myeloid frac =========
    e = obs.loc[epi_mask, ["skey", "prog_share", "prog_raw"]].copy()
    g = e.groupby("skey")
    ps = g.agg(mean_prog_share=("prog_share", "mean"),
               mean_prog_raw=("prog_raw", "mean"),
               n_epi_cells=("prog_share", "size")).reset_index()
    ps["cancer_type"] = c
    ps = ps.merge(comp[comp.cancer_type == c].drop(columns="cancer_type"),
                  on="skey", how="left")
    ps["frac_epi_high_decile"] = np.nan
    per_sample_all.append(ps)

    # correlations for this cancer (require >=2 samples w/ myeloid_frac)
    d = ps.dropna(subset=["myeloid_frac", "mean_prog_share"])
    d = d[d.n_epi_cells >= 10]                    # ignore samples w/ <10 epi cells
    def corr_block(x, y, label):
        if len(x) < 4:
            return dict(comparison=label, n=len(x), spearman_r=np.nan, spearman_p=np.nan,
                        pearson_r=np.nan, pearson_p=np.nan)
        sr, sp_ = spearmanr(x, y); pr, pp = pearsonr(x, y)
        return dict(comparison=label, n=len(x), spearman_r=sr, spearman_p=sp_,
                    pearson_r=pr, pearson_p=pp)
    for yname in ["myeloid_frac", "immune_frac", "lymphoid_frac", "bcell_frac", "tcell_frac"]:
        r = corr_block(d["mean_prog_share"].values, d[yname].values, f"{c}:usage~{yname}")
        r["cancer_type"] = c
        codetail_all.append(r)
    print(f"  Task1 n_samples usable = {len(d)}")

    # ================= TASK 2 : marker co-expression high vs low ========
    ei = np.where(epi_mask)[0]
    sh = prog_share[ei]
    order = np.argsort(sh)
    ndec = max(1, len(ei) // 10)
    low_idx  = ei[order[:ndec]]
    high_idx = ei[order[-ndec:]]
    hi_lo = {"high": high_idx, "low": low_idx}
    row = {"cancer_type": c, "program": pk, "n_epi": n_epi, "n_decile": ndec,
           "epi_def": epi_def,
           "high_share_min": float(prog_share[high_idx].min()),
           "high_share_med": float(np.median(prog_share[high_idx])),
           "low_share_max": float(prog_share[low_idx].max())}
    for grp, idx in hi_lo.items():
        bc = obs.index[idx]
        sub = cntdf.loc[bc]
        for g in MARKERS:
            if g in sub:
                row[f"{grp}_{g}_pos"]  = float((sub[g] > 0).mean())
                row[f"{grp}_{g}_ge2"]  = float((sub[g] >= 2).mean())
        anymk = np.zeros(len(sub), bool)
        for g in MARKERS:
            if g in sub: anymk |= (sub[g].values > 0)
        row[f"{grp}_ANY_pos"] = float(anymk.mean())
        allmk = np.ones(len(sub), bool)
        for g in MARKERS:
            if g in sub: allmk &= (sub[g].values > 0)
        row[f"{grp}_ALL_pos"] = float(allmk.mean())
        for g in EPI_CTRL:
            if g in sub:
                row[f"{grp}_{g}_pos"] = float((sub[g] > 0).mean())
    task2_rows.append(row)

    # ================= TASK 4 : concentration across samples ============
    hb = obs.iloc[high_idx].copy()
    hb["anymk"] = False
    for g in MARKERS:
        if g in cntdf: hb["anymk"] |= (cntdf.loc[hb.index, g].values > 0)
    drivers = hb[hb["anymk"]]
    vc = drivers["skey"].value_counts()
    vc_all = hb["skey"].value_counts()
    conc_rows.append({
        "cancer_type": c,
        "n_high_decile": len(hb),
        "n_high_decile_marker_pos": int(len(drivers)),
        "n_samples_with_high_decile_cells": int((vc_all > 0).sum()),
        "n_samples_hosting_marker_pos": int((vc > 0).sum()),
        "top_sample_share_of_marker_pos": float(vc.iloc[0] / len(drivers)) if len(drivers) else np.nan,
        "top3_sample_share_of_marker_pos": float(vc.iloc[:3].sum() / len(drivers)) if len(drivers) else np.nan,
        "top_sample": vc.index[0] if len(drivers) else None,
    })

per_sample = pd.concat(per_sample_all, ignore_index=True)
per_sample.to_csv(f"{OUT}/task1_per_sample.csv", index=False)
cod = pd.DataFrame(codetail_all)
cod.to_csv(f"{OUT}/task1_correlations_by_cancer.csv", index=False)
t2 = pd.DataFrame(task2_rows)
t2.to_csv(f"{OUT}/task2_coexpression.csv", index=False)
conc = pd.DataFrame(conc_rows)
conc.to_csv(f"{OUT}/task4_concentration.csv", index=False)

# ---- POOLED correlations -------------------------------------------------
pool = per_sample.dropna(subset=["myeloid_frac", "mean_prog_share"]).copy()
pool = pool[pool.n_epi_cells >= 10]
pooled_rows = []
for yname in ["myeloid_frac", "immune_frac", "lymphoid_frac", "bcell_frac", "tcell_frac"]:
    x = pool["mean_prog_share"].values; y = pool[yname].values
    sr, sp_ = spearmanr(x, y); pr, pp = pearsonr(x, y)
    pooled_rows.append(dict(comparison=f"POOLED:usage~{yname}", n=len(pool),
                            spearman_r=sr, spearman_p=sp_, pearson_r=pr, pearson_p=pp))
# pooled, rank-within-cancer (removes cross-cancer offset)
pool2 = pool.copy()
pool2["u_rank"] = pool2.groupby("cancer_type")["mean_prog_share"].rank()
pool2["m_rank"] = pool2.groupby("cancer_type")["myeloid_frac"].rank()
sr, sp_ = spearmanr(pool2["u_rank"], pool2["m_rank"])
pooled_rows.append(dict(comparison="POOLED(within-cancer-ranked):usage~myeloid_frac",
                        n=len(pool2), spearman_r=sr, spearman_p=sp_,
                        pearson_r=np.nan, pearson_p=np.nan))
pooled = pd.DataFrame(pooled_rows)
pooled.to_csv(f"{OUT}/task1_correlations_pooled.csv", index=False)

pd.set_option("display.width", 200); pd.set_option("display.max_columns", 40)
print("\n\n########## TASK 1: correlations by cancer ##########")
print(cod[["cancer_type","comparison","n","spearman_r","spearman_p","pearson_r","pearson_p"]].to_string(index=False))
print("\n########## TASK 1: pooled ##########")
print(pooled.to_string(index=False))
print("\n########## TASK 2: co-expression (fractions) ##########")
show = ["cancer_type","program","n_epi","n_decile",
        "high_PTPRC_pos","low_PTPRC_pos","high_CD68_pos","low_CD68_pos",
        "high_LYZ_pos","low_LYZ_pos","high_ANY_pos","low_ANY_pos",
        "high_ALL_pos","low_ALL_pos","high_EPCAM_pos","low_EPCAM_pos"]
print(t2[[x for x in show if x in t2.columns]].to_string(index=False))
print("\n########## TASK 2: co-expression at count>=2 ##########")
show2 = ["cancer_type","high_PTPRC_ge2","low_PTPRC_ge2","high_CD68_ge2","low_CD68_ge2","high_LYZ_ge2","low_LYZ_ge2"]
print(t2[[x for x in show2 if x in t2.columns]].to_string(index=False))
print("\n########## TASK 4: concentration across samples ##########")
print(conc.to_string(index=False))
print("\nwrote:", OUT)
