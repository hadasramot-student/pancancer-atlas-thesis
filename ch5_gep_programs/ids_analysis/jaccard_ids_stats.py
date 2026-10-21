#!/usr/bin/env python
"""Point 11 sub-issues 11e / 11g / 11f -- recompute numbers from existing Ch.5 output.

11e  Jaccard vs IDS correlation across all cNMF-GeneNMF program pairs (orthogonality test)
11g  Hallmark enrichment: BH-correct within each program (50 tests), recount annotated programs
11f  Kidney GeneNMF MP consolidation recount
"""
import os, json, glob
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ATLAS   = "/groups/ofircohen-group/users/hadasra/temp_objects"
CNMF    = os.path.join(ATLAS, "cnmf_output")
GENENMF = os.path.join(ATLAS, "genenmf_output")
IDSOUT  = os.path.join(ATLAS, "ids_output")
HALLMARK_JSON = "/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/merge/spectra/priors/hallmark.json"
HC85 = "/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/ids_analysis/hc85_ids_pairs.csv"
OUT  = os.path.dirname(os.path.abspath(__file__))

CANCERS_XM = ["Brain","Breast","Colorectal","Kidney","Liver_Biliary","Lung",
              "Ovarian","Pancreas","Prostate","Sarcoma","Skin"]           # have GeneNMF + IDS
CANCERS_HALL = ["Brain","Breast","Colorectal","Head_and_Neck","Kidney","Liver_Biliary",
                "Lung","Ovarian","Pancreas","Prostate","Sarcoma","Skin"]  # 12 -> 240 programs

N_TOP = 50
MIN_OVERLAP = 3
PCUT = 0.05

# ---------------------------------------------------------------------------
def bh_correct(pvals):
    """Benjamini-Hochberg FDR, same routine as meta_program_catalogue.py."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out

def jacc(a, b):
    a, b = set(a), set(b)
    u = len(a | b)
    return len(a & b) / u if u else 0.0

# ===========================================================================
# 11e -- Jaccard vs IDS across ALL cNMF-GeneNMF program pairs
# ===========================================================================
def run_11e():
    rows = []
    for c in CANCERS_XM:
        jf = os.path.join(GENENMF, c, f"{c}_jaccard_MP_vs_cNMF.csv")
        idf = os.path.join(IDSOUT, c, f"{c}_ids_pairs.csv")
        J = pd.read_csv(jf, index_col=0)                       # MP (rows) x cNMF (cols)
        Jl = J.reset_index().melt(id_vars="index", var_name="cnmf", value_name="jaccard")
        Jl = Jl.rename(columns={"index": "mp"})

        ids = pd.read_csv(idf)
        ids = ids[ids["pair_type"] == "GeneNMF_cNMF"].copy()
        def norm(r):
            a, b = r["program_a"], r["program_b"]
            cn = a if a.startswith("cNMF") else b
            mp = b if a.startswith("cNMF") else a
            return pd.Series({"cnmf": cn, "mp": mp, "ids": r["ids"]})
        ids = ids.apply(norm, axis=1)

        m = Jl.merge(ids, on=["cnmf", "mp"], how="inner")
        m["cancer"] = c
        rows.append(m)
    allpairs = pd.concat(rows, ignore_index=True)
    allpairs.to_csv(os.path.join(OUT, "11e_all_cnmf_genenmf_pairs.csv"), index=False)

    def corr_block(df, label):
        x, y = df["jaccard"].values, df["ids"].values
        pr, pp = pearsonr(x, y)
        sr, sp = spearmanr(x, y)
        print(f"  {label:38s} N={len(df):5d}  Pearson r={pr:+.3f} (p={pp:.2e})   "
              f"Spearman rho={sr:+.3f} (p={sp:.2e})")
        return dict(set=label, N=len(df), pearson_r=pr, pearson_p=pp,
                    spearman_rho=sr, spearman_p=sp)

    print("\n===== 11e : Jaccard vs IDS correlation =====")
    summ = []
    summ.append(corr_block(allpairs, "ALL cNMF-GeneNMF pairs (11 cancers)"))
    summ.append(corr_block(allpairs[~allpairs.cancer.isin(["Breast","Sarcoma"])],
                           "excl. Breast & Sarcoma (thesis 5.3.4 set)"))
    # restricted to HC pairs only (expanded 84-pair set)
    hc = pd.read_csv(HC85)
    summ.append(corr_block(hc.rename(columns={"jaccard":"jaccard","IDS":"ids"}),
                           "HC pairs only (Jaccard>0.0538, n=84)"))
    # per-cancer
    print("\n  per-cancer (all pairs):")
    percanc = []
    for c, g in allpairs.groupby("cancer"):
        pr, pp = pearsonr(g.jaccard, g.ids)
        sr, sp = spearmanr(g.jaccard, g.ids)
        print(f"    {c:15s} N={len(g):4d}  Pearson r={pr:+.3f} (p={pp:.1e})  Spearman rho={sr:+.3f} (p={sp:.1e})")
        percanc.append(dict(cancer=c, N=len(g), pearson_r=pr, pearson_p=pp,
                            spearman_rho=sr, spearman_p=sp))
    pd.DataFrame(summ).to_csv(os.path.join(OUT, "11e_correlation_summary.csv"), index=False)
    pd.DataFrame(percanc).to_csv(os.path.join(OUT, "11e_correlation_per_cancer.csv"), index=False)

    # scatter
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    cancers = sorted(allpairs.cancer.unique())
    cmap = plt.cm.tab20(np.linspace(0, 1, len(cancers)))
    for ax, (xs, title) in zip(axes, [("linear", "linear x"), ("log", "log-x (Jaccard>0)")]):
        for col, c in zip(cmap, cancers):
            g = allpairs[allpairs.cancer == c]
            ax.scatter(g.jaccard, g.ids, s=14, color=col, alpha=.65, label=c, edgecolors="none")
        ax.set_xlabel("pairwise Jaccard (top-50 gene lists)")
        ax.set_ylabel("pairwise IDS")
        ax.set_title(f"cNMF-GeneNMF program pairs ({title})")
        if xs == "log":
            ax.set_xscale("log")
        ax.grid(alpha=.25)
    # regression line + stats on linear panel
    pr, pp = pearsonr(allpairs.jaccard, allpairs.ids)
    sr, sp = spearmanr(allpairs.jaccard, allpairs.ids)
    b1, b0 = np.polyfit(allpairs.jaccard, allpairs.ids, 1)
    xx = np.linspace(0, allpairs.jaccard.max(), 50)
    axes[0].plot(xx, b0 + b1 * xx, "k--", lw=1.5)
    axes[0].text(.98, .03,
                 f"N={len(allpairs)}\nPearson r={pr:.3f}\nSpearman ρ={sr:.3f}",
                 transform=axes[0].transAxes, ha="right", va="bottom",
                 bbox=dict(boxstyle="round", fc="white", alpha=.85))
    axes[1].legend(fontsize=7, ncol=2, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "11e_jaccard_vs_ids_scatter.png"), dpi=140)
    print(f"\n  scatter -> {os.path.join(OUT, '11e_jaccard_vs_ids_scatter.png')}")
    return allpairs

# ===========================================================================
# 11g -- Hallmark enrichment BH correction within each program
# ===========================================================================
def run_11g():
    with open(HALLMARK_JSON) as f:
        H = json.load(f)
    H = {k.replace("HALLMARK_", ""): set(v) for k, v in H.items()}
    names = list(H)

    recs = []
    pmat_rows = []
    for c in CANCERS_HALL:
        gs = pd.read_csv(os.path.join(CNMF, c, f"cnmf_{c}",
                        f"cnmf_{c}.gene_spectra_score.k_20.dt_0_1.txt"),
                        sep="\t", index_col=0)
        nbg = gs.shape[1]
        for prog in gs.index:
            top = set(gs.loc[prog].sort_values(ascending=False).head(N_TOP).index)
            pv = []
            for hn in names:
                ov = len(top & H[hn])
                if ov < MIN_OVERLAP:
                    pv.append(1.0)
                else:
                    pv.append(hypergeom.sf(ov - 1, nbg, len(H[hn]), N_TOP))
            pv = np.array(pv)
            qv = bh_correct(pv)
            recs.append(dict(cancer=c, program=f"cNMF_{prog}",
                             min_p=pv.min(), min_q=qv.min(),
                             n_p_lt_05=int((pv < PCUT).sum()),
                             n_q_lt_05=int((qv < PCUT).sum()),
                             best_hallmark=names[int(np.argmin(pv))]))
            pmat_rows.append(dict(cancer=c, program=f"cNMF_{prog}",
                                  **{names[i]: pv[i] for i in range(50)}))
    df = pd.DataFrame(recs)
    pd.DataFrame(pmat_rows).to_csv(os.path.join(OUT, "11g_program_x_hallmark_pvalues.csv"), index=False)
    df.to_csv(os.path.join(OUT, "11g_per_program_summary.csv"), index=False)

    old = (df.min_p < PCUT).sum()
    new = (df.min_q < PCUT).sum()
    n = len(df)
    print("\n===== 11g : Hallmark BH correction =====")
    print(f"  programs total                         : {n}")
    print(f"  >=1 Hallmark at uncorrected p<0.05     : {old}  ({old/n*100:.1f}%)   [thesis: 109 / 45.4%]")
    print(f"  >=1 Hallmark surviving BH q<0.05       : {new}  ({new/n*100:.1f}%)")
    print(f"  lost under BH                          : {old - new}")

    print("\n  distribution: #programs by #Hallmark hits surviving BH (q<0.05):")
    vc = df.n_q_lt_05.value_counts().sort_index()
    for k, v in vc.items():
        print(f"    {k:2d} surviving hit(s) : {v:3d} programs")
    print("\n  (for reference) #programs by #Hallmark hits at uncorrected p<0.05:")
    vc2 = df.n_p_lt_05.value_counts().sort_index()
    for k, v in vc2.items():
        print(f"    {k:2d} hit(s) : {v:3d} programs")

    print("\n  per-cancer (annotated old p<0.05  ->  new BH q<0.05, out of 20):")
    for c, g in df.groupby("cancer"):
        print(f"    {c:15s} {int((g.min_p<PCUT).sum()):2d}  -> {int((g.min_q<PCUT).sum()):2d}")
    return df

# ===========================================================================
# 11f -- Kidney GeneNMF MP consolidation
# ===========================================================================
def run_11f():
    g = pd.read_csv(os.path.join(GENENMF, "Kidney", "Kidney_genenmf_MP_top50_genes.csv"))
    mps = {mp: sub.sort_values("rank").gene.tolist() for mp, sub in g.groupby("meta_program")}
    order = sorted(mps, key=lambda m: int(m[2:]))
    Jm = pd.DataFrame(index=order, columns=order, dtype=float)
    for a in order:
        for b in order:
            Jm.loc[a, b] = 1.0 if a == b else jacc(mps[a], mps[b])
    Jm.to_csv(os.path.join(OUT, "11f_kidney_MP_x_MP_jaccard.csv"))

    iu = [(order[i], order[j], Jm.loc[order[i], order[j]])
          for i in range(len(order)) for j in range(i + 1, len(order))]
    iu.sort(key=lambda t: -t[2])
    print("\n===== 11f : Kidney GeneNMF MP consolidation =====")
    print(f"  MP-MP top-50 Jaccard: max={iu[0][2]:.3f}  ({iu[0][0]}-{iu[0][1]})")
    print("  top 12 MP-MP pairs:")
    for a, b, v in iu[:12]:
        print(f"    {a:5s} {b:5s} {v:.3f}")
    offdiag = np.array([v for _, _, v in iu])
    print(f"  off-diagonal Jaccard quantiles: "
          f"p50={np.percentile(offdiag,50):.3f} p75={np.percentile(offdiag,75):.3f} "
          f"p90={np.percentile(offdiag,90):.3f} max={offdiag.max():.3f}")

    # HC Kidney MPs (expanded threshold, from hc85 file)
    hc = pd.read_csv(HC85)
    kid_hc = sorted(hc[hc.cancer_type == "Kidney"].GeneNMF_program_id.unique(),
                    key=lambda m: int(m[2:]))
    print(f"\n  Kidney HC MPs (Jaccard>0.0538): n={len(kid_hc)} -> {kid_hc}")

    def consolidate(thr, universe):
        parent = {m: m for m in universe}
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        def union(a, b):
            parent[find(a)] = find(b)
        for a, b, v in iu:
            if a in universe and b in universe and v >= thr:
                union(a, b)
        comp = {}
        for m in universe:
            comp.setdefault(find(m), []).append(m)
        return list(comp.values())

    for thr in (0.2, 0.25, 0.3, 0.4):
        for uni_name, uni in [("all 19 MPs", order), (f"{len(kid_hc)} HC MPs", kid_hc)]:
            groups = consolidate(thr, uni)
            multi = [sorted(x, key=lambda m:int(m[2:])) for x in groups if len(x) > 1]
            print(f"\n  thr>={thr}  [{uni_name}]: {len(uni)} MPs -> {len(groups)} groups "
                  f"({len(multi)} multi-MP)")
            for mm in sorted(multi, key=lambda x: -len(x)):
                print(f"      {{{', '.join(mm)}}}")

    # headline recount at chosen thresholds, on HC set
    print("\n  --- Kidney concordance recount (HC set) ---")
    print(f"  before consolidation: {len(kid_hc)} HC Kidney MPs "
          f"(thesis v88 wording: 16/19; expanded-threshold: {len(kid_hc)}/19)")
    for thr in (0.25, 0.3, 0.4):
        groups = consolidate(thr, kid_hc)
        print(f"  after consolidation thr>={thr}: {len(kid_hc)} HC MPs collapse into "
              f"{len(groups)} independent groups")

if __name__ == "__main__":
    run_11e()
    run_11g()
    run_11f()
    print("\nDONE. Outputs in", OUT)
