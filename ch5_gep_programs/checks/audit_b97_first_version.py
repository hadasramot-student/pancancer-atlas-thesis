#!/usr/bin/env python
"""Table B97 contamination audit: reproduce the 22-cMP Gavish-style catalogue
(Jaccard>0.1, recurrence>=2 cancer types), regenerate full consensus top-50
gene lists, score each cMP against lineage-contamination marker panels, and
independently recompute the Gavish-2023 Cancer-MP best-Jaccard for all 22."""
import os, re
import numpy as np
import pandas as pd

CNMF_BASE = "/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
GAVISH_CSV = "/groups/ofircohen-group/users/hadasra/temp_objects/gavish_2023/gavish_MP_gene_lists_long.csv"
CATALOGUE_CSV = "/groups/ofircohen-group/users/hadasra/pancancer-atlas-thesis/ch5_gep_programs/meta_programs/meta_program_catalogue_J_gt_0.10.csv"
OUTDIR = os.path.dirname(os.path.abspath(__file__))

K_PROGRAMS = 20; DT = "dt_0_1"; N_TOP_GENES = 50; THRESHOLD = 0.1; MIN_CANCER_TYPES = 2
CANCERS = ["Brain","Breast","Colorectal","Head_and_Neck","Kidney","Liver_Biliary",
           "Lung","Ovarian","Pancreas","Prostate","Sarcoma","Skin"]

MARKER_PANELS = {
 "myeloid":["PTPRC","CD68","LYZ","C1QB","C1QA","C1QC","LAPTM5","AIF1","TYROBP","FCER1G",
            "CD14","CD163","CSF1R","ITGAM","FCGR3A","MRC1","APOC1","SPP1","MARCO","VSIG4"],
 "myeloid_core":["PTPRC","CD68","LYZ","C1QB","C1QA","LAPTM5","AIF1"],
 "T_NK":["CD3D","CD3E","CD3G","CD2","TRBC1","TRBC2","TRAC","NKG7","GNLY","GZMB","GZMA",
         "GZMK","KLRD1","CCL5","IL7R","CD8A","CD8B","CD7","PRF1","KLRB1"],
 "T_NK_core":["CD3D","CD3E","CD2","TRBC1","NKG7","GNLY"],
 "B_plasma":["MS4A1","CD79A","CD79B","CD19","CD37","MZB1","JCHAIN","DERL3","TNFRSF17","FCRLA",
             "CXCR5","BANK1","VPREB3","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2","IGHM",
             "IGHD","IGKC","IGLC1","IGLC2","IGLC3","IGLL5","IGLL1"],
 "B_plasma_core":["MS4A1","CD79A","CD79B"],
 "fibroblast":["COL1A1","COL1A2","DCN","LUM","PDGFRB","PDGFRA","COL3A1","COL6A1","COL6A2",
               "COL6A3","COL5A1","COL5A2","FN1","SPARC","POSTN","THBS2","FBLN1","MMP2","AEBP1",
               "CTHRC1","FSTL1","COL4A1","COL4A2","TAGLN","ACTA2","FAP","PDPN","COL10A1",
               "COL11A1","MGP","COL8A1"],
 "fibroblast_core":["COL1A1","COL1A2","DCN","LUM","PDGFRB","PDGFRA"],
 "endothelial":["PECAM1","VWF","CDH5","CLDN5","PLVAP","EGFL7","RAMP2","AQP1","CLEC14A","ENG",
                "FLT1","KDR","EMCN","ESAM","ERG","TIE1"],
 "endothelial_core":["PECAM1","VWF","CDH5"],
}
def is_mito(g): return g.startswith("MT-")
def is_ribo(g): return bool(re.match(r"^RP[LS]P?\d", g)) or g in ("RPLP0","RPLP1","RPLP2","RPSA")
def is_ig(g):   return bool(re.match(r"^IG[HKL][VDJACGM]", g)) or g in ("IGKC","IGHM","IGHD","IGHE","JCHAIN")

def load_score(c):
    return pd.read_csv(os.path.join(CNMF_BASE,c,f"cnmf_{c}",
        f"cnmf_{c}.gene_spectra_score.k_{K_PROGRAMS}.{DT}.txt"), sep="\t", index_col=0)

def load_all():
    progs, nbg = {}, None
    for c in CANCERS:
        df = load_score(c); nbg = df.shape[1] if nbg is None else nbg
        assert df.shape[1]==nbg
        for pidx in df.index:
            top = list(df.loc[pidx].sort_values(ascending=False).head(N_TOP_GENES).index)
            progs[f"{c}_cNMF_{pidx}"] = {"cancer":c,"idx":pidx,"top_genes":top,"top_set":set(top)}
    return progs, nbg

def jac(a,b):
    u=len(a|b); return len(a&b)/u if u else 0.0

def build_J(progs):
    keys=list(progs); sets=[progs[k]["top_set"] for k in keys]; n=len(keys)
    J=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            v=jac(sets[i],sets[j]); J[i,j]=J[j,i]=v
    return keys,J

def greedy(keys,J,t):
    n=len(keys); pool=set(range(n)); groups=[]
    while pool:
        pl=list(pool); bi,bc,bp=None,-1,None
        for i in pl:
            partners=[j for j in pl if j!=i and J[i,j]>t]
            if len(partners)>bc: bi,bc,bp=i,len(partners),partners
        if bc<=0: break
        g=[bi]+bp; groups.append(g); pool-=set(g)
    return groups, sorted(pool)

def consensus(progs,keys,gidx):
    freq,rs={}, {}
    for idx in gidx:
        for rank,g in enumerate(progs[keys[idx]]["top_genes"]):
            freq[g]=freq.get(g,0)+1; rs[g]=rs.get(g,0)+rank
    ordered=sorted(freq, key=lambda g:(-freq[g], rs[g]/freq[g]))
    return ordered, freq

def load_gavish():
    df=pd.read_csv(GAVISH_CSV); df=df[df.cell_type=="Cancer"]
    return {mp:set(s.gene) for mp,s in df.groupby("meta_program")}

def gavish_match(ts,gmps):
    best_mp,best_j,allj=None,0.0,{}
    for mp,genes in gmps.items():
        j=jac(ts,genes); allj[mp]=j
        if j>best_j: best_mp,best_j=mp,j
    return best_mp,best_j,allj

def main():
    progs,nbg=load_all(); keys,J=build_J(progs)
    print(f"{len(keys)} cNMF programs, background {nbg} genes")
    groups,singles=greedy(keys,J,THRESHOLD)
    valid=[]
    for g in groups:
        cts=sorted(set(progs[keys[i]]["cancer"] for i in g))
        if len(cts)>=MIN_CANCER_TYPES: valid.append((g,cts))
    valid.sort(key=lambda x:-len(x[0]))
    gmps=load_gavish(); cat=pd.read_csv(CATALOGUE_CSV)

    rows=[]; detail={}
    for num,(g,cts) in enumerate(valid,start=1):
        mp=f"cMP{num}"; members=[keys[i] for i in g]
        ordered,freq=consensus(progs,keys,g); top=ordered[:N_TOP_GENES]; ts=set(top)
        hits={}; hg={}
        for cn in ["myeloid","T_NK","B_plasma","fibroblast","endothelial",
                   "myeloid_core","T_NK_core","B_plasma_core","fibroblast_core","endothelial_core"]:
            panel=set(MARKER_PANELS[cn]); x=[y for y in top if y in panel]
            hits[cn]=len(x); hg[cn]=x
        for nm,fn in [("mito",is_mito),("ribo",is_ribo),("ig",is_ig)]:
            x=[y for y in top if fn(y)]; hits[nm]=len(x); hg[nm]=x
        bmp,bj,allj=gavish_match(ts,gmps)
        cr=cat[cat.meta_program_id==mp].iloc[0]
        rows.append(dict(cMP=mp,n_member_programs=len(members),n_cancer_types=len(cts),
            cancer_types=";".join(cts),myeloid=hits["myeloid"],T_NK=hits["T_NK"],
            B_plasma=hits["B_plasma"],fibroblast=hits["fibroblast"],endothelial=hits["endothelial"],
            mito=hits["mito"],ribo=hits["ribo"],ig=hits["ig"],
            myeloid_core=hits["myeloid_core"],T_NK_core=hits["T_NK_core"],
            B_plasma_core=hits["B_plasma_core"],fibroblast_core=hits["fibroblast_core"],
            endothelial_core=hits["endothelial_core"],
            best_hallmark=cr["best_hallmark_geneset"],hallmark_q=cr["hallmark_BH_qval"],
            gavish_match_recomputed=bmp,gavish_jaccard_recomputed=round(bj,4),
            gavish_jaccard_in_csv=cr["gavish_jaccard"],consensus_top15=";".join(ordered[:15])))
        detail[mp]=dict(members=members,cts=cts,ordered=ordered,freq=freq,hits=hits,hit_genes=hg,
            best_mp=bmp,best_j=bj,allj=allj,csv_row=cr.to_dict())

    out=pd.DataFrame(rows); out.to_csv(f"{OUTDIR}/b97_marker_hits.csv",index=False)
    pd.set_option("display.width",260); pd.set_option("display.max_columns",50)

    print("\n===== reproduced catalogue vs CSV (sanity) =====")
    for r in rows:
        c=cat[cat.meta_program_id==r["cMP"]].iloc[0]
        ok=(r["n_member_programs"]==c["n_member_cnmf_programs"]) and (r["n_cancer_types"]==c["n_distinct_cancer_types"])
        print(f'{r["cMP"]:6s} members {r["n_member_programs"]:2d}(csv {c["n_member_cnmf_programs"]:2d}) '
              f'cts {r["n_cancer_types"]}(csv {c["n_distinct_cancer_types"]}) {"OK" if ok else "*** MISMATCH ***"}')

    print("\n===== MARKER-CATEGORY HITS in consensus top-50 =====")
    print(out[["cMP","n_member_programs","n_cancer_types","myeloid","T_NK","B_plasma",
               "fibroblast","endothelial","mito","ribo","ig"]].to_string(index=False))
    print("\n===== core-marker hits (tight panels) =====")
    print(out[["cMP","myeloid_core","T_NK_core","B_plasma_core","fibroblast_core","endothelial_core"]].to_string(index=False))

    print("\n===== GAVISH sub-threshold confirmation =====")
    for mp in ["cMP6","cMP10","cMP14","cMP20","cMP21"]:
        d=detail[mp]; t3=sorted(d["allj"].items(),key=lambda kv:-kv[1])[:3]
        print(f'{mp}: recomputed best Jaccard = {d["best_j"]:.4f} vs {d["best_mp"]}  '
              f'(CSV {d["csv_row"]["gavish_jaccard"]}, {d["csv_row"]["gavish_2023_MP_match"]})')
        print("      top-3: "+", ".join(f"{m}={j:.4f}" for m,j in t3))

    print("\n===== cMP8 DETAIL =====")
    d=detail["cMP8"]
    print("members:",d["members"]); print("cancer types:",d["cts"])
    print("consensus top-50:\n  "+", ".join(d["ordered"][:50]))
    print("fibroblast hits (top-50):",d["hit_genes"]["fibroblast"])
    print("fibroblast_core hits:",d["hit_genes"]["fibroblast_core"])
    print("myeloid:",d["hit_genes"]["myeloid"]," T_NK:",d["hit_genes"]["T_NK"],
          " B/plasma:",d["hit_genes"]["B_plasma"]," endo:",d["hit_genes"]["endothelial"])

    with open(f"{OUTDIR}/b97_consensus_lists.txt","w") as f:
        for mp in [f"cMP{i}" for i in range(1,len(rows)+1)]:
            d=detail[mp]
            f.write(f'### {mp}  members={len(d["members"])}  cancer_types={len(d["cts"])} ({";".join(d["cts"])})\n')
            f.write("members: "+";".join(d["members"])+"\n")
            f.write("consensus_top50: "+", ".join(d["ordered"][:50])+"\n")
            for cn in ["myeloid","T_NK","B_plasma","fibroblast","endothelial"]:
                f.write(f'  {cn}: n={d["hits"][cn]} genes={d["hit_genes"][cn]} (core={d["hit_genes"][cn+"_core"]})\n')
            f.write(f'  mito: n={d["hits"]["mito"]} {d["hit_genes"]["mito"]}\n')
            f.write(f'  ribo: n={d["hits"]["ribo"]} {d["hit_genes"]["ribo"]}\n')
            f.write(f'  ig:   n={d["hits"]["ig"]} {d["hit_genes"]["ig"]}\n')
            f.write(f'  gavish: best={d["best_mp"]} J={d["best_j"]:.4f}\n\n')
    print(f'\nwrote {OUTDIR}/b97_marker_hits.csv and b97_consensus_lists.txt')

if __name__=="__main__":
    main()
