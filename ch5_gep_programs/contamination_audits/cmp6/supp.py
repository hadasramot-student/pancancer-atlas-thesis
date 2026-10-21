#!/usr/bin/env python
"""Supplementary: (a) mean MHC-II program usage by consensus meta-type *within*
each cNMF epithelial object; (b) extra lineage markers for the co-expression
test (B/T markers), high vs low decile, restricted to consensus-Epithelial."""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

BASE="/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output"
OUT="/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/cmp6-myeloid-contamination-83948d/cmp6_check/out"
PROG={"Brain":"10","Head_and_Neck":"8","Liver_Biliary":"7","Ovarian":"4","Prostate":"8"}
XMARK=["PTPRC","CD68","LYZ","CD3D","CD3E","CD79A","MS4A1","CD19","C1QB","AIF1","EPCAM","KRT8","KRT18"]

byct=[]; deciles=[]
for c,pk in PROG.items():
    A=ad.read_h5ad(f"{BASE}/{c}/cnmf_{c}.h5ad")
    obs=A.obs.copy()
    U=pd.read_csv(f"{BASE}/{c}/cnmf_{c}/cnmf_{c}.usages.k_20.dt_0_1.consensus.txt",sep="\t",index_col=0).loc[obs.index]
    share=(U.div(U.sum(1),axis=0))[pk].values
    obs["share"]=share
    ct = obs["Consensus_Cell_Type"].astype(str)
    if c=="Brain":  # no consensus-epi; show author cell_type instead
        ct = ("author:"+obs["cell_type"].astype(str))
    g=pd.DataFrame({"ct":ct.values,"share":share}).groupby("ct")
    tab=g.agg(n=("share","size"),mean_share=("share","mean"),median_share=("share","median")).reset_index()
    tab.insert(0,"cancer_type",c); tab.insert(1,"program",pk)
    byct.append(tab)

    # decile co-expression on consensus-Epithelial only (Brain: author Malignant)
    if c=="Brain":
        m=obs["cell_type"].astype(str).eq("Malignant").values
    else:
        m=obs["Consensus_Cell_Type"].astype(str).eq("Epithelial").values
    ei=np.where(m)[0]; sh=share[ei]; order=np.argsort(sh); nd=max(1,len(ei)//10)
    lo=ei[order[:nd]]; hi=ei[order[-nd:]]
    Xc=A.X.tocsc() if sp.issparse(A.X) else A.X
    idx={g:A.var_names.get_loc(g) for g in XMARK if g in A.var_names}
    def frac(rows,gene,thr):
        col=Xc[rows,idx[gene]]
        v=np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
        return float((v>=thr).mean())
    r={"cancer_type":c,"program":pk,"n_epi":int(m.sum()),"n_decile":nd,
       "high_share_med":float(np.median(share[hi])),"low_share_med":float(np.median(share[lo]))}
    for g in idx:
        r[f"hi_{g}_pos"]=frac(hi,g,1); r[f"lo_{g}_pos"]=frac(lo,g,1)
    deciles.append(r)

bc=pd.concat(byct,ignore_index=True); bc.to_csv(f"{OUT}/supp_usage_by_metatype.csv",index=False)
dd=pd.DataFrame(deciles); dd.to_csv(f"{OUT}/supp_decile_markers.csv",index=False)
pd.set_option("display.width",250);pd.set_option("display.max_columns",60);pd.set_option("display.max_rows",200)
print("###### mean MHC-II program usage share by consensus meta-type (within the cNMF epithelial object) ######")
print(bc.to_string(index=False))
print("\n###### extra decile markers (fraction with count>0), consensus-Epithelial cells only ######")
cols=["cancer_type","program","n_epi","n_decile","high_share_med","low_share_med",
      "hi_PTPRC_pos","lo_PTPRC_pos","hi_CD79A_pos","lo_CD79A_pos","hi_MS4A1_pos","lo_MS4A1_pos",
      "hi_CD3E_pos","lo_CD3E_pos","hi_C1QB_pos","lo_C1QB_pos","hi_EPCAM_pos","lo_EPCAM_pos"]
print(dd[[x for x in cols if x in dd.columns]].to_string(index=False))
