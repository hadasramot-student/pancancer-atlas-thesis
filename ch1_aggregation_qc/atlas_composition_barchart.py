# atlas_composition_barchart.py
#
# Ch.1 (Aggregation & QC) - atlas composition figure: per-cancer-type
# meta-type (consensus cell-type) composition as a stacked horizontal bar
# chart.
#
# Inputs: none (composition counts are the published per-cancer-type x
# meta-type totals for this atlas; update `data` below if the atlas
# composition changes, or wire it up to Table3_metatype_composition.csv for
# a data-driven version).
# Outputs: <ATLAS_FIGURES_DIR>/Fig_metatype_per_cancer.{png,pdf}
#
# Usage:
#   ATLAS_FIGURES_DIR=/path/to/plots_and_outputs python atlas_composition_barchart.py

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

data = [
    ("Brain","Endothelial",3252),("Brain","Fibroblasts",15998),("Brain","Hematopoietic_Stem_and_Progenitor_cells",5478),("Brain","Myeloid",20678),("Brain","T-cells",4262),("Brain","other",125715),
    ("Breast","B-cell",81397),("Breast","Endothelial",74867),("Breast","Epithelial",1164224),("Breast","Fibroblasts",102742),("Breast","Hematopoietic_Stem_and_Progenitor_cells",6547),("Breast","Myeloid",142742),("Breast","NK-cells",5821),("Breast","T-cells",368672),("Breast","other",984910),
    ("Colorectal","B-cell",85312),("Colorectal","Endothelial",14918),("Colorectal","Epithelial",397587),("Colorectal","Fibroblasts",14963),("Colorectal","Hematopoietic_Stem_and_Progenitor_cells",1484),("Colorectal","Myeloid",95072),("Colorectal","NK-cells",9541),("Colorectal","T-cells",235316),("Colorectal","other",128284),
    ("Head_and_Neck","B-cell",20819),("Head_and_Neck","Endothelial",65),("Head_and_Neck","Epithelial",5978),("Head_and_Neck","Myeloid",4909),("Head_and_Neck","NK-cells",1028),("Head_and_Neck","T-cells",35304),("Head_and_Neck","other",3793),
    ("Kidney","B-cell",6069),("Kidney","Endothelial",21798),("Kidney","Epithelial",99920),("Kidney","Fibroblasts",1244),("Kidney","Hematopoietic_Stem_and_Progenitor_cells",102),("Kidney","Myeloid",65401),("Kidney","NK-cells",11093),("Kidney","T-cells",112260),("Kidney","other",53717),
    ("Liver_Biliary","B-cell",7035),("Liver_Biliary","Endothelial",11649),("Liver_Biliary","Epithelial",23817),("Liver_Biliary","Fibroblasts",354),("Liver_Biliary","Myeloid",19750),("Liver_Biliary","NK-cells",12430),("Liver_Biliary","T-cells",93666),("Liver_Biliary","other",19517),
    ("Lung","B-cell",50469),("Lung","Endothelial",62138),("Lung","Epithelial",203074),("Lung","Fibroblasts",11940),("Lung","Hematopoietic_Stem_and_Progenitor_cells",210),("Lung","Myeloid",345764),("Lung","NK-cells",29568),("Lung","T-cells",242879),("Lung","other",112243),
    ("Neuroendocrine","B-cell",957),("Neuroendocrine","Endothelial",1200),("Neuroendocrine","Epithelial",2184),("Neuroendocrine","Fibroblasts",1481),("Neuroendocrine","Hematopoietic_Stem_and_Progenitor_cells",8),("Neuroendocrine","Myeloid",4920),("Neuroendocrine","NK-cells",793),("Neuroendocrine","T-cells",2229),("Neuroendocrine","other",115787),
    ("Other_Models","B-cell",32395),("Other_Models","Endothelial",14547),("Other_Models","Epithelial",102548),("Other_Models","Fibroblasts",6040),("Other_Models","Myeloid",25369),("Other_Models","NK-cells",84),("Other_Models","T-cells",63375),("Other_Models","other",67093),
    ("Ovarian","B-cell",8960),("Ovarian","Endothelial",7442),("Ovarian","Epithelial",119078),("Ovarian","Fibroblasts",31641),("Ovarian","Hematopoietic_Stem_and_Progenitor_cells",23),("Ovarian","Myeloid",67617),("Ovarian","NK-cells",6411),("Ovarian","T-cells",53300),("Ovarian","other",88079),
    ("Pancreas","B-cell",4701),("Pancreas","Endothelial",23305),("Pancreas","Epithelial",134687),("Pancreas","Fibroblasts",60623),("Pancreas","Hematopoietic_Stem_and_Progenitor_cells",16),("Pancreas","Myeloid",41197),("Pancreas","NK-cells",564),("Pancreas","T-cells",41900),("Pancreas","other",77684),
    ("Prostate","B-cell",1451),("Prostate","Endothelial",6680),("Prostate","Epithelial",104740),("Prostate","Fibroblasts",1904),("Prostate","Hematopoietic_Stem_and_Progenitor_cells",27),("Prostate","Myeloid",12703),("Prostate","NK-cells",2006),("Prostate","T-cells",15834),("Prostate","other",65371),
    ("Sarcoma","B-cell",102),("Sarcoma","Endothelial",1808),("Sarcoma","Epithelial",1558),("Sarcoma","Fibroblasts",6041),("Sarcoma","Hematopoietic_Stem_and_Progenitor_cells",1785),("Sarcoma","Myeloid",17112),("Sarcoma","NK-cells",138),("Sarcoma","T-cells",2345),("Sarcoma","other",49793),
    ("Skin","B-cell",13320),("Skin","Endothelial",1804),("Skin","Epithelial",55675),("Skin","Fibroblasts",5296),("Skin","Hematopoietic_Stem_and_Progenitor_cells",1326),("Skin","Myeloid",58539),("Skin","NK-cells",1531),("Skin","T-cells",107986),("Skin","other",123735),
    ("unknown","B-cell",57611),("unknown","Endothelial",15815),("unknown","Epithelial",99870),("unknown","Fibroblasts",10121),("unknown","Hematopoietic_Stem_and_Progenitor_cells",15429),("unknown","Myeloid",57678),("unknown","NK-cells",7337),("unknown","T-cells",143272),("unknown","other",155385),
]

df = pd.DataFrame(data, columns=["cancer_type","cell_type","n"])
label_map = {
    "Epithelial":"Epithelial/Cancer","T-cells":"T-cells","other":"Other/Unknown",
    "Myeloid":"Myeloid","B-cell":"B-cells","Fibroblasts":"Fibroblasts/Stroma",
    "Endothelial":"Endothelial/Vascular","Hematopoietic_Stem_and_Progenitor_cells":"HSC/Stem-like","NK-cells":"NK-cells",
}
df["cell_type"] = df["cell_type"].map(label_map)

pivot = df.pivot_table(index="cancer_type", columns="cell_type", values="n", aggfunc="sum").fillna(0)
pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100

order = ["Epithelial/Cancer","T-cells","Myeloid","B-cells","NK-cells",
         "Fibroblasts/Stroma","Endothelial/Vascular","HSC/Stem-like","Other/Unknown"]
pivot_pct = pivot_pct[order]
pivot_pct = pivot_pct.sort_values("Epithelial/Cancer", ascending=True)

colors = ["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#B09C85"]

fig, ax = plt.subplots(figsize=(13, 7))
bottoms = np.zeros(len(pivot_pct))
for i, col in enumerate(order):
    vals = pivot_pct[col].values
    ax.barh(range(len(pivot_pct)), vals, left=bottoms, color=colors[i], label=col, height=0.75)
    bottoms += vals

ax.set_yticks(range(len(pivot_pct)))
ax.set_yticklabels(pivot_pct.index, fontsize=11)
ax.set_xlabel("% of cells", fontsize=12)
ax.set_xlim(0, 100)
ax.set_title("Cell type composition per cancer type", fontsize=14, fontweight="bold", pad=12)
ax.spines[['top','right']].set_visible(False)
ax.legend(handles=[mpatches.Patch(color=colors[i], label=order[i]) for i in range(len(order))],
          bbox_to_anchor=(1.01,1), loc='upper left', fontsize=10, frameon=False,
          title="Meta-type", title_fontsize=11)

plt.tight_layout()
out_dir = os.environ.get("ATLAS_FIGURES_DIR", "outputs/figures")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "Fig_metatype_per_cancer.png")
plt.savefig(out, dpi=300, bbox_inches="tight")
plt.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
print("Saved:", out)
