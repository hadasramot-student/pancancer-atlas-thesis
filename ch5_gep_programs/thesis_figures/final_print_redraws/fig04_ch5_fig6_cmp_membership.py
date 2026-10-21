"""ch5_fig6 (thesis Figure 6.6): consensus meta-program (cMP) membership across cancer types.

Data: ch5_gep_programs/meta_programs/meta_program_catalogue_J_gt_0.10.csv (22 cMPs built from
the 240 cNMF programs, Jaccard > 0.10; written by meta_program_catalogue.py).
A dot = the cancer type contributes >= 1 cNMF program to that cMP; dot size = number of
contributed programs; colour = cMP class (pan-cancer, >= 7 cancer types, vs. cancer-type-specific).
The original plotting code was not saved on the cluster; this redraw uses the same table.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from printstyle import save, WIDTH, plain

CAT = "/groups/ofircohen-group/users/hadasra/pancancer-atlas-thesis/ch5_gep_programs/meta_programs/meta_program_catalogue_J_gt_0.10.csv"
cat = pd.read_csv(CAT)
CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney", "Liver_Biliary",
           "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"]
COL = {"pan-cancer": "#b2182b", "cancer-type-specific": "#666666"}   # grey at #666666 for print

rows = []
for i, r in cat.iterrows():
    for prog in r.member_cnmf_programs.split(";"):
        rows.append((i, r.meta_program_id, prog.rsplit("_cNMF_", 1)[0], r.classification))
mem = pd.DataFrame(rows, columns=["i", "cmp", "cancer", "cls"]).groupby(["i", "cmp", "cancer", "cls"]).size().reset_index(name="n")
print(mem.groupby("cmp", sort=False).n.sum().to_string())

def size(n):
    return 14 + 26 * (n - 1)

fig, ax = plt.subplots(figsize=(WIDTH, 5.3))
fig.subplots_adjust(left=0.145, right=0.75, top=0.98, bottom=0.215)
ax.scatter([CANCERS.index(c) for c in mem.cancer], mem.i, s=[size(n) for n in mem.n],
           c=[COL[c] for c in mem.cls], edgecolors="black", linewidths=0.4, zorder=3)
ax.set_xticks(range(len(CANCERS)))
ax.set_xticklabels([plain(c) for c in CANCERS], rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(cat))); ax.set_yticklabels(cat.meta_program_id)
ax.set_xlim(-0.6, len(CANCERS) - 0.4); ax.set_ylim(len(cat) - 0.5, -0.5)
ax.grid(color="#c8c8c8", linewidth=0.4, zorder=0)
ax.set_xlabel("Cancer type"); ax.set_ylabel("Consensus meta-program")

nmax = int(mem.n.max())
h_size = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="black",
                 markersize=size(n) ** 0.5, label=str(n)) for n in range(1, nmax + 1)]
h_cls = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=COL[k], markeredgecolor="black",
                markersize=6, label=lab) for k, lab in [("pan-cancer", "Pan-cancer\n(≥ 7 types)"),
                                                         ("cancer-type-specific", "Cancer-type-\nspecific")]]
l1 = fig.legend(handles=h_size, title="cNMF programs", loc="upper left", bbox_to_anchor=(0.765, 0.98),
                frameon=False, labelspacing=0.9, alignment="left")
fig.legend(handles=h_cls, title="cMP class", loc="upper left", bbox_to_anchor=(0.765, 0.62),
           frameon=False, labelspacing=0.8, alignment="left")
save(fig, "ch5_fig6")
