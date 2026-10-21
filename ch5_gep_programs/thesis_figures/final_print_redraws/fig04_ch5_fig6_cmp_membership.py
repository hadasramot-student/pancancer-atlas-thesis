"""ch5_fig6 (thesis Figure 6.6): consensus meta-program (cMP) membership across cancer types.

Original encoding of the thesis figure:
- one dot per (cMP, cancer type) where the cMP has >= 1 member cNMF program from that cancer type;
- all dots the same size;
- colour = annotation category of the cMP, from the catalogue's own columns:
  Gavish-matched  : gavish_2023_MP_match != "Unmatched" (consensus top-50 vs Gavish MP Jaccard > 0.05)
  Hallmark only   : not Gavish-matched, best_hallmark_geneset != "Unannotated" (hypergeometric p < 0.05)
  Neither         : both unmatched.
Data: ch5_gep_programs/meta_programs/meta_program_catalogue_J_gt_0.10.csv (22 cMPs, written by
meta_program_catalogue.py). The previous rebuild (dot size = programs, colour = pan-cancer) is
fig04_ch5_fig6_alt_cmp_membership.py -> ch5_fig6_alt.png.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from printstyle import save, WIDTH, plain

CAT = "/groups/ofircohen-group/users/hadasra/pancancer-atlas-thesis/ch5_gep_programs/meta_programs/meta_program_catalogue_J_gt_0.10.csv"
cat = pd.read_csv(CAT)
CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney", "Liver_Biliary",
           "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"]
CATS = [("gavish", "#1f77b4", "Gavish-matched\n(Jaccard > 0.05)"),
        ("hallmark", "#ff7f0e", "Hallmark-annotated\nonly"),
        ("neither", "#808080", "Neither")]
COL = {k: c for k, c, _ in CATS}

cat["category"] = ["gavish" if g != "Unmatched" else ("hallmark" if h != "Unannotated" else "neither")
                   for g, h in zip(cat.gavish_2023_MP_match, cat.best_hallmark_geneset)]
rows = []
for i, r in cat.iterrows():
    for c in sorted({p.rsplit("_cNMF_", 1)[0] for p in r.member_cnmf_programs.split(";")}, key=CANCERS.index):
        rows.append((i, c, r.category))
mem = pd.DataFrame(rows, columns=["i", "cancer", "category"])
for k, _, lab in CATS:
    print(f"{lab.replace(chr(10), ' ')}: {', '.join(cat.meta_program_id[cat.category == k])}")
for i, r in cat.iterrows():
    print(r.meta_program_id, ":", ", ".join(plain(c) for c in mem.cancer[mem.i == i]))
print("dots:", len(mem))

fig, ax = plt.subplots(figsize=(WIDTH, 5.3))
fig.subplots_adjust(left=0.145, right=0.70, top=0.98, bottom=0.215)
ax.scatter([CANCERS.index(c) for c in mem.cancer], mem.i, s=36, c=[COL[k] for k in mem.category],
           edgecolors="black", linewidths=0.4, zorder=3)
ax.set_xticks(range(len(CANCERS)))
ax.set_xticklabels([plain(c) for c in CANCERS], rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(cat))); ax.set_yticklabels(cat.meta_program_id)
ax.set_xlim(-0.6, len(CANCERS) - 0.4); ax.set_ylim(len(cat) - 0.5, -0.5)
ax.grid(color="#c8c8c8", linewidth=0.4, zorder=0)
ax.set_xlabel("Cancer type"); ax.set_ylabel("Consensus meta-program")
handles = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=c, markeredgecolor="black",
                  markeredgewidth=0.4, markersize=6, label=lab) for _, c, lab in CATS]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.715, 0.98), frameon=False,
           labelspacing=0.8, handletextpad=0.3)
save(fig, "ch5_fig6")
