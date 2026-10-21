#!/usr/bin/env python
# make_ch5_fig5_gavish_dotplot.py  (comment 15/16, 2026-10-21)
#
# Print-size redraw of the Gavish et al. 2023 meta-program dot plot (thesis Fig 6.5).
# Derived from pancancer-atlas-thesis/ch5_gep_programs/hallmark_enrichment/make_gavish_mp_dotplot.py
# Same data, same filtering (matched_0.1), same dot size (n programs, area ratio 1:2)
# and dot colour (max Jaccard per cancer x MP, viridis, min-max normalised).
# Changes: no title; 5.4 in wide, all text >= 8 pt; all row labels black, narrative
# category shown as a coloured square left of the label; "Liver/Biliary";
# legends + colorbar below the plot.

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

IN_CSV = "/groups/ofircohen-group/users/hadasra/temp_objects/gavish_2023/combined_all12_gavish_jaccard.csv"
OUT = "/groups/ofircohen-group/users/hadasra/thesis_figures/comment15_16_2026-10-21/ch5_fig5"

CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney",
           "Liver_Biliary", "Lung", "Ovarian", "Pancreas", "Prostate",
           "Sarcoma", "Skin"]
CANCER_LABEL = {"Head_and_Neck": "Head and Neck", "Liver_Biliary": "Liver/Biliary"}

CATEGORY_COLOR = {"Cell Cycle": "#D55E00", "EMT": "#0072B2", "MHC-II/IFN-gamma": "#009E73"}
CATEGORY_LABEL = {"Cell Cycle": "Cell cycle", "EMT": "EMT", "MHC-II/IFN-gamma": "MHC-II/IFN-γ"}
CATEGORY_ORDER = ["Cell Cycle", "EMT", "MHC-II/IFN-gamma"]
MP_CATEGORY = {
    "MP1  Cell Cycle - G2/M": "Cell Cycle", "MP3  Cell Cylce HMG-rich": "Cell Cycle",
    "MP12 EMT-I": "EMT", "MP13 EMT-II": "EMT", "MP14 EMT-III": "EMT", "MP15 EMT IV": "EMT",
    "MP17 Interferon/MHC-II (I)": "MHC-II/IFN-gamma", "MP18 Interferon/MHC-II (II)": "MHC-II/IFN-gamma",
}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42,
                     "text.color": "black", "axes.labelcolor": "black",
                     "xtick.color": "black", "ytick.color": "black"})


def mp_sort_key(mp):
    return int(mp.strip().split()[0].replace("MP", ""))


def main():
    df = pd.read_csv(IN_CSV)
    m = df[df["matched_0.1"]].copy()
    assert len(m) == 49, len(m)
    dots = (m.groupby(["cancer", "best_gavish_mp"])
              .agg(n_programs=("program", "count"), jaccard=("jaccard", "max"))
              .reset_index())
    mps = sorted(dots["best_gavish_mp"].unique(), key=mp_sort_key)
    assert len(mps) == 17, len(mps)

    def group_rank(mp):
        cat = MP_CATEGORY.get(mp)
        return (CATEGORY_ORDER.index(cat) if cat in CATEGORY_ORDER else 3, mp_sort_key(mp))
    mps = sorted(mps, key=group_rank)
    y_pos = {mp: i for i, mp in enumerate(mps)}
    x_pos = {c: i for i, c in enumerate(CANCERS)}
    dots["x"] = dots["cancer"].map(x_pos)
    dots["y"] = dots["best_gavish_mp"].map(y_pos)

    norm = matplotlib.colors.Normalize(vmin=dots["jaccard"].min(), vmax=dots["jaccard"].max())
    SIZE_1, SIZE_2 = 48, 96   # same 1:2 area ratio as the original (380, 760)

    fig = plt.figure(figsize=(5.4, 5.3))
    ax = fig.add_axes([0.40, 0.30, 0.585, 0.685])
    sc = ax.scatter(dots["x"], dots["y"], c=dots["jaccard"], cmap="viridis", norm=norm,
                    s=dots["n_programs"].map({1: SIZE_1, 2: SIZE_2}).values,
                    edgecolor="black", linewidth=0.5, zorder=3)
    ax.set_xlim(-0.6, len(CANCERS) - 0.4)
    ax.set_ylim(len(mps) - 0.4, -0.6)  # categories at top
    ax.set_xticks(range(len(CANCERS)))
    ax.set_xticklabels([CANCER_LABEL.get(c, c) for c in CANCERS], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8)
    ax.set_yticks(range(len(mps)))
    ax.set_yticklabels([" ".join(mp.split()) for mp in mps], fontsize=8)
    ax.tick_params(axis="y", pad=11, length=2.5, width=0.6)
    ax.tick_params(axis="x", length=2.5, width=0.6)
    ax.grid(True, color="0.85", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_linewidth(0.6)

    # narrative-category squares left of the row labels (no square = other)
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    for mp in mps:
        cat = MP_CATEGORY.get(mp)
        if cat:
            ax.scatter([-0.022], [y_pos[mp]], marker="s", s=22, color=CATEGORY_COLOR[cat],
                       transform=tr, clip_on=False, zorder=4, linewidth=0)

    # colorbar below, left
    cax = fig.add_axes([0.06, 0.075, 0.28, 0.022])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("Jaccard similarity", fontsize=9, labelpad=2)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=8, length=2, width=0.6)
    cb.outline.set_linewidth(0.6)

    size_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=0.5, markersize=SIZE_1 ** 0.5, label="1"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=0.5, markersize=SIZE_2 ** 0.5, label="2")]
    leg1 = fig.legend(handles=size_handles, title="cNMF programs\nmatching the MP",
                      loc="lower left", bbox_to_anchor=(0.40, 0.02), ncol=2, fontsize=8,
                      title_fontsize=8, frameon=False, handletextpad=0.3, columnspacing=1.0)
    cat_handles = [Line2D([0], [0], marker="s", color="none", markerfacecolor=CATEGORY_COLOR[c],
                          markeredgewidth=0, markersize=6, label=CATEGORY_LABEL[c])
                   for c in CATEGORY_ORDER]
    fig.legend(handles=cat_handles, title="Narrative category", loc="lower left",
               bbox_to_anchor=(0.66, 0.0), fontsize=8, title_fontsize=8, frameon=False,
               handletextpad=0.3, labelspacing=0.3, borderaxespad=0.2)
    fig.add_artist(leg1)

    fig.savefig(OUT + ".png", dpi=400)
    fig.savefig(OUT + ".pdf")
    print("rows (meta-programs):", len(mps))
    print("columns (cancer types):", len(CANCERS))
    print("program-MP matches (matched_0.1 rows):", len(m))
    print("dots (cancer x MP pairs):", len(dots))
    print("dots with 2 programs:\n", dots[dots.n_programs > 1].to_string(index=False))
    print("Jaccard range: %.3f-%.3f" % (dots.jaccard.min(), dots.jaccard.max()))
    print("row order:", [" ".join(x.split()) for x in mps])


if __name__ == "__main__":
    main()
