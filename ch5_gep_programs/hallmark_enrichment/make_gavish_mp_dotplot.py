#!/usr/bin/env python
# make_gavish_mp_dotplot.py
#
# Ch.5 -- cNMF programs (epithelial/malignant compartment) vs. Gavish et al.
# 2023 (Tirosh lab) 41 pan-cancer meta-programs, Jaccard > 0.1.
#
# y = the 17 Gavish MPs (of 41) with >=1 matched_0.1==True hit anywhere in
# the atlas; x = all 12 cancer types. One dot per (cancer, MP) pair with at
# least one matching cNMF program. Dot size = number of cNMF programs
# matching that MP in that cancer type (1, except Prostate/MP12 EMT-I = 2).
# Dot color = Jaccard score (continuous).
#
# y-axis labels color-coded by narrative category (Cell Cycle / EMT /
# MHC-II-IFNg vs. neutral gray for everything else), grouped at the top of
# the axis rather than alphabetical.
#
# Data is unchanged from combined_all12_gavish_jaccard.csv (no Jaccard
# recomputation) -- this script only filters to matched_0.1 and plots.
#
# Input:
#   temp_objects/gavish_2023/combined_all12_gavish_jaccard.csv
# Output:
#   plots_and_outputs/Fig_Gavish_MP_dotplot_Ch5.png

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = "/mnt/new_groups/ofircohen-group/users/hadasra"
IN_CSV = f"{BASE}/temp_objects/gavish_2023/combined_all12_gavish_jaccard.csv"
OUT_PNG = f"{BASE}/create_sc_object_workflow/plots_and_outputs/Fig_Gavish_MP_dotplot_Ch5.png"

CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney",
           "Liver_Biliary", "Lung", "Ovarian", "Pancreas", "Prostate",
           "Sarcoma", "Skin"]

# narrative categories -> Okabe-Ito colorblind-safe palette
CATEGORY_COLOR = {
    "Cell Cycle": "#D55E00",       # vermillion
    "EMT": "#0072B2",              # blue
    "MHC-II/IFN-gamma": "#009E73", # green
}
NEUTRAL = "0.35"

MP_CATEGORY = {
    "MP1  Cell Cycle - G2/M": "Cell Cycle",
    "MP3  Cell Cylce HMG-rich": "Cell Cycle",
    "MP12 EMT-I": "EMT",
    "MP13 EMT-II": "EMT",
    "MP14 EMT-III": "EMT",
    "MP15 EMT IV": "EMT",
    "MP17 Interferon/MHC-II (I)": "MHC-II/IFN-gamma",
    "MP18 Interferon/MHC-II (II)": "MHC-II/IFN-gamma",
}
CATEGORY_ORDER = ["Cell Cycle", "EMT", "MHC-II/IFN-gamma"]


def mp_sort_key(mp):
    # "MP17 Interferon/MHC-II (I)" -> 17 ; "MP1  Cell Cycle - G2/M" -> 1
    return int(mp.strip().split()[0].replace("MP", ""))


def main():
    df = pd.read_csv(IN_CSV)
    m = df[df["matched_0.1"]].copy()
    assert len(m) == 49, f"expected 49 matched_0.1 rows, got {len(m)}"

    # one dot per (cancer, MP): size = program count, color = max jaccard for that pair
    dots = (m.groupby(["cancer", "best_gavish_mp"])
              .agg(n_programs=("program", "count"), jaccard=("jaccard", "max"))
              .reset_index())

    mps_present = sorted(dots["best_gavish_mp"].unique(), key=mp_sort_key)
    assert len(mps_present) == 17, f"expected 17 MPs, got {len(mps_present)}"

    # order y-axis: the 3 narrative categories first (grouped, numeric within
    # group), then everything else (neutral), numeric order
    def group_rank(mp):
        cat = MP_CATEGORY.get(mp)
        return (CATEGORY_ORDER.index(cat) if cat in CATEGORY_ORDER else len(CATEGORY_ORDER),
                mp_sort_key(mp))

    mps_ordered = sorted(mps_present, key=group_rank)
    y_pos = {mp: i for i, mp in enumerate(mps_ordered)}
    x_pos = {c: i for i, c in enumerate(CANCERS)}

    dots["x"] = dots["cancer"].map(x_pos)
    dots["y"] = dots["best_gavish_mp"].map(y_pos)

    cmap = plt.get_cmap("viridis")
    norm = matplotlib.colors.Normalize(vmin=dots["jaccard"].min(), vmax=dots["jaccard"].max())

    SIZE_1, SIZE_2 = 380, 760

    fig, ax = plt.subplots(figsize=(13, 10), dpi=300)

    sizes = dots["n_programs"].map({1: SIZE_1, 2: SIZE_2}).values
    sc = ax.scatter(dots["x"], dots["y"], c=dots["jaccard"], cmap=cmap, norm=norm,
                     s=sizes, edgecolor="0.2", linewidth=1.0, zorder=3)

    ax.set_xlim(-0.6, len(CANCERS) - 0.4)
    ax.set_ylim(-0.7, len(mps_ordered) - 0.3)
    ax.set_xticks(range(len(CANCERS)))
    ax.set_xticklabels([c.replace("_", " ") for c in CANCERS], rotation=40, ha="right", fontsize=15)
    ax.set_yticks(range(len(mps_ordered)))
    ax.set_yticklabels(mps_ordered, fontsize=15)

    # color-code y tick labels by category
    for tick, mp in zip(ax.get_yticklabels(), mps_ordered):
        cat = MP_CATEGORY.get(mp)
        tick.set_color(CATEGORY_COLOR.get(cat, NEUTRAL))
        if cat is not None:
            tick.set_fontweight("bold")

    ax.invert_yaxis()  # categories at top
    ax.grid(True, color="0.9", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_edgecolor("0.4")

    ax.set_title("cNMF programs (Ch.5, epithelial/malignant compartment) vs.\n"
                 "Gavish et al. 2023 (Tirosh lab) 41 meta-programs, Jaccard > 0.1",
                 fontsize=19, pad=16)

    # colorbar (continuous jaccard)
    cbar = fig.colorbar(sc, ax=ax, orientation="vertical", fraction=0.035, pad=0.03)
    cbar.set_label("Jaccard similarity", fontsize=16)
    cbar.ax.tick_params(labelsize=13)

    # size legend (program count)
    size_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="0.6",
               markeredgecolor="0.2", markersize=(SIZE_1 ** 0.5) * 0.85, label="1 program"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="0.6",
               markeredgecolor="0.2", markersize=(SIZE_2 ** 0.5) * 0.85, label="2 programs"),
    ]
    leg = ax.legend(handles=size_handles, title="cNMF programs\nmatching MP",
                     loc="upper left", bbox_to_anchor=(1.16, 1.0),
                     fontsize=14, title_fontsize=15, frameon=True, labelspacing=1.4,
                     borderpad=1.1)
    ax.add_artist(leg)

    # category legend
    cat_handles = [Line2D([0], [0], marker="s", color="none", markerfacecolor=CATEGORY_COLOR[c],
                          markersize=13, label=c) for c in CATEGORY_ORDER]
    cat_handles.append(Line2D([0], [0], marker="s", color="none", markerfacecolor=NEUTRAL,
                              markersize=13, label="Other"))
    ax.legend(handles=cat_handles, title="Narrative category",
              loc="upper left", bbox_to_anchor=(1.16, 0.55),
              fontsize=14, title_fontsize=15, frameon=True, labelspacing=1.0,
              borderpad=1.1)

    fig.subplots_adjust(right=0.78, bottom=0.2, left=0.28, top=0.86)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")
    print(f"Total matched_0.1 rows (input): {len(m)}")
    print(f"Unique (cancer, MP) dot positions: {len(dots)}")
    print(f"MPs on y-axis: {len(mps_ordered)}")
    print(f"Cancer types on x-axis: {len(CANCERS)}")
    print("\nDots with n_programs > 1:")
    print(dots[dots['n_programs'] > 1].to_string(index=False))


if __name__ == "__main__":
    main()
