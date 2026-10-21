#!/usr/bin/env python
# make_ch5_fig5.py  (comment 15/16, 2026-10-21)
#
# Redraw of pancancer-atlas-thesis/ch5_gep_programs/hallmark_enrichment/make_gavish_mp_dotplot.py
# at final print size (5.4 in wide). Same data, same filter (matched_0.1), same dot
# encoding (size = number of matching cNMF programs, 1 or 2, area ratio 1:2; color =
# max Jaccard per (cancer, MP), viridis, min-max normalised). Changes: no title, black
# row labels with a colored category square to their left, "Liver/Biliary",
# legends + colorbar below the plot, all text >= 8 pt.
#
# Input : temp_objects/gavish_2023/combined_all12_gavish_jaccard.csv
# Output: thesis_figures/comment15_16_2026-10-21/ch5_fig5.png (400 dpi) + .pdf

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

BASE = "/groups/ofircohen-group/users/hadasra"
IN_CSV = f"{BASE}/temp_objects/gavish_2023/combined_all12_gavish_jaccard.csv"
OUT = f"{BASE}/thesis_figures/comment15_16_2026-10-21/ch5_fig5"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "pdf.fonttype": 42, "axes.linewidth": 0.6})

CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney",
           "Liver_Biliary", "Lung", "Ovarian", "Pancreas", "Prostate",
           "Sarcoma", "Skin"]
CANCER_LABEL = {"Head_and_Neck": "Head and Neck", "Liver_Biliary": "Liver/Biliary"}

CATEGORY_COLOR = {"Cell Cycle": "#D55E00", "EMT": "#0072B2", "MHC-II/IFN-gamma": "#009E73"}
CATEGORY_LABEL = {"Cell Cycle": "Cell cycle", "EMT": "EMT", "MHC-II/IFN-gamma": "MHC-II/IFN-γ"}
MP_CATEGORY = {
    "MP1  Cell Cycle - G2/M": "Cell Cycle", "MP3  Cell Cylce HMG-rich": "Cell Cycle",
    "MP12 EMT-I": "EMT", "MP13 EMT-II": "EMT", "MP14 EMT-III": "EMT", "MP15 EMT IV": "EMT",
    "MP17 Interferon/MHC-II (I)": "MHC-II/IFN-gamma",
    "MP18 Interferon/MHC-II (II)": "MHC-II/IFN-gamma",
}
CATEGORY_ORDER = ["Cell Cycle", "EMT", "MHC-II/IFN-gamma"]


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
    dots["x"] = dots["cancer"].map({c: i for i, c in enumerate(CANCERS)})
    dots["y"] = dots["best_gavish_mp"].map(y_pos)

    norm = matplotlib.colors.Normalize(vmin=dots["jaccard"].min(), vmax=dots["jaccard"].max())
    SIZE_1, SIZE_2 = 55, 110  # same 1:2 area ratio as the original (380:760)

    W, H = 5.4, 5.2
    fig = plt.figure(figsize=(W, H))
    # main panel (inches -> figure fractions)
    L, R, B, T = 2.05, 0.15, 1.72, 0.08
    ax = fig.add_axes([L / W, B / H, 1 - (L + R) / W, 1 - (B + T) / H])

    sc = ax.scatter(dots["x"], dots["y"], c=dots["jaccard"], cmap="viridis", norm=norm,
                    s=dots["n_programs"].map({1: SIZE_1, 2: SIZE_2}).values,
                    edgecolor="black", linewidth=0.5, zorder=3)
    ax.set_xlim(-0.6, len(CANCERS) - 0.4)
    ax.set_ylim(len(mps) - 0.4, -0.6)
    ax.set_xticks(range(len(CANCERS)))
    ax.set_xticklabels([CANCER_LABEL.get(c, c) for c in CANCERS], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8, color="black")
    ax.set_yticks(range(len(mps)))
    ax.set_yticklabels([" ".join(mp.split()) for mp in mps], fontsize=8, color="black")
    ax.tick_params(length=2.5, width=0.6, colors="black", labelcolor="black")
    ax.grid(True, color="#BFBFBF", linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_edgecolor("black")

    # category squares to the left of each (black) row label
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = ax.transAxes.inverted()
    for tick, mp in zip(ax.get_yticklabels(), mps):
        cat = MP_CATEGORY.get(mp)
        if cat is None:
            continue
        bb = tick.get_window_extent(rend)
        x0, _ = inv.transform((bb.x0, 0))
        y_disp = ax.transData.transform((0, y_pos[mp]))[1]
        _, yc = inv.transform((0, y_disp))
        sq_w = 0.065 / (W - L - R)            # 0.065 in square
        sq_h = 0.065 / (H - B - T)
        gap = 0.04 / (W - L - R)
        ax.add_patch(Rectangle((x0 - gap - sq_w, yc - sq_h / 2), sq_w, sq_h,
                               transform=ax.transAxes, clip_on=False,
                               facecolor=CATEGORY_COLOR[cat], edgecolor="none"))

    # bottom strip: colorbar | size legend | category legend
    cax = fig.add_axes([0.06, 0.06, 0.30, 0.025])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("Jaccard similarity", fontsize=9, color="black", labelpad=3)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=8, length=2, width=0.6, colors="black")
    cb.set_ticks([0.2, 0.3, 0.4, 0.5, 0.6])
    cb.outline.set_linewidth(0.6)

    size_handles = [
        Line2D([0], [0], marker="o", ls="none", markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=0.5, markersize=SIZE_1 ** 0.5, label="1"),
        Line2D([0], [0], marker="o", ls="none", markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=0.5, markersize=SIZE_2 ** 0.5, label="2"),
    ]
    fig.legend(handles=size_handles, title="Matching cNMF\nprograms", loc="lower left",
               bbox_to_anchor=(0.43, 0.025), ncol=2, fontsize=8, title_fontsize=8,
               frameon=False, handletextpad=0.3, columnspacing=1.0, borderaxespad=0)

    cat_handles = [Line2D([0], [0], marker="s", ls="none", markerfacecolor=CATEGORY_COLOR[c],
                          markeredgecolor="none", markersize=6, label=CATEGORY_LABEL[c])
                   for c in CATEGORY_ORDER]
    fig.legend(handles=cat_handles, title="Narrative category", loc="lower left",
               bbox_to_anchor=(0.68, 0.025), fontsize=8, title_fontsize=8, frameon=False,
               handletextpad=0.3, labelspacing=0.3, borderaxespad=0, alignment="left")

    for ext, kw in (("png", {"dpi": 400}), ("pdf", {})):
        fig.savefig(f"{OUT}.{ext}", **kw)
        print(f"wrote {OUT}.{ext}")
    print(f"rows (Gavish MPs): {len(mps)}")
    print(f"matched program-MP rows (Jaccard > 0.1): {len(m)}")
    print(f"dots (cancer x MP pairs): {len(dots)}")
    print(f"cancer types with >=1 match: {dots.cancer.nunique()}")
    print("pairs with 2 programs:", dots[dots.n_programs > 1][["cancer", "best_gavish_mp"]].values.tolist())
    print(f"Jaccard range of dots: {dots.jaccard.min():.3f}-{dots.jaccard.max():.3f}")
    print("row order:", [" ".join(x.split()) for x in mps])


if __name__ == "__main__":
    main()
