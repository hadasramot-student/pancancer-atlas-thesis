#!/usr/bin/env python
# 02_fig6_1_hallmark_dotplot.py
#
# Figure 6.1 (ch5_fig1) -- cNMF program Hallmark enrichment by pathway, all 12
# cancer types. Reads hallmark_top_per_program_bh.csv from script 01.
#   - each cNMF program is assigned its top Hallmark pathway (smallest BH q);
#     "enriched" = q < 0.05 (BH within the program's 50 tests) and >= 3 shared genes
#   - cell (pathway x cancer) is drawn if >= 1 enriched program falls in it
#   - filled dot (dark edge) = >= 1 HC program (cNMF-GeneNMF Jaccard > 0.0538)
#     open circle = only single-method programs
#   - dot size = -log10(q) of the best program in the cell, capped at 10
#   - "x n" label = n >= 2 HC programs in the cell
#   - last row "HC, no enriched pathway" = HC programs whose best q >= 0.05
# Output: ch5_fig1.png (400 dpi) / ch5_fig1.pdf, fig6_1_cells.csv

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = "/groups/ofircohen-group/users/hadasra/thesis_figures/comment12_2026-10-21"
CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney", "Liver_Biliary",
           "Lung", "Ovarian", "Pancreas", "Prostate", "Sarcoma", "Skin"]
LABEL = {"Head_and_Neck": "Head and Neck", "Liver_Biliary": "Liver/Biliary"}
CAP = 10.0                      # -log10(q) cap for the size scale (q <= 1e-10)
S_MIN, S_MAX = 10.0, 75.0       # marker area (pt^2) at -log10(0.05) and at the cap
UNANNOT = "HC, no enriched pathway"

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#B8860B", "#8B0000",
           "#4B0082", "#006400", "#E69F00", "#1F4E79", "#800080", "#2E8B57"]

PRETTY = {
    "TNFA_SIGNALING_VIA_NFKB": "TNF-α signaling via NF-κB",
    "EPITHELIAL_MESENCHYMAL_TRANSITION": "Epithelial–mesenchymal trans.",
    "INTERFERON_GAMMA_RESPONSE": "Interferon-γ response",
    "INTERFERON_ALPHA_RESPONSE": "Interferon-α response",
    "IL2_STAT5_SIGNALING": "IL2–STAT5 signaling",
    "IL6_JAK_STAT3_SIGNALING": "IL6–JAK–STAT3 signaling",
    "PI3K_AKT_MTOR_SIGNALING": "PI3K–AKT–mTOR signaling",
    "MTORC1_SIGNALING": "mTORC1 signaling",
    "TGF_BETA_SIGNALING": "TGF-β signaling",
    "WNT_BETA_CATENIN_SIGNALING": "WNT–β-catenin signaling",
    "REACTIVE_OXYGEN_SPECIES_PATHWAY": "Reactive oxygen species",
    "UNFOLDED_PROTEIN_RESPONSE": "Unfolded protein response",
    "KRAS_SIGNALING_UP": "KRAS signaling (up)",
    "KRAS_SIGNALING_DN": "KRAS signaling (down)",
    "UV_RESPONSE_UP": "UV response (up)",
    "UV_RESPONSE_DN": "UV response (down)",
    "ESTROGEN_RESPONSE_EARLY": "Estrogen response (early)",
    "ESTROGEN_RESPONSE_LATE": "Estrogen response (late)",
    "MYC_TARGETS_V1": "MYC targets v1",
    "MYC_TARGETS_V2": "MYC targets v2",
    "E2F_TARGETS": "E2F targets",
    "G2M_CHECKPOINT": "G2M checkpoint",
    "P53_PATHWAY": "p53 pathway",
    "DNA_REPAIR": "DNA repair",
    "HEDGEHOG_SIGNALING": "Hedgehog signaling",
    "NOTCH_SIGNALING": "NOTCH signaling",
}


def pretty(h):
    return PRETTY.get(h, h.replace("_", " ").capitalize())


def size(v):
    v = np.clip(v, -np.log10(0.05), CAP)
    return S_MIN + (v + np.log10(0.05)) / (CAP + np.log10(0.05)) * (S_MAX - S_MIN)


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    b = pd.read_csv(f"{OUT}/hallmark_top_per_program_bh.csv")
    b["HC"] = b["HC_0.0538"].astype(bool)
    e = b[b.top_hallmark_bh != "Unannotated"].copy()
    e["mlq"] = -np.log10(e.qval)
    cells = (e.groupby(["top_hallmark_bh", "cancer"])
               .agg(n_programs=("program", "size"), n_hc=("HC", "sum"), best_q=("qval", "min"),
                    programs=("program", lambda s: ";".join(s)))
               .reset_index())
    cells["mlq"] = -np.log10(cells.best_q)

    # row order: cancers with HC dot (desc), then cancers with any dot (desc), then name
    rk = cells.groupby("top_hallmark_bh").agg(n_hc_c=("n_hc", lambda s: (s > 0).sum()),
                                              n_c=("cancer", "size"))
    rows = rk.sort_values(["n_hc_c", "n_c"], ascending=False).index.tolist()
    rows = sorted(rows, key=lambda r: (-rk.loc[r, "n_hc_c"], -rk.loc[r, "n_c"], pretty(r)))

    un = b[(b.top_hallmark_bh == "Unannotated") & b.HC].groupby("cancer").size()
    show_un = len(un) > 0
    ylabels = [pretty(r) for r in rows] + ([UNANNOT] if show_un else [])
    nrow = len(ylabels)
    color = {r: PALETTE[i % len(PALETTE)] for i, r in enumerate(rows)}
    cx = {c: i for i, c in enumerate(CANCERS)}

    top_in, bot_in, left_in, right_in = 0.10, 2.15, 1.92, 0.08
    fig_h = top_in + 0.150 * nrow + bot_in
    fig = plt.figure(figsize=(5.4, fig_h))
    ax = fig.add_axes([left_in / 5.4, bot_in / fig_h,
                       1 - (left_in + right_in) / 5.4, 1 - (top_in + bot_in) / fig_h])

    for _, r in cells.iterrows():
        y = rows.index(r.top_hallmark_bh)
        x = cx[r.cancer]
        col = color[r.top_hallmark_bh]
        if r.n_hc > 0:
            ax.scatter(x, y, s=size(r.mlq), facecolor=col, edgecolor="black", linewidth=0.6, zorder=3)
        else:
            ax.scatter(x, y, s=size(r.mlq), facecolor="none", edgecolor=col, linewidth=1.0, zorder=3)
        if r.n_hc >= 2:
            ax.text(x + 0.2, y - 0.08, f"×{int(r.n_hc)}", fontsize=8, ha="left", va="center",
                    color="black", zorder=4)
    if show_un:
        y = len(rows)
        for c, n in un.items():
            ax.scatter(cx[c], y, s=size(-np.log10(0.05)) + 12, facecolor="black",
                       edgecolor="black", marker="s", linewidth=0.6, zorder=3)
            if n >= 2:
                ax.text(cx[c] + 0.2, y - 0.08, f"×{n}", fontsize=8, ha="left", va="center", zorder=4)
        ax.axhline(len(rows) - 0.5, color="black", linewidth=0.5)

    ax.set_xlim(-0.6, len(CANCERS) - 0.2)
    ax.set_ylim(nrow - 0.4, -0.6)
    ax.set_xticks(range(len(CANCERS)))
    ax.set_xticklabels([LABEL.get(c, c) for c in CANCERS], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8.5, color="black")
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(ylabels, fontsize=8, color="black")
    ax.tick_params(length=2, width=0.6, colors="black", pad=2)
    ax.set_axisbelow(True)
    ax.grid(color="#D9D9D9", linewidth=0.4)
    for s in ax.spines.values():
        s.set_linewidth(0.6)

    # ---- legend below the plot
    dark = "#1F4E79"
    h1 = [Line2D([], [], marker="o", ls="", markersize=6.5, markerfacecolor=dark,
                 markeredgecolor="black", markeredgewidth=0.6,
                 label="High-confidence (cNMF and GeneNMF, Jaccard > 0.0538)"),
          Line2D([], [], marker="o", ls="", markersize=6.5, markerfacecolor="none",
                 markeredgecolor=dark, markeredgewidth=1.0, label="Single-method only")]
    if show_un:
        h1.append(Line2D([], [], marker="s", ls="", markersize=4.5, markerfacecolor="black",
                         markeredgecolor="black", label="HC program, no enriched pathway (fixed size)"))
    h1.append(Line2D([], [], ls="", label="×n = number of HC programs in the cell"))
    fig.legend(handles=h1, loc="upper left", bbox_to_anchor=(0.02, 1.38 / fig_h),
                      frameon=False, fontsize=8, handletextpad=0.4, borderaxespad=0, labelspacing=0.35)
    hs = [Line2D([], [], marker="o", ls="", markersize=np.sqrt(size(v)), markerfacecolor="black",
                 markeredgecolor="black", label=lab)
          for v, lab in [(-np.log10(0.05), "q = 0.05"), (5, "q = 1e-5"), (10, "q ≤ 1e-10")]]
    fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.02, 0.62 / fig_h),
               frameon=False, fontsize=8, ncol=3, title="−log10 q (Hallmark, BH-adjusted)",
               title_fontsize=8, alignment="left", columnspacing=1.4, handletextpad=0.3, borderaxespad=0)

    fig.savefig(f"{OUT}/ch5_fig1.png", dpi=400)
    fig.savefig(f"{OUT}/ch5_fig1.pdf")

    cells["filled_HC"] = cells.n_hc > 0
    cells.to_csv(f"{OUT}/fig6_1_cells.csv", index=False)
    print(f"rows (pathways): {len(rows)}; Unannotated row shown: {show_un}")
    print(f"cells: {len(cells)}; filled (HC): {(cells.n_hc > 0).sum()}; open: {(cells.n_hc == 0).sum()}; "
          f"x>=2: {(cells.n_hc >= 2).sum()}; capped: {(cells.mlq >= CAP).sum()}")
    print("HC programs without enriched pathway:", un.to_dict())
    print(b[(b.top_hallmark_bh == "Unannotated") & b.HC][["cancer", "program", "pval", "qval", "overlap",
                                                           "top_hallmark_old"]].to_string(index=False))


if __name__ == "__main__":
    main()
