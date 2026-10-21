#!/usr/bin/env python3
"""
05_figures.py -- build the Figure 1 / Figure 2 master table and render both figures.

Inputs  (--root):
  drawn_cells.csv            row_uid, cell_id, project_name, cancer_type,
                             Consensus_Cell_Type, maj3_cell, maj2_cell,
                             maj2_cluster, source_repo, agreement
  fig12_200k_coords.csv      row_uid, PCA_UMAP1/2, scANVI_UMAP1/2, PC1..30, scANVI_1..30
  fig12_200k_scanvi_pred.csv row_uid, scanvi_predicted_label  (optional)

Outputs (--root/figures):
  fig12_master_table.csv                  restyle-without-rerun table
  Figure1_atlas_umap.{png,pdf,svg}        2x2 scANVI-UMAP, 4 colourings
  Figure2_batch_correction.{png,pdf,svg}  raw PCA UMAP vs scANVI UMAP, by repo
"""
import os, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RNG = np.random.default_rng(42)

META_ORDER = ["Epithelial", "T-cells", "Myeloid", "B-cell", "Fibroblasts",
              "Endothelial", "NK-cells", "Hematopoietic_Stem_and_Progenitor_cells", "other"]
META_COL = dict(zip(META_ORDER, ["#3b7dd8", "#e4572e", "#f2b134", "#5bb867", "#9b5de5",
                                 "#00bbf9", "#ff85a1", "#8d99ae", "#c9cad1"]))
DISPLAY = {"Hematopoietic_Stem_and_Progenitor_cells": "HSPC"}
REPO_COL = {"3CA": "#e4572e", "ScBaseCount": "#3b7dd8", "Breastscape": "#5bb867"}
AGREE_ORDER = ["Full", "Partial", "Other/Unknown"]
AGREE_COL = {"Full": "#2a9d8f", "Partial": "#e9c46a", "Other/Unknown": "#c9cad1"}


def cat_palette(levels):
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return {lv: base[i % len(base)] for i, lv in enumerate(levels)}


def scatter_panel(ax, xy, cats, order, colmap, title, size=1.4):
    idx = RNG.permutation(len(cats))
    xy, cats = xy[idx], np.asarray(cats)[idx]
    for lv in order:
        m = cats == lv
        if not m.any():
            continue
        ax.scatter(xy[m, 0], xy[m, 1], s=size, c=[colmap[lv]], linewidths=0,
                   rasterized=True, label=f"{lv} ({m.sum():,})")
    ax.set_title(title, fontsize=12, pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)


def legend_outside(ax, order, colmap, counts, title, ncol=1):
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                      markerfacecolor=colmap[lv], markeredgewidth=0,
                      label=f"{DISPLAY.get(lv, lv)}  ({counts.get(lv, 0):,})")
               for lv in order if counts.get(lv, 0) > 0]
    ax.legend(handles=handles, title=title, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, fontsize=8, title_fontsize=9, ncol=ncol, handletextpad=0.4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi")
    args = ap.parse_args()
    R = args.root
    figdir = os.path.join(R, "figures")
    os.makedirs(figdir, exist_ok=True)

    drawn = pd.read_csv(os.path.join(R, "drawn_cells.csv"))
    coords = pd.read_csv(os.path.join(R, "fig12_200k_coords.csv"))
    keep = ["row_uid", "PCA_UMAP1", "PCA_UMAP2", "scANVI_UMAP1", "scANVI_UMAP2"]
    df = drawn.merge(coords[keep], on="row_uid", how="inner")
    pred_f = os.path.join(R, "fig12_200k_scanvi_pred.csv")
    if os.path.exists(pred_f):
        df = df.merge(pd.read_csv(pred_f), on="row_uid", how="left")
    print(f"master table: {len(df):,} cells (drawn {len(drawn):,}, coords {len(coords):,})")

    mt = df[["row_uid", "cell_id", "project_name", "cancer_type", "Consensus_Cell_Type",
             "agreement", "source_repo", "PCA_UMAP1", "PCA_UMAP2",
             "scANVI_UMAP1", "scANVI_UMAP2"]].copy()
    mt.columns = ["row_uid", "cell_barcode", "project_name", "cancer_type",
                  "Consensus_Cell_Type", "agreement_level", "source_repository",
                  "rawPCA_UMAP1", "rawPCA_UMAP2", "scANVI_UMAP1", "scANVI_UMAP2"]
    mt.to_csv(os.path.join(figdir, "fig12_master_table.csv"), index=False)
    print("wrote fig12_master_table.csv")

    S = df["scANVI_UMAP1"].values, df["scANVI_UMAP2"].values
    Sxy = np.column_stack(S)
    Pxy = df[["PCA_UMAP1", "PCA_UMAP2"]].values

    cancers = sorted(df["cancer_type"].unique(),
                     key=lambda c: -(df["cancer_type"] == c).sum())
    CAN_COL = cat_palette(cancers)

    # ---------------- Figure 1 : 2x2 scANVI UMAP -----------------------------
    fig, axs = plt.subplots(2, 2, figsize=(13.5, 12))
    panels = [
        ("Consensus_Cell_Type", META_ORDER, META_COL, "a  Consensus meta-type", 1),
        ("cancer_type", cancers, CAN_COL, "b  Cancer type", 1),
        ("source_repo", ["3CA", "ScBaseCount", "Breastscape"], REPO_COL, "c  Source repository", 1),
        ("agreement", AGREE_ORDER, AGREE_COL, "d  Consensus agreement level", 1),
    ]
    for ax, (col, order, cmap, title, nc) in zip(axs.flat, panels):
        cnt = df[col].value_counts().to_dict()
        scatter_panel(ax, Sxy, df[col].values, order, cmap, title)
        legend_outside(ax, order, cmap, cnt, col, ncol=nc)
    fig.suptitle("Figure 1  |  Pan-cancer atlas scANVI embedding  (200,000-cell subsample, seed 42)",
                 fontsize=14, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(figdir, f"Figure1_atlas_umap.{ext}"), dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print("wrote Figure1_atlas_umap.{png,pdf,svg}")

    # ---------------- Figure 2 : before / after, by repo -------------------
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.2))
    order = ["3CA", "ScBaseCount", "Breastscape"]
    cnt = df["source_repo"].value_counts().to_dict()
    scatter_panel(axs[0], Pxy, df["source_repo"].values, order, REPO_COL,
                  "a  Before  —  raw PCA UMAP", size=1.4)
    scatter_panel(axs[1], Sxy, df["source_repo"].values, order, REPO_COL,
                  "b  After  —  scANVI UMAP", size=1.4)
    legend_outside(axs[1], order, REPO_COL, cnt, "source repository")
    fig.suptitle("Figure 2  |  Batch correction: raw PCA vs scANVI  (same 200,000 cells, coloured by source repository)",
                 fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(figdir, f"Figure2_batch_correction.{ext}"), dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print("wrote Figure2_batch_correction.{png,pdf,svg}")


if __name__ == "__main__":
    main()
