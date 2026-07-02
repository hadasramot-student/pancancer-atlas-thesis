#!/usr/bin/env python
# plot_ids_crossmethod.py
#
# Ch.5 (GEP programs) — cross-method comparison figures (Jaccard + IDS) for
# the EPITHELIAL / TUMOR CELL COMPARTMENT.
#
# IMPORTANT — compartment scope:
#   cNMF and GeneNMF were both run on EPITHELIAL / MALIGNANT cells only, so
#   every figure/output here is the *epithelial (tumor) cell compartment*
#   cross-method comparison: cNMF vs GeneNMF, validated by IDS.
#   The immune-lineage comparison (Spectra vs GeneNMF) is a SEPARATE analysis
#   on a different compartment and is not included here.
#
# Inputs : <OUT_BASE>/all_cancers_ids_pairs.csv
# Outputs (-> <OUT_BASE>/figures/):
#   1. epithelial_median_ids_by_pairtype_heatmap.{png,csv}
#        median IDS per (cancer x pair-type)
#   2. epithelial_crossmethod_pairs_ranked.csv
#        all cNMF-GeneNMF pairs, sorted by IDS desc
#   3. epithelial_top_crossmethod_pairs.png
#        bar chart of the strongest cNMF-GeneNMF couplings

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

OUT_BASE = os.path.join(os.environ.get("ATLAS_DATA_DIR", "data/temp_objects"), "ids_output")
FIG_DIR = os.path.join(OUT_BASE, "figures")
COMPARTMENT = "Epithelial / tumor cell compartment"

# friendly labels + fixed column order for the pair-type axis
PT_LABEL = {"cNMF_cNMF": "cNMF-cNMF",
            "GeneNMF_cNMF": "cNMF-GeneNMF",
            "GeneNMF_GeneNMF": "GeneNMF-GeneNMF"}
PT_ORDER = ["cNMF-cNMF", "cNMF-GeneNMF", "GeneNMF-GeneNMF"]
N_TOP = 30  # bars in the top cross-method figure


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    df = pd.read_csv(os.path.join(OUT_BASE, "all_cancers_ids_pairs.csv"))
    df["pair_label"] = df["pair_type"].map(PT_LABEL)

    # ---- 1. median IDS by pair type x cancer -----------------------------
    med = (df.groupby(["cancer", "pair_label"])["ids"].median()
             .unstack("pair_label")
             .reindex(columns=PT_ORDER))
    med = med.sort_values("cNMF-GeneNMF", ascending=False)
    med.to_csv(os.path.join(FIG_DIR, "epithelial_median_ids_by_pairtype_heatmap.csv"))

    fig, ax = plt.subplots(figsize=(7, 0.5 * len(med) + 2.2))
    sns.heatmap(med, annot=True, fmt=".2f", cmap="viridis", vmin=0,
                cbar_kws={"label": "median IDS-max"}, linewidths=0.5,
                linecolor="white", ax=ax)
    ax.set_title(f"Median inter-program IDS by method pair\n{COMPARTMENT} — cNMF vs GeneNMF",
                 fontsize=12, pad=12)
    ax.set_xlabel("Program pair type"); ax.set_ylabel("Cancer type")
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    fig.text(0.01, 0.01,
             "cNMF & GeneNMF run on epithelial/malignant cells only. "
             "IDS-max (Radhakrishnan et al., PNAS 2025), per-cell program activity.",
             fontsize=7, style="italic", color="0.4")
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(os.path.join(FIG_DIR, "epithelial_median_ids_by_pairtype_heatmap.png"), dpi=200)
    plt.close(fig)

    # ---- 2. ranked cross-method (cNMF-GeneNMF) pairs ---------------------
    cross = (df[df["pair_label"] == "cNMF-GeneNMF"]
             .copy().sort_values("ids", ascending=False))
    # normalize column orientation: cNMF program first, GeneNMF (MP) second
    def orient(r):
        if r["method_a"] == "cNMF":
            return pd.Series([r["program_a"], r["program_b"]])
        return pd.Series([r["program_b"], r["program_a"]])
    cross[["cNMF_program", "GeneNMF_MP"]] = cross.apply(orient, axis=1)
    cross["compartment"] = COMPARTMENT
    ranked = cross[["compartment", "cancer", "cNMF_program", "GeneNMF_MP", "ids"]]
    ranked = ranked.reset_index(drop=True)
    ranked.to_csv(os.path.join(FIG_DIR, "epithelial_crossmethod_pairs_ranked.csv"), index=False)

    # ---- 3. top cross-method pairs bar chart -----------------------------
    top = ranked.head(N_TOP).iloc[::-1]
    labels = [f"{c}: {a}~{b}" for c, a, b in
              zip(top["cancer"], top["cNMF_program"], top["GeneNMF_MP"])]
    fig, ax = plt.subplots(figsize=(8, 0.34 * len(top) + 1.8))
    colors = sns.color_palette("rocket_r", len(top))
    ax.barh(range(len(top)), top["ids"], color=colors)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("IDS-max"); ax.set_xlim(0, 1)
    ax.set_title(f"Top {N_TOP} cross-method program couplings (cNMF ~ GeneNMF)\n{COMPARTMENT}",
                 fontsize=12, pad=10)
    for i, v in enumerate(top["ids"]):
        ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "epithelial_top_crossmethod_pairs.png"), dpi=200)
    plt.close(fig)

    # ---- console summary -------------------------------------------------
    print(f"Compartment: {COMPARTMENT}")
    print(f"Cross-method (cNMF-GeneNMF) pairs: {len(ranked)}")
    print("\nMedian IDS by pair type x cancer:")
    print(med.round(3).to_string())
    print("\nTop 15 cross-method couplings:")
    print(ranked.head(15).to_string(index=False))
    print(f"\nFigures + tables -> {FIG_DIR}")


if __name__ == "__main__":
    main()
