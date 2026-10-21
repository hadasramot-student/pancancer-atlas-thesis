"""Shared print settings for the final thesis figure redraws (2026-10-21).

Every figure is drawn at its final printed size (5.4 in wide, <= 8.5 in tall),
so the font sizes below are the sizes on paper.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/groups/ofircohen-group/users/hadasra/thesis_figures/final_figure_fixes_2026-10-21"
WIDTH = 5.4
MAX_H = 8.5

plt.rcParams.update({
    "font.family": ["Nimbus Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.title_fontsize": 8,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})

# Plain print labels
PLAIN = {
    "Liver_Biliary": "Liver/Biliary", "Liver/\nBiliary": "Liver/Biliary",
    "Head_and_Neck": "Head and Neck", "Other_Models": "Other models",
    "unknown": "Unknown", "Unknown": "Unknown", "other": "Other",
    "T-cells": "T cells", "B-cell": "B cells", "B-cells": "B cells",
    "NK-cells": "NK cells", "HSPC": "HSPCs",
    "Hematopoietic_Stem_and_Progenitor_cells": "HSPCs",
}


def plain(x):
    return PLAIN.get(x, x)


def save(fig, name, **kw):
    """Save <name>.png (400 dpi) and <name>.pdf; the figure size is kept as drawn."""
    w, h = fig.get_size_inches()
    assert abs(w - WIDTH) < 1e-6 and h <= MAX_H + 1e-6, (w, h)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=400, **kw)
    plt.close(fig)
    print(f"saved {name}.png/.pdf  ({w:.2f} x {h:.2f} in)")
