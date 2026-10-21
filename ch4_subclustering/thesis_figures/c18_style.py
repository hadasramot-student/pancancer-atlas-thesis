"""Shared print style for the comment-18 redraws (5.4 in wide, text >= 8 pt)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG_W = 5.4          # final printed width (in)
MAX_H = 8.5          # one full page (in)
DPI = 400

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
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
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})


def save(fig, out_base):
    """PNG at 400 dpi + PDF (scatter layers are rasterized at 400 dpi)."""
    fig.savefig(out_base + ".png", dpi=DPI)
    fig.savefig(out_base + ".pdf", dpi=DPI)
    w, h = fig.get_size_inches()
    print(f"saved {out_base}.png/.pdf  ({w:.2f} x {h:.2f} in)")
