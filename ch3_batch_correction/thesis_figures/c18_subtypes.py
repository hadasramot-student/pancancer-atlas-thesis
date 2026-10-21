"""Comment 18: batch-correction UMAPs coloured by Scimilarity subtype (3 + 2 panels, legend in slot 6).
Matches the cluster's appendix_a_fig10/12/14/54 layout. Usage: python c18_subtypes.py <umapdata.csv.gz> <out_base> [n]"""
import sys, re, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # c18_style.py sits next to this script
import c18_style  # noqa: rcParams
from c18_style import save, plt
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]

PAL = ["#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "#984ea3", "#a65628",
       "#f781bf", "#66c2a5", "#fc8d62", "#8da0cb"]
OTHER = "#bbbbbb"
RENAME = {
    "CD14-positive, CD16-positive monocyte": "CD14+ CD16+ monocyte",
    "CD8-positive, alpha-beta cytotoxic T cell": "CD8+ cytotoxic T cell",
    "CD8-positive, alpha-beta T cell": "CD8+ T cell",
    "CD16-negative, CD56-bright natural killer cell, human": "CD16- CD56bright NK cell",
    "CD16-positive, CD56-dim natural killer cell, human": "CD16+ CD56dim NK cell",
    "mature NK T cell": "Mature NKT cell",
    "natural killer cell": "NK cell",
    "T cell": "T cell (unspecified)",
    "B cell": "B cell (unspecified)",
}
def nice(s):
    s = RENAME.get(s, s)
    s = s[0].upper() + s[1:]
    import textwrap
    return "\n".join(textwrap.wrap(s, 25))

src, out = sys.argv[1], sys.argv[2]
d = pd.read_csv(src)
top = d["scimilarity_celltype"].value_counts()
top = top[top.index != "other"]
keep = sorted(top.index[:10], key=lambda k: nice(k).lower())                     # top 10 by count, legend alphabetical (as cluster)
lab = d["scimilarity_celltype"].where(d["scimilarity_celltype"].isin(keep), "Other")
cmap = {k: PAL[i] for i, k in enumerate(keep)}; cmap["Other"] = OTHER
rng = np.random.default_rng(0); order = rng.permutation(len(d))
n = len(d)
s = 0.5 if n > 40000 else 0.9 if n > 10000 else 2.5 if n > 3000 else 6
s_other = s

W, H = 5.4, 3.65
fig = plt.figure(figsize=(W, H))
pw, ph = 0.315, 0.429             # panel size (fig fraction), from the cluster figures
xs = [0.0093, 0.3429, 0.6757]
y1, y2 = 1 - 0.0603 - ph, 1 - 0.5577 - ph
pos = [(xs[0], y1), (xs[1], y1), (xs[2], y1), (xs[0], y2), (xs[1], y2)]
for (x, y), m in zip(pos, ["Uncorrected", "MNN", "ComBat", "Harmony", "scANVI"]):
    ax = fig.add_axes([x, y, pw, ph])
    xy = d[[f"{m}_1", f"{m}_2"]].to_numpy()[order]; L = lab.to_numpy()[order]
    oth = L == "Other"
    ax.scatter(xy[oth, 0], xy[oth, 1], s=s_other, c=OTHER, lw=0, alpha=0.6, rasterized=True)
    ax.scatter(xy[~oth, 0], xy[~oth, 1], s=s, c=[cmap[v] for v in L[~oth]], lw=0, alpha=0.7, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(m, pad=2)
    for sp in ax.spines.values(): sp.set_color("#555555"); sp.set_linewidth(0.6)
handles = [plt.Line2D([], [], ls="", marker="o", ms=5, mec="none", mfc=cmap[k], label=nice(k)) for k in keep]
handles.append(plt.Line2D([], [], ls="", marker="o", ms=5, mec="none", mfc=OTHER, label="Other"))
wrapped = any("\n" in h.get_label() for h in handles)
if wrapped:   # extra lines: anchor at the bottom of the slot and tighten spacing so "Other" stays on the page
    fig.legend(handles=handles, title="Cell type", loc="lower left", bbox_to_anchor=(xs[2] - 0.004, 0.004),
               frameon=False, alignment="left", handletextpad=0.1, labelspacing=0.12, borderaxespad=0)
else:
    fig.legend(handles=handles, title="Cell type", loc="upper left", bbox_to_anchor=(xs[2] - 0.004, y2 + ph + 0.035),
               frameon=False, alignment="left", handletextpad=0.1, labelspacing=0.25, borderaxespad=0)
save(fig, out)
print("n cells", n, "| shown:", [nice(k).replace(chr(10)," ") for k in keep])
