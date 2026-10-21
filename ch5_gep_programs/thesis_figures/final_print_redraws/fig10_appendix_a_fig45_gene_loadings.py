"""appendix_a_fig45a / appendix_a_fig45b (thesis Figure A.35): gene loadings of the top 8 HC
cNMF-GeneNMF program pairs, split over two pages (4 pairs each) so gene names are 8 pt.

Same 8 pairs and same values as the original figure (top 20 genes per program):
- pairs: hc_pairs_0.0538.csv (thesis_figures/comment12_2026-10-21), top 8 by Jaccard;
- cNMF bars: gene_spectra_score (k = 20, dt 0.1) of that program, divided by its maximum;
- GeneNMF bars: rank weight 1 - (rank - 1)/50 from <cancer>_genenmf_MP_top50_genes.csv.
"""
import pandas as pd
import matplotlib.pyplot as plt
from printstyle import save, WIDTH, plain

PAIRS = "/groups/ofircohen-group/users/hadasra/thesis_figures/comment12_2026-10-21/hc_pairs_0.0538.csv"
CNMF = "/groups/ofircohen-group/users/hadasra/temp_objects/cnmf_output/{c}/cnmf_{c}/cnmf_{c}.gene_spectra_score.k_20.dt_0_1.txt"
GNMF = "/groups/ofircohen-group/users/hadasra/temp_objects/genenmf_output/{c}/{c}_genenmf_MP_top50_genes.csv"
N_GENES = 20
C_CNMF, C_MP = "#4a98b9", "#e8606e"

pairs = pd.read_csv(PAIRS).sort_values("jaccard", ascending=False, kind="stable").head(8)
# keep the original panel order for the Kidney tie at J = 0.2821 (cNMF_6/MP13 before cNMF_7/MP15)
pairs["o"] = [float(i) for i in range(len(pairs))]
pairs.loc[(pairs.cancer == "Kidney") & (pairs.genenmf_mp == "MP13"), "o"] = 1.5
pairs = pairs.sort_values("o").reset_index(drop=True)
print(pairs[["cancer", "cnmf_program", "genenmf_mp", "jaccard"]].to_string())

spectra = {}


def cnmf_top(cancer, prog):
    if cancer not in spectra:
        spectra[cancer] = pd.read_csv(CNMF.format(c=cancer), sep="\t", index_col=0)
    s = spectra[cancer].loc[int(prog.split("_")[1])].sort_values(ascending=False).head(N_GENES)
    return s / s.max()


def mp_top(cancer, mp):
    d = pd.read_csv(GNMF.format(c=cancer))
    d = d[d.meta_program == mp].sort_values("rank").head(N_GENES)
    return pd.Series(1 - (d["rank"].values - 1) / 50, index=d.gene.values)


def bars(ax, s, color, xlabel, title):
    ax.barh(range(len(s)), s.values, color=color, height=0.75)
    ax.set_yticks(range(len(s))); ax.set_yticklabels(s.index, fontsize=8)
    ax.set_ylim(len(s) - 0.5, -0.6); ax.set_xlim(0, 1.05)
    ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.set_xlabel(xlabel, fontsize=8, labelpad=1)
    ax.set_title(title, fontsize=8, fontweight="bold", color="black", pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


H = 8.3
for page, name in enumerate(["appendix_a_fig45a", "appendix_a_fig45b"]):
    fig = plt.figure(figsize=(WIDTH, H))
    for k in range(4):
        r = pairs.iloc[page * 4 + k]
        row, col = divmod(k, 2)
        x0 = 0.0 + col * 2.72                # pair block left edge (in)
        y_top = H - 0.08 - row * 4.12         # pair block top (in)
        letter = "abcdefgh"[page * 4 + k]
        num = r.cnmf_program.split("_")[1]
        fig.text((x0 + 0.03) / WIDTH, y_top / H, f"{letter}  {plain(r.cancer)}: cNMF {num} vs {r.genenmf_mp} (J = {r.jaccard:.2f})",
                 fontsize=9, fontweight="bold", va="top")
        pw, ph = 0.6, 3.15                    # bar-chart box (in)
        for j, (s, colr, xl, t) in enumerate([
                (cnmf_top(r.cancer, r.cnmf_program), C_CNMF, "Loading", f"cNMF {num}"),
                (mp_top(r.cancer, r.genenmf_mp), C_MP, "Rank weight", r.genenmf_mp)]):
            ax = fig.add_axes([(x0 + 0.66 + j * 1.34) / WIDTH, (y_top - 0.5 - ph) / H, pw / WIDTH, ph / H])
            bars(ax, s, colr, xl, t)
    save(fig, name)
