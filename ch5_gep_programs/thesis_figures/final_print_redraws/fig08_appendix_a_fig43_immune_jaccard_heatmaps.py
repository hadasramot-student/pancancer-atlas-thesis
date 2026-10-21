"""appendix_a_fig43 (thesis Figure A.33): Spectra vs GeneNMF Jaccard per immune lineage, print version.

Same matrices as the original 5-panel figure: rows = the 19 GeneNMF meta-programs of the lineage,
columns = the Spectra programs that are the best match of at least one meta-program in that
lineage (S<n> = Spectra program n), colour = top-50-gene Jaccard (0-0.4, shared colour bar).
Red row labels = HC meta-programs (best-match Jaccard > 0.1). Layout 2 + 2 + 1 panels.
"""
import matplotlib.pyplot as plt
from printstyle import save, WIDTH
from immune_common import best, full, LINEAGES, LABEL, MPS, HC_RED

fig = plt.figure(figsize=(WIDTH, 8.4))
W, Hh = 2.2, 2.05                         # heatmap box (in)
cols_x = [0.42, 3.12]
rows_y = [5.95, 3.15, 0.35]
pos = [(cols_x[0], rows_y[0]), (cols_x[1], rows_y[0]), (cols_x[0], rows_y[1]),
       (cols_x[1], rows_y[1]), (cols_x[0], rows_y[2])]
for lin, (x0, y0) in zip(LINEAGES, pos):
    sp = best[best.lineage == lin].best_spectra_program.unique()
    sp = sorted(sp, key=lambda s: int(s.split("_")[1]))
    m = (full[full.lineage == lin].pivot(index="meta_program", columns="spectra_program", values="jaccard")
         .reindex(index=MPS, columns=sp))
    hc = set(best[(best.lineage == lin) & best.hc].meta_program)
    ax = fig.add_axes([x0 / WIDTH, y0 / 8.4, W / WIDTH, Hh / 8.4])
    im = ax.imshow(m.values, cmap="Blues", vmin=0, vmax=0.4, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(sp))); ax.set_xticklabels([s.replace("Spectra_", "S") for s in sp], rotation=90)
    ax.set_yticks(range(len(MPS))); ax.set_yticklabels(MPS)
    for t, mp in zip(ax.get_yticklabels(), MPS):
        if mp in hc:
            t.set_color(HC_RED); t.set_fontweight("bold")
    ax.tick_params(length=2, pad=1)
    ax.set_title(LABEL[lin], fontsize=9, fontweight="bold", pad=3)
    print(lin, "MPs", m.shape[0], "Spectra cols", m.shape[1], "HC", len(hc))

cax = fig.add_axes([(cols_x[1] + 0.2) / WIDTH, (rows_y[2] + 0.9) / 8.4, 1.8 / WIDTH, 0.12 / 8.4])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("Jaccard similarity (top-50 genes)", fontsize=8)
cb.set_ticks([0, 0.1, 0.2, 0.3, 0.4])
save(fig, "appendix_a_fig43")
