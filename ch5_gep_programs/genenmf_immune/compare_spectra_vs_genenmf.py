#!/usr/bin/env python
"""
compare_spectra_vs_genenmf.py

Ch.5 (GEP programs) — immune-compartment cross-method comparison:
Spectra programs vs GeneNMF meta-programs (MPs) on the SAME immune cells,
per lineage, via top-gene Jaccard overlap. Combines ALL lineages into single
output files (NOT per-lineage subdirectories).

Inputs (per lineage):
  spectra_out/<lin>_spectra_gene_factors.csv   programs x genes (weights)
  genenmf_immune_output/<lin>/<lin>_genenmf_MP_top50_genes.csv  long: meta_program,rank,gene

Outputs (top level of --genenmf-dir, combined across lineages):
  Spectra_vs_GeneNMF_jaccard_full_all_lineages.csv    lineage,meta_program,spectra_program,jaccard,n_shared
  Spectra_vs_GeneNMF_jaccard_best_match_all_lineages.csv  best Spectra match per MP

Usage:
  python compare_spectra_vs_genenmf.py \
    --spectra-dir /…/clustered_lineages/spectra_out \
    --genenmf-dir /…/genenmf_immune_output \
    --top-n 50
"""
import argparse
import os
import glob
import pandas as pd


def log(m):
    import datetime
    print(f"[{datetime.datetime.now()}] {m}", flush=True)


def top_genes_per_spectra(gene_factors_csv, top_n):
    """Spectra gene_factors: rows=programs, cols=genes, values=weights."""
    df = pd.read_csv(gene_factors_csv, index_col=0)
    out = {}
    for prog, row in df.iterrows():
        out[str(prog)] = set(row.sort_values(ascending=False).head(top_n).index)
    return out


def genes_per_mp(mp_csv, top_n):
    """GeneNMF MP_top50_genes: long format meta_program,rank,gene."""
    df = pd.read_csv(mp_csv)
    out = {}
    for mp, sub in df.groupby("meta_program"):
        genes = sub.sort_values("rank")["gene"].head(top_n).tolist()
        out[str(mp)] = set(genes)
    return out


def jaccard(a, b):
    u = len(a | b)
    return (len(a & b) / u) if u else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectra-dir", required=True)
    ap.add_argument("--genenmf-dir", required=True)
    ap.add_argument("--lineages", nargs="*", default=[
        "B-cell", "NK-cells", "Hematopoietic_Stem_and_Progenitor_cells",
        "T-cells", "Myeloid"])
    ap.add_argument("--top-n", type=int, default=50)
    args = ap.parse_args()

    full_rows = []
    best_rows = []
    included, skipped = [], []

    for lin in args.lineages:
        sp_csv = os.path.join(args.spectra_dir, f"{lin}_spectra_gene_factors.csv")
        mp_csv = os.path.join(args.genenmf_dir, lin, f"{lin}_genenmf_MP_top50_genes.csv")
        if not (os.path.exists(sp_csv) and os.path.exists(mp_csv)):
            skipped.append(lin)
            continue
        spectra = top_genes_per_spectra(sp_csv, args.top_n)
        mps = genes_per_mp(mp_csv, args.top_n)
        log(f"{lin}: {len(mps)} GeneNMF MPs x {len(spectra)} Spectra programs")

        for mp, mp_genes in mps.items():
            best_prog, best_j, best_shared = None, -1.0, 0
            for prog, prog_genes in spectra.items():
                j = jaccard(mp_genes, prog_genes)
                shared = len(mp_genes & prog_genes)
                full_rows.append({"lineage": lin, "meta_program": mp,
                                  "spectra_program": prog,
                                  "jaccard": round(j, 4), "n_shared": shared})
                if j > best_j:
                    best_prog, best_j, best_shared = prog, j, shared
            best_rows.append({"lineage": lin, "meta_program": mp,
                              "best_spectra_program": best_prog,
                              "best_jaccard": round(best_j, 4),
                              "n_shared_genes": best_shared})
        included.append(lin)

    full = pd.DataFrame(full_rows)
    best = pd.DataFrame(best_rows)
    full_path = os.path.join(args.genenmf_dir,
                             "Spectra_vs_GeneNMF_jaccard_full_all_lineages.csv")
    best_path = os.path.join(args.genenmf_dir,
                             "Spectra_vs_GeneNMF_jaccard_best_match_all_lineages.csv")
    full.to_csv(full_path, index=False)
    best.to_csv(best_path, index=False)

    log(f"Included lineages: {included}")
    if skipped:
        log(f"Skipped (missing Spectra or GeneNMF output): {skipped}")
    log(f"Wrote {full_path}  ({len(full)} pairs)")
    log(f"Wrote {best_path}  ({len(best)} MPs)")
    if len(best):
        log("Best-match Jaccard summary (median per lineage):")
        print(best.groupby("lineage")["best_jaccard"].median().to_string(), flush=True)


if __name__ == "__main__":
    main()
