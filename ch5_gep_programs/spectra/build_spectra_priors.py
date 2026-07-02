#!/usr/bin/env python
"""
build_spectra_priors.py

Ch.5 (GEP programs) — Spectra pipeline step D. Assemble the Spectra nested
gene_set_dictionary for one immune lineage:

    {
      "global":   {SET: [genes], ...},          # applied to ALL cells
      "<celltype_A>": {SET: [genes], ...},       # compartment sets, per cell type
      "<celltype_B>": {SET: [genes], ...},
      ...
    }

- "global"  = MSigDB Hallmark (from export_hallmark_gmt.R) + TME 'global' sets.
- per-type   = the lineage's compartment sets (T / Myeloid / B / NK / HSPC),
               attached to every cell-type label present in the data so Spectra
               can learn where each program is active.

All gene sets are intersected with the data's gene vocabulary and pruned to
>= --min-genes, exactly the form Spectra's check_gene_set_dictionary expects.

Usage:
  python build_spectra_priors.py \
    --h5ad        /path/to/clustered_lineages/spectra_out/T-cells_immune.h5ad \
    --cell-type-key scimilarity_celltype \
    --compartment T \
    --hallmark-json priors/hallmark.json \
    --tme-json      tme_gene_sets.json \
    --out-json      priors/T-cells_gene_set_dictionary.json
"""
import argparse
import json
import os
import anndata as ad


def log(msg):
    import datetime
    print(f"[{datetime.datetime.now()}] {msg}", flush=True)


# Immune-relevant MSigDB Hallmark sets (drops developmental/tissue sets like
# spermatogenesis, myogenesis, pancreas-beta, adipogenesis that are noise for
# immune cells). Curating these also shrinks the Spectra vocab a lot (the
# memory-dominant genes^2 graph term), which is needed for the many-cell-type
# lineages (T-cells, Myeloid).
IMMUNE_HALLMARK = {
    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_APOPTOSIS",
    "HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
}


def prune(sets, vocab, min_genes):
    """Intersect each gene set with vocab; drop sets below min_genes."""
    out = {}
    for name, genes in sets.items():
        kept = sorted(set(genes) & vocab)
        if len(kept) >= min_genes:
            out[name] = kept
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cell-type-key", required=True)
    p.add_argument("--compartment", required=True,
                   choices=["T", "Myeloid", "B", "NK", "HSPC"],
                   help="Which TME compartment block applies to this lineage")
    p.add_argument("--hallmark-json", required=True)
    p.add_argument("--tme-json", required=True)
    p.add_argument("--out-json", required=True)
    p.add_argument("--min-genes", type=int, default=3)
    p.add_argument("--no-hallmark", action="store_true",
                   help="Exclude Hallmark; use curated TME sets only")
    p.add_argument("--hallmark-subset", choices=["all", "immune"], default="all",
                   help="'immune' keeps only immune-relevant Hallmark sets "
                        "(smaller vocab -> fits many-cell-type lineages on GPU)")
    args = p.parse_args()

    # Read only structure (var_names + the label column), not the full X.
    log(f"Reading {args.h5ad} (backed)")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    vocab = set(map(str, adata.var_names))
    if args.cell_type_key not in adata.obs.columns:
        raise KeyError(f"cell_type_key '{args.cell_type_key}' not in obs")
    labels = sorted(adata.obs[args.cell_type_key].astype(str).unique().tolist())
    log(f"{adata.n_vars} genes in vocab; {len(labels)} cell types: {labels}")

    with open(args.hallmark_json) as fh:
        hallmark = json.load(fh)
    with open(args.tme_json) as fh:
        tme = json.load(fh)
    tme = {k: v for k, v in tme.items() if not k.startswith("_")}

    # ---- global block: Hallmark + TME global ----
    global_raw = {}
    if not args.no_hallmark:
        hk = {f"HALLMARK_{k}" if not k.startswith("HALLMARK_") else k: v
              for k, v in hallmark.items()}
        if args.hallmark_subset == "immune":
            before = len(hk)
            hk = {k: v for k, v in hk.items() if k in IMMUNE_HALLMARK}
            log(f"Hallmark subset 'immune': {len(hk)}/{before} sets kept")
        global_raw.update(hk)
    global_raw.update(tme.get("global", {}))
    global_sets = prune(global_raw, vocab, args.min_genes)
    log(f"global block: {len(global_sets)}/{len(global_raw)} sets pass "
        f"(>= {args.min_genes} genes in vocab)")

    # ---- per-cell-type block: this lineage's compartment sets ----
    comp_raw = tme.get(args.compartment, {})
    comp_sets = prune(comp_raw, vocab, args.min_genes)
    log(f"compartment '{args.compartment}': {len(comp_sets)}/{len(comp_raw)} sets pass")
    if not comp_sets:
        raise RuntimeError(f"No compartment sets survived pruning for "
                           f"'{args.compartment}'. Check gene symbols vs vocab.")

    gsd = {"global": global_sets}
    for lab in labels:
        # fresh copy per label so Spectra treats them independently
        gsd[lab] = {name: list(genes) for name, genes in comp_sets.items()}

    n_total = len(global_sets) + len(comp_sets) * len(labels)
    log(f"Final dictionary: 'global' ({len(global_sets)} sets) + "
        f"{len(labels)} cell types x {len(comp_sets)} sets = {n_total} factor priors")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(gsd, fh, indent=1)
    log(f"Wrote {args.out_json}")


if __name__ == "__main__":
    main()
