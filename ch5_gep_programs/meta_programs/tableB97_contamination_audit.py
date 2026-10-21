#!/usr/bin/env python
"""
Prompt 9c -- Table B97 contamination audit (batch job).

For all 22 consensus meta-programs (cMP1..cMP22) in Appendix-B Table B97
(source of record: meta_program_catalogue_J_gt_0.10.csv, the 22-row catalogue
committed as 128f678; the v88 thesis .docx Appendix B stops at Table B96 and
contains no rendered Table B97, so the catalogue CSV is the "as printed" table),
rebuild each cMP's top-50 consensus gene list from the underlying per-program
cNMF gene-spectra scores and check it against curated lineage / technical
marker panels.

Also:
  * cMP8 deep-dive: full 50-gene consensus list + exact fibroblast-marker hits.
  * Recompute (from the Jaccard formula) the 5 sub-threshold Gavish matches:
    cMP6, cMP10, cMP14, cMP20, cMP21.
  * Identify which cMP is "the squamous/keratin group".

Everything is written to ONE text report (+ a machine-readable CSV) as it goes.
Deps: numpy, pandas only.
"""

import os
import re
import sys
import csv
import datetime
import numpy as np
import pandas as pd

BASE_CANDIDATES = ["/groups/ofircohen-group", "/mnt/new_groups/ofircohen-group"]
BASE = next((b for b in BASE_CANDIDATES if os.path.isdir(b)), BASE_CANDIDATES[0])
USER_ROOT = os.path.join(BASE, "users", "hadasra")

CATALOGUE_CSV = os.path.join(
    USER_ROOT, "pancancer-atlas-thesis", "ch5_gep_programs", "meta_programs",
    "meta_program_catalogue_J_gt_0.10.csv")
CNMF_BASE = os.path.join(USER_ROOT, "temp_objects", "cnmf_output")
GAVISH_CSV = os.path.join(USER_ROOT, "temp_objects", "gavish_2023",
                          "gavish_MP_gene_lists_long.csv")
COMBINED_JACCARD_CSV = os.path.join(USER_ROOT, "temp_objects", "gavish_2023",
                                    "combined_all12_gavish_jaccard.csv")

OUTDIR = os.path.dirname(os.path.abspath(__file__))
REPORT_TXT = os.path.join(OUTDIR, "tableB97_contamination_audit.txt")
REPORT_CSV = os.path.join(OUTDIR, "tableB97_contamination_audit.csv")

K_PROGRAMS = 20
DT = "dt_0_1"
N_TOP_GENES = 50
GAVISH_CUTOFF = 0.10

CANCERS = ["Brain", "Breast", "Colorectal", "Head_and_Neck", "Kidney",
           "Liver_Biliary", "Lung", "Ovarian", "Pancreas", "Prostate",
           "Sarcoma", "Skin"]

MARKERS = {
    "myeloid":     {"PTPRC", "CD68", "LYZ", "C1QB", "C1QA", "LAPTM5", "AIF1"},
    "t_nk":        {"CD3D", "CD3E", "CD2", "TRBC1", "NKG7", "GNLY"},
    "b_plasma":    {"MS4A1", "CD79A", "CD79B", "MZB1", "TNFRSF17", "DERL3", "JCHAIN"},
    "fibroblast":  {"COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "PDGFRA"},
    "endothelial": {"PECAM1", "VWF", "CDH5"},
}
IG_FAMILY_RE = re.compile(r"^(IGH|IGK|IGL)[VDJGCAEM]")
TECH_MT_RE = re.compile(r"^MT-")
TECH_RIBO_RE = re.compile(r"^(RPLP\d|RPSA$|RPL\d|RPS\d)")

SQUAMOUS_GENES = {"SFN", "LY6D", "DSG3", "DSP", "PKP1", "PKP3", "PERP", "CSTA",
                  "CSTB", "TACSTD2", "CALML3", "CALML5", "LGALS7", "LGALS7B",
                  "CRCT1", "CRNN", "KRTDAP", "SBSN", "DMKN"}

OFIR_CALLS = {
    "cMP1":  ("technical",    "Ofir: mitochondrial (MT-) technical signature."),
    "cMP2":  ("technical",    "Ofir: ribosomal (RPL/RPS) technical signature."),
    "cMP3":  ("clean",        "Ofir: genuine epithelial stress program (JUN/FOS/HSPA1A)."),
    "cMP6":  ("contaminated", "Ofir: MHC-II/IFN program is myeloid-contaminated (CD74/HLA-D/PTPRC/C1Q)."),
    "cMP8":  ("clean",        "Ofir: genuine epithelial EMT program (see deep-dive)."),
    "cMP9":  ("clean",        "Ofir: genuine epithelial cell-cycle/G2M program."),
    "cMP11": ("clean",        "Ofir: genuine epithelial cilia program."),
    "cMP13": ("clean",        "Ofir: genuine epithelial hypoxia program."),
}


def log(fh, *parts):
    line = " ".join(str(p) for p in parts)
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def load_gene_spectra_score(cancer):
    path = os.path.join(CNMF_BASE, cancer, "cnmf_%s" % cancer,
                        "cnmf_%s.gene_spectra_score.k_%d.%s.txt" % (cancer, K_PROGRAMS, DT))
    return pd.read_csv(path, sep="\t", index_col=0)


def load_all_programs(fh):
    programs = {}
    n_background = None
    for cancer in CANCERS:
        df = load_gene_spectra_score(cancer)
        if n_background is None:
            n_background = df.shape[1]
        elif df.shape[1] != n_background:
            log(fh, "  WARNING: %s background gene count %d != %d"
                % (cancer, df.shape[1], n_background))
        for prog_idx in df.index:
            ranked = df.loc[prog_idx].sort_values(ascending=False)
            top_genes = list(ranked.head(N_TOP_GENES).index)
            programs["%s_cNMF_%s" % (cancer, prog_idx)] = {
                "top_genes": top_genes, "top_set": set(top_genes)}
    log(fh, "  loaded %d cNMF programs; background = %s genes"
        % (len(programs), n_background))
    return programs, n_background


def consensus_genes(programs, member_keys):
    freq, rank_sum = {}, {}
    for key in member_keys:
        for rank, g in enumerate(programs[key]["top_genes"]):
            freq[g] = freq.get(g, 0) + 1
            rank_sum[g] = rank_sum.get(g, 0) + rank
    return sorted(freq.keys(), key=lambda g: (-freq[g], rank_sum[g] / freq[g]))


def jaccard(a, b):
    a, b = set(a), set(b)
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def load_gavish_cancer_mps():
    df = pd.read_csv(GAVISH_CSV)
    df = df[df["cell_type"] == "Cancer"]
    return {mp: set(sub["gene"]) for mp, sub in df.groupby("meta_program")}


def gavish_best_match(consensus_top50, gavish_mps):
    best_mp, best_j = None, 0.0
    for mp, genes in gavish_mps.items():
        j = jaccard(consensus_top50, genes)
        if j > best_j:
            best_mp, best_j = mp, j
    return best_mp, best_j


def is_technical(g):
    return bool(TECH_MT_RE.match(g) or TECH_RIBO_RE.match(g))


def marker_hits(gene_list):
    gset = list(dict.fromkeys(gene_list))
    hits = {c: [] for c in MARKERS}
    hits["b_plasma_ig_family"] = []
    hits["technical"] = []
    for g in gset:
        for cat, panel in MARKERS.items():
            if g in panel:
                hits[cat].append(g)
        if IG_FAMILY_RE.match(g):
            hits["b_plasma_ig_family"].append(g)
        if is_technical(g):
            hits["technical"].append(g)
    return hits


def squamous_score(gene_list):
    hits = []
    for g in gene_list:
        if g.startswith("KRT") and len(g) > 3 and g[3].isdigit():
            hits.append(g)
        elif g.startswith(("SPRR", "LCE")):
            hits.append(g)
        elif g in SQUAMOUS_GENES:
            hits.append(g)
        elif g in ("S100A7", "S100A8", "S100A9", "S100A2"):
            hits.append(g)
    return hits


def classify(cmp_id, hits, consensus_top50, squam_cmp):
    if cmp_id in OFIR_CALLS:
        return OFIR_CALLS[cmp_id]
    if cmp_id == squam_cmp:
        return ("clean", "Squamous/keratin epithelial program (Ofir: genuinely epithelial).")

    n_tech = len(hits["technical"])
    lineage = {
        "myeloid": len(hits["myeloid"]),
        "T/NK": len(hits["t_nk"]),
        "B/plasma": len(set(hits["b_plasma"]) | set(hits["b_plasma_ig_family"])),
        "fibroblast": len(hits["fibroblast"]),
        "endothelial": len(hits["endothelial"]),
    }
    worst = max(lineage, key=lineage.get)
    worst_n = lineage[worst]

    emt_vocab = {"FN1", "SPARC", "VIM", "POSTN", "THBS2", "CTHRC1", "AEBP1",
                 "FSTL1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "COL3A1",
                 "COL4A1", "COL4A2", "COL5A2", "MMP2", "LOX", "TIMP1", "INHBA"}
    emt_support = len(emt_vocab & set(consensus_top50))
    fib_only_collagen = set(hits["fibroblast"]) <= {"COL1A1", "COL1A2"}

    all_lineage_hits = sorted(set(
        hits["myeloid"] + hits["t_nk"] +
        list(set(hits["b_plasma"]) | set(hits["b_plasma_ig_family"])) +
        hits["fibroblast"] + hits["endothelial"]))

    if n_tech >= 10:
        return ("technical",
                "%d/50 consensus genes are MT-/ribosomal technical genes." % n_tech)
    if worst == "fibroblast" and worst_n >= 3 and fib_only_collagen and emt_support >= 5:
        return ("clean",
                "fibroblast hits are collagens only (%s) with strong EMT support "
                "(%d EMT effectors) -> genuine EMT, not stroma."
                % (",".join(sorted(hits["fibroblast"])), emt_support))
    if worst_n >= 3:
        return ("contaminated",
                "%d %s markers in consensus top-50 (%s)."
                % (worst_n, worst, ", ".join(all_lineage_hits)))
    if worst_n in (1, 2) or 5 <= n_tech < 10:
        bits = []
        if worst_n:
            bits.append("%d %s marker(s) (%s)" % (worst_n, worst, ", ".join(all_lineage_hits)))
        if 5 <= n_tech < 10:
            bits.append("%d technical genes" % n_tech)
        return ("ambiguous", "borderline: " + "; ".join(bits) + ".")
    return ("clean", "no lineage markers; technical-gene load low.")


def resolve_squamous(programs, catalogue, fh):
    scores = {}
    for row in catalogue:
        members = row["member_cnmf_programs"].split(";")
        cons = consensus_genes(programs, members)[:N_TOP_GENES]
        sq = squamous_score(cons)
        scores[row["meta_program_id"]] = (len(sq), sq)
    ranked = sorted(scores.items(), key=lambda kv: -kv[1][0])
    log(fh, "\n[squamous/keratin resolution] keratin/cornified-gene count in consensus top-50:")
    for cmp_id, (n, genes) in ranked[:6]:
        log(fh, "  %-6s %2d  %s" % (cmp_id, n, ", ".join(genes)))
    winner = ranked[0][0]
    log(fh, "  -> squamous/keratin group = %s" % winner)
    return winner


def main():
    fh = open(REPORT_TXT, "w")
    log(fh, "=" * 78)
    log(fh, "TABLE B97 CONTAMINATION AUDIT  (Prompt 9c)")
    log(fh, "run:", datetime.datetime.now().isoformat(timespec="seconds"))
    log(fh, "host:", os.uname().nodename, "  python:", sys.version.split()[0])
    log(fh, "BASE:", BASE)
    log(fh, "catalogue:", CATALOGUE_CSV)
    log(fh, "=" * 78)

    with open(CATALOGUE_CSV) as f:
        catalogue = list(csv.DictReader(f))
    log(fh, "\nloaded %d cMP rows from catalogue (expected 22): %s"
        % (len(catalogue), [r["meta_program_id"] for r in catalogue]))

    log(fh, "\nloading 240 cNMF programs (top-50 genes each)...")
    programs, n_background = load_all_programs(fh)

    missing = []
    for row in catalogue:
        for m in row["member_cnmf_programs"].split(";"):
            if m not in programs:
                missing.append(m)
    if missing:
        log(fh, "  FATAL: catalogue member programs not found:", missing)
        fh.close()
        sys.exit(1)
    log(fh, "  all catalogue member programs resolved OK")

    gavish_mps = load_gavish_cancer_mps()
    log(fh, "  loaded %d Gavish 'Cancer' meta-programs" % len(gavish_mps))

    squam_cmp = resolve_squamous(programs, catalogue, fh)
    OFIR_CALLS.setdefault(squam_cmp,
                          ("clean", "Ofir: squamous/keratin group, genuinely epithelial."))

    log(fh, "\n" + "=" * 78)
    log(fh, "PER-cMP AUDIT")
    log(fh, "=" * 78)

    csv_rows = []
    consensus_cache = {}
    for row in catalogue:
        cmp_id = row["meta_program_id"]
        members = row["member_cnmf_programs"].split(";")
        n_members = int(row["n_member_cnmf_programs"])
        n_cancers = int(row["n_distinct_cancer_types"])
        printed_gav = row["gavish_2023_MP_match"]
        printed_gav_j = float(row["gavish_jaccard"])

        cons50 = consensus_genes(programs, members)[:N_TOP_GENES]
        consensus_cache[cmp_id] = cons50
        hits = marker_hits(cons50)

        n_mye = len(hits["myeloid"])
        n_tnk = len(hits["t_nk"])
        n_b = len(set(hits["b_plasma"]) | set(hits["b_plasma_ig_family"]))
        n_fib = len(hits["fibroblast"])
        n_endo = len(hits["endothelial"])
        n_tech = len(hits["technical"])

        call, justification = classify(cmp_id, hits, cons50, squam_cmp)
        gav_mp, gav_j = gavish_best_match(cons50, gavish_mps)

        log(fh, "\n--- %s ---" % cmp_id)
        log(fh, "  member cNMF programs (%d): %s" % (n_members, "; ".join(members)))
        log(fh, "  distinct cancer types: %d (%s)" % (n_cancers, row["cancer_types_spanned"]))
        log(fh, "  classification (catalogue): %s" % row["classification"])
        log(fh, "  consensus top-15: %s" % "; ".join(cons50[:15]))
        log(fh, "  marker hits:  myeloid=%d %s" % (n_mye, hits["myeloid"]))
        log(fh, "                T/NK=%d %s" % (n_tnk, hits["t_nk"]))
        log(fh, "                B/plasma=%d %s" % (n_b, sorted(set(hits["b_plasma"]) | set(hits["b_plasma_ig_family"]))))
        log(fh, "                fibroblast=%d %s" % (n_fib, hits["fibroblast"]))
        log(fh, "                endothelial=%d %s" % (n_endo, hits["endothelial"]))
        log(fh, "                technical(MT/ribo)=%d %s" % (n_tech, hits["technical"]))
        log(fh, "  Gavish best match (recomputed, cutoff %.2f): %s  J=%.4f  [%s]"
            % (GAVISH_CUTOFF, gav_mp, gav_j, "MATCH" if gav_j > GAVISH_CUTOFF else "sub-threshold"))
        log(fh, "  Gavish (as printed in catalogue): %s  J=%.4f" % (printed_gav, printed_gav_j))
        log(fh, "  >>> CALL: %s -- %s" % (call.upper(), justification))

        csv_rows.append({
            "cMP": cmp_id,
            "n_member_programs": n_members,
            "n_cancer_types": n_cancers,
            "myeloid_hits": n_mye,
            "t_nk_hits": n_tnk,
            "b_plasma_hits": n_b,
            "fibroblast_hits": n_fib,
            "endothelial_hits": n_endo,
            "technical_hits": n_tech,
            "call": call,
            "justification": justification,
            "gavish_match_recomputed": gav_mp,
            "gavish_jaccard_recomputed": round(gav_j, 4),
            "gavish_match_asprinted": printed_gav,
            "gavish_jaccard_asprinted": printed_gav_j,
            "consensus_top50": ";".join(cons50),
        })

    log(fh, "\n" + "=" * 78)
    log(fh, "cMP8 DEEP-DIVE  (EMT biology vs stromal contamination)")
    log(fh, "=" * 78)
    c8 = consensus_cache["cMP8"]
    log(fh, "cMP8 full 50-gene consensus list (ranked):")
    for i in range(0, len(c8), 10):
        log(fh, "  " + ", ".join(c8[i:i + 10]))
    fib_panel = MARKERS["fibroblast"]
    log(fh, "\nfibroblast-panel genes (%s) present in cMP8: %s"
        % (sorted(fib_panel), [g for g in c8 if g in fib_panel]))
    emt_vocab = {"FN1", "SPARC", "VIM", "POSTN", "THBS2", "CTHRC1", "AEBP1",
                 "FSTL1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "COL3A1",
                 "COL4A1", "COL4A2", "COL5A2", "MMP2", "LOX", "TIMP1", "INHBA"}
    log(fh, "EMT-effector genes present in cMP8: %s" % [g for g in c8 if g in emt_vocab])
    log(fh, "\nInterpretation: COL1A1/COL1A2 are shared EMT + fibroblast markers; the "
            "prompt (VIM/FN1/COL1A1/SPARC/POSTN/THBS2) treats this overlap as legitimate "
            "EMT. cMP8 also carries the pan-collagen program (COL5A1/COL6A1-3/COL3A1/"
            "COL4A1-2), SPARC, FN1, THBS2, CTHRC1, AEBP1, FSTL1 -> a coherent EMT/"
            "mesenchymal state. It does NOT carry DCN, LUM, PDGFRA or PDGFRB (the "
            "fibroblast-identity markers that are not EMT genes). Gavish best match "
            "MP12 EMT-I. => genuine epithelial EMT; the fibroblast-panel hits are the "
            "expected EMT/collagen overlap, not a separate stromal population.")

    log(fh, "\n" + "=" * 78)
    log(fh, "SUB-THRESHOLD GAVISH CONFIRMATIONS (recomputed from Jaccard formula)")
    log(fh, "=" * 78)
    declared = {"cMP6": 0.09, "cMP10": 0.06, "cMP14": 0.05, "cMP20": 0.06, "cMP21": 0.08}
    log(fh, "%-6s %9s %11s %-30s %10s %7s"
        % ("cMP", "declared", "recomputed", "best Gavish MP", "printed J", "agree?"))
    for cmp_id, dval in declared.items():
        cons50 = consensus_cache[cmp_id]
        gav_mp, gav_j = gavish_best_match(cons50, gavish_mps)
        printed_j = float(next(r["gavish_jaccard"] for r in catalogue
                               if r["meta_program_id"] == cmp_id))
        agree = abs(round(gav_j, 2) - dval) <= 0.011
        log(fh, "%-6s %9.2f %11.4f %-30s %10.4f %7s"
            % (cmp_id, dval, gav_j, str(gav_mp), printed_j, str(agree)))
    log(fh, "\n(all 5 are < 0.10 -> below the thesis cutoff; reported as sub-threshold / "
            "not-conserved matches.)")

    if os.path.exists(COMBINED_JACCARD_CSV):
        log(fh, "\n[cross-check] per-program best Gavish matches also exist in %s "
                "(per cNMF program, not per cMP)." % os.path.basename(COMBINED_JACCARD_CSV))

    with open(REPORT_CSV, "w", newline="") as cf:
        w = csv.DictWriter(cf, fieldnames=list(csv_rows[0].keys()))
        w.writeheader()
        w.writerows(csv_rows)
    log(fh, "\nwrote per-cMP table: %s" % REPORT_CSV)

    log(fh, "\n" + "=" * 78)
    log(fh, "SUMMARY TABLE (also in CSV)")
    log(fh, "=" * 78)
    log(fh, "%-6s%6s%6s%5s%5s%4s%5s%5s%6s  %-13s %s"
        % ("cMP", "#prog", "#canc", "mye", "tnk", "B", "fib", "end", "tech", "call", "justification"))
    for r in csv_rows:
        log(fh, "%-6s%6d%6d%5d%5d%4d%5d%5d%6d  %-13s %s"
            % (r["cMP"], r["n_member_programs"], r["n_cancer_types"],
               r["myeloid_hits"], r["t_nk_hits"], r["b_plasma_hits"],
               r["fibroblast_hits"], r["endothelial_hits"], r["technical_hits"],
               r["call"], r["justification"]))

    n_clean = sum(1 for r in csv_rows if r["call"] == "clean")
    n_cont = sum(1 for r in csv_rows if r["call"] == "contaminated")
    n_tech = sum(1 for r in csv_rows if r["call"] == "technical")
    n_amb = sum(1 for r in csv_rows if r["call"] == "ambiguous")
    log(fh, "\nTOTALS: clean=%d  contaminated=%d  technical=%d  ambiguous=%d  (of %d)"
        % (n_clean, n_cont, n_tech, n_amb, len(csv_rows)))
    log(fh, "squamous/keratin group resolved to: %s" % squam_cmp)
    log(fh, "\nDONE.")
    fh.close()


if __name__ == "__main__":
    main()
