#!/usr/bin/env python
"""
Re-apply the Table B93 (script: tableB97_contamination_audit.py) marker-panel
rules to cMP10, cMP17, cMP15, cMP3. Imports the audit module unchanged so the
panels, consensus-gene construction and classify() rules are the same code.

Reports, per cMP: consensus top-50, per-panel hit counts + genes, the call the
script assigns (incl. hard-coded OFIR_CALLS overrides) and the call from the
automatic rules alone (override bypassed).
"""
import os, sys, csv, re
import pandas as pd

AUDIT_DIR = "/groups/ofircohen-group/users/hadasra/pancancer-atlas-thesis/ch5_gep_programs/meta_programs"
sys.path.insert(0, AUDIT_DIR)
import tableB97_contamination_audit as A   # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
TARGETS = ["cMP10", "cMP17", "cMP15", "cMP3"]
MTRNR2L_RE = re.compile(r"^MTRNR2L\d+$")   # not matched by the audit's technical regex


class _Sink:
    def write(self, s): pass
    def flush(self): pass


with open(A.CATALOGUE_CSV) as f:
    catalogue = {r["meta_program_id"]: r for r in csv.DictReader(f)}
programs, _ = A.load_all_programs(_Sink())
squam = A.resolve_squamous(programs, list(catalogue.values()), _Sink())

rows = []
for cmp_id in TARGETS:
    row = catalogue[cmp_id]
    members = row["member_cnmf_programs"].split(";")
    cons50 = A.consensus_genes(programs, members)[:A.N_TOP_GENES]
    hits = A.marker_hits(cons50)
    b = sorted(set(hits["b_plasma"]) | set(hits["b_plasma_ig_family"]))
    call_script, just_script = A.classify(cmp_id, hits, cons50, squam)
    saved = A.OFIR_CALLS.pop(cmp_id, None)
    call_rule, just_rule = A.classify(cmp_id, hits, cons50, squam)
    if saved is not None:
        A.OFIR_CALLS[cmp_id] = saved
    mtrnr = [g for g in cons50 if MTRNR2L_RE.match(g)]

    print("=" * 78)
    print(f"{cmp_id}  members={members}  cancers={row['cancer_types_spanned']}")
    print("consensus top-50:", ", ".join(cons50))
    for k, v in [("myeloid", hits["myeloid"]), ("T/NK", hits["t_nk"]), ("B/plasma", b),
                 ("fibroblast", hits["fibroblast"]), ("endothelial", hits["endothelial"]),
                 ("technical (MT-/RPL/RPS)", hits["technical"]),
                 ("[not in rules] MTRNR2L*", mtrnr)]:
        print(f"  {k:26s} {len(v):2d}  {v}")
    print(f"  CALL as script assigns : {call_script}  -- {just_script}")
    print(f"  CALL from rules only   : {call_rule}  -- {just_rule}")
    rows.append(dict(cMP=cmp_id, members=";".join(members),
                     cancers=row["cancer_types_spanned"],
                     myeloid=len(hits["myeloid"]), t_nk=len(hits["t_nk"]), b_plasma=len(b),
                     fibroblast=len(hits["fibroblast"]), endothelial=len(hits["endothelial"]),
                     technical=len(hits["technical"]), technical_genes=";".join(hits["technical"]),
                     mtrnr2l_not_counted=len(mtrnr), mtrnr2l_genes=";".join(mtrnr),
                     call_script=call_script, call_rules_only=call_rule,
                     justification_rules_only=just_rule, consensus_top50=";".join(cons50)))

pd.DataFrame(rows).to_csv(os.path.join(OUT, "task2_panel_overlap_4cmps.csv"), index=False)
print("\nwrote", os.path.join(OUT, "task2_panel_overlap_4cmps.csv"))
