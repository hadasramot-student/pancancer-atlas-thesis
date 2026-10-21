"""
Step 3 — aggregate raw cell counts by cancer type across all three repositories.

Inputs are the verified per-project counts from step 2 and this step's
own investigation (ScBaseCount split children, out-of-scope bucketing,
the ERX10956700 double-count fix). Every number below traces back to a
script run recorded in this conversation; see step3_summary.md for the
full derivation and audit list.

Step-5 correction applied here: SRP434473 (Lung, part of the 159 clean
single-type ScBaseCount studies) was duplicated as two identical rows in
scbasecount_cancer_types.tsv, inflating sc_clean["lung"] by 5,400 cells.
Fixed below (1,236,921 -> 1,231,521); see step5_full_status_audit.md.
Study/project counts are unaffected (it was a duplicated row, not an
extra study) -- sc_clean_projects["lung"] stays 26.
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- 3CA: 90 studies, folder-based cancer type, Hematologic -> Other_Models ----
tca = {
    "brain": 255713, "breast": 850986, "colorectal": 391772,
    "head and neck": 233067, "kidney": 376003, "liver/biliary": 166928,
    "lung": 501008, "neuroendocrine": 142840, "ovarian": 263690,
    "pancreatic": 320096, "prostate": 57716, "sarcoma": 80682, "skin": 364904,
}
tca_projects = {
    "brain": 11, "breast": 13, "colorectal": 2, "head and neck": 4,
    "kidney": 5, "liver/biliary": 5, "lung": 10, "neuroendocrine": 4,
    "ovarian": 10, "pancreatic": 5, "prostate": 2, "sarcoma": 2, "skin": 8,
}
tca_other_models = {"Other_Models": 272478}  # Hematologic, 9 studies - see flag below
tca_other_models_projects = {"Other_Models": 9}

# ---- ScBaseCount: 159 clean single-type studies ----
sc_clean = {
    "brain": 22726, "breast": 2822866, "colorectal": 877073, "kidney": 29868,
    "liver/biliary": 33131, "lung": 1231521, "ovarian": 1658276,
    "pancreatic": 176406, "prostate": 285591, "skin": 87507,
}
sc_clean_projects = {
    "brain": 2, "breast": 53, "colorectal": 26, "kidney": 3,
    "liver/biliary": 3, "lung": 26, "ovarian": 12, "pancreatic": 11,
    "prostate": 15, "skin": 8,
}

# ---- ScBaseCount: 25 multi-type studies with real per-cancer-type splits ----
# NB: 25 parent studies produce 41 child files (some parents split into 3-5
# cancer types), not 39 - corrected count, previously mislabeled. Cell sums
# were always computed from the full 41-row table so are unaffected.
sc_split_real = {
    "brain": 6012, "breast": 1028878, "colorectal": 63475, "head and neck": 3025,
    "kidney": 6413, "liver/biliary": 69379, "lung": 356064, "ovarian": 74083,
    "pancreatic": 15776, "prostate": 105804, "sarcoma": 1154, "skin": 9483,
}
sc_split_projects = {
    "brain": 1, "breast": 8, "colorectal": 7, "head and neck": 1, "kidney": 1,
    "liver/biliary": 2, "lung": 9, "ovarian": 4, "pancreatic": 3,
    "prostate": 3, "sarcoma": 1, "skin": 1,
}
sc_split_unknown = 281085  # cells labeled "unknown" tissue, excluded from every split
# no project count for split_unknown: it's a leftover cell fraction within
# the 25 parent studies, not a distinct file/project of its own

# ---- ScBaseCount: 53 out-of-scope/non-tumor studies ----
sc_53_other_models = 672039  # 24 studies: real cancer type outside the 13 (gastric, cervical, etc.)
sc_53_other_models_projects = 24
sc_53_unknown = 289396       # 29 studies: non-tumor / too vague to classify
sc_53_unknown_projects = 29

# ---- Breastscape: 28 studies (29 minus ERX10956700, which duplicates ScBaseCount's ERP148297) ----
breastscape = {"breast": 756858}
breastscape_projects = {"breast": 28}

# ---- Assemble ----
categories = ["brain", "breast", "colorectal", "head and neck", "kidney",
              "liver/biliary", "lung", "neuroendocrine", "ovarian",
              "pancreatic", "prostate", "sarcoma", "skin"]

totals = {c: 0 for c in categories}
totals["Other_Models"] = 0
totals["Unknown"] = 0
proj = {c: 0 for c in categories}
proj["Other_Models"] = 0
proj["Unknown"] = 0

for src in (tca, sc_clean, sc_split_real, breastscape):
    for k, v in src.items():
        totals[k] += v

for src in (tca_projects, sc_clean_projects, sc_split_projects, breastscape_projects):
    for k, v in src.items():
        proj[k] += v

totals["Other_Models"] += tca_other_models["Other_Models"] + sc_53_other_models
totals["Unknown"] += sc_split_unknown + sc_53_unknown
proj["Other_Models"] += tca_other_models_projects["Other_Models"] + sc_53_other_models_projects
proj["Unknown"] += sc_53_unknown_projects  # split_unknown has no project count of its own

grand_total = sum(totals.values())
grand_total_projects = sum(proj.values())

df = pd.DataFrame([(k, v, proj[k]) for k, v in totals.items()],
                   columns=["cancer_type", "n_cells", "n_projects"])
df["pct"] = df["n_cells"] / grand_total * 100
df = df.sort_values("n_cells", ascending=False).reset_index(drop=True)
df = df[["cancer_type", "n_cells", "pct", "n_projects"]]

df.to_csv("step3_cancer_type_table.csv", index=False)
print(df.to_string(index=False))
print(f"\nGrand total: {grand_total:,} cells")
print(f"Expected (final corrected baseline, post SRP434473 fix): 15,241,772")
print(f"Match: {grand_total == 15241772}")
print(f"\nGrand total project-units: {grand_total_projects}")
print("(354 whole studies + 16 extra units from 25 split parents -> 41 children = 354-25+41)")

# ---- Chart: styled to match atlas_cancer_distribution.py (Fig 1) ----
label_map = {
    "liver/biliary": "Liver/Biliary", "head and neck": "Head & Neck",
    "Other_Models": "Other Models", "Unknown": "Unknown",
}
df["label"] = df["cancer_type"].apply(lambda x: label_map.get(x, x.capitalize()))

# Fixed cancer-type -> hex color mapping, pulled directly from Figure 1's
# pixels. Assigned by name, not sort rank, so colors are stable across every
# figure in the thesis (this one, and Step 4/5 charts to come) regardless of
# how the bars happen to be ranked in any given chart.
COLOR_MAP = {
    "breast": "#fdbf6f", "lung": "#cab2d6", "colorectal": "#a6cee3",
    "Unknown": "#80c680", "pancreatic": "#ffdc91", "ovarian": "#b09c85",
    "kidney": "#7e6148", "skin": "#dc0000", "Other_Models": "#91d1c2",
    "prostate": "#8491b4", "liver/biliary": "#f39b7f", "brain": "#3c5488",
    "neuroendocrine": "#00a087", "sarcoma": "#4dbbd5", "head and neck": "#e64b35",
}
df["color"] = df["cancer_type"].map(COLOR_MAP)
assert df["color"].notna().all(), "missing color for: " + str(df.loc[df["color"].isna(), "cancer_type"].tolist())

df_plot = df.sort_values("n_cells", ascending=True)

fig, ax1 = plt.subplots(figsize=(10, 7))

bars = ax1.barh(df_plot["label"], df_plot["n_cells"] / 1e6, color=df_plot["color"], height=0.7, edgecolor="white")
for bar, (_, row) in zip(bars, df_plot.iterrows()):
    ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
             f"{row['n_cells']/1e6:.2f}M ({row['pct']:.1f}%) — {int(row['n_projects'])} projects",
             va="center", fontsize=9)
ax1.set_xlabel("Number of cells (millions)", fontsize=11)
ax1.set_xlim(0, df_plot["n_cells"].max()/1e6 * 1.6)
ax1.spines[["top","right"]].set_visible(False)
ax1.tick_params(axis="y", labelsize=10)

fig.suptitle(f"Raw, pre-QC baseline: {grand_total/1e6:.2f}M cells across cancer types",
             fontsize=13, fontweight="bold", y=1.02)

plt.tight_layout()
plt.savefig("step3_raw_baseline_by_cancer_type.png", dpi=300, bbox_inches="tight")
plt.savefig("step3_raw_baseline_by_cancer_type.pdf", bbox_inches="tight")
print("\nSaved chart.")

# Save the color map alongside the table so Step 4/5 charts reuse it exactly.
import json
with open("cancer_type_color_map.json", "w") as f:
    json.dump(COLOR_MAP, f, indent=2)
print("Saved cancer_type_color_map.json for reuse in Step 4/5.")
