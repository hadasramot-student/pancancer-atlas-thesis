"""Build the project-count audit tables from project_cancer_pairs_raw.csv (read-only on the atlas)."""
import re
import pandas as pd

raw = pd.read_csv("project_cancer_pairs_raw.csv")
assert raw.n_cells.sum() == 8107665


def repo(p, has_srx):
    if has_srx:
        return "ScBaseCount"
    if re.search(r"_et_al_", p, re.I):
        return "Breastscape"
    return "3CA"


def study_key(p):
    # 3CA splits one study into numbered sub-datasets: Pelka2021_Colorectal_1.._4 -> Pelka2021_Colorectal
    return re.sub(r"_\d+$", "", p) if not re.search(r"_et_al_", p, re.I) else p


def author_year(p):
    m = re.match(r"^([A-Za-z_]+?)(\d{4})_", p)  # 3CA: Author2020_Type
    if m:
        return f"{m.group(1).upper()}{m.group(2)}"
    m = re.match(r"^([A-Z_]+?)_[A-Z_]*ET_AL_(\d{4})_", p, re.I)  # Breastscape: XU_K_ET_AL_2021_P32
    if m:
        return f"{m.group(1).split('_')[0].upper()}{m.group(2)}"
    return None


raw["repository"] = [repo(p, s) for p, s in zip(raw.project_name, raw.has_SRX)]
proj = (raw.groupby("project_name")
        .agg(n_cells=("n_cells", "sum"), repository=("repository", "first"),
             cancer_types=("cancer_type", lambda s: ";".join(sorted(s))),
             n_cancer_types=("cancer_type", "nunique"))
        .reset_index())
proj["study_key_suffix_collapsed"] = proj.project_name.map(study_key)
proj["author_year"] = proj.project_name.map(author_year)
proj = proj.sort_values(["repository", "project_name"])
proj.to_csv("projects_207.csv", index=False)

# --- same-study groups ---
rows = []
for key, g in proj.groupby("study_key_suffix_collapsed"):
    if len(g) > 1:
        rows.append(("numeric suffix (_1/_2/...) of one 3CA study", key, g))
for ay, g in proj.dropna(subset=["author_year"]).groupby("author_year"):
    if g.study_key_suffix_collapsed.nunique() > 1:
        rows.append(("same author+year, different name (check manually)", ay, g))
grp = pd.concat([g.assign(group_type=t, group_key=k)[
    ["group_type", "group_key", "project_name", "repository", "cancer_types", "n_cells"]]
    for t, k, g in rows])
grp.to_csv("same_study_groups.csv", index=False)

# --- (project, cancer_type) pairs and per-cancer summary ---
raw[["project_name", "cancer_type", "n_cells", "repository"]].sort_values(
    ["cancer_type", "project_name"]).to_csv("project_cancer_pairs.csv", index=False)
raw["study"] = raw.project_name.map(study_key)
summ = (raw.groupby("cancer_type")
        .agg(n_cells=("n_cells", "sum"), n_project_names=("project_name", "nunique"),
             n_studies_suffix_collapsed=("study", "nunique"))
        .reset_index().sort_values("n_cells", ascending=False))
summ.loc[len(summ)] = ["TOTAL (sum over rows)", summ.n_cells.sum(), summ.n_project_names.sum(),
                       summ.n_studies_suffix_collapsed.sum()]
summ.to_csv("per_cancer_type_summary.csv", index=False)

# --- console report ---
print("project names:", proj.shape[0], " (project, cancer_type) pairs:", raw.shape[0])
print(proj.repository.value_counts().to_string())
print("multi-cancer projects:\n", proj[proj.n_cancer_types > 1][["project_name", "cancer_types"]].to_string(index=False))
print("\nstudies after collapsing numeric suffixes:", proj.study_key_suffix_collapsed.nunique())
print(proj.assign(s=proj.study_key_suffix_collapsed).groupby(["repository"]).s.nunique().to_string())
print("\n", grp.to_string(index=False))
print("\n", summ.to_string(index=False))
