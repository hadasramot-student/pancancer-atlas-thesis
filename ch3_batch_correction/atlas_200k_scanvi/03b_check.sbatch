#!/bin/bash
#SBATCH --job-name=fig12_03b_check
#SBATCH --chdir=/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi/logs
#SBATCH --output=03b_check.%j.out
#SBATCH --error=03b_check.%j.err
#SBATCH --partition=main
#SBATCH --time=00:20:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2

source /storage/modules/packages/anaconda/etc/profile.d/conda.sh
conda activate scanvi-env
python - <<'EOF'
import scanpy as sc, numpy as np
p="/groups/ofircohen-group/users/hadasra/temp_objects/figure1_2_scanvi/fig12_200k.h5ad"
a=sc.read_h5ad(p, backed="r")
print("shape:", a.shape)
print("obs cols:", list(a.obs.columns))
print("layers:", list(a.layers.keys()))
print("obs_names[:2]:", list(a.obs_names[:2]))
print("n project_name:", a.obs['project_name'].nunique())
for c in ["cancer_type","Consensus_Cell_Type","source_repo","agreement"]:
    print(f"\n{c}:")
    print(a.obs[c].value_counts().to_string())
print("\n'NA' string counts:", {c:int((a.obs[c].astype(str)=='NA').sum()) for c in a.obs.columns})
sub=a.to_memory()[:300]
X=sub.X.toarray() if hasattr(sub.X,'toarray') else np.asarray(sub.X)
C=sub.layers['counts'].toarray() if hasattr(sub.layers['counts'],'toarray') else np.asarray(sub.layers['counts'])
print("X   min/max/allint:", float(X.min()), round(float(X.max()),3), bool(np.allclose(X,np.round(X))))
print("cnt min/max/allint:", float(C.min()), round(float(C.max()),3), bool(np.allclose(C,np.round(C))))
print("OK")
EOF
