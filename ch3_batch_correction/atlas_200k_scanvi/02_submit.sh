#!/bin/bash
# 02_submit.sh [--dry-run]
# Fan 02_assemble_counts.R over the drawn projects. One array task = one project,
# binned into 4 memory classes by _step4.RDS file size. Restartable.
set -euo pipefail
DRY=false; [[ "${1:-}" == "--dry-run" ]] && DRY=true

WT_HOST=/groups/ofircohen-group/users/hadasra/create_sc_object_workflow/.claude/worktrees/figure-1-2-scanvi-embedding-4e0a23/fig12
TO_HOST=/groups/ofircohen-group/users/hadasra/temp_objects
ROOT=$TO_HOST/figure1_2_scanvi
BCDIR=$ROOT/per_project_barcodes
LOG=$ROOT/logs
CL=$ROOT/classlists
mkdir -p "$LOG" "$ROOT/per_project_counts" "$CL"
[[ -d "$BCDIR" ]] || { echo "no $BCDIR -- run 01 first"; exit 1; }

declare -A MEM=(  [S]=64G      [M]=160G     [L]=384G     [XL]=512G )
declare -A TLIM=( [S]=02:00:00 [M]=03:00:00 [L]=06:00:00 [XL]=12:00:00 )
declare -A CONC=( [S]=24       [M]=10       [L]=5        [XL]=1 )
for c in S M L XL; do : > "$CL/$c.txt"; done

for bc in "$BCDIR"/*.txt; do
  p=$(basename "$bc" .txt)
  f="$TO_HOST/${p}_step4.RDS"
  [[ -f "$f" ]] || { echo "!! no step4 for $p"; exit 1; }
  sz=$(stat -c %s "$f")
  if   (( sz > 20000000000 )); then c=XL
  elif (( sz >  8000000000 )); then c=L
  elif (( sz >  2500000000 )); then c=M
  else                              c=S; fi
  echo "$p" >> "$CL/$c.txt"
done

for c in S M L XL; do
  n=$(wc -l < "$CL/$c.txt" | tr -d ' ')
  echo "class $c: $n projects  mem=${MEM[$c]} time=${TLIM[$c]} conc=%${CONC[$c]}"
  [[ "$n" -eq 0 || "$DRY" == true ]] && continue
  sbatch --job-name="fig12_02_$c" --partition=main \
    --mem="${MEM[$c]}" --cpus-per-task=4 --time="${TLIM[$c]}" \
    --array="1-${n}%${CONC[$c]}" \
    --output="$LOG/02_${c}_%A_%a.out" --error="$LOG/02_${c}_%A_%a.err" \
    "$WT_HOST/02_task.sh" "$CL/$c.txt"
done
$DRY || echo "monitor: squeue -u hadasra | grep fig12_02"
