#!/bin/bash
# Robust per-case Vina eval driver.
# Runs eval_split.py independently per case so a single docking crash
# (e.g. a Vina subprocess segfault on one molecule) does not abort the
# whole arm. Each case gets a 300s timeout. Results are appended to
# summary.txt for later aggregation.
#
# Usage: run_vina_percase.sh <arm_name>
#   expects /root/autodl-tmp/vina_eval/<arm_name>/sampled_results/*result_*.pt
set -u
# Ensure docking CLIs resolve even when launched from a non-login shell
# (setsid/nohup do not source ~/.bashrc, so conda base bin is absent from
# PATH). pdb2pqr30 lives in the base env bin; prepare_receptor4 in apo2mol.
export PATH=/root/miniconda3/envs/apo2mol/bin:/root/miniconda3/bin:$PATH
PY=/root/miniconda3/envs/apo2mol/bin/python
REPO=/root/autodl-tmp/apo2mol
cd "$REPO"
NAME="$1"
S=/root/autodl-tmp/vina_eval/$NAME/sampled_results
BASEOUT=/root/autodl-tmp/vina_eval/$NAME/percase
rm -rf "$BASEOUT"; mkdir -p "$BASEOUT"
N=$(ls "$S"/*result_*.pt | wc -l)
echo "ARM=$NAME N=$N"
: > "$BASEOUT/summary.txt"
for i in $(seq 0 $((N-1))); do
  O="$BASEOUT/case_$i"; mkdir -p "$O"
  timeout 300 "$PY" eval_split.py --sample_path "$S" --result_path "$O" \
    --docking_mode vina_score --pocket_type gen \
    --eval_start_index "$i" --eval_end_index "$i" \
    --protein_root ./apo2mol_dataset/data_folder > "$O/run.log" 2>&1
  rc=$?
  v=$(grep -E "Vina Min" "$O/run.log" | tail -1 | sed -E 's/.*Mean:[[:space:]]*([-0-9.naN]+).*/\1/')
  q=$(grep -E "QED:" "$O/run.log" | tail -1 | sed -E 's/.*Mean:[[:space:]]*([-0-9.naN]+).*/\1/')
  sa=$(grep -E "SA:" "$O/run.log" | tail -1 | sed -E 's/.*Mean:[[:space:]]*([-0-9.naN]+).*/\1/')
  em=$(grep -E "evaluated mols" "$O/run.log" | tail -1 | sed -E 's/.*evaluated mols:[[:space:]]*([0-9]+).*/\1/')
  echo "case=$i rc=$rc evaluated=${em:-0} vina=${v:-NA} qed=${q:-NA} sa=${sa:-NA}" | tee -a "$BASEOUT/summary.txt"
done
echo "ARM_DONE=$NAME"
