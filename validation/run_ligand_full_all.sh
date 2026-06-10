#!/bin/bash
# Serial paper-aligned ligand eval over all A/B arms (single-process; 2GB box).
set -u
export PATH=/root/miniconda3/envs/apo2mol/bin:/root/miniconda3/bin:$PATH
PY=/root/miniconda3/envs/apo2mol/bin/python
REPO=/root/autodl-tmp/apo2mol
cd "$REPO"
export PYTHONPATH="$REPO"
BASE=/root/autodl-tmp/vina_eval
REFCACHE=$BASE/ref_vina_cache.json
ARMS="baseline_static5 ours_shell4_w025 random_top4 distance_top4_hard"
for NAME in $ARMS; do
  echo "START $NAME $(date)" >> $BASE/ligand_full_progress.log
  "$PY" validation/eval_ligand_full.py \
    --sample_path "$BASE/$NAME/sampled_results" \
    --out "$BASE/$NAME/ligand_full.json" \
    --protein_root ./apo2mol_dataset/data_folder \
    --ref-cache "$REFCACHE" > "$BASE/$NAME/ligand_full.log" 2>&1
  echo "DONE $NAME $(date)" >> $BASE/ligand_full_progress.log
done
echo "ALL_LIGAND_FULL_DONE" > $BASE/ligand_full_done.flag
