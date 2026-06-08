#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

PY=".venv310/bin/python"
ACTIVE_RUN="validation/ab_runs/hard20_active_set_candidate_repeat_steps1000_n3"
HARD8_BASELINE_RUN="validation/ab_runs/hard8_original_static5_baseline_steps1000_n3"
HARD20_BASELINE_RUN="validation/ab_runs/hard20_original_static5_baseline_steps1000_n3"
ACTIVE_LOG="${ACTIVE_RUN}/hard20_tmux_cpu.log"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

result_count() {
  local run_dir="$1"
  find "$run_dir" -path "*/sampled_results/result_*.pt" -print 2>/dev/null | wc -l | tr -d " "
}

log_age_seconds() {
  local log_file="$1"
  if [[ ! -f "$log_file" ]]; then
    echo 999999
    return
  fi
  local now mtime
  now="$(date +%s)"
  mtime="$(stat -f "%m" "$log_file")"
  echo $((now - mtime))
}

run_active_set_hard20() {
  echo "[$(timestamp)] Ensuring hard20 active-set run is complete."

  while [[ ! -f "${ACTIVE_RUN}/results.json" ]]; do
    local count age
    count="$(result_count "$ACTIVE_RUN")"
    age="$(log_age_seconds "$ACTIVE_LOG")"
    echo "[$(timestamp)] active-set hard20 status: results=${count}/120, log_age=${age}s"

    if [[ "$count" -ge 120 || "$age" -gt 3600 ]]; then
      echo "[$(timestamp)] Running/resuming active-set hard20 through the official runner."
      "$PY" validation/run_new_method_ab.py \
        --run \
        --active-set-candidate-repeat \
        --run-dir "$ACTIVE_RUN" \
        --num-cases 20 \
        --num-samples 3 \
        --batch-size 1 \
        --num-steps 1000 \
        --docking-mode none \
        --device cpu
    else
      sleep 300
    fi
  done

  echo "[$(timestamp)] Refreshing active-set hard20 summaries."
  "$PY" validation/run_new_method_ab.py \
    --summarize-only \
    --active-set-candidate-repeat \
    --run-dir "$ACTIVE_RUN" \
    --num-cases 20 \
    --num-samples 3 \
    --batch-size 1 \
    --num-steps 1000 \
    --docking-mode none \
    --device cpu
  "$PY" validation/analyze_geometry_diagnostics.py --run-dir "$ACTIVE_RUN"
  "$PY" validation/analyze_ligand_quality.py --run-dir "$ACTIVE_RUN"
}

run_static5_baseline() {
  local run_dir="$1"
  local num_cases="$2"
  echo "[$(timestamp)] Running original static5 baseline: ${run_dir}"
  "$PY" validation/run_new_method_ab.py \
    --run \
    --baseline-only \
    --run-dir "$run_dir" \
    --num-cases "$num_cases" \
    --num-samples 3 \
    --batch-size 1 \
    --num-steps 1000 \
    --docking-mode none \
    --device cpu
  "$PY" validation/analyze_geometry_diagnostics.py --run-dir "$run_dir"
  "$PY" validation/analyze_ligand_quality.py --run-dir "$run_dir"
}

main() {
  echo "[$(timestamp)] Full experiment queue started."
  run_active_set_hard20
  run_static5_baseline "$HARD8_BASELINE_RUN" 8
  run_static5_baseline "$HARD20_BASELINE_RUN" 20
  echo "[$(timestamp)] Full experiment queue completed."
}

main "$@"
