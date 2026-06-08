#!/usr/bin/env bash
set -euo pipefail

cd "${REPO_DIR:-/workspace/Apo2Mol}"

export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/.mplconfig}"
export WANDB_MODE="${WANDB_MODE:-disabled}"
mkdir -p "$MPLCONFIGDIR"

python -m pip install --upgrade "pip<26" "setuptools<81" wheel

python -m pip install \
  numpy==1.24.1 \
  omegaconf==2.3.0 \
  hydra-core==1.3.2 \
  pytorch-lightning==2.2.3 \
  tqdm==4.66.5 \
  easydict==1.13 \
  lmdb==1.5.1 \
  biopython \
  kornia \
  einops==0.8.0 \
  numpy-quaternion==2022.4.4 \
  rdkit \
  openbabel-wheel \
  meeko==0.1.dev3 \
  matplotlib==3.9.2 \
  wandb==0.18.3 \
  huggingface_hub

python -m pip install \
  torch-geometric==2.0.4 \
  torch-scatter==2.1.2 \
  torch-sparse==0.6.18 \
  torch-cluster==1.6.3 \
  torch-spline-conv==1.2.2 \
  -f "https://data.pyg.org/whl/torch-2.4.0+cu124.html"

python - <<'PY'
import torch
print("torch", torch.__version__)
print("cuda_available", torch.cuda.is_available())
if torch.cuda.is_available():
    print("cuda_device", torch.cuda.get_device_name(0))
PY

python validation/download_apo2mol_data.py
python validation/run_new_method_ab.py

RUN_DIR="${RUN_DIR:-validation/ab_runs/hf_gpu}"
run_args=(
  --run \
  --run-dir "$RUN_DIR" \
  --num-cases "${NUM_CASES:-8}" \
  --num-samples "${NUM_SAMPLES:-3}" \
  --num-steps "${NUM_STEPS:-100}" \
  --batch-size "${BATCH_SIZE:-3}" \
  --device "${DEVICE:-auto}" \
  --docking-mode "${DOCKING_MODE:-none}"
)

if [[ "${INCLUDE_CONTROLS:-0}" == "1" || "${INCLUDE_CONTROLS:-0}" == "true" ]]; then
  run_args+=(--include-controls)
fi

python validation/run_new_method_ab.py "${run_args[@]}"

tar -czf ab_results.tgz "$RUN_DIR"
echo "Packed results at $PWD/ab_results.tgz"

if [[ -n "${HF_RESULTS_REPO:-}" ]]; then
  remote_name="apo2mol_ab_results_$(date -u +%Y%m%dT%H%M%SZ).tgz"
  hf upload "$HF_RESULTS_REPO" ab_results.tgz "$remote_name" --type dataset
  echo "Uploaded results to dataset repo: $HF_RESULTS_REPO/$remote_name"
fi
