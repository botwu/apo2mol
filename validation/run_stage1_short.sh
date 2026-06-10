#!/bin/bash
# Stage-1 short training: gate-only, freeze backbone, bf16-mixed, single GPU.
# Run on AutoDL instance with WANDB already logged in.
set -euo pipefail

cd /root/autodl-tmp/apo2mol
source /root/miniconda3/etc/profile.d/conda.sh
conda activate apo2mol

# Run identifier with timestamp so wandb / outputs/ don't collide
RUN_NAME="stage1_gateonly_$(date +%m%d_%H%M)"
mkdir -p /tmp/diag_logs

echo "=== Starting stage-1 training: $RUN_NAME ==="
echo "Logs: /tmp/diag_logs/${RUN_NAME}.log"
echo "GPU:" && nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv | head -2
echo

python train_pl.py \
  experiment_name="$RUN_NAME" \
  sys.devices=1 \
  sys.strategy=auto \
  sys.precision=bf16-mixed \
  model.pocket_router_mode=cross_attn_gate \
  model.cross_attn_gate.stage2_two_forward=false \
  train.freeze_backbone=true \
  train.batch_size=32 \
  train.val_batch_size=32 \
  train.num_workers=12 \
  train.val_num_workers=4 \
  train.prefetch_factor=4 \
  train.pin_memory=false \
  train.max_steps=1250 \
  train.val_freq=1 \
  train.limit_val_batches=0.2 \
  train.patience=15 \
  train.train_report_iter=200 \
  train.pretrained_ckpt=./apo2mol_dataset/apo2mol_checkpoint.ckpt \
  wandb.wandb_status=online \
  wandb.wandb_project=apo2mol-stage1 \
  wandb.wandb_task="$RUN_NAME"
