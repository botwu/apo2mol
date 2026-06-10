"""Run hard20 sampling + Vina with the stage-1 learned cross_attn_gate ckpt.

Thin wrapper around run_new_method_ab.py: prepares only the
`pocket_router_cross_attn_gate` arm, points it at the stage-1 checkpoint,
samples on the same hard20 cases as previous arms, and runs the same
eval_ligand_full.py used for the baseline / fixed-rule arms.

Usage on AutoDL:
    cd /root/autodl-tmp/apo2mol
    python validation/run_hard20_learned_gate.py \
        --ckpt outputs/2026-06-10/.../checkpoints/last.ckpt \
        [--num-samples 3] [--num-cases 20] [--device cuda:0]

If --ckpt is omitted, the script picks the most recently modified ckpt under
outputs/.

The arm directory is laid out like:
    validation/ab_runs/hard20_learned_gate_<ts>/
        hard_test_positions.txt
        hard_original_indices.txt
        selected_cases.json
        pocket_router_cross_attn_gate/
            sampling.yaml
            training.yaml
            samples/result_*.pt
        ligand_full.json   (Vina + QED + High Affinity)

After it finishes, run validation/eval_ligand_full.py if --skip-vina was used.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from glob import glob
from pathlib import Path

# Reuse helpers from the existing AB runner so we stay aligned with how the
# baseline/active-set arms were prepared.
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "validation"))
import run_new_method_ab as ab  # noqa: E402


HARD20_BASELINE_DIR = (
    REPO_ROOT / "validation" / "ab_runs"
    / "hard20_original_static5_baseline_steps1000_n3"
)


def find_latest_ckpt() -> Path | None:
    """Find the most-recent .ckpt under outputs/ (lightning convention)."""
    candidates = []
    for path in glob(str(REPO_ROOT / "outputs" / "*" / "*" / "checkpoints" / "*.ckpt")):
        p = Path(path)
        candidates.append((p.stat().st_mtime, p))
    if not candidates:
        return None
    candidates.sort(reverse=True)
    return candidates[0][1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ckpt", type=Path, default=None,
                    help="path to stage-1 lightning ckpt (auto-detect latest if omitted)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="output dir (default: validation/ab_runs/hard20_learned_gate_<ts>)")
    ap.add_argument("--num-cases", type=int, default=20)
    ap.add_argument("--num-samples", type=int, default=3)
    ap.add_argument("--num-steps", type=int, default=1000)
    ap.add_argument("--batch-size", type=int, default=3)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--hard-threshold", type=float, default=2.0,
                    help="apo-holo RMSD threshold for selecting hard cases (matches baseline)")
    ap.add_argument("--reuse-baseline-ids", action="store_true", default=True,
                    help="use the exact hard20 indices from the existing baseline arm so cases align")
    ap.add_argument("--skip-vina", action="store_true",
                    help="only sample, skip eval_ligand_full")
    args = ap.parse_args()

    ckpt = args.ckpt or find_latest_ckpt()
    if ckpt is None or not ckpt.exists():
        sys.exit("No stage-1 ckpt found under outputs/. Pass --ckpt explicitly.")
    print(f"[ckpt] {ckpt}")

    ts = time.strftime("%m%d_%H%M")
    out_dir = args.out_dir or (
        REPO_ROOT / "validation" / "ab_runs" / f"hard20_learned_gate_{ts}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[out_dir] {out_dir}")

    # 1) Reuse hard20 case list from the baseline arm so it's an apples-to-apples
    #    comparison; copy the three selection files into out_dir.
    if args.reuse_baseline_ids and HARD20_BASELINE_DIR.exists():
        for fname in ("hard_test_positions.txt", "hard_original_indices.txt",
                      "selected_cases.json"):
            src = HARD20_BASELINE_DIR / fname
            if src.exists():
                shutil.copy(src, out_dir / fname)
                print(f"[reuse] {fname}")
            else:
                print(f"[reuse] WARN: {fname} not found in baseline dir")

    ids_path = out_dir / "hard_test_positions.txt"
    original_ids_path = out_dir / "hard_original_indices.txt"
    if not ids_path.exists():
        sys.exit(f"hard test positions file missing: {ids_path}")

    # 2) Build the cross_attn_gate arm payload. We reuse the exact arm spec
    #    from run_new_method_ab so the model overrides match what stage-1 was
    #    trained with (router=cross_attn_gate, no other update gates).
    ablation_arms = ab.make_active_set_ablation_arms()
    arm = next(a for a in ablation_arms if a["name"] == "pocket_router_cross_attn_gate")

    sample_overrides = dict(arm["sample"])
    sample_overrides["num_steps"] = args.num_steps
    sample_overrides["checkpoint"] = str(ckpt)

    model_overrides = dict(arm["model"])

    base_sampling = ab.BASE_SAMPLING_CONFIG.read_text(encoding="utf-8")
    base_training = ab.BASE_TRAINING_CONFIG.read_text(encoding="utf-8")
    sampling_yaml = ab.make_sampling_config(base_sampling, sample_overrides)
    training_yaml = ab.make_training_config(base_training, model_overrides)

    arm_dir = out_dir / arm["name"]
    arm_dir.mkdir(parents=True, exist_ok=True)
    sample_cfg_path = arm_dir / "sampling.yaml"
    train_cfg_path = arm_dir / "training.yaml"
    sample_cfg_path.write_text(sampling_yaml, encoding="utf-8")
    train_cfg_path.write_text(training_yaml, encoding="utf-8")
    print(f"[wrote] {sample_cfg_path}")
    print(f"[wrote] {train_cfg_path}")

    sample_path = arm_dir / "samples"
    eval_path = arm_dir / "eval"
    sample_path.mkdir(parents=True, exist_ok=True)
    eval_path.mkdir(parents=True, exist_ok=True)

    # 3) Run sample_split.py with this arm's sampling.yaml + training.yaml
    print("=== sampling ===")
    cmd = [
        sys.executable, "sample_split.py",
        "--config", str(sample_cfg_path),
        "--train_config", str(train_cfg_path),
        "--device", args.device,
        "--num_samples", str(args.num_samples),
        "--batch_size", str(args.batch_size),
        "--result_path", str(sample_path),
        "--data_ids_file", str(ids_path),
    ]
    print("  cmd:", " ".join(cmd))
    extra_env = {
        "APO2MOL_PREPROCESS_INDEX_FILE": str(original_ids_path),
        "APO2MOL_PROCESSED_PATH": "",  # let sample_split derive default
    }
    env = os.environ.copy()
    env.update({k: v for k, v in extra_env.items() if v})
    subprocess.run(cmd, check=True, cwd=REPO_ROOT, env=env)

    # 4) Vina eval (full-sample + High Affinity, same protocol as previous arms)
    if not args.skip_vina:
        print("=== eval_ligand_full ===")
        eval_script = REPO_ROOT / "validation" / "eval_ligand_full.py"
        if not eval_script.exists():
            print(f"WARN: {eval_script} missing; run it manually later")
        else:
            ligand_full_json = arm_dir / "ligand_full.json"
            cmd2 = [
                sys.executable, str(eval_script),
                "--sample_path", str(sample_path),
                "--out", str(ligand_full_json),
                "--pocket_type", "gen",
            ]
            print("  cmd:", " ".join(cmd2))
            subprocess.run(cmd2, check=True, cwd=REPO_ROOT)

    print(f"DONE. arm_dir={arm_dir}")


if __name__ == "__main__":
    main()
