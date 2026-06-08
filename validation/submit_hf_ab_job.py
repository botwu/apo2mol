#!/usr/bin/env python3
"""Submit the Apo2Mol hard-tail A/B experiment to Hugging Face Jobs.

This script assumes the current validation patches have been pushed to a git
branch that HF Jobs can clone. The job needs an authenticated Hugging Face
account and a token with access to the gated Apo2Mol dataset.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def capture(cmd: list[str]) -> str:
    try:
        return subprocess.check_output(cmd, cwd=REPO_ROOT, text=True).strip()
    except subprocess.CalledProcessError:
        return ""


def default_repo_url() -> str:
    return capture(["git", "remote", "get-url", "origin"]) or "https://github.com/AIDD-LiLab/Apo2Mol.git"


def default_ref() -> str:
    return capture(["git", "branch", "--show-current"]) or "main"


def build_remote_command(args: argparse.Namespace) -> str:
    return (
        "set -euo pipefail; "
        "git clone --depth 1 --branch \"$REPO_REF\" \"$REPO_URL\" /workspace/Apo2Mol; "
        "cd /workspace/Apo2Mol; "
        "bash validation/hf_job_entrypoint.sh"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-url", default=default_repo_url())
    parser.add_argument("--repo-ref", default=default_ref())
    parser.add_argument("--flavor", default="l4x1")
    parser.add_argument("--timeout", default="8h")
    parser.add_argument("--image", default="pytorch/pytorch:2.4.0-cuda12.4-cudnn9-runtime")
    parser.add_argument("--num-cases", type=int, default=8)
    parser.add_argument("--num-samples", type=int, default=3)
    parser.add_argument("--num-steps", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=3)
    parser.add_argument("--run-dir", default="validation/ab_runs/hf_gpu")
    parser.add_argument("--docking-mode", default="none")
    parser.add_argument("--include-controls", action="store_true")
    parser.add_argument("--results-repo", default="")
    parser.add_argument("--detach", action="store_true", default=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    dirty = capture(["git", "status", "--short"])
    if dirty:
        print(
            "WARNING: local worktree has uncommitted changes. Push the branch with these "
            "validation patches before submitting, or the remote job will not see them.",
            file=sys.stderr,
        )

    env_pairs = {
        "REPO_URL": args.repo_url,
        "REPO_REF": args.repo_ref,
        "RUN_DIR": args.run_dir,
        "NUM_CASES": str(args.num_cases),
        "NUM_SAMPLES": str(args.num_samples),
        "NUM_STEPS": str(args.num_steps),
        "BATCH_SIZE": str(args.batch_size),
        "DOCKING_MODE": args.docking_mode,
        "DEVICE": "auto",
        "INCLUDE_CONTROLS": "1" if args.include_controls else "0",
    }
    if args.results_repo:
        env_pairs["HF_RESULTS_REPO"] = args.results_repo

    cmd = [
        str(REPO_ROOT / ".venv310" / "bin" / "hf"),
        "jobs",
        "run",
        "--flavor",
        args.flavor,
        "--timeout",
        args.timeout,
        "--secrets",
        "HF_TOKEN",
    ]
    if args.detach:
        cmd.append("--detach")
    for key, value in env_pairs.items():
        cmd.extend(["--env", f"{key}={value}"])
    cmd.extend([args.image, "bash", "-lc", build_remote_command(args)])

    print(shlex.join(cmd))
    if args.dry_run:
        return

    whoami = subprocess.run(
        [str(REPO_ROOT / ".venv310" / "bin" / "hf"), "auth", "whoami"],
        cwd=REPO_ROOT,
        text=True,
    )
    if whoami.returncode != 0:
        raise SystemExit("HF CLI is not logged in. Run `.venv310/bin/hf auth login` first.")

    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


if __name__ == "__main__":
    main()
