#!/usr/bin/env python3
"""Preflight check for running Apo2Mol/LASC on a remote server.

This script intentionally avoids importing project modules. It checks the file
layout, required assets, data samples, optional LMDB cache, and key Python/GPU
dependencies before launching training.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import os
import pickle
import subprocess
import sys
from pathlib import Path


REQUIRED_FILES = [
    "train_pl.py",
    "sample_split.py",
    "configs/config.yaml",
    "configs/training.yaml",
    "configs/sampling.yaml",
    "models/ligand_residue_cross_attn.py",
    "models/molopt_score_model.py",
    "models/pl_model.py",
    "apo2mol_dataset/apo2mol_checkpoint.ckpt",
    "apo2mol_dataset/apo2mol_version/split_druglike.pt",
    "apo2mol_dataset/apo2mol_version/selected_index_apo_druglike.pkl",
]

REQUIRED_MODULES = [
    "torch",
    "torch_geometric",
    "torch_scatter",
    "torch_sparse",
    "torch_cluster",
    "torch_spline_conv",
    "pytorch_lightning",
    "sklearn",
    "rdkit",
    "Bio",
    "lmdb",
    "omegaconf",
    "kornia",
    "wandb",
    "einops",
    "psutil",
    "quaternion",
]


def status(ok: bool) -> str:
    return "OK" if ok else "MISSING"


def print_check(name: str, ok: bool, detail: str = "") -> None:
    suffix = f" - {detail}" if detail else ""
    print(f"[{status(ok):7}] {name}{suffix}")


def import_module(name: str):
    if importlib.util.find_spec(name) is None:
        return False, "not found"
    try:
        module = importlib.import_module(name)
    except Exception as exc:
        return False, repr(exc)
    version = getattr(module, "__version__", None)
    return True, str(version) if version is not None else ""


def run_git_status(repo_root: Path) -> None:
    git_dir = repo_root / ".git"
    if not git_dir.exists():
        print_check("git repository", False, ".git not found")
        return
    try:
        result = subprocess.run(
            ["git", "status", "--short"],
            cwd=repo_root,
            check=False,
            text=True,
            capture_output=True,
        )
    except OSError as exc:
        print_check("git status", False, str(exc))
        return
    dirty_lines = [line for line in result.stdout.splitlines() if line.strip()]
    print_check("git repository", result.returncode == 0, f"{len(dirty_lines)} local change(s)")


def check_torch_cuda(expected_gpus: int) -> bool:
    torch_ok, _ = import_module("torch")
    if not torch_ok:
        print_check("torch cuda", False, "torch is not importable")
        return False
    try:
        import torch
    except Exception as exc:
        print_check("torch import", False, repr(exc))
        return False
    print_check("torch import", True, torch.__version__)
    cuda_ok = torch.cuda.is_available()
    if cuda_ok:
        devices = []
        for idx in range(torch.cuda.device_count()):
            try:
                devices.append(torch.cuda.get_device_name(idx))
            except Exception:
                devices.append(f"cuda:{idx}")
        enough_devices = len(devices) >= expected_gpus
        print_check(
            "cuda devices",
            enough_devices,
            f"{len(devices)} device(s), expected >= {expected_gpus}: {devices}",
        )
        return enough_devices
    else:
        print_check("cuda devices", False, "torch.cuda.is_available() is false")
        return False


def check_required_files(repo_root: Path) -> bool:
    all_ok = True
    for rel in REQUIRED_FILES:
        path = repo_root / rel
        ok = path.exists()
        all_ok = all_ok and ok
        detail = str(path) if ok else f"expected at {path}"
        print_check(rel, ok, detail)
    return all_ok


def check_data_layout(repo_root: Path, data_dir: Path, sample_cases: int) -> bool:
    all_ok = True
    data_ok = data_dir.exists() and data_dir.is_dir()
    all_ok = all_ok and data_ok
    detail = str(data_dir.resolve()) if data_ok else f"expected at {data_dir}"
    print_check("data/data_folder", data_ok, detail)

    link_path = repo_root / "apo2mol_dataset" / "data_folder"
    link_ok = link_path.exists()
    link_detail = ""
    if link_path.is_symlink():
        link_detail = f"symlink -> {os.readlink(link_path)}"
    elif link_ok:
        link_detail = "exists but is not a symlink"
    else:
        link_detail = f"expected symlink or directory at {link_path}"
    print_check("apo2mol_dataset/data_folder", link_ok, link_detail)
    all_ok = all_ok and link_ok

    index_path = repo_root / "apo2mol_dataset" / "apo2mol_version" / "selected_index_apo_druglike.pkl"
    if not index_path.exists():
        print_check("sample file check", False, "index file missing")
        return False

    try:
        with index_path.open("rb") as handle:
            index = pickle.load(handle)
    except Exception as exc:
        print_check("load selected index", False, repr(exc))
        return False

    print_check("load selected index", True, f"{len(index)} entries")
    checked = 0
    missing: list[str] = []
    for item in index:
        if checked >= sample_cases:
            break
        if not item or item[0] is None or item[1] is None:
            continue
        holo_pocket, apo_pocket, ligand = item[:3]
        for rel in (holo_pocket, apo_pocket, ligand):
            if rel and not (data_dir / rel).exists():
                missing.append(rel)
        checked += 1

    ok = checked > 0 and not missing
    all_ok = all_ok and ok
    detail = f"checked {checked} index entries"
    if missing:
        detail += f"; missing examples: {missing[:10]}"
    print_check("sample data files", ok, detail)
    return all_ok


def check_lmdb(repo_root: Path) -> None:
    lmdb_path = repo_root / "apo2mol_dataset" / "data_folder_both_apo2mol_final.lmdb"
    if lmdb_path.exists():
        print_check("processed LMDB cache", True, str(lmdb_path))
    else:
        print_check(
            "processed LMDB cache",
            False,
            "not found; first run will build it, use 1 GPU/process for the first run",
        )


def check_modules(expected_gpus: int) -> bool:
    all_ok = True
    for name in REQUIRED_MODULES:
        ok, detail = import_module(name)
        all_ok = all_ok and ok
        print_check(f"python module {name}", ok, detail)
    cuda_ok = check_torch_cuda(expected_gpus)
    return all_ok and cuda_ok


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="Path to the Apo2Mol repository root.",
    )
    parser.add_argument(
        "--data-dir",
        default=None,
        help="Path to data_folder. Defaults to <repo-root>/data/data_folder.",
    )
    parser.add_argument("--sample-cases", type=int, default=20)
    parser.add_argument(
        "--expected-gpus",
        type=int,
        default=1,
        help="Minimum number of CUDA devices required for this run.",
    )
    parser.add_argument("--skip-modules", action="store_true")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    data_dir = Path(args.data_dir).resolve() if args.data_dir else repo_root / "data" / "data_folder"

    print(f"Repo root: {repo_root}")
    print(f"Data dir : {data_dir}")
    print()

    ok_files = check_required_files(repo_root)
    print()
    ok_data = check_data_layout(repo_root, data_dir, args.sample_cases)
    print()
    check_lmdb(repo_root)
    print()
    run_git_status(repo_root)
    print()
    ok_modules = True
    if not args.skip_modules:
        ok_modules = check_modules(args.expected_gpus)
        print()

    if ok_files and ok_data and ok_modules:
        print("RESULT: OK - remote setup is ready for a smoke training run.")
        sys.exit(0)
    print("RESULT: NOT READY - fix missing items above before training.")
    sys.exit(1)


if __name__ == "__main__":
    main()
