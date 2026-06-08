#!/usr/bin/env python3
"""Download and place the gated Apo2Mol data_folder.

Requirements:
* The Hugging Face account behind HF_TOKEN must have access to
  AIDD-LiLab/Apo2Mol_Dataset.
* Do not pass the token on the command line. Export HF_TOKEN in the shell.

The script downloads the dataset snapshot into apo2mol_dataset/_hf_download,
then places or extracts data_folder at apo2mol_dataset/data_folder.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tarfile
import zipfile
from pathlib import Path

from huggingface_hub import snapshot_download


REPO_ROOT = Path(__file__).resolve().parents[1]
DATASET_ID = "AIDD-LiLab/Apo2Mol_Dataset"
TARGET_DIR = REPO_ROOT / "apo2mol_dataset" / "data_folder"
DOWNLOAD_DIR = REPO_ROOT / "apo2mol_dataset" / "_hf_download"


def has_pocket_files(path: Path) -> bool:
    return any(path.glob("*/receptor_apo_pocket10.pdb")) and any(path.glob("*/*.sdf"))


def copy_data_folder(source: Path, target: Path) -> None:
    if target.exists():
        raise FileExistsError(f"{target} already exists; remove it manually if you want to replace it")
    shutil.copytree(source, target)


def extract_archive(archive: Path, target: Path) -> bool:
    tmp_extract = target.parent / "_data_extract_tmp"
    if tmp_extract.exists():
        shutil.rmtree(tmp_extract)
    tmp_extract.mkdir(parents=True)

    try:
        if archive.suffix == ".zip":
            with zipfile.ZipFile(archive) as zf:
                zf.extractall(tmp_extract)
        elif archive.name.endswith((".tar.gz", ".tgz", ".tar")):
            with tarfile.open(archive) as tf:
                tf.extractall(tmp_extract)
        else:
            return False

        candidates = [tmp_extract / "data_folder"] + [p for p in tmp_extract.iterdir() if p.is_dir()]
        for candidate in candidates:
            if candidate.exists() and has_pocket_files(candidate):
                copy_data_folder(candidate, target)
                return True
        return False
    finally:
        shutil.rmtree(tmp_extract, ignore_errors=True)


def place_data_folder(download_dir: Path, target: Path) -> None:
    direct = download_dir / "data_folder"
    if direct.exists() and has_pocket_files(direct):
        copy_data_folder(direct, target)
        return

    archives = sorted(
        [
            *download_dir.glob("*.tar.gz"),
            *download_dir.glob("*.tgz"),
            *download_dir.glob("*.tar"),
            *download_dir.glob("*.zip"),
        ]
    )
    for archive in archives:
        if extract_archive(archive, target):
            return

    raise FileNotFoundError(
        "Downloaded snapshot did not contain a recognizable data_folder or archive. "
        f"Inspect {download_dir} manually."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-id", default=DATASET_ID)
    parser.add_argument("--download-dir", type=Path, default=DOWNLOAD_DIR)
    parser.add_argument("--target-dir", type=Path, default=TARGET_DIR)
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    if args.target_dir.exists() and has_pocket_files(args.target_dir):
        print(f"data_folder already exists: {args.target_dir}")
        return

    token = os.environ.get("HF_TOKEN")
    if not token:
        raise SystemExit(
            "HF_TOKEN is not set. Create a Hugging Face token with access to "
            f"{args.repo_id}, then run: HF_TOKEN=... {sys.executable} validation/download_apo2mol_data.py"
        )

    args.download_dir.mkdir(parents=True, exist_ok=True)
    snapshot_download(
        repo_id=args.repo_id,
        repo_type="dataset",
        token=token,
        local_dir=args.download_dir,
        allow_patterns=[
            "data_folder/**",
            "data_folder.tar",
            "data_folder.tar.gz",
            "data_folder.tgz",
            "*.tar",
            "*.tar.gz",
            "*.tgz",
            "*.zip",
        ],
        force_download=args.force_download,
    )
    place_data_folder(args.download_dir, args.target_dir)
    print(f"Placed data_folder at {args.target_dir}")


if __name__ == "__main__":
    main()
