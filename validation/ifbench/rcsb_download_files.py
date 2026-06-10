#!/usr/bin/env python3
"""Download RCSB structure and ligand coordinate files for IFBench candidates."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


STRUCTURE_URL = "https://files.rcsb.org/download/{entry_id}.cif"
LIGAND_URL = "https://models.rcsb.org/v1/{entry_id}/ligand"


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def resolve_path(raw_root: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = raw_root / path
    return path


def download(url: str, out_file: Path, timeout: int, retries: int) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                data = response.read()
            if not data:
                raise RuntimeError("empty response")
            out_file.write_bytes(data)
            return
        except (urllib.error.URLError, TimeoutError, RuntimeError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Failed to download {url}: {last_error}")


def ligand_url(row: dict[str, str], encoding: str) -> str:
    query: dict[str, Any] = {
        "label_comp_id": row["ligand_id"],
        "encoding": encoding,
    }
    if row.get("ligand_label_asym_id"):
        query["label_asym_id"] = row["ligand_label_asym_id"]
    elif row.get("ligand_auth_asym_id"):
        query["auth_asym_id"] = row["ligand_auth_asym_id"]
    return f"{LIGAND_URL.format(entry_id=row['holo_pdb_id'].lower())}?{urllib.parse.urlencode(query)}"


def unique_structure_downloads(rows: list[dict[str, str]], raw_root: Path) -> list[tuple[str, Path]]:
    downloads: dict[Path, str] = {}
    for row in rows:
        for id_key, file_key in [("holo_pdb_id", "holo_structure_file"), ("apo_pdb_id", "apo_structure_file")]:
            if row.get(id_key) and row.get(file_key):
                downloads[resolve_path(raw_root, row[file_key])] = STRUCTURE_URL.format(entry_id=row[id_key].lower())
    return [(url, path) for path, url in sorted(downloads.items(), key=lambda item: str(item[0]))]


def unique_ligand_downloads(rows: list[dict[str, str]], raw_root: Path, encoding: str) -> list[tuple[str, Path]]:
    downloads: dict[Path, str] = {}
    for row in rows:
        if row.get("ligand_file"):
            downloads[resolve_path(raw_root, row["ligand_file"])] = ligand_url(row, encoding=encoding)
    return [(url, path) for path, url in sorted(downloads.items(), key=lambda item: str(item[0]))]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--raw-root", default=".", type=Path, help="Root used to resolve candidate file paths.")
    parser.add_argument("--ligand-encoding", default="sdf", choices=["sdf", "mol", "mol2", "cif"])
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=2)
    args = parser.parse_args()

    rows = read_candidates(args.candidates)
    downloads = unique_structure_downloads(rows, args.raw_root)
    downloads.extend(unique_ligand_downloads(rows, args.raw_root, encoding=args.ligand_encoding))
    planned = []
    completed = 0
    skipped = 0
    failed: list[dict[str, str]] = []
    for url, out_file in downloads:
        planned.append({"url": url, "out_file": str(out_file)})
        if args.skip_existing and out_file.exists() and out_file.stat().st_size > 0:
            skipped += 1
            continue
        if args.dry_run:
            continue
        try:
            print(f"download {url} -> {out_file}", file=sys.stderr)
            download(url, out_file, timeout=args.timeout, retries=args.retries)
            completed += 1
        except Exception as exc:
            failed.append({"url": url, "out_file": str(out_file), "error": repr(exc)})

    summary = {
        "candidate_rows": len(rows),
        "planned_downloads": len(downloads),
        "completed": completed,
        "skipped_existing": skipped,
        "failed": failed,
    }
    if args.dry_run:
        summary["dry_run_plan"] = planned[:20]
    print(json.dumps(summary, indent=2, sort_keys=True))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
