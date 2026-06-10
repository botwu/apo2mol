#!/usr/bin/env python3
"""Check IFBench cases against blocklists such as Apo2Mol processed indices."""

from __future__ import annotations

import argparse
import csv
import json
import pickle
import re
import sys
from pathlib import Path
from typing import Any

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import IFBenchCase, QualityIssue, REJECT, read_jsonl, write_jsonl


PDB_RE = re.compile(r"(?<![A-Za-z0-9])([0-9][A-Za-z0-9]{3})(?![A-Za-z0-9])")


def extract_pdb_ids(value: Any) -> set[str]:
    ids: set[str] = set()
    if value is None:
        return ids
    if isinstance(value, (list, tuple, set)):
        for item in value:
            ids.update(extract_pdb_ids(item))
        return ids
    for match in PDB_RE.finditer(str(value)):
        ids.add(match.group(1).lower())
    return ids


def load_apo2mol_index(path: Path) -> set[str]:
    with path.open("rb") as handle:
        index = pickle.load(handle)
    ids: set[str] = set()
    for item in index:
        ids.update(extract_pdb_ids(item))
    return ids


def load_blocklist_csv(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for key, value in row.items():
                if "pdb" in key.lower() or key.lower() in {"id", "case_id", "path"}:
                    ids.update(extract_pdb_ids(value))
    return ids


def check_records(records: list[IFBenchCase], blocked_pdb_ids: set[str]) -> tuple[list[IFBenchCase], list[dict[str, Any]]]:
    hits: list[dict[str, Any]] = []
    for record in records:
        case_ids = {record.holo_pdb_id.lower(), record.apo_pdb_id.lower()} - {""}
        overlap = sorted(case_ids & blocked_pdb_ids)
        if overlap:
            record.accepted = False
            record.issues.append(
                QualityIssue(
                    code="pdb_id_leakage",
                    severity=REJECT,
                    message=f"Case overlaps blocked PDB IDs: {overlap}",
                )
            )
            hits.append({"case_id": record.case_id, "overlap_pdb_ids": overlap})
    return records, hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="IFBench manifest JSONL.")
    parser.add_argument("--apo2mol-index", type=Path, help="Apo2Mol selected_index pickle to block by PDB ID.")
    parser.add_argument("--blocklist-csv", action="append", type=Path, default=[], help="Extra CSV blocklist, repeatable.")
    parser.add_argument("--out", type=Path, help="Optional output manifest with leakage hits marked rejected.")
    parser.add_argument("--report-json", type=Path, help="Leakage report JSON.")
    args = parser.parse_args()

    blocked_pdb_ids: set[str] = set()
    if args.apo2mol_index:
        blocked_pdb_ids.update(load_apo2mol_index(args.apo2mol_index))
    for path in args.blocklist_csv:
        blocked_pdb_ids.update(load_blocklist_csv(path))

    records = read_jsonl(args.manifest)
    records, hits = check_records(records, blocked_pdb_ids)
    if args.out:
        write_jsonl(args.out, records)
    report = {
        "blocked_pdb_ids": len(blocked_pdb_ids),
        "num_cases": len(records),
        "num_leakage_hits": len(hits),
        "hits": hits,
    }
    if args.report_json:
        args.report_json.parent.mkdir(parents=True, exist_ok=True)
        args.report_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
