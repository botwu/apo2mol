#!/usr/bin/env python3
"""Assign flexibility tiers and leakage-aware splits to accepted IFBench cases."""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import IFBenchCase, read_jsonl, write_flat_csv, write_jsonl


def metric(record: IFBenchCase, name: str, default: float = 0.0) -> float:
    value = record.metrics.get(name)
    if value is None:
        return default
    return float(value)


def flexibility_score(record: IFBenchCase) -> float:
    """Composite score centered on induced-fit difficulty, not ligand size."""
    pocket_rmsd = metric(record, "pocket_ca_rmsd")
    p90_disp = metric(record, "pocket_p90_ca_displacement")
    contact_change = metric(record, "contact_change_fraction")
    motion_core = metric(record, "motion_core_fraction")
    return pocket_rmsd + 0.35 * p90_disp + 1.25 * contact_change + 0.75 * motion_core


def assign_tier(record: IFBenchCase) -> str:
    score = flexibility_score(record)
    pocket_rmsd = metric(record, "pocket_ca_rmsd")
    if pocket_rmsd < 0.75 and score < 1.25:
        return "easy"
    if pocket_rmsd < 1.50 and score < 2.25:
        return "medium"
    if pocket_rmsd < 3.00 and score < 4.00:
        return "hard"
    return "extreme"


def group_key(record: IFBenchCase) -> str:
    return record.cluster_key or record.uniprot_id or f"{record.holo_pdb_id}:{record.holo_chain_id}"


def split_groups(records: list[IFBenchCase], seed: int, valid_fraction: float, test_fraction: float) -> dict[str, str]:
    groups: dict[str, list[IFBenchCase]] = defaultdict(list)
    for record in records:
        groups[group_key(record)].append(record)

    group_items = list(groups.items())
    rng = random.Random(seed)
    rng.shuffle(group_items)

    total_cases = len(records)
    target_valid = round(total_cases * valid_fraction)
    target_test = round(total_cases * test_fraction)
    split_counts = {"train": 0, "valid": 0, "test": 0}
    assignment: dict[str, str] = {}

    for key, items in group_items:
        if split_counts["test"] < target_test:
            split = "test"
        elif split_counts["valid"] < target_valid:
            split = "valid"
        else:
            split = "train"
        assignment[key] = split
        split_counts[split] += len(items)
    return assignment


def summarize(records: list[IFBenchCase]) -> dict[str, Any]:
    by_split = Counter(record.split for record in records)
    by_tier = Counter(record.flexibility_tier for record in records)
    by_split_tier = Counter(f"{record.split}/{record.flexibility_tier}" for record in records)
    return {
        "num_cases": len(records),
        "by_split": dict(by_split),
        "by_tier": dict(by_tier),
        "by_split_tier": dict(by_split_tier),
        "mean_flexibility_score": sum(flexibility_score(record) for record in records) / max(1, len(records)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="Accepted IFBench manifest JSONL.")
    parser.add_argument("--out", required=True, type=Path, help="Output JSONL with split and flexibility_tier.")
    parser.add_argument("--flat-csv", type=Path, help="Optional flattened CSV.")
    parser.add_argument("--splits-json", type=Path, help="Optional split case-id lists.")
    parser.add_argument("--summary-json", type=Path, help="Optional summary JSON.")
    parser.add_argument("--seed", type=int, default=20270609)
    parser.add_argument("--valid-fraction", type=float, default=0.10)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    args = parser.parse_args()

    records = [record for record in read_jsonl(args.manifest) if record.accepted]
    for record in records:
        record.flexibility_tier = assign_tier(record)
        record.metrics["flexibility_score"] = flexibility_score(record)

    assignments = split_groups(records, args.seed, args.valid_fraction, args.test_fraction)
    for record in records:
        record.split = assignments[group_key(record)]

    write_jsonl(args.out, records)
    if args.flat_csv:
        write_flat_csv(args.flat_csv, records)
    if args.splits_json:
        splits: dict[str, list[str]] = {"train": [], "valid": [], "test": []}
        for record in records:
            splits[str(record.split)].append(record.case_id)
        args.splits_json.parent.mkdir(parents=True, exist_ok=True)
        args.splits_json.write_text(json.dumps(splits, indent=2, sort_keys=True) + "\n")
    summary = summarize(records)
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
