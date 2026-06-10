#!/usr/bin/env python3
"""Create a manual review queue for IFBench cases.

Manual review is part of the benchmark, not a post-hoc nicety. This script
exports compact CSV/Markdown tables for accepted hard cases, accepted evidence
checks, and representative rejects.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import IFBenchCase, read_jsonl


def load_records(path: Path | None, bucket: str) -> list[tuple[str, IFBenchCase]]:
    if path is None or not path.exists():
        return []
    return [(bucket, record) for record in read_jsonl(path)]


def metric(record: IFBenchCase, key: str, default: float = 0.0) -> float:
    value = record.metrics.get(key)
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def issue_codes(record: IFBenchCase) -> list[str]:
    return [issue.code for issue in record.issues]


def evidence_sources(record: IFBenchCase) -> list[str]:
    ligand_evidence = record.evidence.get("ligand_relevance", {})
    return list(ligand_evidence.get("sources", []))


def review_priority(bucket: str, record: IFBenchCase) -> tuple[int, str]:
    tier = record.flexibility_tier or ""
    issues = set(issue_codes(record))
    if bucket == "final" and tier in {"hard", "extreme"}:
        return 100, "accepted_hard_extreme"
    if bucket == "final" and metric(record, "motion_core_residues") > 0:
        return 90, "accepted_motion_core"
    if bucket == "final" and metric(record, "ligand_evidence_positive_count") == 0:
        return 85, "accepted_missing_ligand_evidence"
    if bucket == "rejected" and issues & {"structure_parse_failed", "ligand_sanitize_failed"}:
        return 75, "parser_or_ligand_failure"
    if bucket == "rejected" and issues & {"too_few_pocket_residues", "apo_pocket_has_hetero_ligand"}:
        return 70, "pocket_definition_failure"
    if bucket == "rejected" and issues & {"low_sequence_identity", "low_mapping_coverage", "too_few_common_residues"}:
        return 60, "mapping_failure"
    if bucket == "rejected":
        return 40, "other_reject"
    return 30, f"{bucket}_spotcheck"


def case_row(bucket: str, record: IFBenchCase) -> dict[str, Any]:
    priority, reason = review_priority(bucket, record)
    ligand_evidence = record.evidence.get("ligand_relevance", {})
    return {
        "review_priority": priority,
        "review_reason": reason,
        "bucket": bucket,
        "case_id": record.case_id,
        "source": record.source,
        "holo_pdb_id": record.holo_pdb_id,
        "apo_pdb_id": record.apo_pdb_id,
        "holo_chain_id": record.holo_chain_id,
        "apo_chain_id": record.apo_chain_id,
        "ligand_id": record.ligand_id,
        "split": record.split or "",
        "flexibility_tier": record.flexibility_tier or "",
        "pocket_ca_rmsd": metric(record, "pocket_ca_rmsd", math.nan),
        "pocket_p90_ca_displacement": metric(record, "pocket_p90_ca_displacement", math.nan),
        "motion_core_residues": metric(record, "motion_core_residues", math.nan),
        "contact_change_residues": metric(record, "contact_change_residues", math.nan),
        "ligand_evidence_positive_count": metric(record, "ligand_evidence_positive_count", 0.0),
        "ligand_evidence_negative_count": metric(record, "ligand_evidence_negative_count", 0.0),
        "ligand_evidence_sources": ";".join(evidence_sources(record)),
        "issue_codes": ";".join(issue_codes(record)),
        "holo_structure_file": record.holo_structure_file,
        "apo_structure_file": record.apo_structure_file,
        "ligand_file": record.ligand_file,
        "review_decision": "",
        "reviewer": "",
        "review_notes": "",
        "ligand_evidence_json": json.dumps(ligand_evidence, sort_keys=True),
    }


def select_rows(rows: list[dict[str, Any]], max_rows: int, include_rejected_per_reason: int) -> list[dict[str, Any]]:
    rows = sorted(rows, key=lambda row: (-int(row["review_priority"]), row["case_id"]))
    selected: list[dict[str, Any]] = []
    rejected_by_reason: Counter[str] = Counter()
    for row in rows:
        if len(selected) >= max_rows:
            break
        bucket = row["bucket"]
        reason = row["review_reason"]
        if bucket == "rejected":
            if rejected_by_reason[reason] >= include_rejected_per_reason:
                continue
            rejected_by_reason[reason] += 1
        selected.append(row)
    return selected


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    if not fieldnames:
        fieldnames = ["review_priority", "review_reason", "bucket", "case_id", "review_decision", "review_notes"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# IFBench Manual Review Queue\n"]
    lines.append("Review decisions should be one of: `accept`, `reject`, `needs_fix`, `uncertain`.\n")
    lines.append("| priority | reason | bucket | case | holo | apo | ligand | tier | issues |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for row in rows:
        holo = f"{row['holo_pdb_id']}:{row['holo_chain_id']}"
        apo = f"{row['apo_pdb_id']}:{row['apo_chain_id']}"
        lines.append(
            "| {priority} | {reason} | {bucket} | {case_id} | {holo} | {apo} | {ligand_id} | {tier} | {issues} |".format(
                priority=row["review_priority"],
                reason=row["review_reason"],
                bucket=row["bucket"],
                case_id=row["case_id"],
                holo=holo,
                apo=apo,
                ligand_id=row["ligand_id"],
                tier=row["flexibility_tier"],
                issues=row["issue_codes"] or "-",
            )
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted", type=Path)
    parser.add_argument("--rejected", type=Path)
    parser.add_argument("--final", type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--out-md", type=Path)
    parser.add_argument("--max-rows", type=int, default=200)
    parser.add_argument("--include-rejected-per-reason", type=int, default=20)
    args = parser.parse_args()

    tagged_records = []
    tagged_records.extend(load_records(args.final, "final"))
    tagged_records.extend(load_records(args.accepted, "accepted"))
    tagged_records.extend(load_records(args.rejected, "rejected"))
    rows = [case_row(bucket, record) for bucket, record in tagged_records]
    selected = select_rows(rows, max_rows=args.max_rows, include_rejected_per_reason=args.include_rejected_per_reason)
    write_csv(args.out_csv, selected)
    if args.out_md:
        write_markdown(args.out_md, selected)
    summary = {
        "input_records": len(rows),
        "review_rows": len(selected),
        "out_csv": str(args.out_csv),
        "out_md": str(args.out_md) if args.out_md else None,
        "by_reason": dict(Counter(row["review_reason"] for row in selected)),
        "by_bucket": dict(Counter(row["bucket"] for row in selected)),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
