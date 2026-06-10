#!/usr/bin/env python3
"""Build a human-readable IFBench construction report."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Any

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import IFBenchCase, read_jsonl


KEY_METRICS = [
    "ligand_heavy_atoms",
    "ligand_mol_wt",
    "sequence_identity",
    "mapping_coverage",
    "global_ca_rmsd_after_align",
    "pocket_residues",
    "pocket_ca_rmsd",
    "pocket_mean_ca_displacement",
    "pocket_p90_ca_displacement",
    "motion_core_residues",
    "motion_core_fraction",
    "holo_contact_residues",
    "contact_change_residues",
    "contact_change_fraction",
    "ligand_evidence_count",
    "ligand_evidence_positive_count",
    "ligand_evidence_negative_count",
    "flexibility_score",
]


def load_json(path: Path | None) -> dict[str, Any] | None:
    if not path or not path.exists():
        return None
    return json.loads(path.read_text())


def safe_read_jsonl(path: Path | None) -> list[IFBenchCase]:
    if not path or not path.exists():
        return []
    return read_jsonl(path)


def summarize_values(values: list[float]) -> dict[str, float | int | None]:
    clean = [float(value) for value in values if math.isfinite(float(value))]
    if not clean:
        return {"n": 0, "mean": None, "median": None, "min": None, "max": None}
    return {
        "n": len(clean),
        "mean": statistics.fmean(clean),
        "median": statistics.median(clean),
        "min": min(clean),
        "max": max(clean),
    }


def metric_summary(records: list[IFBenchCase]) -> dict[str, dict[str, float | int | None]]:
    summaries: dict[str, dict[str, float | int | None]] = {}
    for key in KEY_METRICS:
        values = []
        for record in records:
            value = record.metrics.get(key)
            if isinstance(value, (int, float)) and value is not None:
                values.append(float(value))
        summaries[key] = summarize_values(values)
    return summaries


def format_num(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    if not rows:
        return "_No rows._\n"
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(format_num(value) for value in row) + " |")
    return "\n".join(lines) + "\n"


def issue_counts(records: list[IFBenchCase]) -> Counter:
    counts: Counter = Counter()
    for record in records:
        for issue in record.issues:
            counts[f"{issue.severity}:{issue.code}"] += 1
    return counts


def build_report(
    accepted: list[IFBenchCase],
    rejected: list[IFBenchCase],
    final: list[IFBenchCase],
    qc_summary: dict[str, Any] | None,
    leakage_report: dict[str, Any] | None,
    split_summary: dict[str, Any] | None,
) -> str:
    lines: list[str] = []
    lines.append("# IFBench Construction Report\n")
    lines.append("## Scope\n")
    lines.append("This report audits candidate quality for the independent induced-fit benchmark. Counts alone are not success; reject reasons and tier balance are the primary sanity checks.\n")

    lines.append("## Candidate Flow\n")
    flow_rows = [
        ["accepted_manifest", len(accepted)],
        ["rejected_manifest", len(rejected)],
        ["final_manifest", len(final)],
    ]
    if qc_summary:
        flow_rows.extend(
            [
                ["qc_num_candidates", qc_summary.get("num_candidates")],
                ["qc_num_accepted", qc_summary.get("num_accepted")],
                ["qc_num_rejected", qc_summary.get("num_rejected")],
            ]
        )
    if leakage_report:
        flow_rows.extend(
            [
                ["blocked_pdb_ids", leakage_report.get("blocked_pdb_ids")],
                ["leakage_hits", leakage_report.get("num_leakage_hits")],
            ]
        )
    lines.append(markdown_table(["item", "value"], flow_rows))

    lines.append("## Reject Reasons\n")
    counts = issue_counts(rejected)
    if qc_summary and qc_summary.get("reject_counts"):
        rows = sorted(qc_summary["reject_counts"].items(), key=lambda item: (-item[1], item[0]))
    else:
        rows = sorted(((key, value) for key, value in counts.items()), key=lambda item: (-item[1], item[0]))
    lines.append(markdown_table(["reason", "count"], rows[:30]))

    lines.append("## Split And Tier Balance\n")
    split_tier_counts = Counter(f"{record.split}/{record.flexibility_tier}" for record in final)
    rows = sorted(split_tier_counts.items())
    lines.append(markdown_table(["split/tier", "cases"], rows))
    if split_summary:
        lines.append("Raw split summary:\n")
        lines.append("```json\n" + json.dumps(split_summary, indent=2, sort_keys=True) + "\n```\n")

    lines.append("## Key Metric Summary\n")
    summaries = metric_summary(final or accepted)
    metric_rows = [
        [key, item["n"], item["mean"], item["median"], item["min"], item["max"]]
        for key, item in summaries.items()
        if item["n"]
    ]
    lines.append(markdown_table(["metric", "n", "mean", "median", "min", "max"], metric_rows))

    lines.append("## Quality Interpretation Checklist\n")
    lines.append("- Accept rate should be low enough to show the gates are doing work, but not dominated by one avoidable parsing failure.\n")
    lines.append("- Hard/extreme tiers must exist in the final test split; otherwise the benchmark does not expose induced-fit failure modes.\n")
    lines.append("- Leakage hits must be zero after filtering against Apo2Mol and any additional blocklists.\n")
    lines.append("- Motion-core and contact-change labels should not collapse to all-zero or all-one distributions.\n")
    lines.append("- Manual inspection should cover the top rejected reasons and at least 20 accepted hard/extreme cases before public release.\n")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted", type=Path, help="Accepted manifest JSONL.")
    parser.add_argument("--rejected", type=Path, help="Rejected manifest JSONL.")
    parser.add_argument("--final", type=Path, help="Final manifest JSONL after leakage, stratification, and split.")
    parser.add_argument("--qc-summary", type=Path)
    parser.add_argument("--leakage-report", type=Path)
    parser.add_argument("--split-summary", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    report = build_report(
        accepted=safe_read_jsonl(args.accepted),
        rejected=safe_read_jsonl(args.rejected),
        final=safe_read_jsonl(args.final),
        qc_summary=load_json(args.qc_summary),
        leakage_report=load_json(args.leakage_report),
        split_summary=load_json(args.split_summary),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(report)
    print(json.dumps({"out": str(args.out), "bytes": len(report.encode("utf-8"))}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
