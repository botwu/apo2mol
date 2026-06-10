#!/usr/bin/env python3
"""Gate an IFBench build against release criteria."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import IFBenchCase, read_jsonl


DEFAULT_GATE = {
    "min_final_cases": 200,
    "min_hard_extreme_cases": 40,
    "min_test_cases": 40,
    "min_test_hard_extreme_cases": 20,
    "min_motion_core_positive_cases": 40,
    "min_contact_change_positive_cases": 80,
    "min_ligand_evidence_positive_cases": 200,
    "max_leakage_hits": 0,
    "min_unique_uniprot_or_cluster": 80,
    "min_unique_ligands": 120,
    "required_splits": ["train", "valid", "test"],
    "required_tiers": ["easy", "medium", "hard", "extreme"],
}


@dataclass
class GateCheck:
    name: str
    passed: bool
    observed: Any
    required: Any

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "passed": self.passed,
            "observed": self.observed,
            "required": self.required,
        }


def load_json(path: Path | None) -> dict[str, Any]:
    if path is None or not path.exists():
        return {}
    return json.loads(path.read_text())


def load_spec(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {"release_gate": DEFAULT_GATE}
    try:
        from omegaconf import OmegaConf
    except Exception as exc:
        raise RuntimeError("omegaconf is required to read benchmark_spec.yaml") from exc
    config = OmegaConf.to_container(OmegaConf.load(path), resolve=True)
    if not isinstance(config, dict):
        raise ValueError(f"Invalid spec: {path}")
    return config


def metric(record: IFBenchCase, key: str, default: float = 0.0) -> float:
    value = record.metrics.get(key)
    if value is None:
        return default
    return float(value)


def is_hard_extreme(record: IFBenchCase) -> bool:
    return record.flexibility_tier in {"hard", "extreme"}


def summarize(records: list[IFBenchCase], leakage_report: dict[str, Any]) -> dict[str, Any]:
    split_counts = Counter(str(record.split) for record in records)
    tier_counts = Counter(str(record.flexibility_tier) for record in records)
    split_tier_counts = Counter(f"{record.split}/{record.flexibility_tier}" for record in records)
    hard_extreme = [record for record in records if is_hard_extreme(record)]
    test_records = [record for record in records if record.split == "test"]
    test_hard_extreme = [record for record in test_records if is_hard_extreme(record)]
    motion_core_positive = [record for record in records if metric(record, "motion_core_residues") > 0]
    contact_change_positive = [record for record in records if metric(record, "contact_change_residues") > 0]
    ligand_evidence_positive = [record for record in records if metric(record, "ligand_evidence_positive_count") > 0]
    clusters = {
        record.cluster_key or record.uniprot_id or f"{record.holo_pdb_id}:{record.holo_chain_id}"
        for record in records
    }
    ligands = {record.ligand_id for record in records if record.ligand_id}
    return {
        "num_final_cases": len(records),
        "num_hard_extreme_cases": len(hard_extreme),
        "num_test_cases": len(test_records),
        "num_test_hard_extreme_cases": len(test_hard_extreme),
        "num_motion_core_positive_cases": len(motion_core_positive),
        "num_contact_change_positive_cases": len(contact_change_positive),
        "num_ligand_evidence_positive_cases": len(ligand_evidence_positive),
        "num_unique_uniprot_or_cluster": len(clusters),
        "num_unique_ligands": len(ligands),
        "num_leakage_hits": int(leakage_report.get("num_leakage_hits", 0)),
        "split_counts": dict(split_counts),
        "tier_counts": dict(tier_counts),
        "split_tier_counts": dict(split_tier_counts),
    }


def check_min(name: str, observed: int, required: int) -> GateCheck:
    return GateCheck(name=name, passed=observed >= required, observed=observed, required=f">= {required}")


def check_max(name: str, observed: int, required: int) -> GateCheck:
    return GateCheck(name=name, passed=observed <= required, observed=observed, required=f"<= {required}")


def run_checks(summary: dict[str, Any], gate: dict[str, Any]) -> list[GateCheck]:
    checks = [
        check_min("final_cases", summary["num_final_cases"], int(gate["min_final_cases"])),
        check_min("hard_extreme_cases", summary["num_hard_extreme_cases"], int(gate["min_hard_extreme_cases"])),
        check_min("test_cases", summary["num_test_cases"], int(gate["min_test_cases"])),
        check_min("test_hard_extreme_cases", summary["num_test_hard_extreme_cases"], int(gate["min_test_hard_extreme_cases"])),
        check_min("motion_core_positive_cases", summary["num_motion_core_positive_cases"], int(gate["min_motion_core_positive_cases"])),
        check_min("contact_change_positive_cases", summary["num_contact_change_positive_cases"], int(gate["min_contact_change_positive_cases"])),
        check_min("ligand_evidence_positive_cases", summary["num_ligand_evidence_positive_cases"], int(gate["min_ligand_evidence_positive_cases"])),
        check_max("leakage_hits", summary["num_leakage_hits"], int(gate["max_leakage_hits"])),
        check_min("unique_uniprot_or_cluster", summary["num_unique_uniprot_or_cluster"], int(gate["min_unique_uniprot_or_cluster"])),
        check_min("unique_ligands", summary["num_unique_ligands"], int(gate["min_unique_ligands"])),
    ]
    split_counts = summary["split_counts"]
    tier_counts = summary["tier_counts"]
    for split in gate.get("required_splits", []):
        checks.append(GateCheck(f"has_split_{split}", split_counts.get(split, 0) > 0, split_counts.get(split, 0), "> 0"))
    for tier in gate.get("required_tiers", []):
        checks.append(GateCheck(f"has_tier_{tier}", tier_counts.get(tier, 0) > 0, tier_counts.get(tier, 0), "> 0"))
    return checks


def write_markdown(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# IFBench Gate Report\n"]
    lines.append(f"Status: **{report['status']}**\n")
    lines.append("## Summary\n")
    for key, value in report["summary"].items():
        if isinstance(value, dict):
            continue
        lines.append(f"- `{key}`: {value}")
    lines.append("\n## Checks\n")
    lines.append("| check | observed | required | status |")
    lines.append("| --- | --- | --- | --- |")
    for check in report["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        lines.append(f"| {check['name']} | {check['observed']} | {check['required']} | {status} |")
    lines.append("\n## Counts\n")
    lines.append("```json")
    lines.append(json.dumps(
        {
            "split_counts": report["summary"]["split_counts"],
            "tier_counts": report["summary"]["tier_counts"],
            "split_tier_counts": report["summary"]["split_tier_counts"],
        },
        indent=2,
        sort_keys=True,
    ))
    lines.append("```")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="Final manifest JSONL after leakage and stratification.")
    parser.add_argument("--spec", type=Path, default=Path(__file__).with_name("benchmark_spec.yaml"))
    parser.add_argument("--leakage-report", type=Path)
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--out-md", type=Path)
    parser.add_argument("--allow-fail", action="store_true", help="Exit 0 even if the gate fails.")
    args = parser.parse_args()

    records = [record for record in read_jsonl(args.manifest) if record.accepted]
    spec = load_spec(args.spec)
    gate = dict(DEFAULT_GATE)
    gate.update(spec.get("release_gate", {}))
    leakage_report = load_json(args.leakage_report)
    summary = summarize(records, leakage_report)
    checks = run_checks(summary, gate)
    passed = all(check.passed for check in checks)
    report = {
        "status": "PASS" if passed else "FAIL",
        "summary": summary,
        "checks": [check.to_dict() for check in checks],
        "spec": str(args.spec) if args.spec else None,
    }
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if args.out_md:
        write_markdown(args.out_md, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    if not passed and not args.allow_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
