#!/usr/bin/env python3
"""Shared schema helpers for IFBench manifests.

The public benchmark should be auditable without importing Apo2Mol internals.
Records are stored as JSONL so each case keeps nested metrics, labels, and QC
issues in one place.
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable


REJECT = "reject"
WARNING = "warning"


@dataclass
class QualityIssue:
    code: str
    severity: str
    message: str


@dataclass
class IFBenchCase:
    case_id: str
    source: str
    holo_pdb_id: str
    apo_pdb_id: str
    holo_chain_id: str
    apo_chain_id: str
    ligand_id: str
    holo_structure_file: str
    apo_structure_file: str
    ligand_file: str
    uniprot_id: str | None = None
    cluster_key: str | None = None
    holo_release_date: str | None = None
    apo_release_date: str | None = None
    holo_resolution: float | None = None
    apo_resolution: float | None = None
    holo_method: str | None = None
    apo_method: str | None = None
    accepted: bool = False
    split: str | None = None
    flexibility_tier: str | None = None
    metrics: dict[str, Any] = field(default_factory=dict)
    labels: dict[str, Any] = field(default_factory=dict)
    evidence: dict[str, Any] = field(default_factory=dict)
    issues: list[QualityIssue] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, item: dict[str, Any]) -> "IFBenchCase":
        issues = [
            issue if isinstance(issue, QualityIssue) else QualityIssue(**issue)
            for issue in item.get("issues", [])
        ]
        item = dict(item)
        item["issues"] = issues
        return cls(**item)


def stable_case_id(parts: Iterable[str]) -> str:
    joined = "|".join(str(part).strip().lower() for part in parts)
    return hashlib.sha1(joined.encode("utf-8")).hexdigest()[:12]


def normalize_pdb_id(value: str | None) -> str:
    return (value or "").strip().lower()


def parse_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    value = str(value).strip()
    if not value or value.lower() in {"na", "nan", "none", "null"}:
        return None
    return float(value)


def read_jsonl(path: Path) -> list[IFBenchCase]:
    records: list[IFBenchCase] = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(IFBenchCase.from_dict(json.loads(line)))
            except Exception as exc:
                raise ValueError(f"Invalid JSONL record at {path}:{line_no}: {exc}") from exc
    return records


def write_jsonl(path: Path, records: Iterable[IFBenchCase]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record.to_dict(), sort_keys=True) + "\n")


def flatten_case(record: IFBenchCase) -> dict[str, Any]:
    flat = {
        "case_id": record.case_id,
        "source": record.source,
        "accepted": record.accepted,
        "split": record.split,
        "flexibility_tier": record.flexibility_tier,
        "holo_pdb_id": record.holo_pdb_id,
        "apo_pdb_id": record.apo_pdb_id,
        "holo_chain_id": record.holo_chain_id,
        "apo_chain_id": record.apo_chain_id,
        "ligand_id": record.ligand_id,
        "uniprot_id": record.uniprot_id,
        "cluster_key": record.cluster_key,
        "num_issues": len(record.issues),
        "reject_issues": sum(1 for issue in record.issues if issue.severity == REJECT),
    }
    for key, value in sorted(record.metrics.items()):
        if isinstance(value, (str, int, float, bool)) or value is None:
            flat[f"metric.{key}"] = value
    return flat


def write_flat_csv(path: Path, records: Iterable[IFBenchCase]) -> None:
    rows = [flatten_case(record) for record in records]
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
