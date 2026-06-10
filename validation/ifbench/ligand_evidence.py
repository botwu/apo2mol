#!/usr/bin/env python3
"""Curated ligand relevance evidence for IFBench candidates.

This module intentionally uses a simple CSV adapter. Binding MOAD, BioLiP,
PDBbind, PLINDER, and manual review tables can all be normalized to this schema
without making the benchmark depend on one upstream database layout.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


POSITIVE_DECISIONS = {"positive", "accept", "accepted", "binder", "biological", "relevant", "yes", "true", "1"}
NEGATIVE_DECISIONS = {"negative", "reject", "rejected", "artifact", "additive", "buffer", "irrelevant", "no", "false", "0"}


@dataclass
class LigandEvidence:
    pdb_id: str
    ligand_id: str
    source: str
    evidence_type: str
    confidence: float | None = None
    decision: str | None = None
    chain_id: str | None = None
    label_asym_id: str | None = None
    notes: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def is_negative(self) -> bool:
        return normalize_token(self.decision) in NEGATIVE_DECISIONS

    @property
    def is_positive(self) -> bool:
        decision = normalize_token(self.decision)
        return not decision or decision in POSITIVE_DECISIONS


def normalize_token(value: Any) -> str:
    return str(value or "").strip().lower()


def normalize_pdb_id(value: Any) -> str:
    return normalize_token(value)


def normalize_ligand_id(value: Any) -> str:
    return str(value or "").strip().upper()


def first_present(row: dict[str, Any], keys: list[str]) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def parse_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return float(text)


def evidence_from_row(row: dict[str, Any]) -> LigandEvidence:
    return LigandEvidence(
        pdb_id=normalize_pdb_id(first_present(row, ["pdb_id", "holo_pdb_id", "entry_id", "rcsb_id"])),
        ligand_id=normalize_ligand_id(first_present(row, ["ligand_id", "comp_id", "ccd_id", "het_id"])),
        source=first_present(row, ["source", "database", "curation_source"]) or "curated",
        evidence_type=first_present(row, ["evidence_type", "type", "annotation", "label"]) or "ligand_relevance",
        confidence=parse_float(first_present(row, ["confidence", "score", "probability"])),
        decision=first_present(row, ["decision", "is_biological", "accepted", "label_decision"]) or None,
        chain_id=first_present(row, ["chain_id", "holo_chain_id", "auth_chain_id"]) or None,
        label_asym_id=first_present(row, ["label_asym_id", "ligand_label_asym_id"]) or None,
        notes=first_present(row, ["notes", "comment", "description"]) or None,
    )


class LigandEvidenceIndex:
    def __init__(self, evidence: list[LigandEvidence]) -> None:
        self.evidence = evidence
        self.by_pair: dict[tuple[str, str], list[LigandEvidence]] = {}
        for item in evidence:
            if item.pdb_id and item.ligand_id:
                self.by_pair.setdefault((item.pdb_id, item.ligand_id), []).append(item)

    @classmethod
    def from_csvs(cls, paths: list[Path]) -> "LigandEvidenceIndex":
        evidence: list[LigandEvidence] = []
        for path in paths:
            with path.open(newline="") as handle:
                for row in csv.DictReader(handle):
                    item = evidence_from_row(row)
                    if item.pdb_id and item.ligand_id:
                        evidence.append(item)
        return cls(evidence)

    def match_candidate(self, row: dict[str, Any]) -> list[LigandEvidence]:
        pdb_id = normalize_pdb_id(first_present(row, ["holo_pdb_id", "pdb_id", "entry_id", "rcsb_id"]))
        ligand_id = normalize_ligand_id(first_present(row, ["ligand_id", "comp_id", "ccd_id", "het_id"]))
        chain_id = first_present(row, ["holo_chain_id", "chain_id", "holo_auth_chain_id", "auth_chain_id"])
        label_asym_id = first_present(row, ["ligand_label_asym_id", "label_asym_id"])
        matches = []
        for item in self.by_pair.get((pdb_id, ligand_id), []):
            if item.chain_id and chain_id and item.chain_id != chain_id:
                continue
            if item.label_asym_id and label_asym_id and item.label_asym_id != label_asym_id:
                continue
            matches.append(item)
        return matches


def summarize_evidence(matches: list[LigandEvidence]) -> dict[str, Any]:
    return {
        "count": len(matches),
        "sources": sorted({item.source for item in matches}),
        "positive_count": sum(1 for item in matches if item.is_positive),
        "negative_count": sum(1 for item in matches if item.is_negative),
        "items": [item.to_dict() for item in matches],
    }


def evidence_summary_json(matches: list[LigandEvidence]) -> str:
    return json.dumps(summarize_evidence(matches), sort_keys=True)
