#!/usr/bin/env python3
"""Build an IFBench QC manifest from candidate apo/holo/ligand triples."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import REJECT, write_flat_csv, write_jsonl
from ifbench.ligand_evidence import LigandEvidenceIndex
from ifbench.structure_qc import QualityConfig, evaluate_candidate


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path, help="Candidate CSV with apo/holo/ligand file paths.")
    parser.add_argument("--raw-root", default=".", type=Path, help="Root directory for relative candidate file paths.")
    parser.add_argument("--out", required=True, type=Path, help="Output JSONL manifest for accepted cases.")
    parser.add_argument("--rejects", type=Path, help="Output JSONL manifest for rejected cases.")
    parser.add_argument("--flat-csv", type=Path, help="Optional flattened CSV summary for accepted and rejected cases.")
    parser.add_argument("--summary-json", type=Path, help="Optional QC summary JSON.")
    parser.add_argument("--max-resolution", type=float, default=2.5)
    parser.add_argument("--min-seq-identity", type=float, default=0.95)
    parser.add_argument("--min-mapping-coverage", type=float, default=0.85)
    parser.add_argument("--min-common-residues", type=int, default=40)
    parser.add_argument("--pocket-cutoff", type=float, default=8.0)
    parser.add_argument("--contact-cutoff", type=float, default=4.5)
    parser.add_argument("--motion-core-cutoff", type=float, default=1.5)
    parser.add_argument("--min-pocket-residues", type=int, default=8)
    parser.add_argument("--ligand-evidence-csv", action="append", type=Path, default=[], help="Curated ligand evidence CSV, repeatable.")
    parser.add_argument("--require-ligand-evidence", action="store_true", help="Reject candidates without positive curated ligand evidence.")
    parser.add_argument("--allow-negative-ligand-evidence", action="store_true", help="Do not reject candidates marked negative by curated evidence.")
    args = parser.parse_args()

    config = QualityConfig(
        max_resolution_a=args.max_resolution,
        min_sequence_identity=args.min_seq_identity,
        min_mapping_coverage=args.min_mapping_coverage,
        min_common_residues=args.min_common_residues,
        pocket_cutoff_a=args.pocket_cutoff,
        contact_cutoff_a=args.contact_cutoff,
        motion_core_cutoff_a=args.motion_core_cutoff,
        min_pocket_residues=args.min_pocket_residues,
        require_ligand_evidence=args.require_ligand_evidence,
        reject_negative_ligand_evidence=not args.allow_negative_ligand_evidence,
    )
    ligand_evidence_index = LigandEvidenceIndex.from_csvs(args.ligand_evidence_csv) if args.ligand_evidence_csv else None

    records = [
        evaluate_candidate(row, args.raw_root, config, ligand_evidence_index=ligand_evidence_index)
        for row in read_candidates(args.candidates)
    ]
    accepted = [record for record in records if record.accepted]
    rejected = [record for record in records if not record.accepted]
    write_jsonl(args.out, accepted)
    if args.rejects:
        write_jsonl(args.rejects, rejected)
    if args.flat_csv:
        write_flat_csv(args.flat_csv, records)

    issue_counts = Counter()
    reject_counts = Counter()
    for record in records:
        for issue in record.issues:
            issue_counts[issue.code] += 1
            if issue.severity == REJECT:
                reject_counts[issue.code] += 1
    summary = {
        "num_candidates": len(records),
        "num_accepted": len(accepted),
        "num_rejected": len(rejected),
        "issue_counts": dict(issue_counts),
        "reject_counts": dict(reject_counts),
    }
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
