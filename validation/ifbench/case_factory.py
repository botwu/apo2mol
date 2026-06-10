#!/usr/bin/env python3
"""Run the source-agnostic IFBench case factory.

Input is a candidate CSV plus raw files. Candidate rows can come from RCSB,
Binding MOAD/PDBbind/PLINDER adapters, hand curation, or synthetic generators.
This script applies the same QC, leakage filtering, stratification, reporting,
export, and release gate to all sources.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


THIS_DIR = Path(__file__).resolve().parent


def run_step(name: str, command: list[str]) -> dict[str, Any]:
    print(f"\n[IFBench] {name}", file=sys.stderr)
    print(" ".join(command), file=sys.stderr)
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    return {
        "name": name,
        "command": command,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def require_success(step: dict[str, Any]) -> None:
    if step["returncode"] != 0:
        raise RuntimeError(f"Step failed: {step['name']} (returncode={step['returncode']})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path, help="Candidate CSV.")
    parser.add_argument("--raw-root", default=".", type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--apo2mol-index", type=Path, help="Optional Apo2Mol selected_index blocklist.")
    parser.add_argument("--ligand-evidence-csv", action="append", type=Path, default=[], help="Curated ligand evidence CSV, repeatable.")
    parser.add_argument("--require-ligand-evidence", action="store_true")
    parser.add_argument("--allow-negative-ligand-evidence", action="store_true")
    parser.add_argument("--spec", type=Path, default=THIS_DIR / "benchmark_spec.yaml")
    parser.add_argument("--skip-export", action="store_true")
    parser.add_argument("--allow-gate-fail", action="store_true", help="Keep artifacts even when release gate fails.")
    args = parser.parse_args()

    args.run_dir.mkdir(parents=True, exist_ok=True)
    accepted = args.run_dir / "manifest.accepted.jsonl"
    rejected = args.run_dir / "manifest.rejected.jsonl"
    qc_csv = args.run_dir / "manifest.qc.csv"
    qc_summary = args.run_dir / "qc_summary.json"
    no_leakage = args.run_dir / "manifest.no_leakage.jsonl"
    leakage_report = args.run_dir / "leakage_report.json"
    final = args.run_dir / "manifest.final.jsonl"
    final_csv = args.run_dir / "manifest.final.csv"
    splits_json = args.run_dir / "splits.json"
    split_summary = args.run_dir / "split_summary.json"
    report_md = args.run_dir / "ifbench_report.md"
    review_csv = args.run_dir / "review_queue.csv"
    review_md = args.run_dir / "review_queue.md"
    gate_json = args.run_dir / "gate_report.json"
    gate_md = args.run_dir / "gate_report.md"
    data_folder = args.run_dir / "data_folder"
    index_pkl = args.run_dir / "selected_index_ifbench.pkl"
    split_pt = args.run_dir / "split_ifbench.pt"
    split_export_json = args.run_dir / "split_ifbench.json"

    steps: list[dict[str, Any]] = []
    python = sys.executable

    build_command = [
        python,
        str(THIS_DIR / "build_manifest.py"),
        "--candidates",
        str(args.candidates),
        "--raw-root",
        str(args.raw_root),
        "--out",
        str(accepted),
        "--rejects",
        str(rejected),
        "--flat-csv",
        str(qc_csv),
        "--summary-json",
        str(qc_summary),
    ]
    for evidence_csv in args.ligand_evidence_csv:
        build_command.extend(["--ligand-evidence-csv", str(evidence_csv)])
    if args.require_ligand_evidence:
        build_command.append("--require-ligand-evidence")
    if args.allow_negative_ligand_evidence:
        build_command.append("--allow-negative-ligand-evidence")

    steps.append(run_step(
        "build_manifest",
        build_command,
    ))
    require_success(steps[-1])

    if args.apo2mol_index:
        steps.append(run_step(
            "leakage_check",
            [
                python,
                str(THIS_DIR / "leakage_check.py"),
                "--manifest",
                str(accepted),
                "--apo2mol-index",
                str(args.apo2mol_index),
                "--out",
                str(no_leakage),
                "--report-json",
                str(leakage_report),
            ],
        ))
        require_success(steps[-1])
    else:
        shutil.copy2(accepted, no_leakage)
        leakage_report.write_text(json.dumps({"num_cases": 0, "num_leakage_hits": 0, "hits": []}, indent=2) + "\n")

    steps.append(run_step(
        "stratify_and_split",
        [
            python,
            str(THIS_DIR / "stratify_and_split.py"),
            "--manifest",
            str(no_leakage),
            "--out",
            str(final),
            "--flat-csv",
            str(final_csv),
            "--splits-json",
            str(splits_json),
            "--summary-json",
            str(split_summary),
        ],
    ))
    require_success(steps[-1])

    steps.append(run_step(
        "make_report",
        [
            python,
            str(THIS_DIR / "make_report.py"),
            "--accepted",
            str(accepted),
            "--rejected",
            str(rejected),
            "--final",
            str(final),
            "--qc-summary",
            str(qc_summary),
            "--leakage-report",
            str(leakage_report),
            "--split-summary",
            str(split_summary),
            "--out",
            str(report_md),
        ],
    ))
    require_success(steps[-1])

    steps.append(run_step(
        "make_review_queue",
        [
            python,
            str(THIS_DIR / "make_review_queue.py"),
            "--accepted",
            str(accepted),
            "--rejected",
            str(rejected),
            "--final",
            str(final),
            "--out-csv",
            str(review_csv),
            "--out-md",
            str(review_md),
        ],
    ))
    require_success(steps[-1])

    if not args.skip_export:
        steps.append(run_step(
            "export_apo2mol",
            [
                python,
                str(THIS_DIR / "export_apo2mol.py"),
                "--manifest",
                str(final),
                "--raw-root",
                str(args.raw_root),
                "--out-data-folder",
                str(data_folder),
                "--out-index",
                str(index_pkl),
                "--out-split-pt",
                str(split_pt),
                "--out-split-json",
                str(split_export_json),
            ],
        ))
        require_success(steps[-1])

    gate_command = [
        python,
        str(THIS_DIR / "gate_manifest.py"),
        "--manifest",
        str(final),
        "--spec",
        str(args.spec),
        "--leakage-report",
        str(leakage_report),
        "--out-json",
        str(gate_json),
        "--out-md",
        str(gate_md),
    ]
    if args.allow_gate_fail:
        gate_command.append("--allow-fail")
    steps.append(run_step("gate_manifest", gate_command))
    if not args.allow_gate_fail:
        require_success(steps[-1])

    summary = {
        "run_dir": str(args.run_dir),
        "candidates": str(args.candidates),
        "raw_root": str(args.raw_root),
        "ligand_evidence_csv": [str(path) for path in args.ligand_evidence_csv],
        "require_ligand_evidence": args.require_ligand_evidence,
        "steps": [{"name": step["name"], "returncode": step["returncode"]} for step in steps],
        "artifacts": {
            "accepted": str(accepted),
            "rejected": str(rejected),
            "final": str(final),
            "report_md": str(report_md),
            "review_csv": str(review_csv),
            "review_md": str(review_md),
            "gate_json": str(gate_json),
            "gate_md": str(gate_md),
            "data_folder": str(data_folder) if not args.skip_export else None,
            "selected_index": str(index_pkl) if not args.skip_export else None,
            "split_pt": str(split_pt) if not args.skip_export else None,
        },
    }
    summary_path = args.run_dir / "case_factory_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
