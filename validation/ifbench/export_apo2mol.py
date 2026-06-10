#!/usr/bin/env python3
"""Export accepted IFBench cases into an Apo2Mol-compatible data layout."""

from __future__ import annotations

import argparse
import json
import pickle
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import IFBenchCase, read_jsonl
from ifbench.structure_qc import parse_structure, resolve_path

try:
    from Bio.PDB import PDBIO, Select
except Exception:  # pragma: no cover
    PDBIO = None
    Select = object


class PocketSelect(Select):
    def __init__(self, chain_id: str, residue_numbers: set[int]) -> None:
        self.chain_id = chain_id
        self.residue_numbers = residue_numbers

    def accept_chain(self, chain: Any) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue: Any) -> bool:
        hetflag, resseq, _ = residue.id
        return hetflag == " " and int(resseq) in self.residue_numbers


def safe_name(value: str) -> str:
    keep = []
    for char in value:
        if char.isalnum() or char in {"-", "_", "."}:
            keep.append(char)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "case"


def pocket_residue_numbers(record: IFBenchCase, key: str) -> set[int]:
    values = set()
    for item in record.labels.get("pocket_residues", []):
        if key in item:
            values.add(int(item[key]))
    return values


def save_pocket_pdb(source_file: Path, chain_id: str, residue_numbers: set[int], out_file: Path) -> None:
    if PDBIO is None:
        raise RuntimeError("Biopython is required for Apo2Mol export")
    structure = parse_structure(source_file)
    io = PDBIO()
    io.set_structure(structure)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out_file), PocketSelect(chain_id, residue_numbers))


def copy_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def export_record(record: IFBenchCase, raw_root: Path, out_data_folder: Path) -> tuple[Any, ...]:
    case_dir_name = safe_name(f"{record.case_id}__{record.holo_pdb_id}_{record.holo_chain_id}__{record.apo_pdb_id}_{record.apo_chain_id}__{record.ligand_id}")
    case_dir = out_data_folder / case_dir_name

    holo_path = resolve_path(raw_root, record.holo_structure_file)
    apo_path = resolve_path(raw_root, record.apo_structure_file)
    ligand_path = resolve_path(raw_root, record.ligand_file)

    holo_residues = pocket_residue_numbers(record, "holo_resseq")
    apo_residues = pocket_residue_numbers(record, "apo_resseq")
    if not holo_residues or not apo_residues:
        raise ValueError(f"{record.case_id} has no pocket residue labels; run build_manifest first")

    holo_pocket = case_dir / "receptor_holo_pocket_ifbench.pdb"
    apo_pocket = case_dir / "receptor_apo_pocket_ifbench.pdb"
    ligand_out = case_dir / f"{safe_name(record.ligand_id or 'ligand')}.sdf"
    holo_full = case_dir / f"holo_full{holo_path.suffix or '.pdb'}"
    apo_full = case_dir / f"apo_full{apo_path.suffix or '.pdb'}"

    save_pocket_pdb(holo_path, record.holo_chain_id, holo_residues, holo_pocket)
    save_pocket_pdb(apo_path, record.apo_chain_id, apo_residues, apo_pocket)
    copy_file(ligand_path, ligand_out)
    copy_file(holo_path, holo_full)
    copy_file(apo_path, apo_full)

    return (
        str(holo_pocket.relative_to(out_data_folder)),
        str(apo_pocket.relative_to(out_data_folder)),
        str(ligand_out.relative_to(out_data_folder)),
        str(holo_full.relative_to(out_data_folder)),
        str(apo_full.relative_to(out_data_folder)),
        float(record.metrics.get("sequence_identity", 0.0)) * 100.0,
        float(record.metrics.get("mapping_coverage", 0.0)) * 100.0,
        record.holo_release_date or "",
        float(record.metrics.get("pocket_ca_rmsd") or 0.0),
    )


def write_split_pt(path: Path, splits: dict[str, list[int]]) -> None:
    try:
        import torch
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("torch is required to write Apo2Mol split .pt files") from exc
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(splits, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="Final IFBench manifest JSONL.")
    parser.add_argument("--raw-root", default=".", type=Path, help="Root directory for relative raw paths.")
    parser.add_argument("--out-data-folder", required=True, type=Path, help="Apo2Mol-style data_folder output.")
    parser.add_argument("--out-index", required=True, type=Path, help="Output selected_index_ifbench.pkl.")
    parser.add_argument("--out-split-pt", type=Path, help="Optional torch split file.")
    parser.add_argument("--out-split-json", type=Path, help="Optional JSON split file.")
    parser.add_argument("--default-split", choices=["train", "valid", "test"], default="test")
    args = parser.parse_args()

    records = [record for record in read_jsonl(args.manifest) if record.accepted]
    index: list[tuple[Any, ...]] = []
    splits: dict[str, list[int]] = defaultdict(list)
    for record in records:
        item = export_record(record, args.raw_root, args.out_data_folder)
        position = len(index)
        index.append(item)
        splits[record.split or args.default_split].append(position)

    args.out_index.parent.mkdir(parents=True, exist_ok=True)
    with args.out_index.open("wb") as handle:
        pickle.dump(index, handle)

    split_dict = {key: list(splits.get(key, [])) for key in ["train", "valid", "test"]}
    if args.out_split_pt:
        write_split_pt(args.out_split_pt, split_dict)
    if args.out_split_json:
        args.out_split_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_split_json.write_text(json.dumps(split_dict, indent=2, sort_keys=True) + "\n")

    summary = {
        "num_exported": len(index),
        "out_data_folder": str(args.out_data_folder),
        "out_index": str(args.out_index),
        "splits": {key: len(value) for key, value in split_dict.items()},
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
