#!/usr/bin/env python3
"""Synthetic smoke test for the IFBench QC path.

This does not validate scientific thresholds. It only verifies that the parser,
ligand QC, residue mapping, pocket metrics, and active-set labels can accept a
well-formed apo/holo/ligand triple.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.structure_qc import QualityConfig, evaluate_candidate


def pdb_atom_line(serial: int, atom_name: str, resname: str, chain: str, resseq: int, xyz: np.ndarray, element: str) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:^4s} {resname:>3s} {chain:1s}{resseq:4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          {element:>2s}\n"
    )


def write_synthetic_pdb(path: Path, holo: bool, n_residues: int = 50) -> None:
    serial = 1
    lines: list[str] = []
    for resseq in range(1, n_residues + 1):
        base = np.asarray([resseq * 1.5, 0.0, 0.0], dtype=np.float64)
        if holo and 21 <= resseq <= 27:
            base = base + np.asarray([0.0, 1.8, 0.0], dtype=np.float64)
        atoms = [
            ("N", base + [-0.55, 0.00, 0.00], "N"),
            ("CA", base + [0.00, 0.00, 0.00], "C"),
            ("C", base + [0.55, 0.00, 0.00], "C"),
            ("O", base + [0.85, 0.45, 0.00], "O"),
            ("CB", base + [0.00, 1.05, 0.00], "C"),
        ]
        for atom_name, xyz, element in atoms:
            lines.append(pdb_atom_line(serial, atom_name, "ALA", "A", resseq, np.asarray(xyz), element))
            serial += 1
    lines.append("TER\nEND\n")
    path.write_text("".join(lines))


def write_synthetic_ligand(path: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CCOc1ccccc1"))
    status = AllChem.EmbedMolecule(mol, randomSeed=20270609)
    if status != 0:
        raise RuntimeError("RDKit failed to embed synthetic ligand")
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    coords = np.asarray([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    center = coords.mean(axis=0)
    target = np.asarray([24.0 * 1.5, 2.2, 0.0], dtype=np.float64)
    shift = target - center
    for idx in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(idx)
        conf.SetAtomPosition(idx, (pos.x + shift[0], pos.y + shift[1], pos.z + shift[2]))
    Chem.MolToMolFile(mol, str(path))


def write_candidate_csv(path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "source",
                "holo_pdb_id",
                "apo_pdb_id",
                "holo_chain_id",
                "apo_chain_id",
                "ligand_id",
                "holo_structure_file",
                "apo_structure_file",
                "ligand_file",
                "uniprot_id",
                "cluster_key",
                "holo_resolution",
                "apo_resolution",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "case_id": "synthetic_ifbench_smoke",
                "source": "synthetic",
                "holo_pdb_id": "9zzz",
                "apo_pdb_id": "8yyy",
                "holo_chain_id": "A",
                "apo_chain_id": "A",
                "ligand_id": "SYN",
                "holo_structure_file": "holo.pdb",
                "apo_structure_file": "apo.pdb",
                "ligand_file": "ligand.sdf",
                "uniprot_id": "SYNTHETIC",
                "cluster_key": "SYNTHETIC_A",
                "holo_resolution": "1.8",
                "apo_resolution": "1.8",
            }
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    write_synthetic_pdb(args.work_dir / "apo.pdb", holo=False)
    write_synthetic_pdb(args.work_dir / "holo.pdb", holo=True)
    write_synthetic_ligand(args.work_dir / "ligand.sdf")
    write_candidate_csv(args.work_dir / "candidates.csv")

    row = next(csv.DictReader((args.work_dir / "candidates.csv").open()))
    record = evaluate_candidate(row, args.work_dir, QualityConfig())
    summary = {
        "accepted": record.accepted,
        "case_id": record.case_id,
        "issues": [issue.__dict__ for issue in record.issues],
        "pocket_residues": record.metrics.get("pocket_residues"),
        "pocket_ca_rmsd": record.metrics.get("pocket_ca_rmsd"),
        "motion_core_residues": record.metrics.get("motion_core_residues"),
        "contact_change_residues": record.metrics.get("contact_change_residues"),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not record.accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
