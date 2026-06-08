"""Lightweight ligand-quality diagnostics for Apo2Mol validation runs.

This complements eval_split.py when docking is unavailable. It reconstructs all
generated samples and reports chemistry metrics for complete molecules without
calling qvina/vina.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import torch
from rdkit import Chem
from rdkit import RDLogger
from tqdm.auto import tqdm

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utils import reconstruct, transforms
from utils.evaluation import analyze, scoring_func


def _stats(values: list[float]) -> dict[str, Any]:
    finite = [float(v) for v in values if np.isfinite(v)]
    if not finite:
        return {"n": 0}
    arr = np.asarray(finite, dtype=float)
    return {
        "n": int(arr.size),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def _ring_ratios(rings: list[Counter[int]]) -> dict[str, float | None]:
    if not rings:
        return {str(size): None for size in range(3, 10)}
    return {
        str(size): float(sum(1 for ring_counter in rings if size in ring_counter) / len(rings))
        for size in range(3, 10)
    }


def _to_jsonable(obj: Any) -> Any:
    if isinstance(obj, Counter):
        return {str(k): int(v) for k, v in obj.items()}
    if isinstance(obj, dict):
        return {str(k): _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def analyze_sample_dir(
    sample_dir: Path,
    eval_step: int,
    atom_enc_mode: str,
) -> dict[str, Any]:
    result_files = sorted(
        glob.glob(str(sample_dir / "result_*.pt")),
        key=lambda path: int(Path(path).stem.split("_")[-1]),
    )

    generated = 0
    stable_mols = 0
    stable_atoms = 0
    total_atoms = 0
    recon_success = 0
    complete = 0
    chem_success = 0

    qed_values: list[float] = []
    sa_values: list[float] = []
    logp_values: list[float] = []
    lipinski_values: list[float] = []
    ring_counters: list[Counter[int]] = []
    atom_types = Counter()
    failures = Counter()

    for result_file in tqdm(result_files, desc=f"{sample_dir.parent.name}", leave=False):
        result = torch.load(result_file, map_location="cpu", weights_only=False)
        ligand_pos_traj = result.get("pred_ligand_ligand_pos_traj")
        ligand_v_traj = result.get("pred_ligand_v_traj")
        if ligand_pos_traj is None or ligand_v_traj is None:
            failures["missing_ligand_traj"] += 1
            continue

        for pred_pos_traj, pred_v_traj in zip(ligand_pos_traj, ligand_v_traj):
            generated += 1
            pred_pos = pred_pos_traj[eval_step]
            pred_v = pred_v_traj[eval_step]
            pred_atom_type = transforms.get_atomic_number_from_index(
                pred_v,
                mode=atom_enc_mode,
            )
            atom_types.update(int(v) for v in pred_atom_type)

            mol_stable, atom_stable, num_atoms = analyze.check_stability(pred_pos, pred_atom_type)
            stable_mols += int(mol_stable)
            stable_atoms += int(atom_stable)
            total_atoms += int(num_atoms)

            try:
                pred_aromatic = transforms.is_aromatic_from_index(pred_v, mode=atom_enc_mode)
                mol = reconstruct.reconstruct_from_generated(pred_pos, pred_atom_type, pred_aromatic)
                recon_success += 1
            except reconstruct.MolReconsError:
                failures["reconstruct"] += 1
                continue

            smiles = Chem.MolToSmiles(mol)
            if "." in smiles:
                failures["fragmented_smiles"] += 1
                continue
            complete += 1

            try:
                chem = scoring_func.get_chem(mol)
            except Exception:
                failures["chem"] += 1
                continue

            chem_success += 1
            qed_values.append(float(chem["qed"]))
            sa_values.append(float(chem["sa"]))
            logp_values.append(float(chem["logp"]))
            lipinski_values.append(float(chem["lipinski"]))
            ring_counters.append(chem["ring_size"])

    rates = {
        "mol_stable": float(stable_mols / generated) if generated else 0.0,
        "atm_stable": float(stable_atoms / total_atoms) if total_atoms else 0.0,
        "recon_success": float(recon_success / generated) if generated else 0.0,
        "complete": float(complete / generated) if generated else 0.0,
        "chem_success": float(chem_success / generated) if generated else 0.0,
        "chem_success_over_complete": float(chem_success / complete) if complete else 0.0,
    }

    return {
        "result_files": len(result_files),
        "generated_samples": generated,
        "counts": {
            "stable_mols": stable_mols,
            "stable_atoms": stable_atoms,
            "total_atoms": total_atoms,
            "recon_success": recon_success,
            "complete": complete,
            "chem_success": chem_success,
        },
        "rates": rates,
        "chemistry": {
            "qed": _stats(qed_values),
            "sa": _stats(sa_values),
            "logp": _stats(logp_values),
            "lipinski": _stats(lipinski_values),
            "ring_ratios": _ring_ratios(ring_counters),
        },
        "atom_types": {str(k): int(v) for k, v in atom_types.items()},
        "failures": {str(k): int(v) for k, v in failures.items()},
    }


def discover_arms(run_dir: Path) -> list[Path]:
    plan_path = run_dir / "experiment_plan.json"
    if plan_path.exists():
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        arms = []
        for arm in plan.get("arms", []):
            sample_path = arm.get("sample_path")
            if sample_path:
                path = Path(sample_path)
                if path.exists():
                    arms.append(path)
        if arms:
            return arms

    return sorted(path / "sampled_results" for path in run_dir.iterdir() if (path / "sampled_results").exists())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument("--eval-step", type=int, default=-1)
    parser.add_argument("--atom-enc-mode", default="add_aromatic")
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")

    run_dir = Path(args.run_dir)
    output_path = Path(args.output) if args.output else run_dir / "ligand_quality.json"
    arms = discover_arms(run_dir)
    if not arms:
        raise SystemExit(f"No sampled_results directories found under {run_dir}")

    payload = {
        "run_dir": str(run_dir),
        "eval_step": args.eval_step,
        "atom_enc_mode": args.atom_enc_mode,
        "arms": {},
    }

    for sample_dir in tqdm(arms, desc="arms"):
        arm_name = sample_dir.parent.name
        payload["arms"][arm_name] = analyze_sample_dir(
            sample_dir=sample_dir,
            eval_step=args.eval_step,
            atom_enc_mode=args.atom_enc_mode,
        )

    output_path.write_text(json.dumps(_to_jsonable(payload), indent=2), encoding="utf-8")
    print(f"Wrote {output_path}")
    print("arm\tgenerated\tcomplete\tqed_mean\tsa_mean\tlogp_mean\tlipinski_mean")
    for arm_name, result in payload["arms"].items():
        chemistry = result["chemistry"]
        print(
            "\t".join(
                [
                    arm_name,
                    str(result["generated_samples"]),
                    f"{result['rates']['complete']:.4f}",
                    f"{chemistry['qed'].get('mean', float('nan')):.4f}"
                    if chemistry["qed"].get("n", 0)
                    else "nan",
                    f"{chemistry['sa'].get('mean', float('nan')):.4f}"
                    if chemistry["sa"].get("n", 0)
                    else "nan",
                    f"{chemistry['logp'].get('mean', float('nan')):.4f}"
                    if chemistry["logp"].get("n", 0)
                    else "nan",
                    f"{chemistry['lipinski'].get('mean', float('nan')):.4f}"
                    if chemistry["lipinski"].get("n", 0)
                    else "nan",
                ]
            )
        )


if __name__ == "__main__":
    main()
