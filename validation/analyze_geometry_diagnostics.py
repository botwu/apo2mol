#!/usr/bin/env python3
"""Geometry diagnostics for sampled Apo2Mol result directories.

The standard eval path focuses on molecule reconstruction and aggregate protein
RMSD. This script adds cheap local checks that are useful for sparse pocket
update experiments:

* ligand-protein steric clashes;
* residue-level ligand contact recovery;
* local discontinuity around an approximated active set;
* protein-protein clashes between non-neighbor residues.

It intentionally works from saved ``result_*.pt`` files, so it can be run after
any existing validation run without re-sampling.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import torch


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def as_numpy(value: Any) -> np.ndarray:
    if isinstance(value, np.ndarray):
        return value.astype(np.float64, copy=False)
    if torch.is_tensor(value):
        return value.detach().cpu().numpy().astype(np.float64, copy=False)
    return np.asarray(value, dtype=np.float64)


def final_array(result: dict[str, Any], key: str) -> np.ndarray | None:
    value = result.get(key)
    if value is None:
        return None
    if isinstance(value, list):
        if not value:
            return None
        value = value[-1]
    return as_numpy(value)


def residue_centers(atom_pos: np.ndarray, residue_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unique_ids = np.unique(residue_ids)
    centers = np.zeros((len(unique_ids), 3), dtype=np.float64)
    for i, residue_id in enumerate(unique_ids):
        centers[i] = atom_pos[residue_ids == residue_id].mean(axis=0)
    return unique_ids, centers


def residue_contacts(
    ligand_pos: np.ndarray,
    protein_pos: np.ndarray,
    residue_ids: np.ndarray,
    cutoff: float,
) -> set[int]:
    if ligand_pos.size == 0 or protein_pos.size == 0:
        return set()
    distances = np.linalg.norm(ligand_pos[:, None, :] - protein_pos[None, :, :], axis=-1)
    atom_hits = np.where(distances.min(axis=0) < cutoff)[0]
    return {int(residue_ids[i]) for i in atom_hits}


def contact_metrics(pred_contacts: set[int], ref_contacts: set[int]) -> dict[str, float | int | None]:
    inter = pred_contacts & ref_contacts
    union = pred_contacts | ref_contacts
    precision = len(inter) / len(pred_contacts) if pred_contacts else None
    recall = len(inter) / len(ref_contacts) if ref_contacts else None
    jaccard = len(inter) / len(union) if union else None
    return {
        "pred_contacts": len(pred_contacts),
        "ref_contacts": len(ref_contacts),
        "contact_intersection": len(inter),
        "contact_precision": precision,
        "contact_recall": recall,
        "contact_jaccard": jaccard,
    }


def pairwise_clash_count(
    left_pos: np.ndarray,
    right_pos: np.ndarray,
    cutoff: float,
    exclude_mask: np.ndarray | None = None,
) -> tuple[int, float | None]:
    if left_pos.size == 0 or right_pos.size == 0:
        return 0, None
    distances = np.linalg.norm(left_pos[:, None, :] - right_pos[None, :, :], axis=-1)
    if exclude_mask is not None:
        distances = np.where(exclude_mask, np.inf, distances)
    return int((distances < cutoff).sum()), float(np.min(distances))


def protein_nonlocal_clashes(
    protein_pos: np.ndarray,
    residue_ids: np.ndarray,
    cutoff: float,
    neighbor_gap: int = 1,
) -> tuple[int, float | None]:
    n_atoms = protein_pos.shape[0]
    if n_atoms <= 1:
        return 0, None
    distances = np.linalg.norm(protein_pos[:, None, :] - protein_pos[None, :, :], axis=-1)
    residue_gap = np.abs(residue_ids[:, None] - residue_ids[None, :])
    exclude = np.triu(np.ones((n_atoms, n_atoms), dtype=bool), k=0) | (residue_gap <= neighbor_gap)
    distances = np.where(exclude, np.inf, distances)
    finite = distances[np.isfinite(distances)]
    min_distance = float(finite.min()) if finite.size else None
    return int((distances < cutoff).sum()), min_distance


def approximate_core_shell(
    ligand_pos: np.ndarray,
    centers: np.ndarray,
    topk: int,
    shell_radius: float,
) -> tuple[np.ndarray, np.ndarray]:
    if centers.size == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    distances_to_ligand = np.linalg.norm(centers[:, None, :] - ligand_pos[None, :, :], axis=-1).min(axis=1)
    k = min(topk, len(centers))
    core = np.argsort(distances_to_ligand)[:k]
    if shell_radius <= 0 or len(core) == 0:
        return core, core.copy()
    center_distances = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=-1)
    shell_mask = (center_distances[:, core] <= shell_radius).any(axis=1)
    active = np.where(shell_mask)[0]
    return core, active


def local_discontinuity_metrics(
    apo_centers: np.ndarray,
    pred_centers: np.ndarray,
    holo_centers: np.ndarray,
    core: np.ndarray,
    active: np.ndarray,
    neighbor_cutoff: float,
) -> dict[str, float | int | None]:
    n_res = pred_centers.shape[0]
    if n_res == 0:
        return {
            "approx_core_residues": 0,
            "approx_active_residues": 0,
            "boundary_pairs": 0,
            "boundary_delta_pred_vs_apo_mean": None,
            "boundary_delta_pred_vs_holo_mean": None,
            "neighbor_displacement_jump_mean": None,
            "active_displacement_mean": None,
            "background_displacement_mean": None,
        }

    active_mask = np.zeros(n_res, dtype=bool)
    active_mask[active] = True
    core_mask = np.zeros(n_res, dtype=bool)
    core_mask[core] = True
    background_mask = ~active_mask

    pred_disp = np.linalg.norm(pred_centers - apo_centers, axis=1)
    active_disp = float(pred_disp[active_mask].mean()) if active_mask.any() else None
    background_disp = float(pred_disp[background_mask].mean()) if background_mask.any() else None

    apo_dist = np.linalg.norm(apo_centers[:, None, :] - apo_centers[None, :, :], axis=-1)
    pred_dist = np.linalg.norm(pred_centers[:, None, :] - pred_centers[None, :, :], axis=-1)
    holo_dist = np.linalg.norm(holo_centers[:, None, :] - holo_centers[None, :, :], axis=-1)
    upper = np.triu(np.ones((n_res, n_res), dtype=bool), k=1)
    neighbor = upper & (apo_dist < neighbor_cutoff)
    boundary = neighbor & (active_mask[:, None] != active_mask[None, :])

    if boundary.any():
        boundary_delta_apo = float(np.abs(pred_dist[boundary] - apo_dist[boundary]).mean())
        boundary_delta_holo = float(np.abs(pred_dist[boundary] - holo_dist[boundary]).mean())
        disp_jump = np.abs(pred_disp[:, None] - pred_disp[None, :])
        jump_mean = float(disp_jump[boundary].mean())
        n_boundary = int(boundary.sum())
    else:
        boundary_delta_apo = None
        boundary_delta_holo = None
        jump_mean = None
        n_boundary = 0

    return {
        "approx_core_residues": int(core_mask.sum()),
        "approx_active_residues": int(active_mask.sum()),
        "boundary_pairs": n_boundary,
        "boundary_delta_pred_vs_apo_mean": boundary_delta_apo,
        "boundary_delta_pred_vs_holo_mean": boundary_delta_holo,
        "neighbor_displacement_jump_mean": jump_mean,
        "active_displacement_mean": active_disp,
        "background_displacement_mean": background_disp,
    }


def summarize(values: list[float | int | None]) -> dict[str, float | int | None]:
    clean = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not clean:
        return {"n": 0, "mean": None, "median": None, "min": None, "max": None}
    arr = np.asarray(clean, dtype=np.float64)
    return {
        "n": int(arr.size),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "min": float(arr.min()),
        "max": float(arr.max()),
    }


def infer_shell_radius(arm_name: str, default_shell_radius: float) -> float:
    if "_shell4_" in arm_name:
        return 4.0
    if "_shell6_" in arm_name:
        return 6.0
    return default_shell_radius


def analyze_result_file(
    result_file: Path,
    arm_name: str,
    topk: int,
    default_shell_radius: float,
    contact_cutoff: float,
    clash_cutoff: float,
    neighbor_cutoff: float,
) -> dict[str, Any]:
    result = torch.load(result_file, map_location="cpu")
    data = result["data"]

    pred_protein_pos = final_array(result, "pred_protein_pos")
    pred_ligand_pos = final_array(result, "pred_ligand_pos")
    if pred_protein_pos is None or pred_ligand_pos is None:
        raise ValueError(f"{result_file} is missing final protein or ligand positions")

    apo_protein_pos = as_numpy(data["protein_pos"])
    holo_protein_pos = as_numpy(data["protein_pos_holo"])
    ref_ligand_pos = as_numpy(data["ligand_pos"])
    residue_ids = as_numpy(data["protein_atom_to_aa_group"]).astype(np.int64)

    residue_unique, apo_centers = residue_centers(apo_protein_pos, residue_ids)
    _, pred_centers = residue_centers(pred_protein_pos, residue_ids)
    _, holo_centers = residue_centers(holo_protein_pos, residue_ids)

    shell_radius = infer_shell_radius(arm_name, default_shell_radius)
    core, active = approximate_core_shell(pred_ligand_pos, pred_centers, topk, shell_radius)
    local_metrics = local_discontinuity_metrics(
        apo_centers,
        pred_centers,
        holo_centers,
        core,
        active,
        neighbor_cutoff,
    )

    pred_contacts = residue_contacts(pred_ligand_pos, pred_protein_pos, residue_ids, contact_cutoff)
    ref_contacts = residue_contacts(ref_ligand_pos, holo_protein_pos, residue_ids, contact_cutoff)
    contacts = contact_metrics(pred_contacts, ref_contacts)

    ligand_protein_clashes, ligand_protein_min_dist = pairwise_clash_count(
        pred_ligand_pos,
        pred_protein_pos,
        clash_cutoff,
    )
    protein_clashes, protein_min_dist = protein_nonlocal_clashes(
        pred_protein_pos,
        residue_ids,
        clash_cutoff,
    )

    output = {
        "result_file": str(result_file),
        "sample_id": int(result_file.stem.replace("result_", "")),
        "num_protein_atoms": int(pred_protein_pos.shape[0]),
        "num_ligand_atoms": int(pred_ligand_pos.shape[0]),
        "num_residues": int(len(residue_unique)),
        "protein_rmsd": float(result.get("pred_protein_pos_rmsd", [np.nan])[-1]),
        "protein_tmscore": float(result.get("pred_protein_pos_tmscore", [np.nan])[-1]),
        "ligand_protein_clashes": ligand_protein_clashes,
        "ligand_protein_min_dist": ligand_protein_min_dist,
        "protein_nonlocal_clashes": protein_clashes,
        "protein_nonlocal_min_dist": protein_min_dist,
        **contacts,
        **local_metrics,
    }
    return output


def analyze_run_dir(args: argparse.Namespace) -> dict[str, Any]:
    run_dir = Path(args.run_dir)
    arms: dict[str, Any] = {}
    for arm_dir in sorted(run_dir.iterdir()):
        sample_dir = arm_dir / "sampled_results"
        if not sample_dir.is_dir():
            continue
        result_files = sorted(sample_dir.glob("result_*.pt"))
        if not result_files:
            continue

        samples = [
            analyze_result_file(
                result_file,
                arm_dir.name,
                args.topk,
                args.shell_radius,
                args.contact_cutoff,
                args.clash_cutoff,
                args.neighbor_cutoff,
            )
            for result_file in result_files
        ]

        metric_names = sorted(
            key
            for key, value in samples[0].items()
            if isinstance(value, (int, float)) or value is None
        )
        summary = {name: summarize([sample.get(name) for sample in samples]) for name in metric_names}
        arms[arm_dir.name] = {
            "num_samples": len(samples),
            "summary": summary,
            "samples": samples,
        }

    return {
        "run_dir": str(run_dir),
        "settings": {
            "topk": args.topk,
            "shell_radius": args.shell_radius,
            "contact_cutoff": args.contact_cutoff,
            "clash_cutoff": args.clash_cutoff,
            "neighbor_cutoff": args.neighbor_cutoff,
        },
        "arms": arms,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, help="Validation run directory containing arm subdirs.")
    parser.add_argument("--out", default=None, help="Output JSON path. Defaults to run_dir/geometry_diagnostics.json.")
    parser.add_argument("--topk", type=int, default=4, help="Top-k residues for approximate active core.")
    parser.add_argument("--shell-radius", type=float, default=0.0, help="Default shell radius for arms without shell in name.")
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--clash-cutoff", type=float, default=2.0)
    parser.add_argument("--neighbor-cutoff", type=float, default=8.0)
    args = parser.parse_args()

    output = analyze_run_dir(args)
    out_path = Path(args.out) if args.out else Path(args.run_dir) / "geometry_diagnostics.json"
    out_path.write_text(json.dumps(output, indent=2, sort_keys=True), encoding="utf-8")

    print(f"Wrote {out_path}")
    print("arm\tprotein_rmsd\tlp_clashes\tlp_min_dist\tcontact_recall\tboundary_jump")
    for arm_name, arm in output["arms"].items():
        summary = arm["summary"]

        def mean(name: str) -> Any:
            return summary.get(name, {}).get("mean")

        def fmt(value: Any) -> str:
            if value is None:
                return ""
            if isinstance(value, float):
                return f"{value:.4f}"
            return str(value)

        print(
            "\t".join(
                [
                    arm_name,
                    fmt(mean("protein_rmsd")),
                    fmt(mean("ligand_protein_clashes")),
                    fmt(mean("ligand_protein_min_dist")),
                    fmt(mean("contact_recall")),
                    fmt(mean("neighbor_displacement_jump_mean")),
                ]
            )
        )


if __name__ == "__main__":
    main()
