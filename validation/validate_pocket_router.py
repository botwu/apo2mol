#!/usr/bin/env python3
"""Validate whether sparse pocket routing is worth implementing.

This script is deliberately data-first. It asks a falsifiable question before
training a learned router:

Can a small set of ligand-relevant residues explain most apo-to-holo pocket
adaptation on hard Apo2Mol cases?

If oracle sparse residue selection cannot improve a holo replay upper bound, a
PocketRouter method is unlikely to help. If oracle improves and simple routers
cover useful motion/contact residues, the route is worth model-level validation.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import random
import statistics
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any

import torch
from omegaconf import OmegaConf

REPO_ROOT = Path(__file__).resolve().parents[1]
import sys

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from datasets import get_dataset


INDEX_PATH = REPO_ROOT / "apo2mol_dataset" / "apo2mol_version" / "selected_index_apo_druglike.pkl"
SPLIT_PATH = REPO_ROOT / "apo2mol_dataset" / "apo2mol_version" / "split_druglike.pt"
TRAINING_CONFIG = REPO_ROOT / "configs" / "training.yaml"


def load_index(path: Path) -> list[tuple[Any, ...]]:
    with path.open("rb") as handle:
        return pickle.load(handle)


def load_split(path: Path) -> dict[str, list[int]]:
    with zipfile.ZipFile(path) as zf:
        data = zf.read("split_druglike/data.pkl")
    split = pickle.loads(data)
    return {str(k): list(v) for k, v in split.items()}


def select_hard_cases(records: list[tuple[Any, ...]], split: dict[str, list[int]], threshold: float, num_cases: int) -> list[dict[str, Any]]:
    candidates = []
    for test_position, original_index in enumerate(split["test"]):
        record = records[original_index]
        candidates.append(
            {
                "test_position": test_position,
                "original_index": int(original_index),
                "rmsd": float(record[-1]),
                "holo_pocket": str(record[0]),
                "apo_pocket": str(record[1]),
                "ligand": str(record[2]),
            }
        )
    hard = [item for item in candidates if item["rmsd"] >= threshold]
    hard.sort(key=lambda item: item["rmsd"], reverse=True)
    if len(hard) < num_cases:
        candidates.sort(key=lambda item: item["rmsd"], reverse=True)
        return candidates[:num_cases]
    return hard[:num_cases]


def residue_slices(atom_to_residue: torch.Tensor) -> list[torch.Tensor]:
    groups = []
    for residue_id in range(int(atom_to_residue.max().item()) + 1):
        groups.append((atom_to_residue == residue_id).nonzero(as_tuple=True)[0])
    return groups


def residue_centers(pos: torch.Tensor, groups: list[torch.Tensor]) -> torch.Tensor:
    return torch.stack([pos[idx].mean(dim=0) for idx in groups], dim=0)


def residue_motion(data: Any, groups: list[torch.Tensor]) -> torch.Tensor:
    motions = []
    for idx in groups:
        disp = torch.sqrt(((data.protein_pos_holo[idx] - data.protein_pos[idx]) ** 2).sum(dim=-1))
        motions.append(disp.mean())
    return torch.stack(motions)


def residue_ligand_min_distance(residue_pos: torch.Tensor, ligand_pos: torch.Tensor) -> torch.Tensor:
    return torch.cdist(residue_pos, ligand_pos).min(dim=1).values


def atom_rmsd(pos: torch.Tensor, ref: torch.Tensor) -> float:
    return float(torch.sqrt(((pos - ref) ** 2).sum(dim=-1)).mean().item())


def replay_sparse_holo(data: Any, selected_residues: set[int], groups: list[torch.Tensor]) -> float:
    replay = data.protein_pos.clone()
    for residue_id in selected_residues:
        replay[groups[residue_id]] = data.protein_pos_holo[groups[residue_id]]
    return atom_rmsd(replay, data.protein_pos_holo)


def topk_from_scores(scores: torch.Tensor, k: int, largest: bool) -> set[int]:
    k = max(1, min(k, scores.numel()))
    indices = torch.topk(scores, k=k, largest=largest).indices.tolist()
    return {int(i) for i in indices}


def summarize(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"n": 0}
    return {
        "n": len(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def evaluate_router(
    data: Any,
    topk_values: list[int],
    random_trials: int,
    rng: random.Random,
) -> dict[str, Any]:
    groups = residue_slices(data.protein_atom_to_aa_group)
    n_res = len(groups)
    apo_rmsd = atom_rmsd(data.protein_pos, data.protein_pos_holo)
    motion = residue_motion(data, groups)
    apo_centers = residue_centers(data.protein_pos, groups)
    holo_centers = residue_centers(data.protein_pos_holo, groups)
    apo_dist = residue_ligand_min_distance(apo_centers, data.ligand_pos)
    holo_dist = residue_ligand_min_distance(holo_centers, data.ligand_pos)

    motion_total = float(motion.sum().item()) + 1e-12
    mobile_budget = max(1, math.ceil(0.2 * n_res))
    mobile_set = topk_from_scores(motion, mobile_budget, largest=True)
    holo_contact_set = {int(i) for i, dist in enumerate(holo_dist.tolist()) if dist <= 4.5}
    contact_change = torch.clamp(apo_dist - holo_dist, min=0.0)

    routers: dict[str, dict[str, Any]] = {
        "distance": {"score": -apo_dist, "largest": True},
        "motion_oracle": {"score": motion, "largest": True},
        "contact_oracle": {"score": -holo_dist, "largest": True},
        "contact_change_oracle": {"score": contact_change, "largest": True},
    }

    result: dict[str, Any] = {
        "num_residues": n_res,
        "apo_atom_rmsd": apo_rmsd,
        "holo_contact_residues": len(holo_contact_set),
        "mobile_residues_top20pct": len(mobile_set),
        "topk": {},
    }

    all_residues = list(range(n_res))
    for k in topk_values:
        k = min(k, n_res)
        k_result: dict[str, Any] = {}
        for name, spec in routers.items():
            selected = topk_from_scores(spec["score"], k, largest=spec["largest"])
            selected_motion = float(motion[list(selected)].sum().item()) if selected else 0.0
            k_result[name] = {
                "k": k,
                "fraction": k / n_res,
                "motion_coverage": selected_motion / motion_total,
                "mobile_recall": len(selected & mobile_set) / max(1, len(mobile_set)),
                "holo_contact_recall": len(selected & holo_contact_set) / max(1, len(holo_contact_set)),
                "sparse_holo_replay_rmsd": replay_sparse_holo(data, selected, groups),
            }

        random_metrics = []
        for _ in range(random_trials):
            selected = set(rng.sample(all_residues, k))
            selected_motion = float(motion[list(selected)].sum().item()) if selected else 0.0
            random_metrics.append(
                {
                    "motion_coverage": selected_motion / motion_total,
                    "mobile_recall": len(selected & mobile_set) / max(1, len(mobile_set)),
                    "holo_contact_recall": len(selected & holo_contact_set) / max(1, len(holo_contact_set)),
                    "sparse_holo_replay_rmsd": replay_sparse_holo(data, selected, groups),
                }
            )
        k_result["random"] = {
            "k": k,
            "fraction": k / n_res,
            **{key: summarize([item[key] for item in random_metrics]) for key in random_metrics[0]},
        }
        result["topk"][str(k)] = k_result
    return result


def aggregate_case_results(case_results: list[dict[str, Any]]) -> dict[str, Any]:
    if not case_results:
        return {}
    topk_keys = sorted(case_results[0]["router"].get("topk", {}).keys(), key=lambda item: int(item))
    router_names = list(next(iter(case_results[0]["router"]["topk"].values())).keys())
    aggregate: dict[str, Any] = {
        "apo_atom_rmsd": summarize([case["router"]["apo_atom_rmsd"] for case in case_results]),
        "topk": {},
    }
    for k in topk_keys:
        aggregate["topk"][k] = {}
        for router in router_names:
            router_items = [case["router"]["topk"][k][router] for case in case_results]
            aggregate["topk"][k][router] = {}
            for metric in ["motion_coverage", "mobile_recall", "holo_contact_recall", "sparse_holo_replay_rmsd"]:
                if router == "random":
                    values = [item[metric]["mean"] for item in router_items]
                else:
                    values = [item[metric] for item in router_items]
                aggregate["topk"][k][router][metric] = summarize(values)
    return aggregate


def write_report(output_dir: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# PocketRouter Data Validation",
        "",
        f"Generated: {payload['generated_at']}",
        "",
        "## Question",
        "",
        "Can sparse residue routing explain hard apo-to-holo pocket adaptation before training a learned router?",
        "",
        "## Selected Cases",
        "",
        "| test position | original index | metadata RMSD | ligand |",
        "|---:|---:|---:|---|",
    ]
    for case in payload["cases"]:
        item = case["case"]
        lines.append(
            f"| {item['test_position']} | {item['original_index']} | {item['rmsd']:.4f} | `{item['ligand']}` |"
        )

    lines.extend(["", "## Aggregate", ""])
    aggregate = payload["aggregate"]
    lines.append(f"Mean apo atom RMSD: {aggregate['apo_atom_rmsd']['mean']:.4f} A")
    lines.extend(
        [
            "",
            "| top-k | router | motion coverage | mobile recall | holo-contact recall | replay RMSD |",
            "|---:|---|---:|---:|---:|---:|",
        ]
    )
    for k, routers in aggregate["topk"].items():
        for router, metrics in routers.items():
            lines.append(
                "| {k} | {router} | {mc:.3f} | {mr:.3f} | {hr:.3f} | {rr:.4f} |".format(
                    k=k,
                    router=router,
                    mc=metrics["motion_coverage"]["mean"],
                    mr=metrics["mobile_recall"]["mean"],
                    hr=metrics["holo_contact_recall"]["mean"],
                    rr=metrics["sparse_holo_replay_rmsd"]["mean"],
                )
            )

    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "- Positive route signal: oracle top-k replay RMSD is much lower than apo RMSD, and distance/non-learned routing beats random on motion or contact recall.",
            "- Negative route signal: oracle top-k replay is close to apo RMSD, meaning sparse pocket memory cannot explain the hard motion.",
            "",
        ]
    )
    (output_dir / "pocket_router_validation.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=str(REPO_ROOT / "validation" / "pocket_router" / "hard8"))
    parser.add_argument("--num-cases", type=int, default=8)
    parser.add_argument("--hard-threshold", type=float, default=2.0)
    parser.add_argument("--topk", default="4,8,12,16")
    parser.add_argument("--random-trials", type=int, default=64)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    output_dir = Path(args.run_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    topk_values = [int(value) for value in args.topk.split(",") if value.strip()]

    records = load_index(INDEX_PATH)
    split = load_split(SPLIT_PATH)
    selected = select_hard_cases(records, split, args.hard_threshold, args.num_cases)

    original_ids_path = output_dir / "hard_original_indices.txt"
    original_ids_path.write_text(
        "\n".join(str(item["original_index"]) for item in selected) + "\n",
        encoding="utf-8",
    )
    (output_dir / "selected_cases.json").write_text(
        json.dumps(selected, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    os.environ.setdefault("APO2MOL_PREPROCESS_WORKERS", "1")
    os.environ["APO2MOL_PREPROCESS_INDEX_FILE"] = str(original_ids_path)
    os.environ["APO2MOL_PROCESSED_PATH"] = str(output_dir / "selected_cases_data.lmdb")

    cfg = OmegaConf.load(TRAINING_CONFIG)
    subsets = get_dataset(cfg.data)
    test_set = subsets["test"]
    rng = random.Random(args.seed)

    case_results = []
    for item in selected:
        data = test_set[item["test_position"]]
        router_result = evaluate_router(data, topk_values, args.random_trials, rng)
        case_results.append({"case": item, "router": router_result})

    payload = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "run_dir": str(output_dir),
        "topk": topk_values,
        "random_trials": args.random_trials,
        "cases": case_results,
        "aggregate": aggregate_case_results(case_results),
    }
    (output_dir / "pocket_router_validation.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_report(output_dir, payload)
    print(f"Wrote {output_dir / 'pocket_router_validation.md'}")


if __name__ == "__main__":
    main()
