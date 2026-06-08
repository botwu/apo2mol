#!/usr/bin/env python3
"""Validate whether Apo2Mol exposes a publishable top-conference direction.

The script intentionally avoids importing the project or PyTorch. It validates
the research cut using artifacts that are available in this checkout:

1. Apo/holo RMSD distribution and hard-tail cases from the dataset index.
2. Split-level distribution to check whether hard cases are underrepresented.
3. Static code/config evidence for underused retrieval prompts, benchmark
   crutches, and asynchronous protein/ligand denoising.

It writes a Markdown memo and a machine-readable JSON summary.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pickle
import re
import statistics
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
INDEX_PATH = REPO_ROOT / "apo2mol_dataset" / "apo2mol_version" / "selected_index_apo_druglike.pkl"
SPLIT_PATH = REPO_ROOT / "apo2mol_dataset" / "apo2mol_version" / "split_druglike.pt"

EXTERNAL_EVIDENCE = [
    {
        "name": "ICLR 2026 reviewer guide",
        "url": "https://iclr.cc/Conferences/2026/ReviewerGuide",
        "use": "New, relevant, impactful knowledge backed by convincing evidence.",
    },
    {
        "name": "NeurIPS reviewer guidelines",
        "url": "https://nips.cc/Conferences/2025/ReviewerGuidelines",
        "use": "Novelty, empirical support, clarity, limitations.",
    },
    {
        "name": "ICML reviewer instructions",
        "url": "https://icml.cc/Conferences/2025/ReviewerInstructions",
        "use": "Significance, technical quality, empirical validation.",
    },
    {
        "name": "AAAI-26 main technical track call",
        "url": "https://aaai.org/conference/aaai/aaai-26/main-technical-track-call/",
        "use": "Originality, rigor, reproducibility, impact.",
    },
    {
        "name": "Apo2Mol",
        "url": "https://arxiv.org/abs/2511.14559",
        "use": "Direct baseline and dataset/method context.",
    },
    {
        "name": "DecompDiff, ICML 2023",
        "url": "https://proceedings.mlr.press/v202/guan23a.html",
        "use": "SBDD diffusion precedent with decomposition as a methodological hook.",
    },
    {
        "name": "IPDiff, ICLR 2024",
        "url": "https://openreview.net/forum?id=qH9nrMNTIW",
        "use": "Interaction-prior diffusion precedent.",
    },
    {
        "name": "FlexSBDD, NeurIPS 2024",
        "url": "https://proceedings.neurips.cc/paper_files/paper/2024/hash/60fb8cf8000f0386063fb24ead366330-Abstract-Conference.html",
        "use": "Flexible-receptor SBDD precedent.",
    },
    {
        "name": "DynamicFlow, ICLR 2025",
        "url": "https://arxiv.org/abs/2503.03989",
        "use": "Protein dynamics for molecular generation precedent.",
    },
]


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return ""


def load_index(path: Path) -> list[tuple[Any, ...]]:
    with path.open("rb") as handle:
        records = pickle.load(handle)
    if not isinstance(records, list):
        raise TypeError(f"Expected list in {path}, got {type(records)!r}")
    return records


def load_split(path: Path) -> dict[str, list[int]]:
    with zipfile.ZipFile(path) as zf:
        data = zf.read("split_druglike/data.pkl")
    split = pickle.loads(data)
    if not isinstance(split, dict):
        raise TypeError(f"Expected dict in {path}, got {type(split)!r}")
    return {str(k): list(v) for k, v in split.items()}


def rmsd(record: tuple[Any, ...]) -> float:
    return float(record[-1])


def date_value(record: tuple[Any, ...]) -> str:
    return str(record[7]) if len(record) > 7 else ""


def pdb_id(record: tuple[Any, ...]) -> str:
    first = str(record[0]).split("/", 1)[0]
    return first.split("__", 1)[0].lower()


def quantile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * q
    lower = int(pos)
    upper = min(lower + 1, len(sorted_values) - 1)
    frac = pos - lower
    return sorted_values[lower] * (1.0 - frac) + sorted_values[upper] * frac


def bin_counts(values: Iterable[float]) -> dict[str, int]:
    bins = {
        "0-0.5": 0,
        "0.5-1": 0,
        "1-1.5": 0,
        "1.5-2": 0,
        "2-3": 0,
        "3-5": 0,
        "5-10": 0,
        "10+": 0,
    }
    for value in values:
        if value < 0.5:
            bins["0-0.5"] += 1
        elif value < 1.0:
            bins["0.5-1"] += 1
        elif value < 1.5:
            bins["1-1.5"] += 1
        elif value < 2.0:
            bins["1.5-2"] += 1
        elif value < 3.0:
            bins["2-3"] += 1
        elif value < 5.0:
            bins["3-5"] += 1
        elif value < 10.0:
            bins["5-10"] += 1
        else:
            bins["10+"] += 1
    return bins


def summarize_values(values: list[float]) -> dict[str, Any]:
    sorted_values = sorted(values)
    return {
        "n": len(values),
        "min": min(values) if values else None,
        "p25": quantile(sorted_values, 0.25),
        "median": statistics.median(values) if values else None,
        "mean": statistics.fmean(values) if values else None,
        "p75": quantile(sorted_values, 0.75),
        "p90": quantile(sorted_values, 0.90),
        "p95": quantile(sorted_values, 0.95),
        "p99": quantile(sorted_values, 0.99),
        "max": max(values) if values else None,
        "ge_1_5": sum(v >= 1.5 for v in values),
        "ge_2": sum(v >= 2.0 for v in values),
        "ge_3": sum(v >= 3.0 for v in values),
        "ge_5": sum(v >= 5.0 for v in values),
        "bins": bin_counts(values),
    }


def split_summary(records: list[tuple[Any, ...]], split: dict[str, list[int]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name, indices in split.items():
        values = [rmsd(records[i]) for i in indices if 0 <= i < len(records)]
        stats = summarize_values(values)
        stats["index_count"] = len(indices)
        stats["valid_index_count"] = len(values)
        out[name] = stats
    return out


def examples(records: list[tuple[Any, ...]], n: int = 8) -> list[dict[str, Any]]:
    ranked = sorted(records, key=rmsd, reverse=True)[:n]
    return [
        {
            "rmsd": rmsd(record),
            "holo_pocket": str(record[0]),
            "apo_pocket": str(record[1]),
            "ligand": str(record[2]),
            "date": date_value(record),
        }
        for record in ranked
    ]


def scan_code() -> dict[str, Any]:
    files = {
        "sampling_config": REPO_ROOT / "configs" / "sampling.yaml",
        "training_config": REPO_ROOT / "configs" / "training.yaml",
        "sample_split": REPO_ROOT / "sample_split.py",
        "model": REPO_ROOT / "models" / "molopt_score_model.py",
        "eval_split": REPO_ROOT / "eval_split.py",
    }
    texts = {name: read_text(path) for name, path in files.items()}

    checks = {
        "sampling_uses_ref_atom_count": bool(
            re.search(r"sample_num_atoms\s*:\s*ref\b", texts["sampling_config"])
        ),
        "training_topk_prompt_zero_count": len(
            re.findall(r"topk_prompt\s*:\s*0\b", texts["training_config"])
        ),
        "model_has_topk_prompt_parameter": "self.topk_prompt" in texts["model"],
        "model_prompt_list_initialized_empty": bool(
            re.search(r"prompt_hbap_ligand_batch_all_list\s*=\s*\[\]", texts["model"])
        ),
        "model_prompt_fallback_self_attention": "h_retrieved=hbap_ligand_aug_batch" in texts["model"],
        "model_has_static_protein_update_steps": "self.protein_update_steps" in texts["model"],
        "sampling_uses_holo_center_for_ligand_init": (
            "batch.protein_pos_holo" in texts["sample_split"]
            and "scatter_mean(batch.protein_pos_holo" in texts["sample_split"]
        ),
        "eval_supports_gen_apo_holo_pocket_types": all(
            token in texts["eval_split"] for token in ["pocket_type", "'gen'", "'apo'", "'holo'"]
        ),
        "eval_has_vina_or_qvina": (
            "VinaDockingTask" in texts["eval_split"] or "QVinaDockingTask" in texts["eval_split"]
        ),
        "data_folder_exists": (REPO_ROOT / "data_folder").exists(),
        "torch_importable_in_current_python": importlib.util.find_spec("torch") is not None,
    }

    return {
        "files": {name: str(path.relative_to(REPO_ROOT)) for name, path in files.items()},
        "checks": checks,
    }


def top_years(records: list[tuple[Any, ...]], n: int = 10) -> list[tuple[str, int]]:
    years = []
    for record in records:
        value = date_value(record)
        if len(value) >= 4:
            years.append(value[:4])
    return Counter(years).most_common(n)


def dataset_summary(records: list[tuple[Any, ...]]) -> dict[str, Any]:
    values = [rmsd(record) for record in records]
    dates = [date_value(record) for record in records if date_value(record)]
    pdb_counts = Counter(pdb_id(record) for record in records)
    pli_dirs = Counter(str(record[0]).rsplit("/", 1)[0] for record in records)
    return {
        "record_count": len(records),
        "record_tuple_length": len(records[0]) if records else None,
        "rmsd": summarize_values(values),
        "date_min": min(dates) if dates else None,
        "date_max": max(dates) if dates else None,
        "top_years": top_years(records),
        "unique_holo_pdb_ids": len(pdb_counts),
        "most_common_holo_pdb_ids": pdb_counts.most_common(10),
        "unique_pli_dirs": len(pli_dirs),
        "duplicated_pli_dirs": sum(count > 1 for count in pli_dirs.values()),
        "hard_examples": examples(records),
    }


def bool_word(value: bool) -> str:
    return "yes" if value else "no"


def pct(count: int, total: int) -> str:
    return f"{100.0 * count / total:.2f}%" if total else "n/a"


def fmt_float(value: Any) -> str:
    if value is None:
        return "n/a"
    try:
        return f"{float(value):.4f}"
    except (TypeError, ValueError):
        return str(value)


def make_markdown(summary: dict[str, Any]) -> str:
    data = summary["dataset"]
    rmsd_stats = data["rmsd"]
    code = summary["code"]["checks"]
    splits = summary["splits"]

    hard_2 = int(rmsd_stats["ge_2"])
    hard_3 = int(rmsd_stats["ge_3"])
    hard_5 = int(rmsd_stats["ge_5"])
    total = int(rmsd_stats["n"])

    lines = [
        "# Apo2Mol Top-Conference Direction Validation",
        "",
        f"Generated: {summary['generated_at']}",
        "",
        "## Bottom line",
        "",
        (
            "The local evidence supports a top-conference-worthy direction only if the work is "
            "framed as a new problem/benchmark plus a substantial method, not as a small patch "
            "on Apo2Mol. The strongest framing is conformational-memory or retrieval-augmented "
            "generation for long-tail apo-to-holo pocket motion."
        ),
        "",
        "## Dataset evidence",
        "",
        f"- Records: {data['record_count']}",
        f"- Date range: {data['date_min']} to {data['date_max']}",
        f"- Unique holo PDB IDs: {data['unique_holo_pdb_ids']}",
        f"- Unique PLI directories: {data['unique_pli_dirs']} (duplicated dirs: {data['duplicated_pli_dirs']})",
        f"- Apo/holo RMSD median: {fmt_float(rmsd_stats['median'])}; p90: {fmt_float(rmsd_stats['p90'])}; p95: {fmt_float(rmsd_stats['p95'])}; max: {fmt_float(rmsd_stats['max'])}",
        f"- Hard cases >=2A: {hard_2} ({pct(hard_2, total)}); >=3A: {hard_3} ({pct(hard_3, total)}); >=5A: {hard_5} ({pct(hard_5, total)})",
        "",
        "RMSD bins:",
        "",
        "| bin (A) | count | share |",
        "|---|---:|---:|",
    ]

    for name, count in rmsd_stats["bins"].items():
        lines.append(f"| {name} | {count} | {pct(int(count), total)} |")

    lines.extend(
        [
            "",
            "Split distribution:",
            "",
            "| split | n | median | p90 | p95 | >=2A | >=3A | >=5A |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for split_name in ["train", "valid", "test"]:
        if split_name not in splits:
            continue
        s = splits[split_name]
        lines.append(
            "| {name} | {n} | {median} | {p90} | {p95} | {ge2} | {ge3} | {ge5} |".format(
                name=split_name,
                n=s["n"],
                median=fmt_float(s["median"]),
                p90=fmt_float(s["p90"]),
                p95=fmt_float(s["p95"]),
                ge2=s["ge_2"],
                ge3=s["ge_3"],
                ge5=s["ge_5"],
            )
        )

    lines.extend(
        [
            "",
            "Highest-RMSD examples:",
            "",
            "| RMSD | date | holo pocket | ligand |",
            "|---:|---|---|---|",
        ]
    )
    for item in data["hard_examples"]:
        lines.append(
            f"| {fmt_float(item['rmsd'])} | {item['date']} | `{item['holo_pocket']}` | `{item['ligand']}` |"
        )

    lines.extend(
        [
            "",
            "## Code/config evidence",
            "",
            f"- Sampling uses reference ligand atom count: {bool_word(code['sampling_uses_ref_atom_count'])}",
            f"- Sampling initializes ligand around holo pocket center: {bool_word(code['sampling_uses_holo_center_for_ligand_init'])}",
            f"- Training config sets `topk_prompt: 0`: {code['training_topk_prompt_zero_count']} occurrence(s)",
            f"- Model has a `topk_prompt` parameter path: {bool_word(code['model_has_topk_prompt_parameter'])}",
            f"- Prompt retrieval list is initialized empty in the model: {bool_word(code['model_prompt_list_initialized_empty'])}",
            f"- Prompt path falls back to self-attention over the current ligand features: {bool_word(code['model_prompt_fallback_self_attention'])}",
            f"- Model has static protein update timesteps: {bool_word(code['model_has_static_protein_update_steps'])}",
            f"- Evaluation supports generated/apo/holo pocket types: {bool_word(code['eval_supports_gen_apo_holo_pocket_types'])}",
            f"- Evaluation has Vina/QVina docking hooks: {bool_word(code['eval_has_vina_or_qvina'])}",
            f"- `data_folder` exists locally: {bool_word(code['data_folder_exists'])}",
            f"- PyTorch importable in current Python: {bool_word(code['torch_importable_in_current_python'])}",
            "",
            "Interpretation:",
            "",
            (
                "The code already contains enough structure to support three research cuts: "
                "a harder benchmark split, real retrieval prompts, and adaptive protein-ligand "
                "co-denoising. It also exposes evaluation risks: reference atom count and holo-center "
                "initialization can make results less realistic unless explicitly separated into oracle "
                "and realistic protocols."
            ),
            "",
            "## Top-conference judgment",
            "",
            "**Can this target a top conference?** Yes, but the paper should not be positioned as "
            "`Apo2Mol + retrieval`. The publishable version needs all of the following:",
            "",
            "1. A new problem definition: long-tail apo-to-holo pocket motion under realistic ligand-generation constraints.",
            "2. A benchmark contribution: Apo2Mol-Hard or Apo2Mol-Realistic with explicit RMSD buckets and no oracle ligand-size/holo-center shortcuts in the main protocol.",
            "3. A method contribution: retrieval as conformational memory, or adaptive asynchronous co-denoising, with clear algorithmic novelty.",
            "4. Strong ablations: retrieval source, hard-tail buckets, oracle-vs-realistic sampling, dynamic-vs-static protein updates, and pocket type.",
            "5. External baselines: Apo2Mol, DiffSBDD/DecompDiff/IPDiff-style SBDD methods, and flexible/dynamic protein-generation baselines where feasible.",
            "",
            "A weak submission is likely if the novelty is only enabling `topk_prompt > 0` or adding a simple nearest-neighbor prompt. The top-conference story becomes credible when the benchmark demonstrates a failure mode that existing methods do not address and the method materially improves hard-tail performance without relying on oracle information.",
            "",
            "## Evidence links to cite",
            "",
        ]
    )

    for item in EXTERNAL_EVIDENCE:
        lines.append(f"- [{item['name']}]({item['url']}): {item['use']}")

    lines.extend(
        [
            "",
            "## Next verification experiments",
            "",
            "1. Create an Apo2Mol-Hard test split with RMSD >=2A and >=3A buckets, then report metrics by bucket.",
            "2. Re-run Apo2Mol under two protocols: oracle/reference and realistic/no-reference atom count plus apo-centered initialization.",
            "3. Implement retrieval memory with strict train-only retrieval and leakage checks by PDB cluster/date.",
            "4. Replace static protein update steps with an adaptive scheduler and ablate compute-matched variants.",
            "5. Report validity, docking/Vina score, interaction recovery, generated-pocket RMSD, diversity, and per-bucket degradation.",
            "",
            "## Local limitations",
            "",
            (
                "This validation is static and dataset-index based. It does not run training or sampling because "
                "the current Python cannot import PyTorch and the full `data_folder` used by the configs is not "
                "present in this checkout. The evidence proves the research opportunity and benchmark risk; it "
                "does not prove the final method will improve metrics."
            ),
            "",
        ]
    )
    return "\n".join(lines)


def build_summary() -> dict[str, Any]:
    records = load_index(INDEX_PATH)
    split = load_split(SPLIT_PATH)
    return {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "repo_root": str(REPO_ROOT),
        "index_path": str(INDEX_PATH.relative_to(REPO_ROOT)),
        "split_path": str(SPLIT_PATH.relative_to(REPO_ROOT)),
        "dataset": dataset_summary(records),
        "splits": split_summary(records, split),
        "code": scan_code(),
        "external_evidence": EXTERNAL_EVIDENCE,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json-out",
        type=Path,
        default=REPO_ROOT / "validation" / "apo2mol_validation.json",
        help="Path for the JSON summary.",
    )
    parser.add_argument(
        "--md-out",
        type=Path,
        default=REPO_ROOT / "validation" / "apo2mol_validation_report.md",
        help="Path for the Markdown report.",
    )
    args = parser.parse_args()

    summary = build_summary()
    markdown = make_markdown(summary)

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.md_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    args.md_out.write_text(markdown, encoding="utf-8")

    print(f"Wrote JSON: {args.json_out}")
    print(f"Wrote Markdown: {args.md_out}")
    print()
    print("Bottom line:")
    print(
        "Top-conference target is plausible only as a benchmark/problem/method paper; "
        "a small Apo2Mol retrieval patch is not enough."
    )
    print(
        "Hard-tail cases: >=2A={ge2}, >=3A={ge3}, >=5A={ge5}, max={max_rmsd:.4f}A".format(
            ge2=summary["dataset"]["rmsd"]["ge_2"],
            ge3=summary["dataset"]["rmsd"]["ge_3"],
            ge5=summary["dataset"]["rmsd"]["ge_5"],
            max_rmsd=summary["dataset"]["rmsd"]["max"],
        )
    )


if __name__ == "__main__":
    main()
