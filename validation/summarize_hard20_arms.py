"""5-arm hard20 summary.

Reads /root/autodl-tmp/vina_eval/<arm>/ligand_full.json for each arm
and prints a side-by-side table (Vina min mean/median, QED mean,
High Affinity, n_docked) plus a paired comparison vs baseline_static5.

Output: stdout markdown table + optional --json out.

Run on the AutoDL instance after run_hard20_learned_gate.py finishes:
    python validation/summarize_hard20_arms.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

DEFAULT_ROOT = Path("/root/autodl-tmp/vina_eval")
ARM_ORDER = [
    "baseline_static5",
    "ours_shell4_w025",
    "distance_top4_hard",
    "random_top4",
    "learned_gate",
]


def load_arm(root: Path, name: str) -> dict | None:
    path = root / name / "ligand_full.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def fmt_pair(stats: dict | None, key: str) -> str:
    if stats is None or key not in stats or stats[key] is None:
        return "  -  "
    sub = stats[key]
    if isinstance(sub, dict):
        m = sub.get("mean")
        med = sub.get("median")
        return f"{m:+.3f} / {med:+.3f}" if m is not None else "  -  "
    return f"{sub:+.4f}"


def fmt_scalar(stats: dict | None, key: str, fmt: str = "{:.3f}") -> str:
    if stats is None or key not in stats or stats[key] is None:
        return "  -  "
    return fmt.format(stats[key])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help=f"directory of <arm>/ligand_full.json (default {DEFAULT_ROOT})")
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args()

    rows = {}
    for name in ARM_ORDER:
        d = load_arm(args.root, name)
        rows[name] = (d or {}).get("summary")

    # Print markdown table
    print(f"# hard20 5-arm comparison")
    print(f"source: {args.root}\n")
    print("| arm | Vina mean / median | QED mean / med | SA mean | HighAff | n_docked |")
    print("|---|---|---|---|---|---|")
    for name in ARM_ORDER:
        s = rows[name]
        if s is None:
            print(f"| {name} | (no ligand_full.json) | | | | |")
            continue
        vina = fmt_pair(s, "vina_min")
        qed = fmt_pair(s, "qed")
        sa = fmt_scalar(s.get("sa", {}), "mean") if s.get("sa") else "  -  "
        ha = fmt_scalar(s, "high_affinity")
        nd = s.get("n_docked", "-")
        print(f"| {name} | {vina} | {qed} | {sa} | {ha} | {nd} |")

    print()
    # Paired comparison vs baseline
    base = load_arm(args.root, "baseline_static5") or {}
    base_per = {(p["case"], p["sample"]): p for p in base.get("per_sample", [])}

    print("## Paired vs baseline_static5 (case+sample matched, Vina delta)")
    print("| arm | n_paired | Δ Vina mean | win rate (lower is better) |")
    print("|---|---|---|---|")
    for name in ARM_ORDER:
        if name == "baseline_static5":
            continue
        d = load_arm(args.root, name)
        if d is None:
            continue
        per = d.get("per_sample", [])
        pairs = []
        for p in per:
            key = (p["case"], p["sample"])
            b = base_per.get(key)
            if b is None or b.get("vina") is None or p.get("vina") is None:
                continue
            pairs.append(p["vina"] - b["vina"])
        if not pairs:
            print(f"| {name} | 0 | - | - |")
            continue
        delta_mean = sum(pairs) / len(pairs)
        wins = sum(1 for d in pairs if d < 0)
        print(f"| {name} | {len(pairs)} | {delta_mean:+.3f} | {wins}/{len(pairs)} ({wins/len(pairs):.1%}) |")

    if args.out_json:
        args.out_json.write_text(json.dumps({"arms": rows}, indent=2))
        print(f"\nSaved: {args.out_json}")


if __name__ == "__main__":
    main()
