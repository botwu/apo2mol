#!/usr/bin/env python3
"""Run or prepare a minimal A/B test for Apo2Mol method variants.

This script is intentionally conservative:

* It selects hard test-set cases from the existing Apo2Mol split by apo/holo RMSD.
* It writes two controlled experiment arms:
  - baseline_realistic_static5: apo-centered, prior atom count, original static updates.
  - adaptive_realistic_residual: same protocol, residual-adaptive protein updates.
* With --include-controls, it also adds compute/schedule controls so residual
  adaptivity is not mistaken for a generic "more pocket updates" effect.
* It preflights data and dependency requirements before running.
* With --run, it launches sampling/eval and summarizes protein RMSD results.

The default mode does not require PyTorch; it prepares the exact run directory
and reports what is missing.
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import os
import pickle
import re
import shutil
import statistics
import subprocess
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
INDEX_PATH = REPO_ROOT / "apo2mol_dataset" / "apo2mol_version" / "selected_index_apo_druglike.pkl"
SPLIT_PATH = REPO_ROOT / "apo2mol_dataset" / "apo2mol_version" / "split_druglike.pt"
BASE_SAMPLING_CONFIG = REPO_ROOT / "configs" / "sampling.yaml"
BASE_TRAINING_CONFIG = REPO_ROOT / "configs" / "training.yaml"
NUM_TIMESTEPS = 1000


ARMS = [
    {
        "name": "baseline_realistic_static5",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "static5",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
        },
    },
    {
        "name": "adaptive_realistic_residual",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "residual_adaptive",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
        },
    },
]

CONTROL_ARMS = [
    {
        "name": "control_realistic_late_dense",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
        },
    },
    {
        "name": "control_realistic_uniform10",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "uniform10",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
        },
    },
]

POCKET_ROUTER_ARMS = [
    {
        "name": "baseline_realistic_static5",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "static5",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "control_realistic_late_dense",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "pocket_router_distance_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_motion_oracle_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "motion_oracle",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_contact_change_oracle_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "contact_change_oracle",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_random_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "random",
            "pocket_router_topk": 12,
        },
    },
]


HARDTAIL_CORE_ROUTER_ARMS = [
    {
        "name": "baseline_realistic_static5",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "static5",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "control_realistic_late_dense",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "pocket_router_random_top4",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "random",
            "pocket_router_topk": 4,
        },
    },
    {
        "name": "pocket_router_distance_top4",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 4,
        },
    },
    {
        "name": "pocket_router_distance_top8",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 8,
        },
    },
]


FULL_ROUTER_ARMS = [
    {
        "name": "baseline_realistic_static5",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "static5",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "control_realistic_late_dense",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "control_realistic_uniform10",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "uniform10",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "adaptive_realistic_residual",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "residual_adaptive",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "none",
            "pocket_router_topk": 0,
        },
    },
    {
        "name": "pocket_router_random_top4",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "random",
            "pocket_router_topk": 4,
        },
    },
    {
        "name": "pocket_router_random_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "random",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_motion_oracle_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "motion_oracle",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_contact_oracle_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "contact_oracle",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_contact_change_oracle_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "contact_change_oracle",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_distance_top4",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 4,
        },
    },
    {
        "name": "pocket_router_distance_top8",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 8,
        },
    },
    {
        "name": "pocket_router_distance_top12",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 12,
        },
    },
    {
        "name": "pocket_router_distance_top16",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 16,
        },
    },
    {
        "name": "pocket_router_distance_top24",
        "sample": {
            "num_steps": NUM_TIMESTEPS,
            "sample_num_atoms": "prior",
            "init_center_mode": "apo",
        },
        "model": {
            "protein_update_schedule": "late_dense",
            "protein_update_interval": 50,
            "protein_update_min_t": 10,
            "protein_update_residual_threshold": 0.5,
            "pocket_router_mode": "distance",
            "pocket_router_topk": 24,
        },
    },
]


def make_distance_topk_sweep_arms(topks: list[int]) -> list[dict[str, Any]]:
    arms = []
    for topk in topks:
        arms.append(
            {
                "name": f"pocket_router_distance_top{topk}",
                "sample": {
                    "num_steps": NUM_TIMESTEPS,
                    "sample_num_atoms": "prior",
                    "init_center_mode": "apo",
                },
                "model": {
                    "protein_update_schedule": "late_dense",
                    "protein_update_interval": 50,
                    "protein_update_min_t": 10,
                    "protein_update_residual_threshold": 0.5,
                    "pocket_router_mode": "distance",
                    "pocket_router_topk": topk,
                },
            }
        )
    return arms


def make_active_set_ablation_arms() -> list[dict[str, Any]]:
    base_sample = {
        "num_steps": NUM_TIMESTEPS,
        "sample_num_atoms": "prior",
        "init_center_mode": "apo",
    }
    base_model = {
        "protein_update_schedule": "late_dense",
        "protein_update_interval": 50,
        "protein_update_min_t": 10,
        "protein_update_residual_threshold": 0.5,
    }
    return [
        {
            "name": "control_realistic_late_dense",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "none",
                "pocket_router_topk": 0,
                "pocket_router_shell_radius": 0.0,
                "pocket_router_shell_weight": 0.0,
                "pocket_router_background_weight": 0.0,
            },
        },
        {
            "name": "pocket_router_random_top4",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "random",
                "pocket_router_topk": 4,
                "pocket_router_shell_radius": 0.0,
                "pocket_router_shell_weight": 0.0,
                "pocket_router_background_weight": 0.0,
            },
        },
        {
            "name": "pocket_router_distance_top4_hard",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "distance",
                "pocket_router_topk": 4,
                "pocket_router_shell_radius": 0.0,
                "pocket_router_shell_weight": 0.0,
                "pocket_router_background_weight": 0.0,
            },
        },
        {
            "name": "active_set_distance_top4_shell4_w025",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "distance",
                "pocket_router_topk": 4,
                "pocket_router_shell_radius": 4.0,
                "pocket_router_shell_weight": 0.25,
                "pocket_router_background_weight": 0.0,
            },
        },
        {
            "name": "active_set_distance_top4_shell6_w025",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "distance",
                "pocket_router_topk": 4,
                "pocket_router_shell_radius": 6.0,
                "pocket_router_shell_weight": 0.25,
                "pocket_router_background_weight": 0.0,
            },
        },
        {
            "name": "active_set_distance_top4_shell6_w050",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "distance",
                "pocket_router_topk": 4,
                "pocket_router_shell_radius": 6.0,
                "pocket_router_shell_weight": 0.50,
                "pocket_router_background_weight": 0.0,
            },
        },
        {
            "name": "active_set_random_top4_shell6_w025",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "random",
                "pocket_router_topk": 4,
                "pocket_router_shell_radius": 6.0,
                "pocket_router_shell_weight": 0.25,
                "pocket_router_background_weight": 0.0,
            },
        },
        {
            "name": "pocket_router_cross_attn_gate",
            "sample": dict(base_sample),
            "model": {
                **base_model,
                "pocket_router_mode": "cross_attn_gate",
                "pocket_router_topk": 0,
                "pocket_router_shell_radius": 0.0,
                "pocket_router_shell_weight": 0.0,
                "pocket_router_background_weight": 0.0,
            },
        },
    ]


def make_active_set_shell_only_arms() -> list[dict[str, Any]]:
    return [
        arm
        for arm in make_active_set_ablation_arms()
        if arm["name"].startswith("active_set_")
    ]


def make_active_set_conservative_shell_sweep_arms() -> list[dict[str, Any]]:
    base_sample = {
        "num_steps": NUM_TIMESTEPS,
        "sample_num_atoms": "prior",
        "init_center_mode": "apo",
    }
    base_model = {
        "protein_update_schedule": "late_dense",
        "protein_update_interval": 50,
        "protein_update_min_t": 10,
        "protein_update_residual_threshold": 0.5,
        "pocket_router_mode": "distance",
        "pocket_router_topk": 4,
        "pocket_router_background_weight": 0.0,
    }
    settings = [
        ("active_set_distance_top4_shell3_w010", 3.0, 0.10),
        ("active_set_distance_top4_shell3_w025", 3.0, 0.25),
        ("active_set_distance_top4_shell4_w010", 4.0, 0.10),
        ("active_set_distance_top4_shell5_w025", 5.0, 0.25),
    ]
    arms = []
    for name, shell_radius, shell_weight in settings:
        arms.append(
            {
                "name": name,
                "sample": dict(base_sample),
                "model": {
                    **base_model,
                    "pocket_router_shell_radius": shell_radius,
                    "pocket_router_shell_weight": shell_weight,
                },
            }
        )
    return arms


def make_active_set_candidate_repeat_arms() -> list[dict[str, Any]]:
    base_sample = {
        "num_steps": NUM_TIMESTEPS,
        "sample_num_atoms": "prior",
        "init_center_mode": "apo",
    }
    base_model = {
        "protein_update_schedule": "late_dense",
        "protein_update_interval": 50,
        "protein_update_min_t": 10,
        "protein_update_residual_threshold": 0.5,
        "pocket_router_topk": 4,
        "pocket_router_background_weight": 0.0,
    }
    settings = [
        ("pocket_router_random_top4", "random", 0.0, 0.0),
        ("pocket_router_distance_top4_hard", "distance", 0.0, 0.0),
        ("active_set_distance_top4_shell4_w025", "distance", 4.0, 0.25),
        ("active_set_distance_top4_shell3_w010", "distance", 3.0, 0.10),
        ("active_set_distance_top4_shell3_w025", "distance", 3.0, 0.25),
        ("active_set_distance_top4_shell4_w010", "distance", 4.0, 0.10),
    ]
    arms = []
    for name, router_mode, shell_radius, shell_weight in settings:
        arms.append(
            {
                "name": name,
                "sample": dict(base_sample),
                "model": {
                    **base_model,
                    "pocket_router_mode": router_mode,
                    "pocket_router_shell_radius": shell_radius,
                    "pocket_router_shell_weight": shell_weight,
                },
            }
        )
    return arms


def parse_topks(raw: str) -> list[int]:
    topks: list[int] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        topk = int(item)
        if topk <= 0:
            raise argparse.ArgumentTypeError("router top-k values must be positive")
        topks.append(topk)
    if not topks:
        raise argparse.ArgumentTypeError("at least one router top-k value is required")
    return topks


def load_index(path: Path) -> list[tuple[Any, ...]]:
    with path.open("rb") as handle:
        return pickle.load(handle)


def load_split(path: Path) -> dict[str, list[int]]:
    with zipfile.ZipFile(path) as zf:
        data = zf.read("split_druglike/data.pkl")
    split = pickle.loads(data)
    return {str(k): list(v) for k, v in split.items()}


def rmsd(record: tuple[Any, ...]) -> float:
    return float(record[-1])


def yaml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    return str(value)


def set_two_space_key(text: str, key: str, value: Any, parent: str) -> str:
    replacement = f"  {key}: {yaml_scalar(value)}"
    pattern = re.compile(rf"(^  {re.escape(key)}:\s*).*$", re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(replacement, text)

    parent_match = re.search(rf"^{re.escape(parent)}:\n", text, flags=re.MULTILINE)
    if not parent_match:
        raise ValueError(f"Cannot find `{parent}:` section")
    insert_at = parent_match.end()
    return text[:insert_at] + replacement + "\n" + text[insert_at:]


def make_sampling_config(base: str, overrides: dict[str, Any]) -> str:
    text = base
    for key, value in overrides.items():
        parent = "model" if key == "checkpoint" else "sample"
        text = set_two_space_key(text, key, value, parent=parent)
    return text


def make_training_config(base: str, overrides: dict[str, Any]) -> str:
    text = base
    for key, value in overrides.items():
        text = set_two_space_key(text, key, value, parent="model")
    return text


def current_sampler_time_sequence(num_steps: int, num_timesteps: int = NUM_TIMESTEPS) -> list[int]:
    """Mirror MolOptScoreModel.sample_diffusion's adjacent-step sampler."""
    return list(reversed(range(num_timesteps - num_steps, num_timesteps)))


def configured_protein_update_steps(model_settings: dict[str, Any], num_timesteps: int = NUM_TIMESTEPS) -> set[int]:
    schedule = model_settings["protein_update_schedule"]
    final_step = num_timesteps - 1
    if schedule == "static5":
        return {int(final_step * (1 - i / 5)) for i in range(1, 5)} | {10}
    if schedule == "none":
        return set()
    if schedule == "uniform10":
        return {int(final_step * (1 - i / 10)) for i in range(1, 10)} | {10}
    if schedule in {"late_dense", "residual_adaptive"}:
        interval = max(1, int(model_settings["protein_update_interval"]))
        min_t = int(model_settings["protein_update_min_t"])
        return set(range(min_t, (final_step // 2) + 1, interval)) | {10}
    raise ValueError(f"Unknown protein_update_schedule: {schedule}")


def sampling_diagnostics(arm: dict[str, Any], num_timesteps: int = NUM_TIMESTEPS) -> dict[str, Any]:
    num_steps = int(arm["sample"]["num_steps"])
    time_seq = current_sampler_time_sequence(num_steps, num_timesteps)
    update_steps = configured_protein_update_steps(arm["model"], num_timesteps)
    hit_steps = sorted(set(time_seq) & update_steps, reverse=True)
    terminal_step = min(time_seq) if time_seq else None
    return {
        "sampler": "current_adjacent_tail_chain",
        "num_timesteps": num_timesteps,
        "num_steps": num_steps,
        "starts_at_t": max(time_seq) if time_seq else None,
        "terminal_t": terminal_step,
        "reaches_t0": terminal_step == 0,
        "configured_update_steps": sorted(update_steps, reverse=True),
        "candidate_update_steps_hit": hit_steps,
        "candidate_update_hit_count": len(hit_steps),
    }


def select_hard_test_positions(
    records: list[tuple[Any, ...]],
    split: dict[str, list[int]],
    threshold: float,
    num_cases: int,
) -> list[dict[str, Any]]:
    candidates = []
    for test_pos, original_idx in enumerate(split["test"]):
        record = records[original_idx]
        candidates.append(
            {
                "test_position": test_pos,
                "original_index": original_idx,
                "rmsd": rmsd(record),
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


def preflight(selected_cases: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    modules = [
        "torch",
        "torch_geometric",
        "torch_scatter",
        "torch_sparse",
        "torch_cluster",
        "torch_spline_conv",
        "pytorch_lightning",
        "rdkit",
        "Bio",
        "lmdb",
        "omegaconf",
        "openbabel",
        "openbabel.pybel",
        "kornia",
        "quaternion",
        "wandb",
        "einops",
        "huggingface_hub",
    ]
    module_status = {name: importlib.util.find_spec(name) is not None for name in modules}
    paths = {
        "data_folder": REPO_ROOT / "apo2mol_dataset" / "data_folder",
        "checkpoint": REPO_ROOT / "apo2mol_dataset" / "apo2mol_checkpoint.ckpt",
        "pmnet": REPO_ROOT / "pretrained_models" / "PMINet",
        "index": INDEX_PATH,
        "split": SPLIT_PATH,
    }
    path_status = {name: path.exists() for name, path in paths.items()}
    missing = [name for name, ok in path_status.items() if not ok]
    missing += [name for name, ok in module_status.items() if not ok]
    data_files: dict[str, Any] = {"checked": False, "missing": []}
    data_folder = paths["data_folder"]
    if path_status["data_folder"] and selected_cases:
        required_files = []
        for item in selected_cases:
            required_files.extend([item["holo_pocket"], item["apo_pocket"], item["ligand"]])
        missing_files = [
            filename for filename in required_files if not (data_folder / filename).exists()
        ]
        data_files = {
            "checked": True,
            "required_count": len(required_files),
            "missing": missing_files,
        }
        if missing_files:
            missing.append("selected_case_files")
    return {
        "paths": {name: {"path": str(path), "exists": path_status[name]} for name, path in paths.items()},
        "modules": module_status,
        "data_files": data_files,
        "missing": missing,
    }


def run_command(cmd: list[str], cwd: Path, extra_env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(REPO_ROOT / ".mplconfig"))
    env.setdefault("WANDB_MODE", "disabled")
    env.setdefault("APO2MOL_PREPROCESS_WORKERS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("VECLIB_MAXIMUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")
    env.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    env.setdefault("KMP_INIT_AT_FORK", "FALSE")
    env.setdefault("KMP_WARNINGS", "FALSE")
    if extra_env:
        env.update(extra_env)
    subprocess.run(cmd, cwd=str(cwd), env=env, check=True)


def aggregate_sampling(sample_dir: Path) -> dict[str, Any]:
    try:
        import torch
    except ImportError:
        return {"available": False, "reason": "torch is not importable"}

    rmsds: list[float] = []
    tmscores: list[float] = []
    times: list[float] = []
    router_counts: list[float] = []
    result_files = sorted(glob.glob(str(sample_dir / "result_*.pt")))
    for filename in result_files:
        result = torch.load(filename, map_location="cpu", weights_only=False)
        for value in result.get("pred_protein_pos_rmsd", []):
            rmsds.append(float(value))
        for value in result.get("pred_protein_pos_tmscore", []):
            tmscores.append(float(value))
        for value in result.get("time", []):
            times.append(float(value))
        for value in result.get("router_selected_counts", []):
            router_counts.append(float(value))

    def stats(values: list[float]) -> dict[str, Any]:
        if not values:
            return {"n": 0}
        return {
            "n": len(values),
            "mean": statistics.fmean(values),
            "median": statistics.median(values),
            "min": min(values),
            "max": max(values),
        }

    return {
        "available": True,
        "result_files": len(result_files),
        "protein_rmsd": stats(rmsds),
        "protein_tmscore": stats(tmscores),
        "seconds": stats(times),
        "router_selected_counts": stats(router_counts),
    }


def parse_eval_log(eval_dir: Path) -> dict[str, Any]:
    log_paths = [eval_dir / "eval_v2.log", eval_dir / "log.txt"]
    existing = [path for path in log_paths if path.exists()]
    if not existing:
        return {"available": False, "reason": "eval logs not found"}
    metrics: dict[str, Any] = {}
    number = r"-?(?:\d+(?:\.\d*)?|\.\d+|nan)"

    def parse_number(value: str) -> float | None:
        if value.lower() == "nan":
            return None
        return float(value)

    for log_path in existing:
        for raw_line in log_path.read_text(errors="ignore").splitlines():
            line = raw_line.strip()
            if "::evaluate::INFO]" in line:
                line = line.split("::evaluate::INFO]", 1)[1].strip()

            count_match = re.match(
                r"^Number of reconstructed mols:\s+(\d+), complete mols:\s+(\d+), evaluated mols:\s+(\d+)",
                line,
            )
            if count_match:
                metrics["Number_of_reconstructed_mols"] = float(count_match.group(1))
                metrics["complete_mols"] = float(count_match.group(2))
                metrics["evaluated_mols"] = float(count_match.group(3))
                continue

            summary_match = re.match(
                rf"^(QED|SA):\s+Mean:\s+({number})\s+Median:\s+({number})",
                line,
                flags=re.IGNORECASE,
            )
            if summary_match:
                key = summary_match.group(1).upper()
                metrics[f"{key}_Mean"] = parse_number(summary_match.group(2))
                metrics[f"{key}_Median"] = parse_number(summary_match.group(3))
                continue

            ring_match = re.match(
                rf"^ring size:\s+(\d+)\s+ratio:\s+({number})",
                line,
                flags=re.IGNORECASE,
            )
            if ring_match:
                metrics[f"ring_size_{ring_match.group(1)}_ratio"] = parse_number(ring_match.group(2))
                continue

            match = re.match(
                rf"^([A-Za-z0-9_ /|.-]+):\s+({number})",
                line,
                flags=re.IGNORECASE,
            )
            if match:
                key = match.group(1).strip().replace(" ", "_")
                metrics[key] = parse_number(match.group(2))
    return {"available": True, "metrics": metrics}


def collect_results(payload: dict[str, Any]) -> dict[str, Any]:
    results: dict[str, Any] = {}
    for arm in payload["arms"]:
        sample_path = Path(arm["sample_path"])
        eval_path = Path(arm["eval_path"])
        results[arm["name"]] = {
            "sampling": aggregate_sampling(sample_path),
            "eval": parse_eval_log(eval_path),
        }

    baseline = results.get("baseline_realistic_static5", {}).get("sampling", {})
    adaptive = results.get("adaptive_realistic_residual", {}).get("sampling", {})
    baseline_rmsd = baseline.get("protein_rmsd", {}).get("mean")
    adaptive_rmsd = adaptive.get("protein_rmsd", {}).get("mean")
    if baseline_rmsd is not None and adaptive_rmsd is not None:
        comparison: dict[str, Any] = {
            "protein_rmsd_mean_delta": adaptive_rmsd - baseline_rmsd,
            "adaptive_improves_protein_rmsd": adaptive_rmsd < baseline_rmsd,
        }
        control_rmsds = {}
        for control_name in ["control_realistic_late_dense", "control_realistic_uniform10"]:
            control_rmsd = (
                results.get(control_name, {})
                .get("sampling", {})
                .get("protein_rmsd", {})
                .get("mean")
            )
            if control_rmsd is not None:
                control_rmsds[control_name] = control_rmsd
                comparison[f"adaptive_delta_vs_{control_name}"] = adaptive_rmsd - control_rmsd
        if control_rmsds:
            comparison["adaptive_beats_all_controls"] = all(
                adaptive_rmsd < control_rmsd for control_rmsd in control_rmsds.values()
            )
        results["comparison"] = comparison
    if "pocket_router_motion_oracle_top12" in results:
        router_comparison: dict[str, Any] = {}
        rmsds = {}
        for name, result in results.items():
            if not isinstance(result, dict):
                continue
            value = result.get("sampling", {}).get("protein_rmsd", {}).get("mean")
            if value is not None:
                rmsds[name] = value
        if rmsds:
            best_name = min(rmsds, key=rmsds.get)
            router_comparison["best_arm"] = best_name
            router_comparison["best_protein_rmsd"] = rmsds[best_name]
            baseline_value = rmsds.get("baseline_realistic_static5")
            global_late_value = rmsds.get("control_realistic_late_dense")
            for name, value in rmsds.items():
                if baseline_value is not None:
                    router_comparison[f"{name}_delta_vs_baseline_static5"] = value - baseline_value
                if global_late_value is not None:
                    router_comparison[f"{name}_delta_vs_realistic_late_dense"] = value - global_late_value
        results["router_comparison"] = router_comparison
    distance_topk_rmsds = {}
    for name, result in results.items():
        if not isinstance(result, dict):
            continue
        match = re.match(r"^pocket_router_distance_top(\d+)$", name)
        if not match:
            continue
        value = result.get("sampling", {}).get("protein_rmsd", {}).get("mean")
        if value is not None:
            distance_topk_rmsds[int(match.group(1))] = value
    if len(distance_topk_rmsds) >= 2:
        best_topk = min(distance_topk_rmsds, key=distance_topk_rmsds.get)
        results["distance_topk_comparison"] = {
            "best_topk": best_topk,
            "best_protein_rmsd": distance_topk_rmsds[best_topk],
            "protein_rmsd_by_topk": {
                str(topk): distance_topk_rmsds[topk]
                for topk in sorted(distance_topk_rmsds)
            },
        }
    core_rmsds = {}
    for name, result in results.items():
        if not isinstance(result, dict):
            continue
        value = result.get("sampling", {}).get("protein_rmsd", {}).get("mean")
        if value is not None and (
            name in {"baseline_realistic_static5", "control_realistic_late_dense"}
            or name.startswith("pocket_router_random_top")
            or name.startswith("pocket_router_distance_top")
        ):
            core_rmsds[name] = value
    if {"baseline_realistic_static5", "control_realistic_late_dense"} <= set(core_rmsds) and len(core_rmsds) >= 3:
        best_name = min(core_rmsds, key=core_rmsds.get)
        baseline_value = core_rmsds["baseline_realistic_static5"]
        late_value = core_rmsds["control_realistic_late_dense"]
        results["hardtail_core_comparison"] = {
            "best_arm": best_name,
            "best_protein_rmsd": core_rmsds[best_name],
            "protein_rmsd_by_arm": {
                name: core_rmsds[name] for name in sorted(core_rmsds)
            },
            "delta_vs_baseline_static5": {
                name: value - baseline_value for name, value in sorted(core_rmsds.items())
            },
            "delta_vs_late_dense": {
                name: value - late_value for name, value in sorted(core_rmsds.items())
            },
        }
    return results


def write_report(run_dir: Path, payload: dict[str, Any]) -> None:
    missing = payload["preflight"]["missing"]
    warnings = payload.get("warnings", [])
    lines = [
        "# New Method A/B Verification",
        "",
        f"Generated: {payload['generated_at']}",
        "",
        "## Status",
        "",
    ]
    if missing:
        lines.append("Not runnable yet. Missing requirements:")
        lines.append("")
        for item in missing:
            lines.append(f"- {item}")
        if "data_folder" in missing:
            lines.extend(
                [
                    "",
                    "Download gated data after exporting an authorized Hugging Face token:",
                    "",
                    "```bash",
                    f"HF_TOKEN=... {sys.executable} validation/download_apo2mol_data.py",
                    "```",
                ]
            )
    else:
        lines.append("Runnable: all local prerequisites were found.")
    if warnings:
        lines.extend(["", "Warnings:"])
        for warning in warnings:
            lines.append(f"- {warning}")
    data_files = payload["preflight"].get("data_files", {})
    if data_files.get("checked"):
        lines.extend(
            [
                "",
                f"Selected-case files checked: {data_files.get('required_count', 0)}",
                f"Missing selected-case files: {len(data_files.get('missing', []))}",
            ]
        )

    lines.extend(
        [
            "",
            "## Selected Hard Test Cases",
            "",
            "| test position | original index | apo/holo RMSD | ligand |",
            "|---:|---:|---:|---|",
        ]
    )
    for item in payload["selected_cases"]:
        lines.append(
            f"| {item['test_position']} | {item['original_index']} | {item['rmsd']:.4f} | `{item['ligand']}` |"
        )

    lines.extend(
        [
            "",
            "## Experiment Arms",
            "",
            "| arm | sampling protocol | protein update | pocket router |",
            "|---|---|---|---|",
        ]
    )
    for arm in payload["arms"]:
        sample = arm["sample"]
        model = arm["model"]
        lines.append(
            "| {name} | atoms={atoms}, center={center}, steps={steps} | {schedule} | {router} topk={topk} |".format(
                name=arm["name"],
                atoms=sample["sample_num_atoms"],
                center=sample["init_center_mode"],
                steps=sample["num_steps"],
                schedule=model["protein_update_schedule"],
                router=model.get("pocket_router_mode", "none"),
                topk=model.get("pocket_router_topk", 0),
            )
        )

    lines.extend(
        [
            "",
            "## Sampling Diagnostics",
            "",
            "The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.",
            "",
            "| arm | timestep range | reaches t=0 | candidate update hits |",
            "|---|---:|---:|---|",
        ]
    )
    for arm in payload["arms"]:
        diag = arm["sampling_diagnostics"]
        hits = ", ".join(str(step) for step in diag["candidate_update_steps_hit"]) or "none"
        lines.append(
            "| {name} | {start} to {terminal} | {reaches} | {hits} |".format(
                name=arm["name"],
                start=diag["starts_at_t"],
                terminal=diag["terminal_t"],
                reaches=diag["reaches_t0"],
                hits=hits,
            )
        )

    lines.extend(
        [
            "",
            "## Run Command",
            "",
            "```bash",
            f"{sys.executable} validation/run_new_method_ab.py --run --run-dir {run_dir}",
            "```",
            "",
            "A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.",
            "",
        ]
    )

    if payload.get("results"):
        lines.extend(["## Results", ""])
        lines.append("```json")
        lines.append(json.dumps(payload["results"], indent=2))
        lines.append("```")

    (run_dir / "preflight_report.md").write_text("\n".join(lines), encoding="utf-8")


def prepare_run(args: argparse.Namespace) -> dict[str, Any]:
    if not 1 <= args.num_steps <= NUM_TIMESTEPS:
        raise SystemExit(f"--num-steps must be between 1 and {NUM_TIMESTEPS}")
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    records = load_index(INDEX_PATH)
    split = load_split(SPLIT_PATH)
    selected = select_hard_test_positions(records, split, args.hard_threshold, args.num_cases)

    ids_path = run_dir / "hard_test_positions.txt"
    ids_path.write_text("\n".join(str(item["test_position"]) for item in selected) + "\n", encoding="utf-8")
    original_ids_path = run_dir / "hard_original_indices.txt"
    original_ids_path.write_text(
        "\n".join(str(item["original_index"]) for item in selected) + "\n",
        encoding="utf-8",
    )
    (run_dir / "selected_cases.json").write_text(
        json.dumps(selected, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    base_sampling = BASE_SAMPLING_CONFIG.read_text(encoding="utf-8")
    base_training = BASE_TRAINING_CONFIG.read_text(encoding="utf-8")
    if args.baseline_only:
        arms = [ARMS[0]]
    elif args.router_full_suite:
        arms = FULL_ROUTER_ARMS
    elif args.router_hardtail_core:
        arms = HARDTAIL_CORE_ROUTER_ARMS
    elif args.router_topk_sweep:
        arms = make_distance_topk_sweep_arms(args.router_topks)
    elif args.active_set_candidate_repeat:
        arms = make_active_set_candidate_repeat_arms()
    elif args.active_set_conservative_shell_sweep:
        arms = make_active_set_conservative_shell_sweep_arms()
    elif args.active_set_shell_only:
        arms = make_active_set_shell_only_arms()
    elif args.active_set_ablation:
        arms = make_active_set_ablation_arms()
    elif args.router_validation:
        arms = POCKET_ROUTER_ARMS
    else:
        arms = ARMS + (CONTROL_ARMS if args.include_controls else [])
    arms_out = []
    warnings_out = []
    for arm in arms:
        arm_dir = run_dir / arm["name"]
        arm_dir.mkdir(parents=True, exist_ok=True)
        sample_overrides = dict(arm["sample"])
        sample_overrides["num_steps"] = args.num_steps
        if arm["name"] == "pocket_router_cross_attn_gate":
            if args.cross_attn_gate_checkpoint:
                sample_overrides["checkpoint"] = args.cross_attn_gate_checkpoint
            else:
                warnings_out.append(
                    "pocket_router_cross_attn_gate needs --cross-attn-gate-checkpoint; "
                    "without it, sampling will fail because the original checkpoint has no cross-attn gate weights."
                )
        model_overrides = dict(arm["model"])
        arm_payload = {
            "name": arm["name"],
            "sample": sample_overrides,
            "model": model_overrides,
        }
        arm_payload["sampling_diagnostics"] = sampling_diagnostics(arm_payload)
        sample_config = make_sampling_config(base_sampling, sample_overrides)
        train_config = make_training_config(base_training, model_overrides)
        sample_config_path = arm_dir / "sampling.yaml"
        train_config_path = arm_dir / "training.yaml"
        sample_config_path.write_text(sample_config, encoding="utf-8")
        train_config_path.write_text(train_config, encoding="utf-8")
        arms_out.append(
            {
                **arm_payload,
                "sample_config": str(sample_config_path),
                "train_config": str(train_config_path),
                "sample_path": str(arm_dir / "sampled_results"),
                "eval_path": str(arm_dir / "eval_results"),
            }
        )

    payload = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "repo_root": str(REPO_ROOT),
        "run_dir": str(run_dir),
        "ids_path": str(ids_path),
        "original_ids_path": str(original_ids_path),
        "processed_lmdb_path": str(run_dir / "selected_cases_data.lmdb"),
        "selected_cases": selected,
        "arms": arms_out,
        "preflight": preflight(selected),
        "warnings": warnings_out,
    }
    if args.num_steps < NUM_TIMESTEPS:
        payload["warnings"].append(
            f"--num-steps={args.num_steps} does not reach t=0 in the current sampler; this is a smoke test, not evidence of generation quality."
        )
    if any(
        arm["sampling_diagnostics"]["candidate_update_hit_count"] == 0
        for arm in arms_out
        if arm["model"]["protein_update_schedule"] != "none"
    ):
        payload["warnings"].append(
            "At least one active protein-update schedule has zero candidate update hits under this num_steps setting."
        )
    (run_dir / "experiment_plan.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_report(run_dir, payload)
    return payload


def run_experiment(payload: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    if payload["preflight"]["missing"]:
        raise SystemExit(
            "Cannot run experiment; missing requirements: "
            + ", ".join(payload["preflight"]["missing"])
        )

    data_env = {
        "APO2MOL_PREPROCESS_INDEX_FILE": payload["original_ids_path"],
        "APO2MOL_PROCESSED_PATH": payload["processed_lmdb_path"],
    }
    expected_ids = [int(item["test_position"]) for item in payload["selected_cases"]]
    for arm in payload["arms"]:
        sample_path = Path(arm["sample_path"])
        eval_path = Path(arm["eval_path"])
        if sample_path.exists() and args.clean:
            shutil.rmtree(sample_path)
        if eval_path.exists() and args.clean:
            shutil.rmtree(eval_path)
        sample_path.mkdir(parents=True, exist_ok=True)
        eval_path.mkdir(parents=True, exist_ok=True)

        existing_ids = {
            int(path.stem.replace("result_", ""))
            for path in sample_path.glob("result_*.pt")
            if path.stem.replace("result_", "").isdigit()
        }
        missing_ids = [idx for idx in expected_ids if idx not in existing_ids]
        did_sample = False
        if not missing_ids and not args.clean:
            print(
                f"Skipping sampling for {arm['name']}; all {len(expected_ids)} result_*.pt files found.",
                flush=True,
            )
        else:
            ids_file = payload["ids_path"]
            if missing_ids and not args.clean:
                ids_path = sample_path.parent / "missing_test_positions.txt"
                ids_path.write_text(
                    "\n".join(str(idx) for idx in missing_ids) + "\n",
                    encoding="utf-8",
                )
                ids_file = str(ids_path)
                print(
                    f"Sampling {arm['name']} for missing test positions: {missing_ids}",
                    flush=True,
                )
            run_command(
                [
                    sys.executable,
                    "sample_split.py",
                    "--config",
                    arm["sample_config"],
                    "--train_config",
                    arm["train_config"],
                    "--device",
                    args.device,
                    "--num_samples",
                    str(args.num_samples),
                    "--batch_size",
                    str(args.batch_size),
                    "--result_path",
                    str(sample_path),
                    "--data_ids_file",
                    ids_file,
                ],
                cwd=REPO_ROOT,
                extra_env=data_env,
            )
            did_sample = True

        if (
            ((eval_path / "log.txt").exists() or (eval_path / "eval_v2.log").exists())
            and not args.clean
            and not did_sample
        ):
            print(f"Skipping eval for {arm['name']}; existing eval logs found.", flush=True)
        else:
            run_command(
                [
                    sys.executable,
                    "eval_split.py",
                    "--sample_path",
                    str(sample_path),
                    "--result_path",
                    str(eval_path),
                    "--docking_mode",
                    args.docking_mode,
                    "--eval_start_index",
                    "0",
                    "--eval_end_index",
                    "1000000",
                ],
                cwd=REPO_ROOT,
                extra_env=data_env,
            )
    return collect_results(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=str(REPO_ROOT / "validation" / "ab_runs" / "latest"))
    parser.add_argument("--hard-threshold", type=float, default=2.0)
    parser.add_argument("--num-cases", type=int, default=8)
    parser.add_argument("--num-samples", type=int, default=3)
    parser.add_argument("--batch-size", type=int, default=3)
    parser.add_argument("--num-steps", type=int, default=NUM_TIMESTEPS)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--cross-attn-gate-checkpoint", default=None)
    parser.add_argument("--docking-mode", default="none", choices=["none", "qvina", "vina_score", "vina_dock"])
    parser.add_argument("--include-controls", action="store_true")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--router-validation", action="store_true")
    parser.add_argument("--router-full-suite", action="store_true")
    parser.add_argument("--router-hardtail-core", action="store_true")
    parser.add_argument("--router-topk-sweep", action="store_true")
    parser.add_argument("--active-set-ablation", action="store_true")
    parser.add_argument("--active-set-shell-only", action="store_true")
    parser.add_argument("--active-set-candidate-repeat", action="store_true")
    parser.add_argument("--active-set-conservative-shell-sweep", action="store_true")
    parser.add_argument("--router-topks", type=parse_topks, default=parse_topks("4,8,12,16,24"))
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()

    payload = prepare_run(args)
    if args.run or args.summarize_only:
        results = collect_results(payload) if args.summarize_only else run_experiment(payload, args)
        payload["results"] = results
        run_dir = Path(payload["run_dir"])
        (run_dir / "results.json").write_text(
            json.dumps(results, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        write_report(run_dir, payload)
    else:
        print("Prepared A/B experiment. Use --run after prerequisites are available.")

    print(f"Run directory: {payload['run_dir']}")
    if payload["preflight"]["missing"]:
        print("Missing:", ", ".join(payload["preflight"]["missing"]))
    else:
        print("All prerequisites found.")


if __name__ == "__main__":
    main()
