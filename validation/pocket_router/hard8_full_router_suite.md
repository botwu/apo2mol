# Hard8 Full Router Suite

Date: 2026-05-28

Run directory:

`validation/ab_runs/hard8_full_router_steps1000_n1`

Raw summary:

`validation/ab_runs/hard8_full_router_steps1000_n1/results.json`

## One-Line Result

On the 8 hardest apo-to-holo cases, the best arm is still:

> `pocket_router_distance_top4`

It reaches mean protein RMSD `1.7924 A`, compared with:

- `baseline_realistic_static5`: `2.4365 A`
- `control_realistic_late_dense`: `2.4142 A`
- `pocket_router_random_top4`: `1.8305 A`
- `pocket_router_distance_top8`: `1.8371 A`
- `pocket_router_distance_top12`: `1.8992 A`
- `pocket_router_distance_top16`: `1.9534 A`
- `pocket_router_distance_top24`: `2.0531 A`

Plain-language interpretation:

> The model does best when it only lets the tiny ligand-near local pocket move.
> Updating more residues gradually makes the protein RMSD worse.

## What We Tested

This experiment asks a narrow question:

> In Apo2Mol sampling, should the protein pocket update be dense, random sparse,
> oracle sparse, or current-ligand-distance sparse?

We did not train a new learned router yet. This is still a validation suite for
the mechanism and the baseline difficulty.

Plain-language version:

> Before building a smart learned router, we first tested whether the simple
> rule "only move the residues closest to the current ligand" is already strong.

## Arms

All arms use:

- 8 hard cases
- 1 generated sample per case
- 1000 denoising steps
- realistic apo initialization
- docking disabled: `--docking-mode none`

| arm | meaning | selected residues/update |
|---|---|---:|
| `baseline_realistic_static5` | original-style sparse-time dense pocket update | 52.25 |
| `control_realistic_late_dense` | later 10 dense pocket updates | 52.25 |
| `control_realistic_uniform10` | uniformly spaced 10 dense updates | 53.71 |
| `adaptive_realistic_residual` | global update triggered by residual schedule | 53.13 |
| `pocket_router_random_top4` | random sparse mask, top4 | 4 |
| `pocket_router_random_top12` | random sparse mask, top12 | 12 |
| `pocket_router_motion_oracle_top12` | oracle based on apo-to-holo motion | 12 |
| `pocket_router_contact_oracle_top12` | oracle based on holo contact | 12 |
| `pocket_router_contact_change_oracle_top12` | oracle based on contact change | 12 |
| `pocket_router_distance_top4` | current ligand distance router, top4 | 4 |
| `pocket_router_distance_top8` | current ligand distance router, top8 | 8 |
| `pocket_router_distance_top12` | current ligand distance router, top12 | 12 |
| `pocket_router_distance_top16` | current ligand distance router, top16 | 16 |
| `pocket_router_distance_top24` | current ligand distance router, top24 | 24 |

## Result Table

Lower protein RMSD is better. Higher TM-score is better.

| arm | mean protein RMSD | mean TM-score | complete | evaluated mols |
|---|---:|---:|---:|---:|
| `baseline_realistic_static5` | 2.4365 | 0.8977 | 0.875 | 0 |
| `control_realistic_late_dense` | 2.4142 | 0.9005 | 0.875 | 0 |
| `control_realistic_uniform10` | 3.1475 | 0.8514 | 0.875 | 0 |
| `adaptive_realistic_residual` | 2.0904 | 0.9141 | 0.875 | 0 |
| `pocket_router_random_top4` | 1.8305 | 0.9248 | 0.625 | 0 |
| `pocket_router_random_top12` | 1.9555 | 0.9196 | 1.000 | 0 |
| `pocket_router_motion_oracle_top12` | 2.0511 | 0.9136 | 0.750 | 0 |
| `pocket_router_contact_oracle_top12` | 2.1195 | 0.9065 | 0.750 | 0 |
| `pocket_router_contact_change_oracle_top12` | 2.2410 | 0.9037 | 0.875 | 0 |
| `pocket_router_distance_top4` | 1.7924 | 0.9249 | 0.750 | 0 |
| `pocket_router_distance_top8` | 1.8371 | 0.9229 | 0.750 | 0 |
| `pocket_router_distance_top12` | 1.8992 | 0.9206 | 0.750 | 0 |
| `pocket_router_distance_top16` | 1.9534 | 0.9180 | 0.875 | 0 |
| `pocket_router_distance_top24` | 2.0531 | 0.9140 | 0.750 | 0 |

Important caveat:

> `evaluated_mols` is 0 because docking was disabled. QED/SA and docking-based
> conclusions are not available from this run.

## Distance Top-k Sweep

The distance router has a clean trend:

| distance top-k | mean protein RMSD | mean TM-score |
|---:|---:|---:|
| 4 | 1.7924 | 0.9249 |
| 8 | 1.8371 | 0.9229 |
| 12 | 1.8992 | 0.9206 |
| 16 | 1.9534 | 0.9180 |
| 24 | 2.0531 | 0.9140 |

Plain-language interpretation:

> More pocket residues is not better. The best setting is the smallest tested
> local update budget.

This is the strongest mechanism-level signal in the suite.

## Main Comparisons

`distance_top4` vs dense baselines:

- vs `baseline_realistic_static5`: `-0.6441 A`
- vs `control_realistic_late_dense`: `-0.6218 A`

`distance_top4` vs sparse baselines:

- vs `random_top4`: `-0.0381 A`
- vs `distance_top8`: `-0.0446 A`
- vs `random_top12`: `-0.1631 A`

Plain-language interpretation:

> Sparse masking itself is very strong. Distance top4 wins, but random top4 is
> close. So the learned router cannot just beat dense update; it must beat
> `random_top4` and `distance_top4`.

## What The Oracle Arms Tell Us

The oracle arms did not win:

| oracle arm | mean protein RMSD |
|---|---:|
| `motion_oracle_top12` | 2.0511 |
| `contact_oracle_top12` | 2.1195 |
| `contact_change_oracle_top12` | 2.2410 |

This is important.

It means the story is probably not:

> Just find residues that move a lot from apo to holo.

The better story is:

> The current ligand state should decide a very small local pocket workspace.

In other words, the router should be fragment/contact-conditioned, not just a
static oracle over globally mobile residues.

## What We Changed To Run This Suite

The method-side PocketRouter support was already in place:

- protein update schedules;
- router modes;
- top-k residue masks;
- router selected count logging.

For this full suite, the experiment harness was extended with:

- a `FULL_ROUTER_ARMS` list in `validation/run_new_method_ab.py`;
- a `--router-full-suite` flag;
- configs for the 14 arms above;
- resume logic based on existing `result_*.pt` files.

During the long run, we also added a small scheduling guard in `sample_split.py`:

- if `result_{id}.pt` already exists, skip that case;
- if `result_{id}.pt.lock` exists, wait for an external prefill job to finish;
- `APO2MOL_IGNORE_RESULT_LOCKS=1` lets the external prefill job itself write
  the locked result.

Plain-language version:

> Several slow arms were being filled in parallel. The lock guard prevented the
> main controller and the prefill workers from accidentally generating the same
> result file at the same time.

## Completion Check

The final sample count is:

```text
14 arms x 8 cases = 112 result_*.pt files
```

Observed:

```text
112 result_*.pt files
```

The controller exited successfully and wrote:

`validation/ab_runs/hard8_full_router_steps1000_n1/results.json`

## Current Top-Conference Conclusion

This run makes the idea stronger, but also raises the bar.

What is strong:

- dense pocket update is not the best way to handle hard apo-to-holo cases;
- very sparse local update is consistently strong;
- distance top4 is the best current hand-crafted baseline;
- top-k sweep has a clean monotonic pattern: top4 > top8 > top12 > top16 > top24.

What is still weak:

- random top4 is very close to distance top4;
- no learned router has beaten distance top4 yet;
- n=8 and n=1 sample per case is still pilot-scale;
- docking/QED/SA are not available in this run;
- ligand quality must be checked before claiming full SBDD improvement.

Practical next step:

> Train a learned fragment-conditioned router, but judge it against
> `distance_top4` and `random_top4`, not only against dense Apo2Mol baselines.

