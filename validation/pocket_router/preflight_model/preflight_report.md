# New Method A/B Verification

Generated: 2026-05-26T21:00:30

## Status

Runnable: all local prerequisites were found.

Selected-case files checked: 3
Missing selected-case files: 0

## Selected Hard Test Cases

| test position | original index | apo/holo RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |

## Experiment Arms

| arm | sampling protocol | protein update | pocket router |
|---|---|---|---|
| baseline_realistic_static5 | atoms=prior, center=apo, steps=1000 | static5 | none topk=0 |
| control_global_late_dense | atoms=prior, center=apo, steps=1000 | late_dense | none topk=0 |
| pocket_router_distance_top12 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=12 |
| pocket_router_motion_oracle_top12 | atoms=prior, center=apo, steps=1000 | late_dense | motion_oracle topk=12 |
| pocket_router_random_top12 | atoms=prior, center=apo, steps=1000 | late_dense | random topk=12 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| baseline_realistic_static5 | 999 to 0 | True | 799, 599, 399, 199, 10 |
| control_global_late_dense | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_motion_oracle_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_random_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/pocket_router/preflight_model
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.
