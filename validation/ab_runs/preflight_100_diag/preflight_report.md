# New Method A/B Verification

Generated: 2026-05-26T15:52:27

## Status

Runnable: all local prerequisites were found.

Warnings:
- --num-steps=100 does not reach t=0 in the current sampler; this is a smoke test, not evidence of generation quality.
- At least one active protein-update schedule has zero candidate update hits under this num_steps setting.

Selected-case files checked: 6
Missing selected-case files: 0

## Selected Hard Test Cases

| test position | original index | apo/holo RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |

## Experiment Arms

| arm | sampling protocol | protein update |
|---|---|---|
| baseline_realistic_static5 | atoms=prior, center=apo, steps=100 | static5 |
| adaptive_realistic_residual | atoms=prior, center=apo, steps=100 | residual_adaptive |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| baseline_realistic_static5 | 999 to 900 | False | none |
| adaptive_realistic_residual | 999 to 900 | False | none |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/preflight_100_diag
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.
