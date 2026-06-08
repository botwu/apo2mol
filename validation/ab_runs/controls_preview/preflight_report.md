# New Method A/B Verification

Generated: 2026-05-26T14:58:39

## Status

Runnable: all local prerequisites were found.

Selected-case files checked: 24
Missing selected-case files: 0

## Selected Hard Test Cases

| test position | original index | apo/holo RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |
| 390 | 24413 | 3.3143 | `8u6b__1__1.A__1.C/1.C.sdf` |
| 310 | 24333 | 2.9210 | `8ow3__1__1.A__1.B/1.B.sdf` |
| 377 | 24400 | 2.7470 | `8sfu__2__1.B__1.E/1.E.sdf` |
| 342 | 24365 | 2.6896 | `8qn5__3__1.C__1.L/1.L.sdf` |
| 365 | 24388 | 2.6712 | `8pqh__1__1.A__1.B/1.B.sdf` |
| 347 | 24370 | 2.4535 | `8sbv__2__1.B__1.G/1.G.sdf` |

## Experiment Arms

| arm | sampling protocol | protein update |
|---|---|---|
| baseline_realistic_static5 | atoms=prior, center=apo, steps=100 | static5 |
| adaptive_realistic_residual | atoms=prior, center=apo, steps=100 | residual_adaptive |
| control_realistic_late_dense | atoms=prior, center=apo, steps=100 | late_dense |
| control_realistic_uniform10 | atoms=prior, center=apo, steps=100 | uniform10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/controls_preview
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.
