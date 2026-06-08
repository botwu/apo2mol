# New Method A/B Verification

Generated: 2026-06-05T01:28:59

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

| arm | sampling protocol | protein update | pocket router |
|---|---|---|---|
| baseline_realistic_static5 | atoms=prior, center=apo, steps=1000 | static5 | none topk=0 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| baseline_realistic_static5 | 999 to 0 | True | 799, 599, 399, 199, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard8_original_static5_baseline_steps1000_n3
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "baseline_realistic_static5": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 24,
        "mean": 2.4449720829725266,
        "median": 2.5210769176483154,
        "min": 1.610964059829712,
        "max": 3.1260809898376465
      },
      "protein_tmscore": {
        "n": 24,
        "mean": 0.8961856736323349,
        "median": 0.9101989410735749,
        "min": 0.8148447168935643,
        "max": 0.9335825020926339
      },
      "seconds": {
        "n": 24,
        "mean": 1952.3153098324935,
        "median": 592.8981920480728,
        "min": 409.076936006546,
        "max": 13610.304042100906
      },
      "router_selected_counts": {
        "n": 120,
        "mean": 52.25,
        "median": 52.5,
        "min": 42.0,
        "max": 60.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.5,
        "atm_stable": 0.9449,
        "recon_success": 0.875,
        "eval_success": 0.875,
        "complete": 0.875,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 21.0,
        "complete_mols": 21.0,
        "evaluated_mols": 0.0,
        "QED_Mean": null,
        "QED_Median": null,
        "SA_Mean": null,
        "SA_Median": null,
        "ring_size_3_ratio": null,
        "ring_size_4_ratio": null,
        "ring_size_5_ratio": null,
        "ring_size_6_ratio": null,
        "ring_size_7_ratio": null,
        "ring_size_8_ratio": null,
        "ring_size_9_ratio": null
      }
    }
  }
}
```