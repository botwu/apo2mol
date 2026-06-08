# New Method A/B Verification

Generated: 2026-06-05T14:30:15

## Status

Runnable: all local prerequisites were found.

Selected-case files checked: 60
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
| 457 | 24480 | 2.3482 | `8j5x__1__1.A__1.B/1.B.sdf` |
| 282 | 24305 | 2.3015 | `8pyk__1__1.A__1.B/1.B.sdf` |
| 344 | 24367 | 2.2580 | `8qn5__1__1.A__1.E/1.E.sdf` |
| 276 | 24299 | 2.2213 | `7gsf__1__1.A__1.C/1.C.sdf` |
| 281 | 24304 | 2.1823 | `8pyn__1__1.A__1.B/1.B.sdf` |
| 279 | 24302 | 2.1791 | `8u8j__1__1.A__1.B/1.B.sdf` |
| 278 | 24301 | 2.1698 | `8u8k__1__1.A__1.B/1.B.sdf` |
| 396 | 24419 | 2.1659 | `8tqg__1__1.A__1.B/1.B.sdf` |
| 305 | 24328 | 2.1472 | `8k5y__1__1.A__1.H/1.H.sdf` |
| 301 | 24324 | 2.1410 | `8k5y__2__1.B__1.O/1.O.sdf` |
| 304 | 24327 | 2.1011 | `8k5x__2__1.B__1.O/1.O.sdf` |
| 248 | 24271 | 2.0943 | `8uwp__1__1.A__1.C/1.C.sdf` |

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
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard20_original_static5_baseline_steps1000_n3
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "baseline_realistic_static5": {
    "sampling": {
      "available": true,
      "result_files": 20,
      "protein_rmsd": {
        "n": 60,
        "mean": 2.196637290716171,
        "median": 2.102377414703369,
        "min": 1.4533493518829346,
        "max": 3.431992530822754
      },
      "protein_tmscore": {
        "n": 60,
        "mean": 0.9143092761179007,
        "median": 0.9260759235662988,
        "min": 0.8148447168935643,
        "max": 0.9468459695728902
      },
      "seconds": {
        "n": 60,
        "mean": 1359.0048889636994,
        "median": 650.5102739334106,
        "min": 392.5809962749481,
        "max": 39455.234869003296
      },
      "router_selected_counts": {
        "n": 300,
        "mean": 54.1,
        "median": 52.0,
        "min": 28.0,
        "max": 77.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 20.0,
        "mol_stable": 0.4167,
        "atm_stable": 0.9227,
        "recon_success": 0.9,
        "eval_success": 0.8833,
        "complete": 0.8833,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 54.0,
        "complete_mols": 53.0,
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