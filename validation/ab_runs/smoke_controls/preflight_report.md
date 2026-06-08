# New Method A/B Verification

Generated: 2026-05-26T15:08:54

## Status

Runnable: all local prerequisites were found.

Selected-case files checked: 3
Missing selected-case files: 0

## Selected Hard Test Cases

| test position | original index | apo/holo RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |

## Experiment Arms

| arm | sampling protocol | protein update |
|---|---|---|
| baseline_realistic_static5 | atoms=prior, center=apo, steps=100 | static5 |
| adaptive_realistic_residual | atoms=prior, center=apo, steps=100 | residual_adaptive |
| control_realistic_late_dense | atoms=prior, center=apo, steps=100 | late_dense |
| control_realistic_uniform10 | atoms=prior, center=apo, steps=100 | uniform10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/smoke_controls
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "baseline_realistic_static5": {
    "sampling": {
      "available": true,
      "result_files": 1,
      "protein_rmsd": {
        "n": 1,
        "mean": 1.780302882194519,
        "median": 1.780302882194519,
        "min": 1.780302882194519,
        "max": 1.780302882194519
      },
      "protein_tmscore": {
        "n": 1,
        "mean": 0.8699539234929352,
        "median": 0.8699539234929352,
        "min": 0.8699539234929352,
        "max": 0.8699539234929352
      },
      "seconds": {
        "n": 1,
        "mean": 8.481096982955933,
        "median": 8.481096982955933,
        "min": 8.481096982955933,
        "max": 8.481096982955933
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 1.0
      }
    }
  },
  "adaptive_realistic_residual": {
    "sampling": {
      "available": true,
      "result_files": 1,
      "protein_rmsd": {
        "n": 1,
        "mean": 5.254971504211426,
        "median": 5.254971504211426,
        "min": 5.254971504211426,
        "max": 5.254971504211426
      },
      "protein_tmscore": {
        "n": 1,
        "mean": 0.6514431289320338,
        "median": 0.6514431289320338,
        "min": 0.6514431289320338,
        "max": 0.6514431289320338
      },
      "seconds": {
        "n": 1,
        "mean": 8.371865034103394,
        "median": 8.371865034103394,
        "min": 8.371865034103394,
        "max": 8.371865034103394
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 1.0
      }
    }
  },
  "control_realistic_late_dense": {
    "sampling": {
      "available": true,
      "result_files": 1,
      "protein_rmsd": {
        "n": 1,
        "mean": 1.780302882194519,
        "median": 1.780302882194519,
        "min": 1.780302882194519,
        "max": 1.780302882194519
      },
      "protein_tmscore": {
        "n": 1,
        "mean": 0.8699539234929352,
        "median": 0.8699539234929352,
        "min": 0.8699539234929352,
        "max": 0.8699539234929352
      },
      "seconds": {
        "n": 1,
        "mean": 8.268973112106323,
        "median": 8.268973112106323,
        "min": 8.268973112106323,
        "max": 8.268973112106323
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 1.0
      }
    }
  },
  "control_realistic_uniform10": {
    "sampling": {
      "available": true,
      "result_files": 1,
      "protein_rmsd": {
        "n": 1,
        "mean": 1.780302882194519,
        "median": 1.780302882194519,
        "min": 1.780302882194519,
        "max": 1.780302882194519
      },
      "protein_tmscore": {
        "n": 1,
        "mean": 0.8699539234929352,
        "median": 0.8699539234929352,
        "min": 0.8699539234929352,
        "max": 0.8699539234929352
      },
      "seconds": {
        "n": 1,
        "mean": 8.250765800476074,
        "median": 8.250765800476074,
        "min": 8.250765800476074,
        "max": 8.250765800476074
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 1.0
      }
    }
  },
  "comparison": {
    "protein_rmsd_mean_delta": 3.4746686220169067,
    "adaptive_improves_protein_rmsd": false,
    "adaptive_delta_vs_control_realistic_late_dense": 3.4746686220169067,
    "adaptive_delta_vs_control_realistic_uniform10": 3.4746686220169067,
    "adaptive_beats_all_controls": false
  }
}
```