# New Method A/B Verification

Generated: 2026-05-26T15:24:26

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
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard1_steps100_controls_n1_gated
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
        "mean": 40.85400414466858,
        "median": 40.85400414466858,
        "min": 40.85400414466858,
        "max": 40.85400414466858
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
        "mean": 41.38965678215027,
        "median": 41.38965678215027,
        "min": 41.38965678215027,
        "max": 41.38965678215027
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
        "mean": 41.6462459564209,
        "median": 41.6462459564209,
        "min": 41.6462459564209,
        "max": 41.6462459564209
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
        "mean": 41.669597864151,
        "median": 41.669597864151,
        "min": 41.669597864151,
        "max": 41.669597864151
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
    "protein_rmsd_mean_delta": 0.0,
    "adaptive_improves_protein_rmsd": false,
    "adaptive_delta_vs_control_realistic_late_dense": 0.0,
    "adaptive_delta_vs_control_realistic_uniform10": 0.0,
    "adaptive_beats_all_controls": false
  }
}
```