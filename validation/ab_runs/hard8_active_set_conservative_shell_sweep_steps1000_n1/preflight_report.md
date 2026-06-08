# New Method A/B Verification

Generated: 2026-05-29T17:01:37

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
| active_set_distance_top4_shell3_w010 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell3_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell4_w010 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell5_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| active_set_distance_top4_shell3_w010 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell3_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell4_w010 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell5_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "active_set_distance_top4_shell3_w010": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.85708636790514,
        "median": 1.8690937161445618,
        "min": 0.949084460735321,
        "max": 3.062455892562866
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9214894397419396,
        "median": 0.9350082895990439,
        "min": 0.8575248088773721,
        "max": 0.9591542271205357
      },
      "seconds": {
        "n": 8,
        "mean": 594.881490200758,
        "median": 602.7455314397812,
        "min": 408.1677029132843,
        "max": 731.4424550533295
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 4.325,
        "median": 4.0,
        "min": 4.0,
        "max": 7.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.75,
        "atm_stable": 0.9928,
        "recon_success": 1.0,
        "eval_success": 0.875,
        "complete": 0.875,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 8.0,
        "complete_mols": 7.0,
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
  },
  "active_set_distance_top4_shell3_w025": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.8583604469895363,
        "median": 1.872806966304779,
        "min": 0.949084460735321,
        "max": 3.0652172565460205
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9214154305127654,
        "median": 0.9350057054492833,
        "min": 0.8575260174943276,
        "max": 0.9591542271205357
      },
      "seconds": {
        "n": 8,
        "mean": 582.9044781625271,
        "median": 587.9029669761658,
        "min": 406.6463279724121,
        "max": 720.1467597484589
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 4.3125,
        "median": 4.0,
        "min": 4.0,
        "max": 7.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.75,
        "atm_stable": 0.9928,
        "recon_success": 1.0,
        "eval_success": 0.875,
        "complete": 0.875,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 8.0,
        "complete_mols": 7.0,
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
  },
  "active_set_distance_top4_shell4_w010": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.8512054905295372,
        "median": 1.8814327120780945,
        "min": 0.9490051865577698,
        "max": 3.031834602355957
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.92177708712141,
        "median": 0.9344162936090465,
        "min": 0.8574430257967203,
        "max": 0.9591453334263392
      },
      "seconds": {
        "n": 8,
        "mean": 585.2540013492107,
        "median": 591.2035095691681,
        "min": 401.3900737762451,
        "max": 721.2929439544678
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 6.1,
        "median": 6.0,
        "min": 4.0,
        "max": 9.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.75,
        "atm_stable": 0.9495,
        "recon_success": 0.875,
        "eval_success": 0.75,
        "complete": 0.75,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 7.0,
        "complete_mols": 6.0,
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
  },
  "active_set_distance_top4_shell5_w025": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.8563024401664734,
        "median": 1.8694559335708618,
        "min": 0.9528132677078247,
        "max": 3.049210548400879
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.921907451241187,
        "median": 0.9342271896987946,
        "min": 0.8573468400306827,
        "max": 0.9590774972098214
      },
      "seconds": {
        "n": 8,
        "mean": 584.4164631962776,
        "median": 591.3809169530869,
        "min": 399.88133096694946,
        "max": 719.1918241977692
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 9.975,
        "median": 9.0,
        "min": 7.0,
        "max": 15.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.625,
        "atm_stable": 0.9386,
        "recon_success": 0.875,
        "eval_success": 0.75,
        "complete": 0.75,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 7.0,
        "complete_mols": 6.0,
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