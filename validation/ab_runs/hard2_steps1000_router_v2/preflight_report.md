# New Method A/B Verification

Generated: 2026-05-27T13:52:05

## Status

Runnable: all local prerequisites were found.

Selected-case files checked: 6
Missing selected-case files: 0

## Selected Hard Test Cases

| test position | original index | apo/holo RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |

## Experiment Arms

| arm | sampling protocol | protein update | pocket router |
|---|---|---|---|
| baseline_realistic_static5 | atoms=prior, center=apo, steps=1000 | static5 | none topk=0 |
| control_realistic_late_dense | atoms=prior, center=apo, steps=1000 | late_dense | none topk=0 |
| pocket_router_distance_top12 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=12 |
| pocket_router_motion_oracle_top12 | atoms=prior, center=apo, steps=1000 | late_dense | motion_oracle topk=12 |
| pocket_router_contact_change_oracle_top12 | atoms=prior, center=apo, steps=1000 | late_dense | contact_change_oracle topk=12 |
| pocket_router_random_top12 | atoms=prior, center=apo, steps=1000 | late_dense | random topk=12 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| baseline_realistic_static5 | 999 to 0 | True | 799, 599, 399, 199, 10 |
| control_realistic_late_dense | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_motion_oracle_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_contact_change_oracle_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_random_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard2_steps1000_router_v2
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "baseline_realistic_static5": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.8380913734436035,
        "median": 2.8380913734436035,
        "min": 2.5501017570495605,
        "max": 3.1260809898376465
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8483045800806217,
        "median": 0.8483045800806217,
        "min": 0.8148447168935643,
        "max": 0.881764443267679
      },
      "seconds": {
        "n": 2,
        "mean": 535.1859279870987,
        "median": 535.1859279870987,
        "min": 408.71689105033875,
        "max": 661.6549649238586
      },
      "router_selected_counts": {
        "n": 10,
        "mean": 48.0,
        "median": 48.0,
        "min": 42.0,
        "max": 54.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 1.0,
        "mol_stable": 0.0,
        "atm_stable": 0.9706,
        "recon_success": 1.0,
        "eval_success": 1.0,
        "complete": 1.0,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 1.0,
        "complete_mols": 1.0,
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
  "control_realistic_late_dense": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.786664366722107,
        "median": 2.786664366722107,
        "min": 2.4542131423950195,
        "max": 3.1191155910491943
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8495568264058815,
        "median": 0.8495568264058815,
        "min": 0.8129701519956684,
        "max": 0.8861435008160945
      },
      "seconds": {
        "n": 2,
        "mean": 551.2952848672867,
        "median": 551.2952848672867,
        "min": 423.6177988052368,
        "max": 678.9727709293365
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 48.0,
        "median": 48.0,
        "min": 42.0,
        "max": 54.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 2.0,
        "mol_stable": 0.5,
        "atm_stable": 0.9844,
        "recon_success": 1.0,
        "eval_success": 1.0,
        "complete": 1.0,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 2.0,
        "complete_mols": 2.0,
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
  "pocket_router_distance_top12": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.1731373071670532,
        "median": 2.1731373071670532,
        "min": 2.1204867362976074,
        "max": 2.225787878036499
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8763723256824816,
        "median": 0.8763723256824816,
        "min": 0.8517396631020524,
        "max": 0.9010049882629108
      },
      "seconds": {
        "n": 2,
        "mean": 555.9710590839386,
        "median": 555.9710590839386,
        "min": 424.3759422302246,
        "max": 687.5661759376526
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 12.0,
        "median": 12.0,
        "min": 12.0,
        "max": 12.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 2.0,
        "mol_stable": 1.0,
        "atm_stable": 1.0,
        "recon_success": 1.0,
        "eval_success": 1.0,
        "complete": 1.0,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 2.0,
        "complete_mols": 2.0,
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
  "pocket_router_motion_oracle_top12": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.6050139665603638,
        "median": 2.6050139665603638,
        "min": 2.1547818183898926,
        "max": 3.055246114730835
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8595585175617377,
        "median": 0.8595585175617377,
        "min": 0.8198989012060386,
        "max": 0.8992181339174369
      },
      "seconds": {
        "n": 2,
        "mean": 550.041466832161,
        "median": 550.041466832161,
        "min": 414.3159077167511,
        "max": 685.7670259475708
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 12.0,
        "median": 12.0,
        "min": 12.0,
        "max": 12.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 2.0,
        "mol_stable": 0.0,
        "atm_stable": 0.9688,
        "recon_success": 1.0,
        "eval_success": 1.0,
        "complete": 1.0,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 2.0,
        "complete_mols": 2.0,
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
  "pocket_router_contact_change_oracle_top12": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.8072317838668823,
        "median": 2.8072317838668823,
        "min": 2.21002459526062,
        "max": 3.4044389724731445
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8492519655459507,
        "median": 0.8492519655459507,
        "min": 0.8015632251701733,
        "max": 0.8969407059217283
      },
      "seconds": {
        "n": 2,
        "mean": 556.9048269987106,
        "median": 556.9048269987106,
        "min": 428.06706285476685,
        "max": 685.7425911426544
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 12.0,
        "median": 12.0,
        "min": 12.0,
        "max": 12.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 2.0,
        "mol_stable": 0.0,
        "atm_stable": 0.9219,
        "recon_success": 1.0,
        "eval_success": 1.0,
        "complete": 1.0,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 2.0,
        "complete_mols": 2.0,
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
  "pocket_router_random_top12": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.591334819793701,
        "median": 2.591334819793701,
        "min": 2.1263606548309326,
        "max": 3.0563089847564697
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8573714191431625,
        "median": 0.8573714191431625,
        "min": 0.8159647019389439,
        "max": 0.8987781363473811
      },
      "seconds": {
        "n": 2,
        "mean": 548.7982380390167,
        "median": 548.7982380390167,
        "min": 421.7739200592041,
        "max": 675.8225560188293
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 12.0,
        "median": 12.0,
        "min": 12.0,
        "max": 12.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 2.0,
        "mol_stable": 0.5,
        "atm_stable": 0.9531,
        "recon_success": 1.0,
        "eval_success": 1.0,
        "complete": 1.0,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 2.0,
        "complete_mols": 2.0,
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
  "router_comparison": {
    "best_arm": "pocket_router_distance_top12",
    "best_protein_rmsd": 2.1731373071670532,
    "baseline_realistic_static5_delta_vs_baseline_static5": 0.0,
    "baseline_realistic_static5_delta_vs_realistic_late_dense": 0.05142700672149658,
    "control_realistic_late_dense_delta_vs_baseline_static5": -0.05142700672149658,
    "control_realistic_late_dense_delta_vs_realistic_late_dense": 0.0,
    "pocket_router_distance_top12_delta_vs_baseline_static5": -0.6649540662765503,
    "pocket_router_distance_top12_delta_vs_realistic_late_dense": -0.6135270595550537,
    "pocket_router_motion_oracle_top12_delta_vs_baseline_static5": -0.23307740688323975,
    "pocket_router_motion_oracle_top12_delta_vs_realistic_late_dense": -0.18165040016174316,
    "pocket_router_contact_change_oracle_top12_delta_vs_baseline_static5": -0.03085958957672119,
    "pocket_router_contact_change_oracle_top12_delta_vs_realistic_late_dense": 0.02056741714477539,
    "pocket_router_random_top12_delta_vs_baseline_static5": -0.24675655364990234,
    "pocket_router_random_top12_delta_vs_realistic_late_dense": -0.19532954692840576
  }
}
```