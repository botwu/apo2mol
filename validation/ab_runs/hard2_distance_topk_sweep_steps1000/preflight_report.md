# New Method A/B Verification

Generated: 2026-05-27T18:23:13

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
| pocket_router_distance_top4 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| pocket_router_distance_top8 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=8 |
| pocket_router_distance_top12 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=12 |
| pocket_router_distance_top16 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=16 |
| pocket_router_distance_top24 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=24 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| pocket_router_distance_top4 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top8 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top12 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top16 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top24 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard2_distance_topk_sweep_steps1000
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "pocket_router_distance_top4": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.062439441680908,
        "median": 2.062439441680908,
        "min": 2.0402297973632812,
        "max": 2.084649085998535
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8796419001731648,
        "median": 0.8796419001731648,
        "min": 0.8575020465913779,
        "max": 0.9017817537549516
      },
      "seconds": {
        "n": 2,
        "mean": 525.8132960796356,
        "median": 525.8132960796356,
        "min": 399.6281011104584,
        "max": 651.9984910488129
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 4.0,
        "median": 4.0,
        "min": 4.0,
        "max": 4.0
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
  "pocket_router_distance_top8": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.1049686670303345,
        "median": 2.1049686670303345,
        "min": 2.1046462059020996,
        "max": 2.1052911281585693
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8783828245237527,
        "median": 0.8783828245237527,
        "min": 0.8555235406353136,
        "max": 0.9012421084121919
      },
      "seconds": {
        "n": 2,
        "mean": 529.0106920003891,
        "median": 529.0106920003891,
        "min": 404.1580250263214,
        "max": 653.8633589744568
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 8.0,
        "median": 8.0,
        "min": 8.0,
        "max": 8.0
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
  "pocket_router_distance_top16": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.240739107131958,
        "median": 2.240739107131958,
        "min": 2.1332733631134033,
        "max": 2.3482048511505127
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8744937171664735,
        "median": 0.8744937171664735,
        "min": 0.8482806013755673,
        "max": 0.9007068329573797
      },
      "seconds": {
        "n": 2,
        "mean": 525.3131849765778,
        "median": 525.3131849765778,
        "min": 402.48929810523987,
        "max": 648.1370718479156
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 16.0,
        "median": 16.0,
        "min": 16.0,
        "max": 16.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 2.0,
        "mol_stable": 0.5,
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
  "pocket_router_distance_top24": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.3320716619491577,
        "median": 2.3320716619491577,
        "min": 2.1631128787994385,
        "max": 2.501030445098877
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8695427022130786,
        "median": 0.8695427022130786,
        "min": 0.8391372126714625,
        "max": 0.8999481917546949
      },
      "seconds": {
        "n": 2,
        "mean": 523.1463080644608,
        "median": 523.1463080644608,
        "min": 400.14674615859985,
        "max": 646.1458699703217
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 24.0,
        "median": 24.0,
        "min": 24.0,
        "max": 24.0
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
  "distance_topk_comparison": {
    "best_topk": 4,
    "best_protein_rmsd": 2.062439441680908,
    "protein_rmsd_by_topk": {
      "4": 2.062439441680908,
      "8": 2.1049686670303345,
      "12": 2.1731373071670532,
      "16": 2.240739107131958,
      "24": 2.3320716619491577
    }
  }
}
```