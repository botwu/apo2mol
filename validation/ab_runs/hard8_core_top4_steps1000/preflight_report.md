# New Method A/B Verification

Generated: 2026-05-28T01:11:35

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
| control_realistic_late_dense | atoms=prior, center=apo, steps=1000 | late_dense | none topk=0 |
| pocket_router_random_top4 | atoms=prior, center=apo, steps=1000 | late_dense | random topk=4 |
| pocket_router_distance_top4 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| pocket_router_distance_top8 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=8 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| baseline_realistic_static5 | 999 to 0 | True | 799, 599, 399, 199, 10 |
| control_realistic_late_dense | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_random_top4 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top4 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top8 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard8_core_top4_steps1000
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
        "n": 8,
        "mean": 2.4365141838788986,
        "median": 2.411094069480896,
        "min": 1.8510633707046509,
        "max": 3.1260809898376465
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.8976762109172821,
        "median": 0.9075472955389338,
        "min": 0.8148447168935643,
        "max": 0.9327622440602908
      },
      "seconds": {
        "n": 8,
        "mean": 619.6938094496727,
        "median": 611.3333287239075,
        "min": 408.71689105033875,
        "max": 840.0011219978333
      },
      "router_selected_counts": {
        "n": 40,
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
        "mol_stable": 0.25,
        "atm_stable": 0.875,
        "recon_success": 0.875,
        "eval_success": 0.875,
        "complete": 0.875,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 7.0,
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
  "control_realistic_late_dense": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 2.4142273515462875,
        "median": 2.359957456588745,
        "min": 1.6291425228118896,
        "max": 3.3319222927093506
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9005225378865523,
        "median": 0.9096084765826067,
        "min": 0.8129701519956684,
        "max": 0.9464549086700674
      },
      "seconds": {
        "n": 8,
        "mean": 616.2993023097515,
        "median": 629.5966629981995,
        "min": 423.6177988052368,
        "max": 773.774124622345
      },
      "router_selected_counts": {
        "n": 80,
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
        "mol_stable": 0.25,
        "atm_stable": 0.8583,
        "recon_success": 0.875,
        "eval_success": 0.875,
        "complete": 0.875,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 7.0,
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
  "pocket_router_random_top4": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.8305286318063736,
        "median": 1.9381644129753113,
        "min": 0.9700620174407959,
        "max": 2.882277727127075
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9247954826069478,
        "median": 0.9364156793735742,
        "min": 0.8577423599293523,
        "max": 0.9587002999441965
      },
      "seconds": {
        "n": 8,
        "mean": 407.8374507725239,
        "median": 392.78439462184906,
        "min": 325.75662112236023,
        "max": 497.8095841407776
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 4.0,
        "median": 4.0,
        "min": 4.0,
        "max": 4.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.375,
        "atm_stable": 0.8881,
        "recon_success": 0.875,
        "eval_success": 0.625,
        "complete": 0.625,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 7.0,
        "complete_mols": 5.0,
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
  "pocket_router_distance_top4": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.792432278394699,
        "median": 1.8086740374565125,
        "min": 1.108406901359558,
        "max": 2.9897732734680176
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9249475269868198,
        "median": 0.9401615720233019,
        "min": 0.8575020465913779,
        "max": 0.9597950128785748
      },
      "seconds": {
        "n": 8,
        "mean": 528.5142754912376,
        "median": 501.4253463745117,
        "min": 399.6281011104584,
        "max": 651.9984910488129
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 4.0,
        "median": 4.0,
        "min": 4.0,
        "max": 4.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.375,
        "atm_stable": 0.9292,
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
  "pocket_router_distance_top8": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.8370685875415802,
        "median": 1.8211688995361328,
        "min": 1.1902565956115723,
        "max": 3.0618085861206055
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9229069043429472,
        "median": 0.9384607420764353,
        "min": 0.8555235406353136,
        "max": 0.9587935812982005
      },
      "seconds": {
        "n": 8,
        "mean": 527.920348316431,
        "median": 500.40989100933075,
        "min": 404.1580250263214,
        "max": 653.8633589744568
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 8.0,
        "median": 8.0,
        "min": 8.0,
        "max": 8.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.25,
        "atm_stable": 0.8417,
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
  "distance_topk_comparison": {
    "best_topk": 4,
    "best_protein_rmsd": 1.792432278394699,
    "protein_rmsd_by_topk": {
      "4": 1.792432278394699,
      "8": 1.8370685875415802
    }
  },
  "hardtail_core_comparison": {
    "best_arm": "pocket_router_distance_top4",
    "best_protein_rmsd": 1.792432278394699,
    "protein_rmsd_by_arm": {
      "baseline_realistic_static5": 2.4365141838788986,
      "control_realistic_late_dense": 2.4142273515462875,
      "pocket_router_distance_top4": 1.792432278394699,
      "pocket_router_distance_top8": 1.8370685875415802,
      "pocket_router_random_top4": 1.8305286318063736
    },
    "delta_vs_baseline_static5": {
      "baseline_realistic_static5": 0.0,
      "control_realistic_late_dense": -0.022286832332611084,
      "pocket_router_distance_top4": -0.6440819054841995,
      "pocket_router_distance_top8": -0.5994455963373184,
      "pocket_router_random_top4": -0.605985552072525
    },
    "delta_vs_late_dense": {
      "baseline_realistic_static5": 0.022286832332611084,
      "control_realistic_late_dense": 0.0,
      "pocket_router_distance_top4": -0.6217950731515884,
      "pocket_router_distance_top8": -0.5771587640047073,
      "pocket_router_random_top4": -0.5836987197399139
    }
  }
}
```