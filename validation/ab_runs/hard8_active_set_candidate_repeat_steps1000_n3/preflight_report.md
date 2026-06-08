# New Method A/B Verification

Generated: 2026-05-30T01:10:15

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
| pocket_router_random_top4 | atoms=prior, center=apo, steps=1000 | late_dense | random topk=4 |
| pocket_router_distance_top4_hard | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell4_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell3_w010 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell3_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell4_w010 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| pocket_router_random_top4 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top4_hard | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell4_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell3_w010 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell3_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell4_w010 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "pocket_router_random_top4": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 24,
        "mean": 1.8152538264791171,
        "median": 1.812263011932373,
        "min": 0.9890567660331726,
        "max": 2.9326016902923584
      },
      "protein_tmscore": {
        "n": 24,
        "mean": 0.9252509293625445,
        "median": 0.9395887560733349,
        "min": 0.8577423599293523,
        "max": 0.9586271075172719
      },
      "seconds": {
        "n": 24,
        "mean": 605.3360419074694,
        "median": 603.7232115268707,
        "min": 392.6833908557892,
        "max": 762.0118458271027
      },
      "router_selected_counts": {
        "n": 240,
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
        "mol_stable": 0.25,
        "atm_stable": 0.9094,
        "recon_success": 0.875,
        "eval_success": 0.7083,
        "complete": 0.7083,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 21.0,
        "complete_mols": 17.0,
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
  "pocket_router_distance_top4_hard": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 24,
        "mean": 1.8147885724902153,
        "median": 1.7952041625976562,
        "min": 0.9508947730064392,
        "max": 3.00980544090271
      },
      "protein_tmscore": {
        "n": 24,
        "mean": 0.9244962411216434,
        "median": 0.9361613542171909,
        "min": 0.8575020465913779,
        "max": 0.9587510463169643
      },
      "seconds": {
        "n": 24,
        "mean": 609.9411842425665,
        "median": 600.6643295288086,
        "min": 408.4447190761566,
        "max": 762.3994739055634
      },
      "router_selected_counts": {
        "n": 240,
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
        "mol_stable": 0.5,
        "atm_stable": 0.9108,
        "recon_success": 0.8333,
        "eval_success": 0.6667,
        "complete": 0.6667,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 20.0,
        "complete_mols": 16.0,
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
  "active_set_distance_top4_shell4_w025": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 24,
        "mean": 1.8206664323806763,
        "median": 1.799104928970337,
        "min": 0.9528805017471313,
        "max": 2.994328498840332
      },
      "protein_tmscore": {
        "n": 24,
        "mean": 0.9243073187772737,
        "median": 0.9361849698153409,
        "min": 0.857395083990821,
        "max": 0.9587244524274554
      },
      "seconds": {
        "n": 24,
        "mean": 1203.465994944175,
        "median": 610.5417084693909,
        "min": 387.2918601036072,
        "max": 15088.137361764908
      },
      "router_selected_counts": {
        "n": 240,
        "mean": 6.083333333333333,
        "median": 6.0,
        "min": 4.0,
        "max": 9.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.5,
        "atm_stable": 0.9173,
        "recon_success": 0.875,
        "eval_success": 0.75,
        "complete": 0.75,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 21.0,
        "complete_mols": 18.0,
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
  "active_set_distance_top4_shell3_w010": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 24,
        "mean": 1.8171459510922432,
        "median": 1.8001397848129272,
        "min": 0.9508947730064392,
        "max": 3.012162208557129
      },
      "protein_tmscore": {
        "n": 24,
        "mean": 0.9243609119031205,
        "median": 0.935780358155921,
        "min": 0.8575248088773721,
        "max": 0.9587510463169643
      },
      "seconds": {
        "n": 24,
        "mean": 603.3906157314777,
        "median": 600.5180221796036,
        "min": 406.90928530693054,
        "max": 758.6079869270325
      },
      "router_selected_counts": {
        "n": 240,
        "mean": 4.245833333333334,
        "median": 4.0,
        "min": 4.0,
        "max": 6.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.5,
        "atm_stable": 0.9147,
        "recon_success": 0.8333,
        "eval_success": 0.6667,
        "complete": 0.6667,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 20.0,
        "complete_mols": 16.0,
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
        "n": 24,
        "mean": 1.8144892180959384,
        "median": 1.8033769726753235,
        "min": 0.9508947730064392,
        "max": 2.996370792388916
      },
      "protein_tmscore": {
        "n": 24,
        "mean": 0.9245075731592048,
        "median": 0.9355151394254089,
        "min": 0.8575260174943276,
        "max": 0.9587510463169643
      },
      "seconds": {
        "n": 24,
        "mean": 2284.6595355570316,
        "median": 699.2717323303223,
        "min": 406.0125222206116,
        "max": 12768.094201087952
      },
      "router_selected_counts": {
        "n": 240,
        "mean": 4.2625,
        "median": 4.0,
        "min": 4.0,
        "max": 7.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.4583,
        "atm_stable": 0.9081,
        "recon_success": 0.8333,
        "eval_success": 0.625,
        "complete": 0.625,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 20.0,
        "complete_mols": 15.0,
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
        "n": 24,
        "mean": 1.8168433333436649,
        "median": 1.7945468425750732,
        "min": 0.9517804980278015,
        "max": 2.989100933074951
      },
      "protein_tmscore": {
        "n": 24,
        "mean": 0.9243982187524636,
        "median": 0.9361958302838311,
        "min": 0.8574430257967203,
        "max": 0.9587398856026785
      },
      "seconds": {
        "n": 24,
        "mean": 586.8278813163439,
        "median": 575.8955800533295,
        "min": 389.20290207862854,
        "max": 738.6122798919678
      },
      "router_selected_counts": {
        "n": 240,
        "mean": 6.079166666666667,
        "median": 6.0,
        "min": 4.0,
        "max": 9.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.4583,
        "atm_stable": 0.9094,
        "recon_success": 0.8333,
        "eval_success": 0.75,
        "complete": 0.75,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 20.0,
        "complete_mols": 18.0,
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