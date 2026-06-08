# New Method A/B Verification

Generated: 2026-05-28T19:33:38

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
| control_realistic_late_dense | atoms=prior, center=apo, steps=1000 | late_dense | none topk=0 |
| pocket_router_random_top4 | atoms=prior, center=apo, steps=1000 | late_dense | random topk=4 |
| pocket_router_distance_top4_hard | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell4_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell6_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell6_w050 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_random_top4_shell6_w025 | atoms=prior, center=apo, steps=1000 | late_dense | random topk=4 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| control_realistic_late_dense | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_random_top4 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| pocket_router_distance_top4_hard | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell4_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell6_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell6_w050 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_random_top4_shell6_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard2_active_set_shell_pilot_steps1000_n1
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
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
        "mean": 530.9685069322586,
        "median": 530.9685069322586,
        "min": 402.50037479400635,
        "max": 659.4366390705109
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
  "pocket_router_random_top4": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.1251946687698364,
        "median": 2.1251946687698364,
        "min": 2.09743332862854,
        "max": 2.152956008911133
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8794582063547877,
        "median": 0.8794582063547877,
        "min": 0.8577423599293523,
        "max": 0.901174052780223
      },
      "seconds": {
        "n": 2,
        "mean": 528.3843539953232,
        "median": 528.3843539953232,
        "min": 408.5506739616394,
        "max": 648.218034029007
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
  "pocket_router_distance_top4_hard": {
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
        "mean": 527.7300126552582,
        "median": 527.7300126552582,
        "min": 401.5647852420807,
        "max": 653.8952400684357
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
  "active_set_distance_top4_shell4_w025": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.0641393661499023,
        "median": 2.0641393661499023,
        "min": 2.044558525085449,
        "max": 2.0837202072143555
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8795998450553063,
        "median": 0.8795998450553063,
        "min": 0.857395083990821,
        "max": 0.9018046061197916
      },
      "seconds": {
        "n": 2,
        "mean": 524.1991533041,
        "median": 524.1991533041,
        "min": 401.8937678337097,
        "max": 646.5045387744904
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 5.95,
        "median": 6.0,
        "min": 4.0,
        "max": 7.0
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
  "active_set_distance_top4_shell6_w025": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.0806291103363037,
        "median": 2.0806291103363037,
        "min": 2.0641345977783203,
        "max": 2.097123622894287
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8793813875406096,
        "median": 0.8793813875406096,
        "min": 0.857334552424969,
        "max": 0.90142822265625
      },
      "seconds": {
        "n": 2,
        "mean": 519.5328899621964,
        "median": 519.5328899621964,
        "min": 397.08777499198914,
        "max": 641.9780049324036
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 17.8,
        "median": 18.5,
        "min": 14.0,
        "max": 22.0
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
  "active_set_distance_top4_shell6_w050": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.0950427055358887,
        "median": 2.0950427055358887,
        "min": 2.0742950439453125,
        "max": 2.115790367126465
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8802055454994464,
        "median": 0.8802055454994464,
        "min": 0.8595440049376032,
        "max": 0.9008670860612896
      },
      "seconds": {
        "n": 2,
        "mean": 521.6375340223312,
        "median": 521.6375340223312,
        "min": 402.28830194473267,
        "max": 640.9867660999298
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 17.5,
        "median": 18.0,
        "min": 13.0,
        "max": 22.0
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
  "active_set_random_top4_shell6_w025": {
    "sampling": {
      "available": true,
      "result_files": 2,
      "protein_rmsd": {
        "n": 2,
        "mean": 2.1625306606292725,
        "median": 2.1625306606292725,
        "min": 2.0899176597595215,
        "max": 2.2351436614990234
      },
      "protein_tmscore": {
        "n": 2,
        "mean": 0.8786022048146578,
        "median": 0.8786022048146578,
        "min": 0.8558750467331889,
        "max": 0.9013293628961268
      },
      "seconds": {
        "n": 2,
        "mean": 519.9250574111938,
        "median": 519.9250574111938,
        "min": 397.8040518760681,
        "max": 642.0460629463196
      },
      "router_selected_counts": {
        "n": 20,
        "mean": 15.5,
        "median": 15.5,
        "min": 10.0,
        "max": 21.0
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
  }
}
```