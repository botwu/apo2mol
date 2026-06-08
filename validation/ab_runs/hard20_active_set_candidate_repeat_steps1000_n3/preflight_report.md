# New Method A/B Verification

Generated: 2026-06-05T01:28:26

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
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard20_active_set_candidate_repeat_steps1000_n3
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "pocket_router_random_top4": {
    "sampling": {
      "available": true,
      "result_files": 20,
      "protein_rmsd": {
        "n": 60,
        "mean": 1.588916689157486,
        "median": 1.5364712476730347,
        "min": 0.8534886837005615,
        "max": 2.9326016902923584
      },
      "protein_tmscore": {
        "n": 60,
        "mean": 0.9421994260426072,
        "median": 0.9509965830130487,
        "min": 0.8577423599293523,
        "max": 0.9664204398992703
      },
      "seconds": {
        "n": 60,
        "mean": 990.1269497712453,
        "median": 619.0901839733124,
        "min": 392.5080511569977,
        "max": 10808.329670906067
      },
      "router_selected_counts": {
        "n": 600,
        "mean": 4.0,
        "median": 4.0,
        "min": 4.0,
        "max": 4.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 20.0,
        "mol_stable": 0.2833,
        "atm_stable": 0.9191,
        "recon_success": 0.9,
        "eval_success": 0.7833,
        "complete": 0.7833,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 54.0,
        "complete_mols": 47.0,
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
      "result_files": 20,
      "protein_rmsd": {
        "n": 60,
        "mean": 1.5886415849129358,
        "median": 1.5545235872268677,
        "min": 0.8731919527053833,
        "max": 3.00980544090271
      },
      "protein_tmscore": {
        "n": 60,
        "mean": 0.941856615299649,
        "median": 0.9488639641324748,
        "min": 0.8575020465913779,
        "max": 0.968493311538309
      },
      "seconds": {
        "n": 60,
        "mean": 664.3423658053081,
        "median": 587.5606954097748,
        "min": 388.13611912727356,
        "max": 1164.3925099372864
      },
      "router_selected_counts": {
        "n": 600,
        "mean": 4.0,
        "median": 4.0,
        "min": 4.0,
        "max": 4.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 20.0,
        "mol_stable": 0.4667,
        "atm_stable": 0.9188,
        "recon_success": 0.8833,
        "eval_success": 0.7667,
        "complete": 0.7667,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 53.0,
        "complete_mols": 46.0,
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
      "result_files": 20,
      "protein_rmsd": {
        "n": 60,
        "mean": 1.5935474306344986,
        "median": 1.5626232028007507,
        "min": 0.8745323419570923,
        "max": 2.994328498840332
      },
      "protein_tmscore": {
        "n": 60,
        "mean": 0.9416859973405435,
        "median": 0.9489476638730041,
        "min": 0.857395083990821,
        "max": 0.9674596641269432
      },
      "seconds": {
        "n": 60,
        "mean": 661.8280358393987,
        "median": 588.2395910024643,
        "min": 386.4373559951782,
        "max": 1157.6435799598694
      },
      "router_selected_counts": {
        "n": 600,
        "mean": 5.6466666666666665,
        "median": 5.0,
        "min": 4.0,
        "max": 9.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 20.0,
        "mol_stable": 0.45,
        "atm_stable": 0.9199,
        "recon_success": 0.9167,
        "eval_success": 0.85,
        "complete": 0.85,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 55.0,
        "complete_mols": 51.0,
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
      "result_files": 20,
      "protein_rmsd": {
        "n": 60,
        "mean": 1.5896905839443207,
        "median": 1.5549184679985046,
        "min": 0.8731919527053833,
        "max": 3.012162208557129
      },
      "protein_tmscore": {
        "n": 60,
        "mean": 0.9418008660351667,
        "median": 0.9488644563514768,
        "min": 0.8575248088773721,
        "max": 0.9685272371708439
      },
      "seconds": {
        "n": 60,
        "mean": 1304.0078563014665,
        "median": 641.158271074295,
        "min": 383.6105341911316,
        "max": 19972.088862895966
      },
      "router_selected_counts": {
        "n": 600,
        "mean": 4.145,
        "median": 4.0,
        "min": 4.0,
        "max": 6.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 20.0,
        "mol_stable": 0.4667,
        "atm_stable": 0.9171,
        "recon_success": 0.8833,
        "eval_success": 0.7667,
        "complete": 0.7667,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 53.0,
        "complete_mols": 46.0,
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
      "result_files": 20,
      "protein_rmsd": {
        "n": 60,
        "mean": 1.5891291946172714,
        "median": 1.5422493815422058,
        "min": 0.8731919527053833,
        "max": 2.996370792388916
      },
      "protein_tmscore": {
        "n": 60,
        "mean": 0.9418353206446155,
        "median": 0.9487052850316998,
        "min": 0.8575260174943276,
        "max": 0.9685341307354457
      },
      "seconds": {
        "n": 60,
        "mean": 1363.3411548256875,
        "median": 646.6227629184723,
        "min": 406.38167667388916,
        "max": 32299.407994031906
      },
      "router_selected_counts": {
        "n": 600,
        "mean": 4.155,
        "median": 4.0,
        "min": 4.0,
        "max": 7.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 20.0,
        "mol_stable": 0.4667,
        "atm_stable": 0.916,
        "recon_success": 0.8833,
        "eval_success": 0.75,
        "complete": 0.75,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 53.0,
        "complete_mols": 45.0,
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
      "result_files": 20,
      "protein_rmsd": {
        "n": 60,
        "mean": 1.5896850854158402,
        "median": 1.562321424484253,
        "min": 0.873341977596283,
        "max": 2.989100933074951
      },
      "protein_tmscore": {
        "n": 60,
        "mean": 0.9418343419897905,
        "median": 0.9486315189019665,
        "min": 0.8574430257967203,
        "max": 0.968302925226047
      },
      "seconds": {
        "n": 60,
        "mean": 836.9177420973778,
        "median": 611.5994780063629,
        "min": 398.2226791381836,
        "max": 9641.820111989975
      },
      "router_selected_counts": {
        "n": 600,
        "mean": 5.67,
        "median": 5.0,
        "min": 4.0,
        "max": 9.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 20.0,
        "mol_stable": 0.45,
        "atm_stable": 0.9166,
        "recon_success": 0.8833,
        "eval_success": 0.8167,
        "complete": 0.8167,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 53.0,
        "complete_mols": 49.0,
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