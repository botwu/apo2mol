# New Method A/B Verification

Generated: 2026-05-29T11:14:34

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
| active_set_distance_top4_shell4_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell6_w025 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_distance_top4_shell6_w050 | atoms=prior, center=apo, steps=1000 | late_dense | distance topk=4 |
| active_set_random_top4_shell6_w025 | atoms=prior, center=apo, steps=1000 | late_dense | random topk=4 |

## Sampling Diagnostics

The current sampler is an adjacent reverse chain. If `num_steps < 1000`, it samples only the high-noise tail and does not reach t=0; use those runs only as pipeline smoke tests.

| arm | timestep range | reaches t=0 | candidate update hits |
|---|---:|---:|---|
| active_set_distance_top4_shell4_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell6_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_distance_top4_shell6_w050 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |
| active_set_random_top4_shell6_w025 | 999 to 0 | True | 460, 410, 360, 310, 260, 210, 160, 110, 60, 10 |

## Run Command

```bash
/Users/wujian1/Downloads/a800_molecular/Apo2Mol/.venv310/bin/python validation/run_new_method_ab.py --run --run-dir /Users/wujian1/Downloads/a800_molecular/Apo2Mol/validation/ab_runs/hard8_active_set_shell_only_steps1000_n1
```

A positive result requires the adaptive arm to reduce mean/best protein RMSD and not degrade validity/reconstruction metrics under the same hard cases and realistic protocol. For top-conference evidence, run with --include-controls and require the adaptive arm to beat the late-dense and uniform schedule controls too.

## Results

```json
{
  "active_set_distance_top4_shell4_w025": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.858922652900219,
        "median": 1.8898723125457764,
        "min": 0.9495697617530823,
        "max": 3.056523561477661
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9213097352357439,
        "median": 0.934025079001259,
        "min": 0.857395083990821,
        "max": 0.9591475132533482
      },
      "seconds": {
        "n": 8,
        "mean": 581.8449841141701,
        "median": 589.8475514650345,
        "min": 397.17040491104126,
        "max": 717.7487111091614
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 6.1125,
        "median": 6.0,
        "min": 4.0,
        "max": 9.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.875,
        "atm_stable": 0.9675,
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
  "active_set_distance_top4_shell6_w025": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.8666208386421204,
        "median": 1.883530616760254,
        "min": 0.9549903869628906,
        "max": 3.042263984680176
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9217480835487953,
        "median": 0.9344082546222862,
        "min": 0.857334552424969,
        "max": 0.9590478515625
      },
      "seconds": {
        "n": 8,
        "mean": 576.6254455447197,
        "median": 581.6412144899368,
        "min": 395.4603250026703,
        "max": 711.0433511734009
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 16.375,
        "median": 16.0,
        "min": 9.0,
        "max": 22.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.375,
        "atm_stable": 0.9747,
        "recon_success": 1.0,
        "eval_success": 1.0,
        "complete": 1.0,
        "JSD_CC_2A": null,
        "JSD_All_12A": null,
        "Atom_type_JS": null,
        "Number_of_reconstructed_mols": 8.0,
        "complete_mols": 8.0,
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
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.9300880655646324,
        "median": 1.938381850719452,
        "min": 0.9673857092857361,
        "max": 3.2678823471069336
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9192357719195905,
        "median": 0.9317049492560262,
        "min": 0.8595440049376032,
        "max": 0.9587510463169643
      },
      "seconds": {
        "n": 8,
        "mean": 589.460047096014,
        "median": 597.2239919900894,
        "min": 396.1496031284332,
        "max": 733.472864151001
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 15.5125,
        "median": 15.0,
        "min": 9.0,
        "max": 22.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.625,
        "atm_stable": 0.9892,
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
  "active_set_random_top4_shell6_w025": {
    "sampling": {
      "available": true,
      "result_files": 8,
      "protein_rmsd": {
        "n": 8,
        "mean": 1.8523517921566963,
        "median": 1.9628440737724304,
        "min": 0.9735156893730164,
        "max": 2.9245963096618652
      },
      "protein_tmscore": {
        "n": 8,
        "mean": 0.9248570714955335,
        "median": 0.9375888702681434,
        "min": 0.8558750467331889,
        "max": 0.9587151227678572
      },
      "seconds": {
        "n": 8,
        "mean": 589.1921036839485,
        "median": 592.258754491806,
        "min": 404.97165298461914,
        "max": 726.9904141426086
      },
      "router_selected_counts": {
        "n": 80,
        "mean": 15.75,
        "median": 16.0,
        "min": 8.0,
        "max": 22.0
      }
    },
    "eval": {
      "available": true,
      "metrics": {
        "Number_of_generated_data": 8.0,
        "mol_stable": 0.25,
        "atm_stable": 0.8556,
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
  }
}
```