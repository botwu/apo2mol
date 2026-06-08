# PocketRouter Model Validation: Hard1 Pilot

Date: 2026-05-26

## Setup

- Run dir: `validation/ab_runs/hard1_steps1000_ab_gated`
- Data: downloaded Apo2Mol data under `apo2mol_dataset/data_folder`
- Case: test position 477, original index 24500, metadata apo/holo RMSD 4.0188
- Sampling: realistic apo initialization, `atoms=prior`, 1000 reverse steps
- Protein updates:
  - `static5`: original-style sparse-in-time global pocket update
  - `late_dense`: 10 update opportunities at timesteps 460, 410, 360, 310, 260,
    210, 160, 110, 60, 10
  - router arms: same late-dense schedule, but only top-12 residues receive the
    high-fidelity rigid/chi update

## Results

| arm | protein RMSD | TM-score | validity/recon | router usage |
|---|---:|---:|---|---:|
| baseline_realistic_static5 | 3.1261 | 0.8148 | recon 1.0, complete 1.0 | n/a |
| control_realistic_late_dense | 3.1191 | 0.8130 | recon 1.0, complete 1.0 | n/a |
| pocket_router_distance_top12 | 2.2258 | 0.8517 | recon 1.0, complete 1.0 | 10 x 12 |
| pocket_router_motion_oracle_top12 | 3.0552 | 0.8199 | recon 1.0, complete 1.0 | 10 x 12 |
| pocket_router_random_top12 | 3.0563 | 0.8160 | recon 1.0, complete 1.0 | 10 x 12 |

Best arm: `pocket_router_distance_top12`

Distance-router deltas:

- vs static5: -0.9003 A protein RMSD
- vs late_dense: -0.8933 A protein RMSD

## Interpretation

This is a positive pilot for sparse pocket updates, but it narrows the idea.

What the result supports:

- The top-conference-worthy point should be a fragment/contact-conditioned
  sparse pocket workspace.
- Updating only a small contact-relevant residue set can beat dense global
  pocket updates under the same generation schedule.

What the result does not support yet:

- It does not support a broad claim that "mobile residue routing" is enough.
  Motion-oracle top12 and random top12 are almost tied in this model-level pilot.
- It does not support a top-conference claim by itself because this is n=1.

Working hypothesis:

> During ligand generation, the immediate ligand-contact subpocket is the useful
> sparse memory for stabilizing induced-fit denoising. Large apo-to-holo mobile
> residues may matter in the final structure, but hard-selecting them without
> current-fragment contact context is not the right computation.

## Next Decisive Tests

1. Run the same router A/B on the hard8 set already used for data-level
   validation.
2. Add `contact_change_oracle_top12` to the model-level arms.
3. Add a learned router target:
   current-fragment distance/contact + apo-to-holo contact change + residue
   displacement, with budgeted top-k selection.
4. Require learned routing to beat:
   `late_dense`, `random_top12`, `distance_top12`, and `motion_oracle_top12`.
5. If learned routing cannot beat distance-only, the defensible paper claim
   becomes an engineered sparse induced-fit baseline, not a main top-conference
   method.
