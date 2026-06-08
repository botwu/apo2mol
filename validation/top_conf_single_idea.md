# Single-Idea Top-Conference Direction

Date: 2026-05-26
Updated: 2026-05-29 with hard8 active-set conservative shell sweep

## Core Idea

Name: Fragment-Guided Contact-Sparse Pocket Memory, or PocketRouter.

One sentence:

During ligand generation, treat the apo pocket as a structured memory bank and
let the current fragment state route to only the ligand-contact-relevant residue
chunks, subpockets, and rotamer modes for high-precision ligand-pocket coupling.

This is not "add attention". The contribution is a state-dependent routing and
memory mechanism for fragment-conditioned apo-to-holo adaptation. The current
validation suggests the useful first instantiation should be contact-conditioned,
not simply "update the most mobile residues".

## Why This Is The One Point To Refine

Top-conference projects often win by making one computational principle concrete:

- Longformer/BigBird: long-context modeling does not need dense all-to-all
  attention; local/global sparse patterns are enough when designed carefully.
- Swin Transformer: vision tokens should interact hierarchically through local
  shifted windows, not globally at every layer.
- Deformable DETR: dense image attention can be replaced by sparse reference
  points around the locations that matter.
- Perceiver/TokenLearner: large inputs can be compressed into a small set of
  learned/adaptive latent tokens before expensive reasoning.
- RAG/DPR/RETRO-style models: parametric models improve when they query an
  external memory rather than storing everything in weights.
- Slot Attention: object-centric reasoning emerges from iterative binding
  between a few slots and the current observation.

The transferable principle is:

> In a large state space, do not let every element interact with every other
> element at every step. Route the current state to a sparse, relevant memory,
> then spend computation there.

For Apo2Mol, the direct mapping is:

- full pocket or all residue states = memory bank;
- current fragment / partial ligand graph / diffusion state = query;
- ligand-contact residues, contact-change residues, side-chain rotamers,
  subpockets = sparse memory;
- next fragment generation and pocket adjustment = sparse reasoning.

## Actual Validation So Far

I ran two local validation layers on the downloaded Apo2Mol data.

Data-level feasibility, 8 hardest test cases:

- Mean apo atom RMSD: 1.7027 A.
- Top-12 distance router captures holo contacts well: holo-contact recall 0.897,
  replay RMSD 1.3415 A.
- Top-12 motion oracle captures mobile residues: mobile recall 1.000, replay RMSD
  0.6426 A.
- Top-12 random replay RMSD is 1.3115 A, close to distance replay, so static
  geometric replay alone is not enough evidence for a method claim.

Model-level 1000-step A/B, one hardest case, realistic apo initialization:

| arm | protein RMSD | TM-score | router selected counts |
|---|---:|---:|---:|
| baseline_realistic_static5 | 3.1261 | 0.8148 | n/a |
| control_realistic_late_dense | 3.1191 | 0.8130 | n/a |
| pocket_router_distance_top12 | 2.2258 | 0.8517 | 10 x 12 |
| pocket_router_motion_oracle_top12 | 3.0552 | 0.8199 | 10 x 12 |
| pocket_router_random_top12 | 3.0563 | 0.8160 | 10 x 12 |

Interpretation:

- There is a real positive signal: distance top-12 sparse residue updates beat
  late-dense by 0.8933 A on this hard case.
- The signal is not "select the most mobile residues". Motion-oracle and random
  are almost tied in the model run.
- The refined hypothesis is therefore: the current fragment needs a
  contact-conditioned sparse pocket workspace. Large apo-to-holo mobile regions
  may be important globally, but hard-masking them without ligand contact context
  does not directly improve generation.
- This is not yet top-conference evidence. It is a promising pilot result that
  must be repeated on a hard-tail set and with learned/contact-change routers.

Model-level 1000-step A/B, two hardest cases, realistic apo initialization
follow-up:

| arm | mean protein RMSD | mean TM-score | delta vs static5 | delta vs late-dense |
|---|---:|---:|---:|---:|
| baseline_realistic_static5 | 2.8381 | 0.8483 | 0.0000 | +0.0514 |
| control_realistic_late_dense | 2.7867 | 0.8496 | -0.0514 | 0.0000 |
| pocket_router_distance_top12 | 2.1731 | 0.8764 | -0.6650 | -0.6135 |
| pocket_router_motion_oracle_top12 | 2.6050 | 0.8596 | -0.2331 | -0.1817 |
| pocket_router_contact_change_oracle_top12 | 2.8072 | 0.8493 | -0.0309 | +0.0206 |
| pocket_router_random_top12 | 2.5913 | 0.8574 | -0.2468 | -0.1953 |

Updated interpretation:

- Distance/contact-proximity routing remains the strongest signal after adding
  the second hard case.
- Random sparse routing improves over static5, so sparse restriction itself may
  regularize pocket denoising, but it remains clearly worse than distance
  routing.
- Motion oracle and contact-change oracle are not reliable main mechanisms in
  the current untrained masked sampler.
- The top-conference claim should therefore be tightened: the learned method
  must predict current-fragment contact relevance and must beat distance-only
  routing. Otherwise, distance routing is only a strong engineered baseline, not
  enough for a main-method paper.

Hard2 distance-router top-k sweep:

| distance router | mean protein RMSD | mean TM-score |
|---|---:|---:|
| top4 | 2.0624 | 0.8796 |
| top8 | 2.1050 | 0.8784 |
| top12 | 2.1731 | 0.8764 |
| top16 | 2.2407 | 0.8745 |
| top24 | 2.3321 | 0.8695 |

This strengthens the mechanism: on the current hard2 set, the best result comes
from a very small contact-local update budget. Increasing top-k monotonically
degrades mean protein RMSD, which suggests that extra residues add update noise
rather than useful induced-fit context.

Updated bar for a learned method:

> It must beat the best hand-crafted distance top-k baseline, currently
> `distance_top4`, not merely random routing or `distance_top12`.

Hard8 core top4 validation:

| arm | mean protein RMSD | mean TM-score | selected residues/update |
|---|---:|---:|---:|
| baseline_realistic_static5 | 2.4365 | 0.8977 | 52.25 |
| control_realistic_late_dense | 2.4142 | 0.9005 | 52.25 |
| pocket_router_random_top4 | 1.8305 | 0.9248 | 4 |
| pocket_router_distance_top4 | 1.7924 | 0.9249 | 4 |
| pocket_router_distance_top8 | 1.8371 | 0.9229 | 8 |

Updated interpretation after hard8:

- The core sparse-update effect generalizes from hard2 to the 8 hardest cases:
  distance top4 improves mean protein RMSD by 0.6441 A over static5 and 0.6218 A
  over late-dense.
- Top4 remains better than top8, strengthening the "very local update" story.
- Random top4 is also strong, only 0.0381 A worse than distance top4. This is a
  serious warning: sparse masking itself is a strong regularizer, and the learned
  router must clearly beat both random top4 and distance top4.
- Ligand-quality evidence is still incomplete in this run: docking was disabled,
  `evaluated_mols` is 0, and QED/SA are null. The current evidence supports
  protein adaptation, not a complete SBDD claim.

Current top-conference bar:

> A paper-grade method should introduce a learned fragment-conditioned router
> that beats `distance_top4` and `random_top4` under the same hard-tail protocol,
> while preserving or improving ligand validity, reconstruction, and docking
> metrics.

Hard8 full router suite:

Full record:

`validation/pocket_router/hard8_full_router_suite.md`

Run directory:

`validation/ab_runs/hard8_full_router_steps1000_n1`

This run completed 14 arms x 8 hard cases = 112 sampled results.

| arm | mean protein RMSD | mean TM-score | selected residues/update |
|---|---:|---:|---:|
| baseline_realistic_static5 | 2.4365 | 0.8977 | 52.25 |
| control_realistic_late_dense | 2.4142 | 0.9005 | 52.25 |
| control_realistic_uniform10 | 3.1475 | 0.8514 | 53.71 |
| adaptive_realistic_residual | 2.0904 | 0.9141 | 53.13 |
| pocket_router_random_top4 | 1.8305 | 0.9248 | 4 |
| pocket_router_random_top12 | 1.9555 | 0.9196 | 12 |
| pocket_router_motion_oracle_top12 | 2.0511 | 0.9136 | 12 |
| pocket_router_contact_oracle_top12 | 2.1195 | 0.9065 | 12 |
| pocket_router_contact_change_oracle_top12 | 2.2410 | 0.9037 | 12 |
| pocket_router_distance_top4 | 1.7924 | 0.9249 | 4 |
| pocket_router_distance_top8 | 1.8371 | 0.9229 | 8 |
| pocket_router_distance_top12 | 1.8992 | 0.9206 | 12 |
| pocket_router_distance_top16 | 1.9534 | 0.9180 | 16 |
| pocket_router_distance_top24 | 2.0531 | 0.9140 | 24 |

Updated interpretation after the full suite:

- The main protein-adaptation signal is now stronger: `distance_top4` beats
  static5 by 0.6441 A and late-dense by 0.6218 A on hard8.
- The top-k trend is clean: distance top4 is best, and larger top-k degrades
  RMSD monotonically through top24.
- Dense update frequency alone is not the answer. `uniform10` is worse than
  static5, and residual-adaptive global update is better than dense baselines
  but still worse than sparse distance top4.
- Oracle arms do not win. Motion/contact/contact-change oracles at top12 are
  worse than distance top4, which reinforces that the useful mechanism is
  current-fragment local contact routing, not simply "move globally mobile
  residues".
- Random top4 remains very close to distance top4. This is the biggest warning:
  sparse masking is a powerful regularizer, and a learned router must prove it
  adds routing intelligence beyond sparsity itself.

Updated paper bar after the full suite:

> The learned router must beat `distance_top4` and `random_top4` on the same
> hard-tail protocol. Beating dense update baselines is no longer enough.

Evidence still missing before a top-conference submission:

- multi-sample statistics, not only n=1 per case;
- broader test set, not only hard8;
- ligand validity, docking, QED/SA, and clash metrics;
- comparison to external SBDD / apo-to-holo baselines;
- learned-router ablation showing real routing value over hand-crafted distance.

Active-set shell pilot:

Full record:

`validation/pocket_router/active_set_shell_pilot.md`

Run directory:

`validation/ab_runs/hard2_active_set_shell_pilot_steps1000_n1`

This run tested whether the idea should remain a hard top-k mask or become a
ligand-conditioned active-set optimization with a strong core, weak local shell,
and anchored background.

| arm | mean protein RMSD | mean TM-score | selected residues/update | mol stable |
|---|---:|---:|---:|---:|
| control_realistic_late_dense | 2.7867 | 0.8496 | 48.00 | 0.5000 |
| pocket_router_random_top4 | 2.1252 | 0.8795 | 4.00 | 0.5000 |
| pocket_router_distance_top4_hard | 2.0624 | 0.8796 | 4.00 | 1.0000 |
| active_set_distance_top4_shell4_w025 | 2.0641 | 0.8796 | 5.95 | 1.0000 |
| active_set_distance_top4_shell6_w025 | 2.0806 | 0.8794 | 17.80 | 0.5000 |
| active_set_distance_top4_shell6_w050 | 2.0950 | 0.8802 | 17.50 | 0.5000 |
| active_set_random_top4_shell6_w025 | 2.1625 | 0.8786 | 15.50 | 0.5000 |

Updated interpretation:

- Soft active-set weighting is technically viable; all arms completed 1000-step
  sampling on hard2.
- A small weak shell is safe: shell4 w0.25 is nearly identical to hard distance
  top4, with only +0.0017 A RMSD.
- Larger shells are not automatically better. shell6 worsens RMSD slightly,
  which matches the hard8 lesson that extra updated residues can add denoising
  noise.
- Distance-conditioned shell still beats random shell, so routing semantics
  matter beyond simply allowing more residues to move.

Refined paper idea:

> Treat apo-to-holo pocket adaptation as ligand-conditioned active-set
> optimization: release a sparse core of fragment-relevant pocket degrees of
> freedom, weakly relax a narrow neighborhood, and anchor the background.

This is stronger than "only update selected residues" and directly addresses
the local-continuity concern. The next evidence needed after the hard8 shell
run is multi-sample repeats, ligand-side quality metrics, and learned-router
training.

Hard8 active-set shell suite:

Full record:

`validation/pocket_router/hard8_active_set_shell_suite.md`

Run directory:

`validation/ab_runs/hard8_active_set_shell_only_steps1000_n1`

This run completed 4 active-set arms x 8 hard cases = 32 sampled results, plus
local geometry diagnostics for ligand-protein clashes, contact recovery, and
active/background boundary displacement.

| arm | mean protein RMSD | mean TM-score | selected residues/update | mol stable | ligand-protein clashes | boundary jump |
|---|---:|---:|---:|---:|---:|---:|
| active_set_distance_top4_shell4_w025 | 1.8589 | 0.9213 | 6.11 | 0.875 | 0.25 | 0.7971 |
| active_set_distance_top4_shell6_w025 | 1.8666 | 0.9217 | 16.38 | 0.375 | 0.25 | 0.4775 |
| active_set_distance_top4_shell6_w050 | 1.9301 | 0.9192 | 15.51 | 0.625 | 0.375 | 0.5336 |
| active_set_random_top4_shell6_w025 | 1.8524 | 0.9249 | 15.75 | 0.250 | 9.50 | 0.3965 |

Updated interpretation:

- A small weak shell is a safer physical variant: shell4 w0.25 keeps strong
  protein adaptation and greatly reduces ligand-protein clash relative to hard
  sparse baselines.
- A larger shell smooths the active/background boundary, but ligand stability
  can collapse. Boundary continuity alone is not a sufficient objective.
- Random shell is the decisive negative control: its protein RMSD looks strong,
  but it creates severe ligand-protein clashes. This proves that protein RMSD
  alone can be misleading.
- The paper idea should not be "make the mask soft" or "relax more residues".
  The main claim should be "learn the ligand-conditioned active set of pocket
  degrees of freedom".

Revised paper bar:

> A learned method must beat hard `distance_top4`, random top4, and the safe
> active-set shell4 baseline while preserving ligand stability, contact recovery,
> and low clash counts. A method that only improves protein RMSD is not enough.

Hard8 conservative shell sweep:

Full record:

`validation/pocket_router/hard8_conservative_shell_sweep.md`

Run directory:

`validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1`

This run completed 4 conservative active-set arms x 8 hard cases = 32 sampled
results. It tested smaller shell radii and weaker shell weights after the
shell-only run showed that large shells can smooth the boundary but damage
ligand stability.

| arm | mean protein RMSD | mean TM-score | selected residues/update | mol stable | complete | ligand-protein clashes | contact recall |
|---|---:|---:|---:|---:|---:|---:|---:|
| active_set_distance_top4_shell3_w010 | 1.8571 | 0.9215 | 4.33 | 0.750 | 0.875 | 0.250 | 0.7342 |
| active_set_distance_top4_shell3_w025 | 1.8584 | 0.9214 | 4.31 | 0.750 | 0.875 | 0.250 | 0.7342 |
| active_set_distance_top4_shell4_w010 | 1.8512 | 0.9218 | 6.10 | 0.750 | 0.750 | 0.250 | 0.7366 |
| active_set_distance_top4_shell5_w025 | 1.8563 | 0.9219 | 9.98 | 0.625 | 0.750 | 0.375 | 0.7537 |

Updated interpretation:

- `shell3_w010` and `shell3_w025` are the best balanced conservative variants:
  they keep low clash, complete 0.875, and reconstruction 1.000 while updating
  only about 4.3 residues per protein-update step.
- `shell4_w010` has the best protein RMSD among this sweep, but weaker ligand
  completion. This is another warning that protein RMSD alone is not enough.
- `shell5_w025` improves contact recall but worsens mol stability and clash,
  so larger shell/contact recovery is not free.
- The fixed-rule shell baselines now serve as a design map for the learned
  method: learn core release, shell relaxation weight, and background anchoring
  instead of hand-picking one radius and weight.

Revised training bar:

> A top-conference method should learn the active-set policy itself. It must
> outperform hard `distance_top4`, random top4, old `shell4_w025`, and the new
> conservative shell3/shell4 baselines on multi-sample hard-tail runs while
> improving ligand validity, low clash, contact recovery, and docking quality.

Hard8 active-set candidate repeat, n=3:

Full record:

`validation/pocket_router/hard8_candidate_repeat_n3.md`

Run directory:

`validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3`

This run repeated the strongest candidate arms with 3 samples per hard case:

| arm | protein RMSD | mol stable | complete | clashes | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 1.8153 | 0.2500 | 0.7083 | 9.7500 | 0.6821 | 0.5024 | 0.6222 |
| `pocket_router_distance_top4_hard` | 1.8148 | 0.5000 | 0.6667 | 1.7500 | 0.7131 | 0.5349 | 0.6143 |
| `active_set_distance_top4_shell4_w025` | 1.8207 | 0.5000 | 0.7500 | 1.3750 | 0.7557 | 0.5540 | 0.5167 |
| `active_set_distance_top4_shell3_w010` | 1.8171 | 0.5000 | 0.6667 | 1.7500 | 0.7155 | 0.5414 | 0.6124 |
| `active_set_distance_top4_shell3_w025` | 1.8145 | 0.4583 | 0.6250 | 1.7500 | 0.7210 | 0.5437 | 0.6097 |
| `active_set_distance_top4_shell4_w010` | 1.8168 | 0.4583 | 0.7500 | 1.8750 | 0.7346 | 0.5487 | 0.6023 |

Updated interpretation:

- The protein-RMSD gap among distance-based arms is tiny under n=3 repeat, so
  RMSD alone is no longer enough to select the method.
- Random top4 remains a critical negative control: it has competitive protein
  RMSD but severe ligand-protein clashes, so sparse masking alone is not a
  publishable method.
- `active_set_distance_top4_shell4_w025` is the best safety candidate in this
  repeat: it has the lowest clashes, highest contact precision/Jaccard, and
  lowest boundary jump, at a small protein-RMSD cost.
- A no-docking ligand-quality pass also favors `shell4_w025`: it has the highest
  QED mean in this repeat and tied-best complete/chem-success rates. True
  docking is still missing because local `vina/qvina` dependencies are not
  available.
- The core paper idea should therefore stay focused on learned
  ligand-conditioned active-set optimization, not fixed top-k masking.

Revised bar after n=3 repeat:

> The learned active-set policy must beat hard distance top4 and random top4 not
> only on protein RMSD, but also on clash, contact precision/Jaccard, ligand
> validity, docking, and background stability. A method that wins RMSD while
> increasing ligand clashes should be treated as a failure.

## Method Sketch

At each selected denoising step:

1. Encode current ligand fragment state:
   - ligand atom features;
   - partial graph / fragment identity;
   - denoising time;
   - current ligand-pocket residual or clash proxy.

2. Build pocket memory slots:
   - residue-level slots from current apo pocket;
   - subpocket slots from geometric clustering;
   - optional rotamer-mode slots from train-set apo-holo transitions.

3. Router selects top-k memory slots:
   - geometry score: distance/contact to current fragment;
   - chemistry score: residue-ligand compatibility;
   - dynamics score: predicted residue mobility / apo-holo uncertainty;
   - generation score: whether this slot affects the next fragment.

4. Only selected slots receive high-precision coupling:
   - full ligand-pocket equivariant update for selected residues;
   - cheap global context update for the rest;
   - optional local side-chain rotamer update only inside selected chunks.

5. Train with both generation and routing losses:
   - ligand diffusion loss;
   - pocket displacement / rotamer reconstruction loss;
   - sparse router supervision from apo-holo labels:
     high displacement residues, contact-change residues, ligand-near residues;
   - budget loss to keep top-k sparse.

## Why This Beats The Current Residual Schedule As A Paper Idea

The current residual-adaptive schedule asks:

> When should the whole pocket be updated?

PocketRouter asks a stronger question:

> Given the current fragment, which subpocket memory should be activated, and
> which residues deserve high-fidelity conformational update?

That is closer to a main-method contribution because it changes the structure of
reasoning, not just the update timetable.

## Fit To Current Apo2Mol Code

The repo already has useful hooks:

- `topk_prompt` exists but is 0 by default.
- `protein_prompt_list` currently only contains the current `hbap_protein_aug`.
- `prompt_hbap_ligand_batch_all_list` is empty, so retrieval/memory prompting is
  structurally present but unused.
- `RetAugmentationLinearAttention` already exists.
- protein update schedules are now configurable, so static/update-frequency
  controls already exist.

Minimal implementation path:

1. Start with the validated non-learned router:
   distance/contact top-k as the primary route; motion-only as an ablation, not
   the main mechanism.
2. Add a learned router that predicts contact-change / fragment relevance from
   current ligand state, residue identity, residue geometry, and denoising time.
3. Replace empty prompt lists with selected residue/subpocket memory features.
4. Add controls:
   random router, distance-only router, late-dense, uniform10, residual-adaptive.
5. Only claim a learned method if it beats both random and the best
   distance-only top-k routing baseline.

## Required Experiments

Must-have controls:

- original Apo2Mol static5;
- late-dense and uniform10 schedules;
- residual-adaptive global update;
- random top-k memory;
- distance-only top-k memory;
- distance top-k sweep;
- oracle top-k from holo contact-change labels.

Must-have metrics:

- hard-tail protein RMSD on apo-holo RMSD buckets;
- coverage of truly mobile/contact-changing residues by router top-k;
- ligand validity/reconstruction/QED/SA;
- clash/contact/PoseBusters-style checks;
- wall-clock time and number of high-precision residue updates.

The killer ablation:

> If contact/distance sparse memory improves but random and motion-only routing
> do not, the idea is specifically fragment-contact routing. If learned routing
> cannot beat distance-only, the paper should be framed as a strong engineered
> sparse induced-fit baseline or dropped for top-conference submission.

## Submission Framing

Do not frame it as:

> We add sparse attention to Apo2Mol.

Frame it as:

> We introduce fragment-guided contact-sparse pocket memory for realistic
> apo-only SBDD. The method routes each partial ligand state to a sparse set of
> contact-relevant pocket memories and performs high-fidelity induced-fit updates
> only where the current fragment demands it.

This gives one clean main innovation:

> state-conditioned sparse memory routing for fragment-conditioned apo-to-holo
> pocket adaptation.

## Reference Patterns

- Longformer: https://github.com/allenai/longformer
- BigBird: https://github.com/google-research/bigbird
- Swin Transformer: https://github.com/microsoft/Swin-Transformer
- Deformable DETR: https://github.com/fundamentalvision/Deformable-DETR
- Perceiver: https://github.com/google-deepmind/deepmind-research/tree/master/perceiver
- TokenLearner: https://github.com/google-research/scenic/tree/main/scenic/projects/token_learner
- DPR: https://github.com/facebookresearch/DPR
- RAG: https://github.com/facebookresearch/fairseq/tree/main/examples/rag
- Slot Attention: https://github.com/google-research/google-research/tree/master/slot_attention
