# Apo2Mol Top-Conference Style Survey

Date: 2026-05-26

This note answers one practical question: can a method built on top of Apo2Mol,
especially an adaptive protein-update/sampling change, plausibly support a
top-conference submission?

Short answer: not as a standalone heuristic. It can become paper-worthy if it is
framed as a rigorous realistic-apo/flexible-pocket problem, supported by hard-tail
benchmarks, compute-matched controls, physical validity checks, and broad
baselines.

## What Recent Top Papers Look Like

Top SBDD papers in 2022-2026 usually have five ingredients:

1. A clear modeling gap, not only an implementation tweak.
2. A method that changes the generative objective, path, representation,
   conditioning, guidance, or evaluation protocol.
3. Strong baselines that include both older autoregressive/flow/diffusion models
   and the closest newer method.
4. Chemistry/physics sanity checks beyond one scalar score.
5. Reproducible code, data splits, trained checkpoints, and ablations.

Relevant public examples:

- Apo2Mol, AAAI 2026, defines the apo-conditioned flexible-pocket generation
  problem and contributes a 24k apo-holo paired dataset plus a joint ligand/pocket
  diffusion model.
  Source: https://ojs.aaai.org/index.php/AAAI/article/view/37138
  Code: https://github.com/AIDD-LiLab/Apo2Mol

- FlexSBDD, NeurIPS 2024, targets the same rigid-pocket limitation but uses flow
  matching, flexible protein-ligand modeling, and augmentation via relaxation and
  side-chain repacking. Its abstract explicitly claims fewer clashes and better
  protein-ligand interactions, not only lower RMSD.
  Source: https://arxiv.org/abs/2409.19645
  Code placeholder: https://github.com/zaixizhang/FlexSBDD

- MolCRAFT/MolPilot, ICML 2024/2025, shows that schedule/path design itself is
  already a competitive top-conference topic only when it is tied to an optimized
  probabilistic objective and strong physical-validity results. MolPilot reports
  VLB-optimal scheduling and a 95.9% PoseBusters-valid rate in its repo summary.
  Source: https://github.com/GenSI-THUAIR/MolCRAFT
  ICML page: https://icml.cc/virtual/2025/poster/44761

- TargetDiff, ICLR 2023, and Pocket2Mol, ICML 2022, are standard de novo SBDD
  baselines. Their public repos also show the expected reproducibility pattern:
  test-set sampling, docking modes, downloaded meta results, and CrossDocked
  evaluation.
  TargetDiff: https://github.com/guanjq/targetdiff
  Pocket2Mol: https://github.com/pengxingang/Pocket2Mol

- DecompDiff, ICML 2023, is a useful baseline because it changes priors and
  reports geometry distributions such as bond distance/angle JSD, not only Vina.
  Source: https://github.com/bytedance/DecompDiff

- DiffSBDD is a strong diffusion baseline with de novo design, inpainting,
  optimization, CrossDocked, Binding MOAD, and sampled molecule evaluation paths.
  Source: https://github.com/arneschneuing/DiffSBDD

- DynamicBind and FlowDock are not direct de novo ligand-generation baselines, but
  they matter because they push the field toward flexible docking/generative
  complex prediction from unbound or apo-like proteins.
  DynamicBind: https://github.com/luwei0917/DynamicBind
  FlowDock: https://github.com/BioinfoMachineLearning/FlowDock

- PoseBusters is important for the evaluation standard. Its key message is that
  RMSD alone is not sufficient; generated/docked poses need steric and energetic
  plausibility checks.
  Source: https://arxiv.org/abs/2308.05777
  Code: https://github.com/maabuu/posebusters

## What This Means For The Current Idea

The current local method candidate is:

- realistic sampling: use apo center instead of holo center;
- realistic atom-count prior instead of reference ligand atom count;
- residual-adaptive protein updates during reverse diffusion;
- hard-test selection from apo-holo RMSD tail.

This is a good experimental starting point, but it is too small for a main
top-conference claim if presented as "we changed when to update the protein." The
reason is that schedule design is already a recognized contribution area, and
papers like MolPilot make that contribution with an explicit objective and broad
evidence. A residual threshold alone will be viewed as a heuristic unless it beats
schedule-matched controls and has a convincing mechanism.

The strongest framing is therefore not:

> Adaptive update schedule improves Apo2Mol.

The stronger framing is:

> Existing apo-conditioned SBDD evaluation hides the hard long-tail of induced-fit
> pocket motion and still uses benchmark-only information. We introduce a
> realistic apo-only hard-tail protocol and an adaptive asynchronous co-denoising
> controller that updates pocket degrees of freedom only when ligand-pocket
> evidence requires it.

## Required Evidence Before Claiming "It Works"

The experiment must separate three effects:

- original Apo2Mol static pocket update;
- more frequent/different scheduled pocket updates;
- truly adaptive residual-triggered updates.

That is why `validation/run_new_method_ab.py` now supports `--include-controls`.
The minimum serious comparison is:

- `baseline_realistic_static5`
- `control_realistic_late_dense`
- `control_realistic_uniform10`
- `adaptive_realistic_residual`

The adaptive arm only looks real if it beats all three under the same cases,
sample budget, number of diffusion steps, atom-count protocol, and apo-centered
initialization.

Minimum metrics:

- pocket: generated-vs-holo pocket RMSD, TM-score, bucket-wise RMSD on `>=2A`
  and `>=3A` cases;
- ligand: validity, reconstruction, QED, SA, diversity, bond-length/angle JSD;
- complex: Vina score/dock where available, steric clashes, H-bonds/contact
  fingerprints, PoseBusters-style checks;
- efficiency: wall-clock time, samples/sec, number of protein updates;
- statistics: paired comparison across the same test cases and at least 3 seeds
  or bootstrap confidence intervals.

The current local `docking_mode=none` path is useful for smoke tests, but not
enough for a submission. A paper result needs full docking/pose validity on a
Linux/conda/GPU environment where Vina and AutoDockTools are stable.

## Best Paper Directions From This Repo

### Direction 1: Apo2Mol-Hard + Realistic Apo Protocol

This is the most defensible direction. The official Hugging Face dataset has
24,601 apo-holo entries and exposes `holo_apo_pocket_rmsd`. Local metadata checks
already show that the official test split has only 21 cases with RMSD >= 2A and
3 cases with RMSD >= 3A, so average test metrics are likely dominated by easy
small-motion pockets.

Contribution:

- define hard-tail buckets: `0-0.5`, `0.5-1`, `1-1.5`, `1.5-2`, `2-3`, `>=3A`;
- remove reference crutches: no holo center, no reference ligand atom count;
- report where Apo2Mol and baselines fail;
- add a robust baseline such as hard-tail reweighting, atom-count prior, or
  adaptive co-denoising.

Top-conference chance: plausible if the benchmark exposes a real failure mode and
the proposed robust method materially improves hard cases without sacrificing
chemistry.

### Direction 2: Adaptive Asynchronous Co-Denoising

This is the natural method direction from the code. Apo2Mol already updates the
ligand every step and pocket only at a few fixed time steps. That is an implicit
temporal separation. A paper-worthy method would turn it into a learned or
principled controller:

- update protein residues when ligand residual, clash, interaction mismatch, or
  uncertainty crosses a threshold;
- optionally update only local residues instead of the whole pocket;
- compare against static, uniform, late-dense, and compute-matched schedules;
- prove the method improves hard induced-fit cases or reduces compute at equal
  quality.

Top-conference chance: weak as a hand-tuned threshold; moderate if upgraded to a
learned controller with strong controls and speed/quality tradeoff.

### Direction 3: Retrieval-Augmented Apo-Holo Transition Memory

Apo2Mol already has prompt-related config/code paths, but `topk_prompt` is 0 by
default. A strong extension is to retrieve train-only apo-holo transition
examples and use their conformational deltas as prompts.

Critical controls:

- no test leakage;
- exclude same PDB/PLI/protein-family leakage where necessary;
- random retrieval negative control;
- sequence-only vs geometry-only retrieval;
- hard-tail bucket evaluation.

Top-conference chance: moderate to high if retrieval clearly helps large
apo-holo motion and the leakage story is airtight.

### Direction 4: Physical Validity/Reranking For Dynamic Pockets

Apo2Mol jointly generates ligand and holo-like pocket conformation, so it can
produce pocket states that look close by RMSD but are physically poor. A useful
extension is a pocket-aware verifier/reranker:

- ligand PoseBusters checks;
- protein side-chain clash/rotamer checks;
- ligand strain and contact fingerprint checks;
- Vina plus physics-composite selection.

Top-conference chance: moderate alone, stronger when combined with Direction 1
or Direction 2.

### Direction 5: Stochastic Bridge Or Flow For Apo-To-Holo Transitions

Instead of treating the apo-holo transition as a fixed diffusion path, model it
as a stochastic bridge or flow over conformational states. This is more ambitious
and closer to a main-method paper, but engineering cost is high.

Top-conference chance: high only if implemented cleanly and evaluated broadly.

## Submission Judgment

Current residual-adaptive schedule alone:

- AAAI/ICLR/NeurIPS/ICML main track: unlikely.
- Strong workshop or short paper: possible if results are clean.
- Main-track potential: only after reframing as a hard-tail realistic-apo
  benchmark plus a stronger adaptive controller or retrieval/physics extension.

The most pragmatic route is:

1. Run the current A/B as soon as `data_folder` is available.
2. Run `--include-controls` before believing any gain.
3. If adaptive only beats `static5` but loses to `late_dense`, drop the current
   method and pivot to a principled controller or benchmark paper.
4. If adaptive beats all controls on `>=2A` and `>=3A` cases, expand to full
   metrics, more seeds, full Vina/PoseBusters, and at least two external baselines.
5. Build the paper around "realistic apo-only hard-tail SBDD", not around a small
   implementation modification.

## Local Validation Update

After the gated residual trigger fix, a full 1000-step local A/B with controls
was run on one hard test case:

- run directory: `validation/ab_runs/hard1_steps1000_ab_gated`
- test position 477, original index 24500, apo/holo RMSD 4.0188 A
- `sample_num_atoms=prior`, `init_center_mode=apo`, `docking_mode=none`

Protein RMSD on this single case:

- `baseline_realistic_static5`: 3.1261 A
- `adaptive_realistic_residual`: 2.9662 A
- `control_realistic_late_dense`: 3.1191 A
- `control_realistic_uniform10`: 4.0791 A

This is a weak positive signal: adaptive beats the original baseline and the two
schedule controls on n=1. It is not yet a top-conference claim. The result says
the idea is worth scaling to hard8/full-seed/full-metric validation, not that the
method is proven.

Important correction: previous shortened 20/100-step runs should be treated only
as pipeline smoke tests. The current sampler runs an adjacent reverse chain; with
`num_steps=100`, it samples t=999 down to 900 and does not reach t=0, so it can
miss protein-update schedules entirely.

## Data Source Notes

The Apo2Mol Hugging Face dataset is gated and contains `data_folder.tar.gz` plus
`split_druglike_dict.pkl`. The dataset card reports 24,601 entries and a total
file size of 2.88 GB.

Source: https://huggingface.co/datasets/AIDD-LiLab/Apo2Mol_Dataset
