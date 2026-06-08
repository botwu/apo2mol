# Apo2Mol New-Method A/B Verification Runbook

This runbook is for proving or falsifying whether the current new-method
candidate improves over the Apo2Mol baseline on hard apo-to-holo cases.

## Current Local Status

Verified locally:

- Python environment: `.venv310`
- Torch: `2.4.0`
- PyG: `torch-geometric==2.0.4`
- PyG native extensions: `torch-scatter`, `torch-sparse`, `torch-cluster`, `torch-spline-conv`
- Apo2Mol checkpoint load: OK
- PMINet load: OK
- `sample_split.py --help`: OK
- `eval_split.py --help`: OK with `docking_mode=none`
- Apo2Mol data: OK via `apo2mol_dataset/data_folder`

Important sampler caveat:

- The current Apo2Mol sampler is an adjacent reverse chain. `num_steps=1000`
  reaches t=0 and is the only setting here that should count as generation
  evidence.
- `num_steps<1000` runs only the high-noise tail, for example 100 steps runs
  999 down to 900. Use shortened runs only to test wiring, dependency health,
  and crash behavior.

The Hugging Face dataset is gated, so the token must belong to an account that
has access to `AIDD-LiLab/Apo2Mol_Dataset`.

## Local Data Download

Do not paste the token into logs or commit it. Export it only in the shell used
for the download:

```bash
HF_TOKEN=... .venv310/bin/python validation/download_apo2mol_data.py
```

The downloader places the data at:

```text
apo2mol_dataset/data_folder
```

After download, run the preflight again:

```bash
.venv310/bin/python validation/run_new_method_ab.py
```

The preflight must say that all requirements are available and that selected
hard-case files are present.

## Local A/B Run

Default hard-tail experiment:

```bash
.venv310/bin/python validation/run_new_method_ab.py --run --run-dir validation/ab_runs/latest
```

For a faster smoke run after data arrives:

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --run-dir validation/ab_runs/smoke \
  --num-cases 2 \
  --num-samples 1 \
  --num-steps 20 \
  --docking-mode none
```

The smoke report should warn that the run does not reach t=0. That warning is
expected for smoke tests.

For the smallest valid local A/B run:

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --run-dir validation/ab_runs/hard1_steps1000_ab \
  --num-cases 1 \
  --num-samples 1 \
  --num-steps 1000 \
  --docking-mode none
```

For a stronger local run:

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --run-dir validation/ab_runs/hard8_steps1000 \
  --num-cases 8 \
  --num-samples 3 \
  --num-steps 1000 \
  --docking-mode none
```

For a claim that is closer to top-conference evidence, include schedule
controls. This tests whether `residual_adaptive` is genuinely useful instead of
merely benefiting from more protein updates:

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --include-controls \
  --run-dir validation/ab_runs/hard8_steps1000_controls \
  --num-cases 8 \
  --num-samples 3 \
  --num-steps 1000 \
  --docking-mode none
```

## What Counts As Evidence Of Gain

The result must be in `results.json` under the run directory. The adaptive arm
only counts as a positive result if it satisfies both:

- Lower mean protein RMSD than `baseline_realistic_static5`.
- No clear degradation in reconstruction/validity metrics from `eval_split.py`.
- For a strong paper claim, lower mean protein RMSD than
  `control_realistic_late_dense` and `control_realistic_uniform10` as well.

If the result is mixed, the correct conclusion is not "the method works"; it is
that the current adaptive update rule is insufficient and needs redesign.

## Current Local Evidence

Completed on 2026-05-26:

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --include-controls \
  --run-dir validation/ab_runs/hard1_steps1000_ab_gated \
  --num-cases 1 \
  --num-samples 1 \
  --num-steps 1000 \
  --docking-mode none
```

This is a valid full-chain smoke A/B with controls, not a publishable result.
It uses one hard test case: test position 477, original index 24500, apo/holo
RMSD 4.0188 A.

Protein RMSD:

- `baseline_realistic_static5`: 3.1261 A
- `adaptive_realistic_residual`: 2.9662 A
- `control_realistic_late_dense`: 3.1191 A
- `control_realistic_uniform10`: 4.0791 A

Interpretation: the gated residual-adaptive rule improved this one hard case and
beat both schedule controls. This is only a weak positive signal because n=1 and
`docking_mode=none` leaves QED/SA/docking unavailable in the local result. The
next credible step is an 8-case full-chain run with controls, then more seeds
and full chemistry/pose metrics on a Linux/GPU environment.

## Hugging Face Jobs Path

Use this if local CPU/MPS is too slow. Requirements:

- Hugging Face account with Jobs access.
- Authenticated `hf` CLI.
- `HF_TOKEN` secret with access to the gated Apo2Mol dataset.
- A code snapshot containing the current validation patches.

First verify login:

```bash
.venv310/bin/hf auth whoami
```

Then run on a GPU-capable environment. A practical option is to push this
working tree to a branch, then submit:

```bash
.venv310/bin/python validation/submit_hf_ab_job.py \
  --repo-url <your-fork-or-branch-url> \
  --repo-ref <branch-with-validation-patches> \
  --include-controls \
  --flavor l4x1 \
  --timeout 8h \
  --results-repo <optional-dataset-repo-for-artifacts>
```

To inspect the exact command without submitting:

```bash
.venv310/bin/python validation/submit_hf_ab_job.py --dry-run
```

Do not use a remote run as evidence unless the job artifacts include:

- `experiment_plan.json`
- `selected_cases.json`
- both arms' sampled result directories
- `results.json`
- the exact git commit or patch used
