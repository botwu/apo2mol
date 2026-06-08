# Apo2Mol Top-Conference Direction Validation

Generated: 2026-05-26T01:15:00

## Bottom line

The local evidence supports a top-conference-worthy direction only if the work is framed as a new problem/benchmark plus a substantial method, not as a small patch on Apo2Mol. The strongest framing is conformational-memory or retrieval-augmented generation for long-tail apo-to-holo pocket motion.

## Dataset evidence

- Records: 24501
- Date range: 1982-06-25 to 2024-04-29
- Unique holo PDB IDs: 12444
- Unique PLI directories: 24501 (duplicated dirs: 0)
- Apo/holo RMSD median: 0.8457; p90: 1.6976; p95: 2.0914; max: 9.9858
- Hard cases >=2A: 1435 (5.86%); >=3A: 362 (1.48%); >=5A: 90 (0.37%)

RMSD bins:

| bin (A) | count | share |
|---|---:|---:|
| 0-0.5 | 3778 | 15.42% |
| 0.5-1 | 11451 | 46.74% |
| 1-1.5 | 5723 | 23.36% |
| 1.5-2 | 2114 | 8.63% |
| 2-3 | 1073 | 4.38% |
| 3-5 | 272 | 1.11% |
| 5-10 | 90 | 0.37% |
| 10+ | 0 | 0.00% |

Split distribution:

| split | n | median | p90 | p95 | >=2A | >=3A | >=5A |
|---|---:|---:|---:|---:|---:|---:|---:|
| train | 22953 | 0.8543 | 1.7104 | 2.1053 | 1376 | 343 | 87 |
| valid | 1070 | 0.7844 | 1.4922 | 1.7310 | 38 | 16 | 3 |
| test | 478 | 0.6846 | 1.4599 | 1.9453 | 21 | 3 | 0 |

Highest-RMSD examples:

| RMSD | date | holo pocket | ligand |
|---:|---|---|---|
| 9.9858 | 2018-02-09 | `6fpj__1__1.C__1.V/receptor_holo_pocket10.pdb` | `6fpj__1__1.C__1.V/1.V.sdf` |
| 9.9477 | 2009-10-31 | `3ki2__1__1.A__1.C/receptor_holo_pocket10.pdb` | `3ki2__1__1.A__1.C/1.C.sdf` |
| 9.9346 | 2009-11-26 | `3ku6__1__2.A__2.D/receptor_holo_pocket10.pdb` | `3ku6__1__2.A__2.D/2.D.sdf` |
| 9.8435 | 2009-10-31 | `3ki1__1__1.A__1.C/receptor_holo_pocket10.pdb` | `3ki1__1__1.A__1.C/1.C.sdf` |
| 9.3707 | 2018-06-17 | `6a3w__3__1.I__1.M/receptor_holo_pocket10.pdb` | `6a3w__3__1.I__1.M/1.M.sdf` |
| 8.3444 | 2002-04-09 | `1leg__1__1.A__1.E/receptor_holo_pocket10.pdb` | `1leg__1__1.A__1.E/1.E.sdf` |
| 8.2269 | 2011-05-06 | `3rv5__4__1.D__1.GA/receptor_holo_pocket10.pdb` | `3rv5__4__1.D__1.GA/1.GA.sdf` |
| 8.1924 | 2017-04-25 | `5vlo__1__1.A__1.C/receptor_holo_pocket10.pdb` | `5vlo__1__1.A__1.C/1.C.sdf` |

## Code/config evidence

- Sampling uses reference ligand atom count: yes
- Sampling initializes ligand around holo pocket center: yes
- Training config sets `topk_prompt: 0`: 2 occurrence(s)
- Model has a `topk_prompt` parameter path: yes
- Prompt retrieval list is initialized empty in the model: yes
- Prompt path falls back to self-attention over the current ligand features: yes
- Model has static protein update timesteps: yes
- Evaluation supports generated/apo/holo pocket types: yes
- Evaluation has Vina/QVina docking hooks: yes
- `data_folder` exists locally: no
- PyTorch importable in current Python: no

Interpretation:

The code already contains enough structure to support three research cuts: a harder benchmark split, real retrieval prompts, and adaptive protein-ligand co-denoising. It also exposes evaluation risks: reference atom count and holo-center initialization can make results less realistic unless explicitly separated into oracle and realistic protocols.

## Top-conference judgment

**Can this target a top conference?** Yes, but the paper should not be positioned as `Apo2Mol + retrieval`. The publishable version needs all of the following:

1. A new problem definition: long-tail apo-to-holo pocket motion under realistic ligand-generation constraints.
2. A benchmark contribution: Apo2Mol-Hard or Apo2Mol-Realistic with explicit RMSD buckets and no oracle ligand-size/holo-center shortcuts in the main protocol.
3. A method contribution: retrieval as conformational memory, or adaptive asynchronous co-denoising, with clear algorithmic novelty.
4. Strong ablations: retrieval source, hard-tail buckets, oracle-vs-realistic sampling, dynamic-vs-static protein updates, and pocket type.
5. External baselines: Apo2Mol, DiffSBDD/DecompDiff/IPDiff-style SBDD methods, and flexible/dynamic protein-generation baselines where feasible.

A weak submission is likely if the novelty is only enabling `topk_prompt > 0` or adding a simple nearest-neighbor prompt. The top-conference story becomes credible when the benchmark demonstrates a failure mode that existing methods do not address and the method materially improves hard-tail performance without relying on oracle information.

## Evidence links to cite

- [ICLR 2026 reviewer guide](https://iclr.cc/Conferences/2026/ReviewerGuide): New, relevant, impactful knowledge backed by convincing evidence.
- [NeurIPS reviewer guidelines](https://nips.cc/Conferences/2025/ReviewerGuidelines): Novelty, empirical support, clarity, limitations.
- [ICML reviewer instructions](https://icml.cc/Conferences/2025/ReviewerInstructions): Significance, technical quality, empirical validation.
- [AAAI-26 main technical track call](https://aaai.org/conference/aaai/aaai-26/main-technical-track-call/): Originality, rigor, reproducibility, impact.
- [Apo2Mol](https://arxiv.org/abs/2511.14559): Direct baseline and dataset/method context.
- [DecompDiff, ICML 2023](https://proceedings.mlr.press/v202/guan23a.html): SBDD diffusion precedent with decomposition as a methodological hook.
- [IPDiff, ICLR 2024](https://openreview.net/forum?id=qH9nrMNTIW): Interaction-prior diffusion precedent.
- [FlexSBDD, NeurIPS 2024](https://proceedings.neurips.cc/paper_files/paper/2024/hash/60fb8cf8000f0386063fb24ead366330-Abstract-Conference.html): Flexible-receptor SBDD precedent.
- [DynamicFlow, ICLR 2025](https://arxiv.org/abs/2503.03989): Protein dynamics for molecular generation precedent.

## Next verification experiments

1. Create an Apo2Mol-Hard test split with RMSD >=2A and >=3A buckets, then report metrics by bucket.
2. Re-run Apo2Mol under two protocols: oracle/reference and realistic/no-reference atom count plus apo-centered initialization.
3. Implement retrieval memory with strict train-only retrieval and leakage checks by PDB cluster/date.
4. Replace static protein update steps with an adaptive scheduler and ablate compute-matched variants.
5. Report validity, docking/Vina score, interaction recovery, generated-pocket RMSD, diversity, and per-bucket degradation.

## Local limitations

This validation is static and dataset-index based. It does not run training or sampling because the current Python cannot import PyTorch and the full `data_folder` used by the configs is not present in this checkout. The evidence proves the research opportunity and benchmark risk; it does not prove the final method will improve metrics.
