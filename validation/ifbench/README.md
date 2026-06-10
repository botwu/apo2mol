# IFBench Data Pipeline

This directory is the start of an independent induced-fit benchmark pipeline.
It should not depend on Apo2Mol processed indices, splits, or checkpoints.

The curation rules are grounded in prior dataset/benchmark methodology. See
`METHODOLOGY_FROM_PRIOR_WORK.md` before changing case filters or split rules.

## Central Principle

The benchmark case is:

`input = apo pocket`, `target = ligand + locally adapted holo-like pocket`.

A case is useful only if it tests the actual optimization problem:

- ligand geometry and chemistry are valid;
- apo and holo are the same protein chain, not unrelated homologs;
- the ligand-binding site is well-defined;
- induced-fit motion is local and measurable;
- background residues are expected to remain stable;
- train/valid/test and existing Apo2Mol assets cannot leak into the benchmark.

If a candidate cannot pass these checks, it should not enter the public dataset,
even if it increases dataset size.

## Pipeline

Prepare a candidate CSV. Each row is one apo/holo/ligand triple:

```bash
cp validation/ifbench/candidate_template.csv validation/ifbench/candidates.csv
```

Or collect candidates from independent RCSB metadata:

```bash
python validation/ifbench/rcsb_collect_candidates.py \
  --out-candidates validation/ifbench/out/candidates.rcsb.csv \
  --query-json validation/ifbench/out/rcsb_query.json \
  --metadata-cache validation/ifbench/out/rcsb_metadata.json \
  --max-holo-entries 200 \
  --max-apo-entries 2000 \
  --max-apo-per-uniprot 50 \
  --holo-release-after 2024-05-01 \
  --apo-search-mode matched-uniprot \
  --min-ligand-formula-weight 120
```

Download the referenced RCSB structure and ligand files:

```bash
python validation/ifbench/rcsb_download_files.py \
  --candidates validation/ifbench/out/candidates.rcsb.csv \
  --raw-root . \
  --skip-existing
```

Required columns:

- `holo_structure_file`, `apo_structure_file`, `ligand_file`
- `holo_chain_id`, `apo_chain_id`
- `holo_pdb_id`, `apo_pdb_id`, `ligand_id`

Curated ligand evidence can be provided separately:

```bash
cp validation/ifbench/ligand_evidence_template.csv validation/ifbench/ligand_evidence.csv
```

Evidence CSV rows are keyed by `pdb_id + ligand_id`, with optional chain/asym
columns. For release candidates, use `--require-ligand-evidence` so tiny
additives or unverified hetero compounds cannot enter silently.

Optional but strongly recommended:

- `uniprot_id`, `cluster_key`
- `holo_release_date`, `apo_release_date`
- `holo_resolution`, `apo_resolution`
- `holo_method`, `apo_method`

Run structure and ligand QC:

```bash
python validation/ifbench/build_manifest.py \
  --candidates validation/ifbench/candidates.csv \
  --raw-root . \
  --ligand-evidence-csv validation/ifbench/ligand_evidence.csv \
  --out validation/ifbench/out/manifest.accepted.jsonl \
  --rejects validation/ifbench/out/manifest.rejected.jsonl \
  --flat-csv validation/ifbench/out/manifest.qc.csv \
  --summary-json validation/ifbench/out/qc_summary.json
```

Or run the source-agnostic case factory end to end:

```bash
python validation/ifbench/case_factory.py \
  --candidates validation/ifbench/out/candidates.rcsb.csv \
  --raw-root . \
  --run-dir validation/ifbench/out/case_factory_run \
  --apo2mol-index apo2mol_dataset/apo2mol_version/selected_index_apo_druglike.pkl \
  --ligand-evidence-csv validation/ifbench/ligand_evidence.csv \
  --allow-gate-fail
```

`--allow-gate-fail` is appropriate for pilots. A release candidate should pass
`validation/ifbench/benchmark_spec.yaml` without this flag.

Run a synthetic smoke test for the accepted QC path:

```bash
python validation/ifbench/smoke_test.py \
  --work-dir /tmp/ifbench_smoke
```

Check leakage against Apo2Mol and other blocklists:

```bash
python validation/ifbench/leakage_check.py \
  --manifest validation/ifbench/out/manifest.accepted.jsonl \
  --apo2mol-index apo2mol_dataset/apo2mol_version/selected_index_apo_druglike.pkl \
  --out validation/ifbench/out/manifest.no_apo2mol_leakage.jsonl \
  --report-json validation/ifbench/out/leakage_report.json
```

Assign induced-fit tiers and group-safe splits:

```bash
python validation/ifbench/stratify_and_split.py \
  --manifest validation/ifbench/out/manifest.no_apo2mol_leakage.jsonl \
  --out validation/ifbench/out/manifest.final.jsonl \
  --flat-csv validation/ifbench/out/manifest.final.csv \
  --splits-json validation/ifbench/out/splits.json \
  --summary-json validation/ifbench/out/split_summary.json
```

Export final cases into the Apo2Mol data layout:

```bash
python validation/ifbench/export_apo2mol.py \
  --manifest validation/ifbench/out/manifest.final.jsonl \
  --raw-root . \
  --out-data-folder validation/ifbench/out/data_folder \
  --out-index validation/ifbench/out/selected_index_ifbench.pkl \
  --out-split-pt validation/ifbench/out/split_ifbench.pt \
  --out-split-json validation/ifbench/out/split_ifbench.json
```

Build a construction report:

```bash
python validation/ifbench/make_report.py \
  --accepted validation/ifbench/out/manifest.accepted.jsonl \
  --rejected validation/ifbench/out/manifest.rejected.jsonl \
  --final validation/ifbench/out/manifest.final.jsonl \
  --qc-summary validation/ifbench/out/qc_summary.json \
  --leakage-report validation/ifbench/out/leakage_report.json \
  --split-summary validation/ifbench/out/split_summary.json \
  --out validation/ifbench/out/ifbench_report.md
```

Create a manual review queue:

```bash
python validation/ifbench/make_review_queue.py \
  --accepted validation/ifbench/out/manifest.accepted.jsonl \
  --rejected validation/ifbench/out/manifest.rejected.jsonl \
  --final validation/ifbench/out/manifest.final.jsonl \
  --out-csv validation/ifbench/out/review_queue.csv \
  --out-md validation/ifbench/out/review_queue.md
```

## Quality Gates

Current hard rejections:

- missing or unparsable apo/holo/ligand files;
- RDKit ligand sanitization failure;
- disconnected ligand, unsupported elements, extreme formal charge;
- ligand heavy atoms outside the configured range;
- apo/holo chain sequence identity below `0.95`;
- mapping coverage below `0.85`;
- fewer than `40` common residues;
- fewer than `8` pocket residues within `8 A` of the holo ligand;
- resolution worse than `2.5 A` when resolution metadata is provided;
- hetero ligand atoms in the apo pocket within `5 A` of the holo ligand site;
- PDB ID overlap with Apo2Mol or explicit blocklists.

These thresholds are deliberately strict for the first public benchmark pass.
They can be relaxed only after inspecting reject reports and documenting why.

## Metrics And Labels

`build_manifest.py` computes:

- ligand chemistry metrics: heavy atoms, fragments, molecular weight, logP,
  rotatable bonds, rings, formal charge, elements;
- apo/holo chain metrics: sequence identity, mapping coverage, common residues;
- global aligned CA RMSD;
- pocket metrics: pocket CA RMSD, mean/P90/max CA displacement;
- active-set labels: motion core residues, holo contacts, contact-change
  residues, per-residue distances and displacements.

`stratify_and_split.py` assigns `easy`, `medium`, `hard`, and `extreme` tiers
from a composite induced-fit score:

```text
pocket_ca_rmsd
+ 0.35 * pocket_p90_ca_displacement
+ 1.25 * contact_change_fraction
+ 0.75 * motion_core_fraction
```

This score is intentionally about pocket adaptation, not ligand size or model
performance.

## Source Strategy

Candidate discovery should use independent public sources, then pass through
the same QC manifest:

- RCSB PDB as the structure authority;
- Binding MOAD / PDBbind / PLINDER as candidate seeds or metadata helpers;
- manual curation for difficult induced-fit examples.

The benchmark release should publish the candidate-generation script, accepted
manifest, reject summary, leakage report, split file, manual review queue, and
data license notes.

The first real RCSB pilot exposed two practical rules:

- broad "latest apo" scanning is ineffective; collect holo seeds first, then
  search apo candidates by the holo UniProt accession;
- candidate-stage ligand filtering should remove obvious tiny additives, but
  final ligand quality is still decided by `build_manifest.py`.

## Next Engineering Tasks

- Add MOAD/PDBbind/PLINDER seed importers.
- Add automatic manual-inspection bundles for accepted hard/extreme cases.
- Calibrate tier thresholds on the first real RCSB candidate batch.
