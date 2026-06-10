# IFBench Methodology From Prior Work

This benchmark should be literature-driven. The goal is not to invent filters in
isolation, but to adapt proven curation ideas to apo-to-holo induced-fit SBDD.

## Target Problem

Every case must test the same question:

```text
Given an apo pocket, can a model generate a chemically valid ligand and a
localized, sparse, ligand-relevant holo-like pocket adaptation?
```

This makes the benchmark different from rigid holo-pocket generation, affinity
prediction, and standard docking.

## What To Borrow

| Work family | What they do well | What we should borrow | IFBench implementation |
| --- | --- | --- | --- |
| Binding MOAD | Curates biologically relevant ligands and separates valid binders from crystallization artifacts. | Ligand biological relevance is not optional; do not treat every hetero compound as a drug-like ligand. | Add ligand provenance tags, additive/cofactor blacklist, optional MOAD/BioLiP validity cross-check. |
| PDBbind refined/core/CASF | Uses quality tiers, affinity/structure filters, and family-balanced subsets. | Public benchmark must expose quality tiers and not just a raw scrape. | Keep `accepted/rejected/final` manifests, reject summaries, tier balance, and family/group split. |
| PLINDER | Designs splits to reduce protein-ligand interaction leakage, not just random train/test. | Leakage must be pocket/interaction-level, not only PDB ID. | Add protein cluster, ligand similarity, and pocket-contact fingerprint leakage reports. |
| CrossDocked2020 | Generates derived cases from cross-docking across similar pockets, but keeps purpose clear: ML training/evaluation data, not raw experimental truth. | Synthetic or postprocessed cases are valid only if their provenance and objective are explicit. | Split `real_if`, `postprocessed_if`, and `synthetic_if`; never mix them silently. |
| DockGen / PoseBench | Focuses on generalization to novel pockets and interaction diversity. | Benchmark should measure whether hard cases are genuinely novel and diverse. | Add ECOD/UniProt/pocket-cluster diversity and interaction-type histograms. |
| PoseBusters | Evaluates physical validity beyond ligand RMSD. | Generated and reference cases need chemical/geometric sanity checks. | Add RDKit sanitization, intramolecular strain, protein-ligand clash, and contact plausibility checks. |
| AH-DB / AHoJ / apo-holo resources | Systematically match apo/holo structures and define binding residues by mapped pockets. | Apo/holo matching must be explicit and auditable. | Keep residue mapping coverage, sequence identity, pocket residue coverage, and apo/holo alignment metrics. |
| Apo2Mol | Builds apo-holo-ligand triplets with strict sequence identity and resolution filters. | Their high-level triplet formulation is useful, but our benchmark needs independent cases and stronger induced-fit stratification. | Use Apo2Mol only as a blocklist/baseline, not as the public benchmark source. |

## Case Types

### 1. Real-IF

Official structures only. A valid case has:

- experimentally resolved apo structure;
- experimentally resolved holo structure;
- same protein chain or very high-confidence matched chain;
- biologically meaningful ligand;
- local pocket motion and contact-change labels.

Use this for the public benchmark.

### 2. Postprocessed-IF

Official structures with additional computational cleanup or pairing:

- structure cleanup;
- protonation/tautomer normalization;
- pocket cropping;
- apo/holo alignment;
- optional side-chain repair when the repair is explicitly logged.

Use this for training and ablation, and only for benchmark if provenance is
clear.

### 3. Synthetic-IF

Controlled cases generated from official structures:

- perturb contact residues to construct an apo-like state;
- preserve background residues;
- keep ligand/holo state as reference;
- produce known active-set labels.

Use this for model training, diagnostics, stress tests, and sanity checks. It
should not be the main public benchmark claim.

## Mandatory Gates

### Source And Ligand Gates

- protein-ligand complex must come from official structural resources;
- ligand must pass RDKit sanitization;
- reject ions, waters, tiny additives, crystallization buffers, peptides, DNA/RNA ligands, and metal clusters unless a specific task explicitly includes them;
- prefer biological ligand labels from curated sources such as Binding MOAD/BioLiP when available;
- record ligand provenance and rejection reason.

### Apo/Holo Matching Gates

- same UniProt or explicit chain mapping;
- high sequence identity;
- high mapping coverage;
- no missing pocket backbone atoms;
- pocket residues must map between apo and holo;
- apo pocket must not contain a ligand in the same site unless the task is explicitly multi-state.

### Induced-Fit Gates

- compute pocket CA/backbone RMSD after global alignment;
- compute side-chain/contact changes;
- label `motion_core`, `contact_changed`, `shell`, and `stable_background`;
- reject cases with no measurable pocket or no reliable pocket mapping;
- stratify into easy/medium/hard/extreme based on induced-fit difficulty.

### Leakage Gates

- PDB ID blocklist;
- UniProt/family cluster split;
- ligand similarity split;
- pocket-contact fingerprint split;
- Apo2Mol overlap blocklist.

### Physical Validity Gates

- ligand chemistry validity;
- ligand internal geometry sanity;
- protein-ligand clash check;
- protein self-clash check near pocket;
- contact plausibility and pocket completeness.

## Near-Term Implementation Order

1. Convert the current scripts into a source-agnostic `CaseFactory`.
2. Add curated source adapters: RCSB, Binding MOAD/BioLiP labels, PDBbind/PLINDER metadata.
3. Add synthetic IF generation only after real-case QC is stable.
4. Add leakage reports beyond PDB ID: UniProt, ligand fingerprint, pocket-contact fingerprint.
5. Calibrate hard/extreme tiers on real accepted cases, then lock thresholds.

## Non-Negotiable Release Artifacts

- candidate manifest;
- accepted/rejected/final manifest;
- reject reason summary;
- leakage report;
- split file;
- tier balance report;
- provenance labels;
- exact scripts and configs used to build the release.
