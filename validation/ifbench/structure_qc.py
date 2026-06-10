#!/usr/bin/env python3
"""Structure and ligand quality checks for IFBench cases."""

from __future__ import annotations

import difflib
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.Polypeptide import is_aa
except Exception:  # pragma: no cover - dependency may be absent on docs-only machines
    MMCIFParser = None
    PDBParser = None
    is_aa = None

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
except Exception:  # pragma: no cover
    Chem = None
    Descriptors = None
    Lipinski = None
    rdMolDescriptors = None

from ifbench.ligand_evidence import LigandEvidenceIndex, summarize_evidence
from ifbench.schema import IFBenchCase, QualityIssue, REJECT, WARNING, normalize_pdb_id, parse_float, stable_case_id


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}
COMMON_ION_NAMES = {
    "NA",
    "K",
    "CL",
    "CA",
    "MG",
    "MN",
    "ZN",
    "FE",
    "CU",
    "CO",
    "NI",
    "CD",
    "HG",
}
ALLOWED_LIGAND_ELEMENTS = {"B", "C", "N", "O", "F", "P", "S", "CL", "BR", "I"}


@dataclass
class ResidueRecord:
    key: str
    resname: str
    resseq: int
    icode: str
    atoms: dict[str, np.ndarray]

    @property
    def aa(self) -> str:
        return AA3_TO_1.get(self.resname.upper(), "X")


@dataclass
class ChainRecord:
    residues: list[ResidueRecord]
    hetero_atoms: np.ndarray
    hetero_names: list[str]
    warnings: list[str]

    @property
    def sequence(self) -> str:
        return "".join(residue.aa for residue in self.residues)


@dataclass
class LigandRecord:
    coords: np.ndarray
    metrics: dict[str, Any]
    issues: list[QualityIssue]


@dataclass
class QualityConfig:
    max_resolution_a: float = 2.5
    min_sequence_identity: float = 0.95
    min_mapping_coverage: float = 0.85
    min_common_residues: int = 40
    min_pocket_residues: int = 8
    pocket_cutoff_a: float = 8.0
    contact_cutoff_a: float = 4.5
    motion_core_cutoff_a: float = 1.5
    contact_change_cutoff_a: float = 1.0
    apo_pocket_hetero_cutoff_a: float = 5.0
    min_ligand_heavy_atoms: int = 8
    max_ligand_heavy_atoms: int = 80
    max_abs_formal_charge: int = 3
    require_ligand_evidence: bool = False
    reject_negative_ligand_evidence: bool = True


def reject(code: str, message: str) -> QualityIssue:
    return QualityIssue(code=code, severity=REJECT, message=message)


def warn(code: str, message: str) -> QualityIssue:
    return QualityIssue(code=code, severity=WARNING, message=message)


def resolve_path(raw_root: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = raw_root / path
    return path


def parse_structure(path: Path):
    if PDBParser is None or MMCIFParser is None:
        raise RuntimeError("Biopython is required for structure QC")
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def _atom_coord(atom: Any) -> np.ndarray:
    if atom.is_disordered():
        atom = atom.selected_child
    return np.asarray(atom.get_coord(), dtype=np.float64)


def _atom_element(atom: Any) -> str:
    element = getattr(atom, "element", "") or atom.get_name()[0]
    return str(element).upper().strip()


def extract_chain(path: Path, chain_id: str) -> ChainRecord:
    structure = parse_structure(path)
    models = list(structure.get_models())
    warnings: list[str] = []
    if not models:
        raise ValueError(f"{path} has no models")
    if len(models) > 1:
        warnings.append(f"{path.name} has {len(models)} models; using model 0")
    model = models[0]
    if chain_id not in model:
        raise ValueError(f"{path.name} does not contain chain {chain_id}")
    chain = model[chain_id]

    residues: list[ResidueRecord] = []
    hetero_coords: list[np.ndarray] = []
    hetero_names: list[str] = []
    for residue in chain.get_residues():
        hetflag, resseq, icode = residue.id
        resname = residue.get_resname().strip().upper()
        is_standard = hetflag == " " and (is_aa(residue, standard=False) if is_aa is not None else resname in AA3_TO_1)
        if is_standard and resname in AA3_TO_1:
            atoms: dict[str, np.ndarray] = {}
            for atom in residue.get_atoms():
                if _atom_element(atom) == "H":
                    continue
                atoms[atom.get_name().strip()] = _atom_coord(atom)
            if "CA" not in atoms:
                warnings.append(f"{path.name}:{chain_id}:{resseq}{icode} {resname} missing CA")
                continue
            key = f"{resname}:{resseq}:{icode.strip() or '-'}"
            residues.append(ResidueRecord(key=key, resname=resname, resseq=int(resseq), icode=icode.strip(), atoms=atoms))
        elif resname not in WATER_NAMES and resname not in COMMON_ION_NAMES:
            for atom in residue.get_atoms():
                if _atom_element(atom) != "H":
                    hetero_coords.append(_atom_coord(atom))
                    hetero_names.append(resname)

    hetero_atoms = np.asarray(hetero_coords, dtype=np.float64).reshape((-1, 3)) if hetero_coords else np.zeros((0, 3))
    return ChainRecord(residues=residues, hetero_atoms=hetero_atoms, hetero_names=hetero_names, warnings=warnings)


def load_ligand(path: Path, config: QualityConfig) -> LigandRecord:
    issues: list[QualityIssue] = []
    if Chem is None:
        return LigandRecord(np.zeros((0, 3)), {}, [reject("rdkit_missing", "RDKit is required for ligand QC")])
    if not path.exists():
        return LigandRecord(np.zeros((0, 3)), {}, [reject("ligand_missing", f"Missing ligand file: {path}")])

    mol = Chem.MolFromMolFile(str(path), sanitize=True, removeHs=False)
    if mol is None:
        return LigandRecord(np.zeros((0, 3)), {}, [reject("ligand_sanitize_failed", f"RDKit failed to sanitize {path}")])
    if mol.GetNumConformers() == 0:
        issues.append(reject("ligand_no_conformer", f"{path} has no 3D conformer"))
        return LigandRecord(np.zeros((0, 3)), {}, issues)

    frags = Chem.GetMolFrags(mol)
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    elements = sorted({_atom.GetSymbol().upper() for _atom in heavy_atoms})
    bad_elements = sorted(set(elements) - ALLOWED_LIGAND_ELEMENTS)
    formal_charge = int(Chem.GetFormalCharge(mol))

    if len(frags) != 1:
        issues.append(reject("ligand_multiple_fragments", f"Ligand has {len(frags)} disconnected fragments"))
    if len(heavy_atoms) < config.min_ligand_heavy_atoms:
        issues.append(reject("ligand_too_small", f"Ligand has {len(heavy_atoms)} heavy atoms"))
    if len(heavy_atoms) > config.max_ligand_heavy_atoms:
        issues.append(reject("ligand_too_large", f"Ligand has {len(heavy_atoms)} heavy atoms"))
    if bad_elements:
        issues.append(reject("ligand_bad_elements", f"Ligand has unsupported elements: {bad_elements}"))
    if abs(formal_charge) > config.max_abs_formal_charge:
        issues.append(reject("ligand_high_charge", f"Ligand formal charge is {formal_charge}"))

    conf = mol.GetConformer()
    coords = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append([pos.x, pos.y, pos.z])

    metrics = {
        "ligand_heavy_atoms": len(heavy_atoms),
        "ligand_fragments": len(frags),
        "ligand_mol_wt": float(Descriptors.MolWt(mol)) if Descriptors is not None else None,
        "ligand_logp": float(Descriptors.MolLogP(mol)) if Descriptors is not None else None,
        "ligand_rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)) if Lipinski is not None else None,
        "ligand_rings": int(rdMolDescriptors.CalcNumRings(mol)) if rdMolDescriptors is not None else None,
        "ligand_formal_charge": formal_charge,
        "ligand_elements": elements,
    }
    return LigandRecord(np.asarray(coords, dtype=np.float64), metrics, issues)


def residue_mapping(apo: ChainRecord, holo: ChainRecord) -> tuple[list[tuple[int, int]], dict[str, float]]:
    apo_by_number = {(r.resseq, r.icode, r.resname): i for i, r in enumerate(apo.residues)}
    exact: list[tuple[int, int]] = []
    for holo_idx, holo_res in enumerate(holo.residues):
        apo_idx = apo_by_number.get((holo_res.resseq, holo_res.icode, holo_res.resname))
        if apo_idx is not None:
            exact.append((apo_idx, holo_idx))
    if len(exact) >= min(20, int(0.75 * min(len(apo.residues), len(holo.residues)))):
        coverage = min(len(exact) / max(1, len(apo.residues)), len(exact) / max(1, len(holo.residues)))
        return exact, {"sequence_identity": 1.0, "mapping_coverage": coverage, "mapping_mode": "residue_number"}

    matcher = difflib.SequenceMatcher(a=apo.sequence, b=holo.sequence, autojunk=False)
    mapping = []
    matches = 0
    for block in matcher.get_matching_blocks():
        for offset in range(block.size):
            apo_idx = block.a + offset
            holo_idx = block.b + offset
            if apo.residues[apo_idx].aa == holo.residues[holo_idx].aa:
                mapping.append((apo_idx, holo_idx))
                matches += 1
    coverage = min(matches / max(1, len(apo.residues)), matches / max(1, len(holo.residues)))
    identity = matches / max(1, max(len(apo.residues), len(holo.residues)))
    return mapping, {"sequence_identity": identity, "mapping_coverage": coverage, "mapping_mode": "sequence"}


def kabsch_align(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if mobile.shape != target.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("Kabsch inputs must both be (N, 3)")
    mobile_centroid = mobile.mean(axis=0)
    target_centroid = target.mean(axis=0)
    mobile_centered = mobile - mobile_centroid
    target_centered = target - target_centroid
    cov = mobile_centered.T @ target_centered
    u, _, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    aligned = mobile_centered @ rot.T + target_centroid
    trans = target_centroid - mobile_centroid @ rot.T
    return aligned, rot, trans


def apply_transform(coords: np.ndarray, rot: np.ndarray, trans: np.ndarray) -> np.ndarray:
    if coords.size == 0:
        return coords.reshape((-1, 3))
    return coords @ rot.T + trans


def min_distance(left: np.ndarray, right: np.ndarray) -> float:
    if left.size == 0 or right.size == 0:
        return math.inf
    diff = left[:, None, :] - right[None, :, :]
    return float(np.sqrt((diff * diff).sum(axis=-1)).min())


def percentile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    return float(np.percentile(np.asarray(values, dtype=np.float64), q))


def evaluate_candidate(
    row: dict[str, Any],
    raw_root: Path,
    config: QualityConfig,
    ligand_evidence_index: LigandEvidenceIndex | None = None,
) -> IFBenchCase:
    source = (row.get("source") or "candidate").strip()
    holo_pdb_id = normalize_pdb_id(row.get("holo_pdb_id"))
    apo_pdb_id = normalize_pdb_id(row.get("apo_pdb_id"))
    holo_chain_id = (row.get("holo_chain_id") or "").strip()
    apo_chain_id = (row.get("apo_chain_id") or "").strip()
    ligand_id = (row.get("ligand_id") or row.get("ligand_resname") or "").strip().upper()
    holo_structure_file = (row.get("holo_structure_file") or "").strip()
    apo_structure_file = (row.get("apo_structure_file") or "").strip()
    ligand_file = (row.get("ligand_file") or "").strip()
    case_id = (row.get("case_id") or "").strip() or stable_case_id(
        [holo_pdb_id, apo_pdb_id, holo_chain_id, apo_chain_id, ligand_id, ligand_file]
    )

    record = IFBenchCase(
        case_id=case_id,
        source=source,
        holo_pdb_id=holo_pdb_id,
        apo_pdb_id=apo_pdb_id,
        holo_chain_id=holo_chain_id,
        apo_chain_id=apo_chain_id,
        ligand_id=ligand_id,
        holo_structure_file=holo_structure_file,
        apo_structure_file=apo_structure_file,
        ligand_file=ligand_file,
        uniprot_id=(row.get("uniprot_id") or None),
        cluster_key=(row.get("cluster_key") or row.get("uniprot_id") or None),
        holo_release_date=(row.get("holo_release_date") or None),
        apo_release_date=(row.get("apo_release_date") or None),
        holo_resolution=parse_float(row.get("holo_resolution")),
        apo_resolution=parse_float(row.get("apo_resolution")),
        holo_method=(row.get("holo_method") or None),
        apo_method=(row.get("apo_method") or None),
    )

    issues: list[QualityIssue] = []
    if ligand_evidence_index is not None:
        ligand_evidence = ligand_evidence_index.match_candidate(row)
        evidence_summary = summarize_evidence(ligand_evidence)
        record.evidence["ligand_relevance"] = evidence_summary
        record.metrics["ligand_evidence_count"] = evidence_summary["count"]
        record.metrics["ligand_evidence_positive_count"] = evidence_summary["positive_count"]
        record.metrics["ligand_evidence_negative_count"] = evidence_summary["negative_count"]
        if config.reject_negative_ligand_evidence and evidence_summary["negative_count"] > 0:
            issues.append(reject("negative_ligand_evidence", "Curated evidence marks this ligand as not benchmark-suitable"))
        if config.require_ligand_evidence and evidence_summary["positive_count"] == 0:
            issues.append(reject("missing_ligand_evidence", "No curated positive evidence for ligand relevance"))
    elif config.require_ligand_evidence:
        issues.append(reject("missing_ligand_evidence", "No ligand evidence index was provided"))

    for field_name, value in {
        "holo_structure_file": holo_structure_file,
        "apo_structure_file": apo_structure_file,
        "ligand_file": ligand_file,
        "holo_chain_id": holo_chain_id,
        "apo_chain_id": apo_chain_id,
    }.items():
        if not value:
            issues.append(reject(f"missing_{field_name}", f"Candidate is missing {field_name}"))

    for label, resolution in [("holo", record.holo_resolution), ("apo", record.apo_resolution)]:
        if resolution is not None and resolution > config.max_resolution_a:
            issues.append(reject(f"{label}_low_resolution", f"{label} resolution {resolution:.2f} A exceeds {config.max_resolution_a:.2f} A"))

    if issues:
        record.issues = issues
        record.accepted = False
        return record

    holo_path = resolve_path(raw_root, holo_structure_file)
    apo_path = resolve_path(raw_root, apo_structure_file)
    ligand_path = resolve_path(raw_root, ligand_file)
    for label, path in [("holo_structure", holo_path), ("apo_structure", apo_path), ("ligand", ligand_path)]:
        if not path.exists():
            issues.append(reject(f"{label}_missing", f"Missing {label}: {path}"))
    if issues:
        record.issues = issues
        return record

    try:
        holo = extract_chain(holo_path, holo_chain_id)
        apo = extract_chain(apo_path, apo_chain_id)
        for msg in holo.warnings:
            issues.append(warn("holo_structure_warning", msg))
        for msg in apo.warnings:
            issues.append(warn("apo_structure_warning", msg))
    except Exception as exc:
        record.issues = issues + [reject("structure_parse_failed", repr(exc))]
        return record

    ligand = load_ligand(ligand_path, config)
    issues.extend(ligand.issues)
    record.metrics.update(ligand.metrics)

    if not holo.residues:
        issues.append(reject("holo_chain_empty", "Holo chain has no standard amino-acid residues"))
    if not apo.residues:
        issues.append(reject("apo_chain_empty", "Apo chain has no standard amino-acid residues"))
    if ligand.coords.size == 0:
        issues.append(reject("ligand_coords_empty", "Ligand has no usable heavy-atom coordinates"))
    if any(issue.severity == REJECT for issue in issues):
        record.issues = issues
        return record

    mapping, mapping_metrics = residue_mapping(apo, holo)
    record.metrics.update(mapping_metrics)
    record.metrics["apo_residues"] = len(apo.residues)
    record.metrics["holo_residues"] = len(holo.residues)
    record.metrics["common_residues"] = len(mapping)

    if len(mapping) < config.min_common_residues:
        issues.append(reject("too_few_common_residues", f"Only {len(mapping)} common residues"))
    if mapping_metrics["sequence_identity"] < config.min_sequence_identity:
        issues.append(reject("low_sequence_identity", f"Sequence identity {mapping_metrics['sequence_identity']:.3f}"))
    if mapping_metrics["mapping_coverage"] < config.min_mapping_coverage:
        issues.append(reject("low_mapping_coverage", f"Mapping coverage {mapping_metrics['mapping_coverage']:.3f}"))
    if any(issue.severity == REJECT for issue in issues):
        record.issues = issues
        return record

    apo_ca = np.asarray([apo.residues[i].atoms["CA"] for i, _ in mapping], dtype=np.float64)
    holo_ca = np.asarray([holo.residues[j].atoms["CA"] for _, j in mapping], dtype=np.float64)
    aligned_apo_ca, rot, trans = kabsch_align(apo_ca, holo_ca)
    ca_disp = np.sqrt(((aligned_apo_ca - holo_ca) ** 2).sum(axis=1))
    record.metrics["global_ca_rmsd_after_align"] = float(np.sqrt(np.mean(ca_disp * ca_disp)))

    pocket_labels = []
    pocket_displacements: list[float] = []
    contact_change_count = 0
    holo_contact_count = 0
    motion_core_count = 0
    for pair_idx, (apo_idx, holo_idx) in enumerate(mapping):
        apo_res = apo.residues[apo_idx]
        holo_res = holo.residues[holo_idx]
        aligned_apo_atoms = apply_transform(np.asarray(list(apo_res.atoms.values()), dtype=np.float64), rot, trans)
        holo_atoms = np.asarray(list(holo_res.atoms.values()), dtype=np.float64)
        apo_dist = min_distance(aligned_apo_atoms, ligand.coords)
        holo_dist = min_distance(holo_atoms, ligand.coords)
        in_pocket = holo_dist <= config.pocket_cutoff_a
        if not in_pocket:
            continue
        disp = float(ca_disp[pair_idx])
        holo_contact = holo_dist <= config.contact_cutoff_a
        apo_contact = apo_dist <= config.contact_cutoff_a
        contact_changed = (holo_contact != apo_contact) or abs(apo_dist - holo_dist) >= config.contact_change_cutoff_a
        motion_core = disp >= config.motion_core_cutoff_a
        pocket_displacements.append(disp)
        holo_contact_count += int(holo_contact)
        contact_change_count += int(contact_changed)
        motion_core_count += int(motion_core)
        pocket_labels.append(
            {
                "residue_key": holo_res.key,
                "resname": holo_res.resname,
                "holo_resseq": holo_res.resseq,
                "apo_resseq": apo_res.resseq,
                "ca_displacement": disp,
                "apo_ligand_min_dist": apo_dist,
                "holo_ligand_min_dist": holo_dist,
                "motion_core": motion_core,
                "holo_contact": holo_contact,
                "contact_changed": contact_changed,
            }
        )

    n_pocket = len(pocket_displacements)
    record.metrics.update(
        {
            "pocket_residues": n_pocket,
            "pocket_ca_rmsd": float(np.sqrt(np.mean(np.square(pocket_displacements)))) if pocket_displacements else None,
            "pocket_mean_ca_displacement": float(np.mean(pocket_displacements)) if pocket_displacements else None,
            "pocket_p90_ca_displacement": percentile(pocket_displacements, 90),
            "pocket_max_ca_displacement": max(pocket_displacements) if pocket_displacements else None,
            "motion_core_residues": motion_core_count,
            "motion_core_fraction": motion_core_count / n_pocket if n_pocket else None,
            "holo_contact_residues": holo_contact_count,
            "holo_contact_fraction": holo_contact_count / n_pocket if n_pocket else None,
            "contact_change_residues": contact_change_count,
            "contact_change_fraction": contact_change_count / n_pocket if n_pocket else None,
        }
    )
    record.labels["pocket_residues"] = pocket_labels

    if n_pocket < config.min_pocket_residues:
        issues.append(reject("too_few_pocket_residues", f"Only {n_pocket} pocket residues within {config.pocket_cutoff_a:.1f} A"))

    if apo.hetero_atoms.size:
        aligned_apo_hetero = apply_transform(apo.hetero_atoms, rot, trans)
        hetero_dist = min_distance(aligned_apo_hetero, ligand.coords)
        record.metrics["apo_nearest_hetero_to_holo_ligand"] = hetero_dist if math.isfinite(hetero_dist) else None
        if hetero_dist <= config.apo_pocket_hetero_cutoff_a:
            names = sorted(set(apo.hetero_names))[:8]
            issues.append(reject("apo_pocket_has_hetero_ligand", f"Apo chain has hetero atoms {names} within {hetero_dist:.2f} A of holo ligand site"))

    record.issues = issues
    record.accepted = not any(issue.severity == REJECT for issue in issues)
    return record
