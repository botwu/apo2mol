#!/usr/bin/env python
"""Paper-aligned ligand evaluation for one A/B arm.

Unlike eval_split.py (which keeps only the best-Vina molecule per pocket),
this script keeps EVERY generated sample so the reported Vina/QED/SA match
the paper's "report mean & median over all generated ligands" protocol, and
it additionally computes High Affinity = fraction of generated ligands whose
Vina min beats the reference ligand docked into its native holo pocket.

Design notes:
- Strictly single-process; the eval box has a ~2GB cgroup memory cap, so we
  never parallelize docking.
- Per-sample try/except: one bad molecule never aborts the arm.
- Reference Vina is cached per case (shared across arms if you reuse the
  --ref-cache file).

Usage:
  PYTHONPATH=<repo> python validation/eval_ligand_full.py \
      --sample_path .../sampled_results \
      --out .../ligand_full.json \
      --protein_root ./apo2mol_dataset/data_folder \
      --ref-cache .../ref_vina_cache.json
"""
import argparse
import json
import os
import statistics as st
import traceback
from glob import glob

import numpy as np
import torch
from rdkit import Chem
from rdkit import RDLogger
from Bio.PDB import PDBParser, PDBIO

from utils import reconstruct, transforms
from utils.evaluation import scoring_func
from utils.evaluation.docking_vina import VinaDockingTask

RDLogger.DisableLog("rdApp.*")


def reconstruct_mol(pred_pos, pred_v, atom_mode):
    atom_type = transforms.get_atomic_number_from_index(pred_v, mode=atom_mode)
    aromatic = transforms.is_aromatic_from_index(pred_v, mode=atom_mode)
    mol = reconstruct.reconstruct_from_generated(pred_pos, atom_type, aromatic)
    return mol


def write_gen_pocket_pdb(protein_root, ligand_filename, pred_protein_pos, out_pdb):
    protein_fn = os.path.join(os.path.dirname(ligand_filename), "receptor_apo_pocket10.pdb")
    protein_path = os.path.join(protein_root, protein_fn)
    structure = PDBParser(QUIET=True).get_structure("p", protein_path)
    for i, atom in enumerate(structure.get_atoms()):
        if i < len(pred_protein_pos):
            atom.set_coord(pred_protein_pos[i])
        else:
            break
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb)


def dock_minimize(mol, ligand_filename, protein_root, pocket_type, pred_protein_path=None):
    task = VinaDockingTask.from_generated_mol(
        mol, ligand_filename, pred_protein_path,
        protein_root=protein_root, pocket_type=pocket_type)
    res = task.run(mode="minimize", exhaustiveness=8)
    return float(res[0]["affinity"])


def reference_vina(protein_root, ligand_filename, work_dir, ref_cache):
    if ligand_filename in ref_cache:
        return ref_cache[ligand_filename]
    val = None
    try:
        ref_path = os.path.join(protein_root, ligand_filename)
        ref_mol = next(iter(Chem.SDMolSupplier(ref_path, sanitize=True)))
        if ref_mol is not None:
            val = dock_minimize(ref_mol, ligand_filename, protein_root, pocket_type="holo")
    except Exception:
        traceback.print_exc()
        val = None
    ref_cache[ligand_filename] = val
    return val


def summarize(values):
    vals = [v for v in values if v is not None and not (isinstance(v, float) and np.isnan(v))]
    if not vals:
        return {"n": 0, "mean": None, "median": None}
    return {"n": len(vals), "mean": st.mean(vals), "median": st.median(vals)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_path", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--protein_root", default="./apo2mol_dataset/data_folder")
    ap.add_argument("--atom_enc_mode", default="add_aromatic")
    ap.add_argument("--pocket_type", default="gen", choices=["gen", "holo", "apo"])
    ap.add_argument("--ref-cache", default=None)
    ap.add_argument("--work-dir", default="./tmp_ligand_full")
    args = ap.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)
    ref_cache = {}
    if args.ref_cache and os.path.exists(args.ref_cache):
        with open(args.ref_cache) as f:
            ref_cache = json.load(f)

    files = sorted(glob(os.path.join(args.sample_path, "*result_*.pt")),
                   key=lambda x: int(os.path.basename(x)[:-3].split("_")[-1]))

    per_sample = []  # dict per generated sample
    n_recon = 0
    n_complete = 0
    n_total = 0

    for fpath in files:
        try:
            r = torch.load(fpath, weights_only=False)
        except Exception:
            traceback.print_exc()
            continue
        data = r["data"]
        ligand_filename = data.ligand_filename
        case_id = ligand_filename.replace(".sdf", "").split("/")[0]

        lig_pos_list = r["pred_ligand_ligand_pos_traj"]
        lig_v_list = r["pred_ligand_v_traj"]
        prot_pos_list = r["pred_protein_pos_traj"]

        ref_v = reference_vina(args.protein_root, ligand_filename, args.work_dir, ref_cache)

        for s_idx in range(len(lig_pos_list)):
            n_total += 1
            rec = {"case": case_id, "sample": s_idx, "ref_vina": ref_v,
                   "vina": None, "qed": None, "sa": None, "logp": None, "lipinski": None}
            try:
                pred_pos = lig_pos_list[s_idx][-1]
                pred_v = lig_v_list[s_idx][-1]
                pred_protein_pos = prot_pos_list[s_idx][-1]
                mol = reconstruct_mol(pred_pos, pred_v, args.atom_enc_mode)
                smiles = Chem.MolToSmiles(mol)
                n_recon += 1
                if "." in smiles:
                    per_sample.append(rec)
                    continue
                n_complete += 1
                chem = scoring_func.get_chem(mol)
                rec["qed"] = float(chem["qed"])
                rec["sa"] = float(chem["sa"])
                rec["logp"] = float(chem["logp"])
                rec["lipinski"] = int(chem["lipinski"])

                pred_pdb = os.path.join(args.work_dir, f"{case_id}_{s_idx}_genpocket.pdb")
                write_gen_pocket_pdb(args.protein_root, ligand_filename, pred_protein_pos, pred_pdb)
                rec["vina"] = dock_minimize(mol, ligand_filename, args.protein_root,
                                            pocket_type=args.pocket_type, pred_protein_path=pred_pdb)
            except Exception:
                traceback.print_exc()
            per_sample.append(rec)
        print(f"[{case_id}] ref_vina={ref_v} samples={len(lig_pos_list)}", flush=True)

    if args.ref_cache:
        with open(args.ref_cache, "w") as f:
            json.dump(ref_cache, f, indent=2)

    vina_vals = [d["vina"] for d in per_sample]
    qed_vals = [d["qed"] for d in per_sample]
    sa_vals = [d["sa"] for d in per_sample]
    logp_vals = [d["logp"] for d in per_sample]
    lip_vals = [d["lipinski"] for d in per_sample]

    # High Affinity over samples that have both vina and ref_vina
    pairs = [(d["vina"], d["ref_vina"]) for d in per_sample
             if d["vina"] is not None and d["ref_vina"] is not None]
    high_aff = (sum(1 for v, rv in pairs if v < rv) / len(pairs)) if pairs else None

    summary = {
        "n_results": len(files),
        "n_total_samples": n_total,
        "n_recon": n_recon,
        "n_complete": n_complete,
        "n_docked": summarize(vina_vals)["n"],
        "vina_min": summarize(vina_vals),
        "qed": summarize(qed_vals),
        "sa": summarize(sa_vals),
        "logp": summarize(logp_vals),
        "lipinski": summarize(lip_vals),
        "high_affinity": high_aff,
        "high_affinity_n_pairs": len(pairs),
    }
    out = {"summary": summary, "per_sample": per_sample}
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
