"""Probe validation/loss=nan: load ckpt, run get_diffusion_loss on one val batch
across the same 10 t's the validation_step uses, print each loss component and
flag the first nan/inf source.

Usage (single-GPU, run on the AutoDL instance):

    cd /root/autodl-tmp/apo2mol
    python validation/diagnose_val_nan.py precision=32-true
    python validation/diagnose_val_nan.py precision=bf16-mixed

Both runs will be done; bf16 vs fp32 comparison localises the failure mode.
"""

import argparse
import sys
import os
import numpy as np
import torch
from omegaconf import OmegaConf

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils.transforms as trans
from datasets import get_dataset
from datasets.pl_data import FOLLOW_BATCH
from models.molopt_score_model import ScorePosNet3D
from models.pl_model import MoleculeTrainer
from graphbap.bapnet import BAPNet
from torch_geometric.loader import DataLoader
from torch_geometric.transforms import Compose


def is_bad(x):
    if x is None:
        return False
    if torch.is_tensor(x):
        return bool(torch.isnan(x).any() or torch.isinf(x).any())
    return not np.isfinite(float(x))


def fmt(x):
    if x is None:
        return "None"
    if torch.is_tensor(x):
        v = x.detach().float().item()
    else:
        v = float(x)
    return f"{v:.6g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precision", choices=["fp32", "bf16"], default="bf16",
                    help="dtype context for the forward (mixed precision uses autocast bf16)")
    ap.add_argument("--ckpt", default="./apo2mol_dataset/apo2mol_checkpoint.ckpt")
    ap.add_argument("--config", default="configs/training.yaml")
    ap.add_argument("--n_batches", type=int, default=2, help="val batches to sweep")
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.set_float32_matmul_precision("high")

    cfg = OmegaConf.load(args.config)
    # tiny val loader for fast probing
    cfg.train.val_batch_size = 2

    protein_featurizer = trans.FeaturizeProteinAtom()
    ligand_featurizer = trans.FeaturizeLigandAtom(cfg.data.transform.ligand_atom_mode)
    transform = Compose([
        protein_featurizer,
        ligand_featurizer,
        trans.FeaturizeLigandBond(),
    ])
    subsets = get_dataset(config=cfg.data, transform=transform)
    val_set = subsets["valid"]
    print(f"val set size: {len(val_set)}")

    val_loader = DataLoader(
        val_set,
        batch_size=cfg.train.val_batch_size,
        shuffle=False,
        follow_batch=FOLLOW_BATCH,
        exclude_keys=["ligand_nbh_list"],
        num_workers=0,
    )

    net_cond = BAPNet(ckpt_path=cfg.net_cond.ckpt_path, hidden_nf=cfg.net_cond.hidden_dim)
    net_cond.freeze_the_model()
    net_cond.to(device).eval()

    pl_model = MoleculeTrainer(cfg, protein_featurizer, ligand_featurizer, net_cond)
    ckpt = torch.load(args.ckpt, map_location="cpu", weights_only=False)
    state_dict = ckpt.get("state_dict", ckpt)
    missing, unexpected = pl_model.load_state_dict(state_dict, strict=False)
    print(f"loaded ckpt: missing={len(missing)} unexpected={len(unexpected)}")
    pl_model.to(device).eval()

    model = pl_model.model

    # 10 ts identical to validation_step
    ts = np.linspace(0, model.num_timesteps - 1, 10).astype(int)
    print(f"sweeping ts (10 points): {ts.tolist()}")
    print(f"precision: {args.precision}")
    print()

    autocast_ctx = (
        torch.amp.autocast(device_type="cuda", dtype=torch.bfloat16)
        if args.precision == "bf16"
        else torch.amp.autocast(device_type="cuda", enabled=False)
    )

    for batch_idx, batch in enumerate(val_loader):
        if batch_idx >= args.n_batches:
            break
        batch = batch.to(device)
        bs = batch.num_graphs
        print(f"=== batch {batch_idx}  (graphs={bs}) ===")
        header = f"{'t':>5} | {'loss':>11} {'lig_pos':>11} {'loss_v':>11} {'prot_tr':>11} {'prot_rot':>11} {'prot_chi':>11} {'gate':>11} | bad_in"
        print(header)
        print("-" * len(header))
        for t in ts:
            t_step = torch.tensor([int(t)] * bs, device=device)
            with torch.no_grad():
                with autocast_ctx:
                    res = model.get_diffusion_loss(
                        net_cond=pl_model.net_cond,
                        data=batch,
                        protein_pos_apo=batch.protein_pos,
                        protein_pos_holo=batch.protein_pos_holo,
                        protein_v=batch.protein_atom_feature.float(),
                        batch_protein=batch.protein_element_batch,
                        ligand_pos=batch.ligand_pos,
                        ligand_v=batch.ligand_atom_feature_full,
                        batch_ligand=batch.ligand_element_batch,
                        time_step=t_step,
                    )
            comps = {
                "loss":     res["loss"],
                "lig_pos":  res["loss_ligang_pos"],
                "loss_v":   res["loss_v"],
                "prot_tr":  res["loss_protein_tr"],
                "prot_rot": res["loss_protein_rot"],
                "prot_chi": res["loss_protein_chi"],
                "gate":     res.get("loss_gate", None),
            }
            bad = [k for k, v in comps.items() if is_bad(v)]
            print(f"{int(t):>5} | "
                  f"{fmt(comps['loss']):>11} {fmt(comps['lig_pos']):>11} {fmt(comps['loss_v']):>11} "
                  f"{fmt(comps['prot_tr']):>11} {fmt(comps['prot_rot']):>11} {fmt(comps['prot_chi']):>11} "
                  f"{fmt(comps['gate']):>11} | {','.join(bad) if bad else '-'}")
        print()


if __name__ == "__main__":
    main()
