"""Probe stage-1 (cross_attn_gate + freeze_backbone) val path; dumps each loss
component across the 10 t's that validation_step sweeps. Run on the AutoDL
instance from /root/autodl-tmp/apo2mol."""

import os
import sys
sys.path.insert(0, ".")
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import torch
from omegaconf import OmegaConf
import utils.transforms as trans
from datasets import get_dataset
from datasets.pl_data import FOLLOW_BATCH
from models.pl_model import MoleculeTrainer
from graphbap.bapnet import BAPNet
from torch_geometric.loader import DataLoader
from torch_geometric.transforms import Compose


def fmt(x):
    if x is None:
        return "  None     "
    return f"{x.detach().float().item():11.4g}"


def main():
    torch.set_float32_matmul_precision("high")
    device = torch.device("cuda")
    cfg = OmegaConf.load("configs/training.yaml")

    cfg.model.pocket_router_mode = "cross_attn_gate"
    cfg.model.cross_attn_gate.stage2_two_forward = False
    cfg.train.freeze_backbone = True
    cfg.train.val_batch_size = 2
    cfg.train.batch_size = 2

    protein_featurizer = trans.FeaturizeProteinAtom()
    ligand_featurizer = trans.FeaturizeLigandAtom(cfg.data.transform.ligand_atom_mode)
    transform = Compose([protein_featurizer, ligand_featurizer, trans.FeaturizeLigandBond()])
    subsets = get_dataset(config=cfg.data, transform=transform)
    val_set = subsets["valid"]
    val_loader = DataLoader(
        val_set, batch_size=2, shuffle=False, follow_batch=FOLLOW_BATCH,
        exclude_keys=["ligand_nbh_list"], num_workers=0,
    )

    net_cond = BAPNet(ckpt_path=cfg.net_cond.ckpt_path, hidden_nf=cfg.net_cond.hidden_dim)
    net_cond.freeze_the_model()
    net_cond.to(device).eval()

    pl_model = MoleculeTrainer(cfg, protein_featurizer, ligand_featurizer, net_cond)
    ckpt = torch.load("./apo2mol_dataset/apo2mol_checkpoint.ckpt", map_location="cpu", weights_only=False)
    state_dict = ckpt.get("state_dict", ckpt)
    missing, unexpected = pl_model.load_state_dict(state_dict, strict=False)
    print(f"loaded ckpt: missing={len(missing)} examples={missing[:3]}, unexpected={len(unexpected)}")
    pl_model.to(device).eval()
    model = pl_model.model

    ga = model.cross_attn_gate
    print("gate q.weight norm:", ga.to_q.weight.float().norm().item())
    print("gate k.weight norm:", ga.to_k.weight.float().norm().item())
    print("gate bias        :", ga.bias.float().item())

    ts = np.linspace(0, model.num_timesteps - 1, 10).astype(int)

    for batch_idx, batch in enumerate(val_loader):
        if batch_idx >= 3:
            break
        batch = batch.to(device)
        bs = batch.num_graphs
        print(f"=== batch {batch_idx} (graphs={bs}) ===")
        print(f"{'t':>4} {'loss':>11} {'lig_pos':>11} {'loss_v':>11} {'p_tr':>11} {'p_rot':>11} {'p_chi':>11} {'gate':>11}  bad")
        for t in ts:
            ts_t = torch.tensor([int(t)] * bs, device=device)
            with torch.no_grad(), torch.amp.autocast(device_type="cuda", dtype=torch.bfloat16):
                res = model.get_diffusion_loss(
                    net_cond=pl_model.net_cond, data=batch,
                    protein_pos_apo=batch.protein_pos,
                    protein_pos_holo=batch.protein_pos_holo,
                    protein_v=batch.protein_atom_feature.float(),
                    batch_protein=batch.protein_element_batch,
                    ligand_pos=batch.ligand_pos,
                    ligand_v=batch.ligand_atom_feature_full,
                    batch_ligand=batch.ligand_element_batch,
                    time_step=ts_t,
                )
            comps = {
                "loss": res["loss"], "lig_pos": res["loss_ligang_pos"],
                "loss_v": res["loss_v"],
                "p_tr": res["loss_protein_tr"], "p_rot": res["loss_protein_rot"],
                "p_chi": res["loss_protein_chi"],
                "gate": res.get("loss_gate", None),
            }
            bad = [k for k, v in comps.items()
                   if v is not None and (torch.isnan(v).any() or torch.isinf(v).any())]
            print(f"{int(t):4d} {fmt(comps['loss'])} {fmt(comps['lig_pos'])} "
                  f"{fmt(comps['loss_v'])} {fmt(comps['p_tr'])} {fmt(comps['p_rot'])} "
                  f"{fmt(comps['p_chi'])} {fmt(comps['gate'])}  {bad}")


if __name__ == "__main__":
    main()
