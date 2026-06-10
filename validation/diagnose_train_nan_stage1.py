"""Run a real training_step under stage-1 (cross_attn_gate + freeze) for a few
batches with detect_anomaly + per-step parameter health checks. Goal: locate
which optimizer step pushes a parameter to nan, and which sub-module was
responsible (forward output, gradient, or update)."""

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


def health(model, label):
    bad_params, bad_grads = [], []
    for n, p in model.named_parameters():
        if not p.requires_grad:
            continue
        if torch.isnan(p).any() or torch.isinf(p).any():
            bad_params.append(n)
        if p.grad is not None and (torch.isnan(p.grad).any() or torch.isinf(p.grad).any()):
            bad_grads.append(n)
    print(f"  [{label}] bad params: {bad_params[:5]}  bad grads: {bad_grads[:5]}")
    return bad_params, bad_grads


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
    train_set = subsets["train"]

    train_loader = DataLoader(
        train_set, batch_size=2, shuffle=True, follow_batch=FOLLOW_BATCH,
        exclude_keys=["ligand_nbh_list"], num_workers=0,
    )

    net_cond = BAPNet(ckpt_path=cfg.net_cond.ckpt_path, hidden_nf=cfg.net_cond.hidden_dim)
    net_cond.freeze_the_model()
    net_cond.to(device).eval()

    pl_model = MoleculeTrainer(cfg, protein_featurizer, ligand_featurizer, net_cond)
    ckpt = torch.load("./apo2mol_dataset/apo2mol_checkpoint.ckpt", map_location="cpu", weights_only=False)
    state_dict = ckpt.get("state_dict", ckpt)
    missing, unexpected = pl_model.load_state_dict(state_dict, strict=False)
    print(f"loaded ckpt: missing={len(missing)}; unexpected={len(unexpected)}")
    pl_model.to(device).train()
    # _apply_freeze_policy already ran in __init__
    trainable = [p for p in pl_model.parameters() if p.requires_grad]
    print(f"trainable param tensors: {len(trainable)}; total elements: {sum(p.numel() for p in trainable)}")
    for n, p in pl_model.named_parameters():
        if p.requires_grad:
            print(f"  trainable: {n} ({tuple(p.shape)})")

    opt = torch.optim.Adam(trainable, lr=cfg.train.optimizer.lr,
                           betas=(cfg.train.optimizer.beta1, cfg.train.optimizer.beta2),
                           weight_decay=cfg.train.optimizer.weight_decay)

    n_steps = 5
    for step, batch in enumerate(train_loader):
        if step >= n_steps:
            break
        batch = batch.to(device)
        opt.zero_grad(set_to_none=True)
        with torch.amp.autocast(device_type="cuda", dtype=torch.bfloat16):
            res = pl_model.model.get_diffusion_loss(
                net_cond=pl_model.net_cond, data=batch,
                protein_pos_apo=batch.protein_pos,
                protein_pos_holo=batch.protein_pos_holo,
                protein_v=batch.protein_atom_feature.float(),
                batch_protein=batch.protein_element_batch,
                ligand_pos=batch.ligand_pos,
                ligand_v=batch.ligand_atom_feature_full,
                batch_ligand=batch.ligand_element_batch,
            )
        loss = res["loss"]
        gate = res.get("loss_gate", None)
        comps = {k: res[k].detach().float().item() for k in
                 ["loss", "loss_ligang_pos", "loss_v",
                  "loss_protein_tr", "loss_protein_rot", "loss_protein_chi"]
                 if k in res}
        if gate is not None:
            comps["loss_gate"] = gate.detach().float().item()
        finite_loss = torch.isfinite(loss).all().item()
        print(f"step {step}: loss={loss.float().item():.4g} finite={finite_loss} "
              f"comps={ {k: f'{v:.3g}' for k, v in comps.items()} }")

        if not finite_loss:
            print("  -> loss already nan/inf in forward; aborting")
            break

        loss.backward()
        # grad-norm
        total_grad_norm = 0.0
        for p in trainable:
            if p.grad is not None:
                total_grad_norm += p.grad.detach().float().pow(2).sum().item()
        total_grad_norm = total_grad_norm ** 0.5
        print(f"  grad_norm = {total_grad_norm:.4g}")
        bad_p, bad_g = health(pl_model, "post-backward")
        if bad_p or bad_g:
            print("  -> bad params/grads; aborting")
            break

        torch.nn.utils.clip_grad_norm_(trainable, max_norm=cfg.train.max_grad_norm)
        opt.step()
        bad_p2, _ = health(pl_model, "post-step")
        if bad_p2:
            print("  -> opt.step() pushed params to nan; aborting")
            break


if __name__ == "__main__":
    main()
