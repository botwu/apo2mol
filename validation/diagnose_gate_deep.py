"""Deep probe: track gate logits, target, BCE input range, loss_gate, and gate
param norms across each train step. fp32 mode (autocast off) so we see the
underlying instability without bf16 noise."""

import os, sys
sys.path.insert(0, ".")
import warnings; warnings.filterwarnings("ignore")
import torch
import torch.nn.functional as F
from omegaconf import OmegaConf
import utils.transforms as trans
from datasets import get_dataset
from datasets.pl_data import FOLLOW_BATCH
from models.pl_model import MoleculeTrainer
from graphbap.bapnet import BAPNet
from torch_geometric.loader import DataLoader
from torch_geometric.transforms import Compose


def main():
    torch.set_float32_matmul_precision("high")
    torch.manual_seed(2021)
    device = torch.device("cuda")
    cfg = OmegaConf.load("configs/training.yaml")
    cfg.model.pocket_router_mode = "cross_attn_gate"
    cfg.model.cross_attn_gate.stage2_two_forward = False
    cfg.train.freeze_backbone = True
    cfg.train.batch_size = 2

    protein_featurizer = trans.FeaturizeProteinAtom()
    ligand_featurizer = trans.FeaturizeLigandAtom(cfg.data.transform.ligand_atom_mode)
    transform = Compose([protein_featurizer, ligand_featurizer, trans.FeaturizeLigandBond()])
    subsets = get_dataset(config=cfg.data, transform=transform)
    loader = DataLoader(subsets["train"], batch_size=2, shuffle=True,
                        follow_batch=FOLLOW_BATCH, exclude_keys=["ligand_nbh_list"], num_workers=0)

    net_cond = BAPNet(ckpt_path=cfg.net_cond.ckpt_path, hidden_nf=cfg.net_cond.hidden_dim)
    net_cond.freeze_the_model()
    net_cond.to(device).eval()
    pl_model = MoleculeTrainer(cfg, protein_featurizer, ligand_featurizer, net_cond)
    ckpt = torch.load("./apo2mol_dataset/apo2mol_checkpoint.ckpt", map_location="cpu", weights_only=False)
    pl_model.load_state_dict(ckpt["state_dict"], strict=False)
    pl_model.to(device).train()
    trainable = [p for p in pl_model.parameters() if p.requires_grad]
    opt = torch.optim.Adam(trainable, lr=5e-4)
    model = pl_model.model

    def gate_norms():
        ga = model.cross_attn_gate
        return (
            ga.bias.float().item(),
            ga.to_q.weight.float().norm().item(),
            ga.to_k.weight.float().norm().item(),
        )

    for step, batch in enumerate(loader):
        if step >= 6:
            break
        batch = batch.to(device)

        # === forward up to gate (no autocast) ===
        opt.zero_grad(set_to_none=True)

        # Inspect raw input features feeding the gate
        # We piggyback on get_diffusion_loss but also peel out gate intermediates
        # by monkey-patching the gate's forward to capture logits
        captured = {}
        orig_forward = model.cross_attn_gate.forward
        def capturing_forward(*a, **kw):
            out = orig_forward(*a, **kw)
            captured["logits"] = out
            return out
        model.cross_attn_gate.forward = capturing_forward
        try:
            res = model.get_diffusion_loss(
                net_cond=pl_model.net_cond, data=batch,
                protein_pos_apo=batch.protein_pos,
                protein_pos_holo=batch.protein_pos_holo,
                protein_v=batch.protein_atom_feature.float(),
                batch_protein=batch.protein_element_batch,
                ligand_pos=batch.ligand_pos,
                ligand_v=batch.ligand_atom_feature_full,
                batch_ligand=batch.ligand_element_batch,
            )
        finally:
            model.cross_attn_gate.forward = orig_forward

        loss = res["loss"]
        lg = res["loss_gate"].item()
        logits = captured["logits"]
        b, q, k = gate_norms()

        print(f"\n=== step {step} ===")
        print(f"  pre-step gate norms: bias={b:+.4g} q.w={q:.4g} k.w={k:.4g}")
        print(f"  logits: shape={tuple(logits.shape)} min={logits.min().item():.4g} "
              f"max={logits.max().item():.4g} mean={logits.mean().item():.4g} "
              f"|abs_max|={logits.abs().max().item():.4g} dtype={logits.dtype}")
        print(f"  loss={loss.item():.4g} loss_gate={lg:.4g} "
              f"loss_ligang_pos={res['loss_ligang_pos'].item():.4g} "
              f"loss_protein_tr={res['loss_protein_tr'].item():.4g} "
              f"loss_protein_rot={res['loss_protein_rot'].item():.4g}")

        
        
        loss.backward()
        # gate-only grads
        for n, p in pl_model.named_parameters():
            if not p.requires_grad:
                continue
            g = p.grad
            gnan = torch.isnan(g).any().item() if g is not None else False
            ginf = torch.isinf(g).any().item() if g is not None else False
            gnorm = g.float().norm().item() if g is not None else float("nan")
            print(f"  grad[{n}] norm={gnorm:.4g} nan={gnan} inf={ginf}")
        if any(torch.isnan(p.grad).any().item() for p in trainable if p.grad is not None):
            break
        opt.step()


if __name__ == "__main__":
    main()
