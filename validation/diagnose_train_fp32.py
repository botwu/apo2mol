"""Same as diagnose_train_nan_stage1.py but force fp32 (no autocast). Verifies
the bf16-is-the-culprit hypothesis."""

import os, sys
sys.path.insert(0, ".")
import warnings; warnings.filterwarnings("ignore")
import torch
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
    print(f"trainable: {[n for n,p in pl_model.named_parameters() if p.requires_grad]}")

    for step, batch in enumerate(loader):
        if step >= 5:
            break
        batch = batch.to(device)
        opt.zero_grad(set_to_none=True)
        # NO autocast — pure fp32
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
        loss.backward()
        gn = sum(p.grad.float().pow(2).sum().item() for p in trainable if p.grad is not None) ** 0.5
        bad = [n for n, p in pl_model.named_parameters()
               if p.grad is not None and (torch.isnan(p.grad).any() or torch.isinf(p.grad).any())]
        print(f"fp32 step {step}: loss={loss.item():.4g} grad_norm={gn:.4g} bad_grads={bad[:3]}")
        if bad:
            break
        opt.step()


if __name__ == "__main__":
    main()
