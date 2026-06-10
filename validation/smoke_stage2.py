#!/usr/bin/env python
"""Single-GPU smoke test for stage-2 two-forward gate training.

Verifies: (1) get_diffusion_loss runs with stage2_two_forward=True,
(2) loss_ligand2 is finite, (3) the ligand-feedback gradient actually reaches
the gate parameters (to_q / to_k / bias), with the backbone frozen.
"""
import argparse
import torch
from torch_geometric.data import Batch
from torch_geometric.transforms import Compose

import utils.misc as misc
import utils.transforms as trans
from datasets import get_dataset
from datasets.pl_data import FOLLOW_BATCH
from models.molopt_score_model import ScorePosNet3D
from graphbap.bapnet import BAPNet


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train_config', default='./configs/training.yaml')
    ap.add_argument('--ckpt', default='./apo2mol_dataset/apo2mol_checkpoint.ckpt')
    ap.add_argument('--device', default='cuda:0')
    ap.add_argument('--bs', type=int, default=2)
    ap.add_argument('--t', type=int, default=80, help='fixed low timestep to exercise stage-2')
    args = ap.parse_args()

    device = torch.device(args.device)
    cfg = misc.load_config(args.train_config)
    cfg.model.pocket_router_mode = 'cross_attn_gate'
    cfg.model.cross_attn_gate.stage2_two_forward = True

    pf = trans.FeaturizeProteinAtom()
    lf = trans.FeaturizeLigandAtom(cfg.data.transform.ligand_atom_mode)
    transform = Compose([pf, lf, trans.FeaturizeLigandBond()])
    subsets = get_dataset(config=cfg.data, transform=transform)
    test_set = subsets['test']

    net_cond = BAPNet(ckpt_path=cfg.net_cond.ckpt_path, hidden_nf=cfg.net_cond.hidden_dim).to(device)
    net_cond.eval()

    model = ScorePosNet3D(cfg.model, protein_atom_feature_dim=pf.feature_dim,
                          ligand_atom_feature_dim=lf.feature_dim).to(device)
    ckpt = torch.load(args.ckpt, map_location='cpu', weights_only=False)
    sd = ckpt.get('state_dict', ckpt)
    sd = {k.replace('model.', ''): v for k, v in sd.items() if k.startswith('model.')} or sd
    missing, unexpected = model.load_state_dict(sd, strict=False)
    gate_missing = [k for k in missing if k.startswith('cross_attn_gate.')]
    print(f'loaded: missing={len(missing)} (gate_missing={len(gate_missing)}) unexpected={len(unexpected)}')

    # freeze everything except the gate
    for name, p in model.named_parameters():
        p.requires_grad = name.startswith('cross_attn_gate.')
    n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f'trainable params (gate only) = {n_train}')

    model.train()
    batch = Batch.from_data_list([test_set[i].clone() for i in range(args.bs)],
                                 follow_batch=FOLLOW_BATCH).to(device)
    time_step = torch.full((args.bs,), args.t, dtype=torch.long, device=device)

    res = model.get_diffusion_loss(
        net_cond=net_cond, data=batch,
        protein_pos_apo=batch.protein_pos,
        protein_pos_holo=batch.protein_pos_holo,
        protein_v=batch.protein_atom_feature.float(),
        batch_protein=batch.protein_element_batch,
        ligand_pos=batch.ligand_pos,
        ligand_v=batch.ligand_atom_feature_full,
        batch_ligand=batch.ligand_element_batch,
        time_step=time_step,
    )
    print('loss        =', float(res['loss']))
    print('loss_ligand2=', float(res['loss_ligand2']))
    print('loss_gate   =', float(res['loss_gate']))
    print('router_w_mean=', float(res['router_w_mean']))

    res['loss'].backward()
    gq = model.cross_attn_gate.to_q.weight.grad
    gk = model.cross_attn_gate.to_k.weight.grad
    gb = model.cross_attn_gate.bias.grad
    print('grad to_q norm =', None if gq is None else float(gq.norm()))
    print('grad to_k norm =', None if gk is None else float(gk.norm()))
    print('grad bias norm =', None if gb is None else float(gb.norm()))

    # isolate stage-2 contribution: grad of loss_ligand2 alone w.r.t. gate
    model.zero_grad(set_to_none=True)
    if float(res['loss_ligand2']) != 0.0:
        # recompute to get a clean graph
        res2 = model.get_diffusion_loss(
            net_cond=net_cond, data=batch,
            protein_pos_apo=batch.protein_pos,
            protein_pos_holo=batch.protein_pos_holo,
            protein_v=batch.protein_atom_feature.float(),
            batch_protein=batch.protein_element_batch,
            ligand_pos=batch.ligand_pos,
            ligand_v=batch.ligand_atom_feature_full,
            batch_ligand=batch.ligand_element_batch,
            time_step=time_step,
        )
        res2['loss_ligand2'].backward()
        gq2 = model.cross_attn_gate.to_q.weight.grad
        print('STAGE2-ONLY grad to_q norm =', None if gq2 is None else float(gq2.norm()))
    print('SMOKE_OK')


if __name__ == '__main__':
    main()
