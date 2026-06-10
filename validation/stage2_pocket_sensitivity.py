#!/usr/bin/env python
"""Stage-2 go/no-go diagnostic: is the frozen ligand head sensitive to the pocket?

Stage-2 (two-forward) only helps the ligand if, holding everything else fixed,
feeding the frozen backbone a BETTER pocket produces a BETTER ligand x0
prediction. This script measures exactly that, with zero training:

  for each (case, timestep t, noise draw):
    perturb the ligand to level t (same as training)
    compute the interaction-prior features (hbap) ONCE from GT (shared)
    run forward TWICE, changing only protein_pos:
        cond "apo"  : protein_pos = apo pocket   (worst)
        cond "holo" : protein_pos = holo pocket  (oracle best)
    record ligand x0 reconstruction error vs GT clean ligand.

If err(holo) < err(apo) systematically (paired), the ligand head responds to
pocket quality -> stage-2 has headroom (GO). If err(holo) ~= err(apo), the
frozen head ignores the pocket -> stage-2 cannot lift ligand metrics (NO-GO).

Single-process, eval-only. Reuses sample_split.build_worker_state for an
identical model/net_cond load.
"""
import argparse
import json
import os
import statistics as st

import numpy as np
import torch
from torch_geometric.data import Batch
from torch_geometric.transforms import Compose

import utils.misc as misc
import utils.transforms as trans
from datasets import get_dataset
from datasets.pl_data import FOLLOW_BATCH
from models.molopt_score_model import center_pos, index_to_log_onehot
from sample_split import build_worker_state


@torch.no_grad()
def ligand_x0_error(model, net_cond, data_batch, protein_pos_cond, ligand_pos_perturbed,
                    ligand_v_perturbed, hbap_protein, hbap_ligand, time_step,
                    gt_ligand_pos, batch_ligand):
    preds = model.forward(
        protein_pos=protein_pos_cond,
        protein_v=data_batch.protein_atom_feature.float(),
        batch_protein=data_batch.protein_element_batch,
        init_ligand_pos=ligand_pos_perturbed,
        init_ligand_v=ligand_v_perturbed,
        batch_ligand=batch_ligand,
        protein_atom_to_aa_group=data_batch.protein_atom_to_aa_group,
        time_step=time_step,
        hbap_protein=hbap_protein,
        hbap_ligand=hbap_ligand,
    )
    pred = preds['pred_ligand_pos']
    if model.model_mean_type == 'C0':
        pred_x0 = pred
    else:  # 'noise'
        eps = pred - ligand_pos_perturbed
        pred_x0 = model._predict_x0_from_eps(xt=ligand_pos_perturbed, eps=eps,
                                             t=time_step, batch=batch_ligand)
    # per-atom L2 distance to GT clean ligand, averaged over atoms
    per_atom = torch.sqrt(((pred_x0 - gt_ligand_pos) ** 2).sum(-1) + 1e-12)
    return float(per_atom.mean().item())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='./configs/sampling.yaml')
    ap.add_argument('--train_config', default='./configs/training.yaml')
    ap.add_argument('--device', default='cuda:0')
    ap.add_argument('--ids_file', default=None,
                    help='test-set positions, one per line; default = first --num_cases')
    ap.add_argument('--num_cases', type=int, default=30)
    ap.add_argument('--timesteps', default='100,300,500,700,900')
    ap.add_argument('--n_noise', type=int, default=3)
    ap.add_argument('--seed', type=int, default=2021)
    ap.add_argument('--out', default='./validation/stage2_sensitivity.json')
    args = ap.parse_args()

    misc.seed_all(args.seed)
    config = misc.load_config(args.config)
    train_config = misc.load_config(args.train_config)
    # Use the base pretrained backbone (no gate) for this diagnostic.
    train_config.model.pocket_router_mode = 'none'

    protein_featurizer = trans.FeaturizeProteinAtom()
    ligand_atom_mode = train_config.data.transform.ligand_atom_mode
    ligand_featurizer = trans.FeaturizeLigandAtom(ligand_atom_mode)
    transform = Compose([protein_featurizer, ligand_featurizer, trans.FeaturizeLigandBond()])
    subsets = get_dataset(config=train_config.data, transform=transform)
    test_set = subsets['test']
    num_test = len(test_set)
    print(f'test set size {num_test}', flush=True)

    if args.ids_file and os.path.exists(args.ids_file):
        with open(args.ids_file) as f:
            ids = [int(l.strip()) for l in f if l.strip() and not l.startswith('#')]
    else:
        ids = list(range(min(args.num_cases, num_test)))
    ids = ids[:args.num_cases]
    timesteps = [int(x) for x in args.timesteps.split(',')]
    print(f'cases={len(ids)} timesteps={timesteps} n_noise={args.n_noise}', flush=True)

    worker_cfg = {'train_config': train_config, 'config': config, 'args': args, 'result_path': '.'}
    state = build_worker_state(args.device, worker_cfg)
    model = state['model']
    net_cond = state['net_cond']
    device = state['device']
    model.eval()
    net_cond.eval()

    records = []  # one per (case,t,noise): {case,t,err_apo,err_holo}
    for ci, data_id in enumerate(ids):
        data = test_set[data_id]
        batch = Batch.from_data_list([data.clone()], follow_batch=FOLLOW_BATCH).to(device)
        batch_protein = batch.protein_element_batch
        batch_ligand = batch.ligand_element_batch

        apo, holo, ligand_pos_c, offset = center_pos(
            batch.protein_pos, batch.protein_pos_holo, batch.ligand_pos,
            batch_protein, batch_ligand, mode=config.sample.center_pos_mode)

        protein_v = batch.protein_atom_feature.float()
        ligand_v = batch.ligand_atom_feature_full
        # interaction-prior features from GT (holo pocket + clean ligand), shared
        gt_p_a = torch.argmax(protein_v[:, :6], dim=1)
        gt_p_r = torch.argmax(protein_v[:, 6:26], dim=1)
        hbap_ligand0, hbap_protein0 = net_cond.extract_features(
            ligand_pos_c, holo, ligand_v, gt_p_a, gt_p_r, batch_ligand, batch_protein)

        for t_val in timesteps:
            time_step = torch.full((1,), t_val, dtype=torch.long, device=device)
            a = model.alphas_cumprod.index_select(0, time_step)
            a_pos = a[batch_ligand].unsqueeze(-1)
            log_v0 = index_to_log_onehot(ligand_v, model.num_classes)
            for _ in range(args.n_noise):
                noise = torch.randn_like(ligand_pos_c)
                lig_pert = a_pos.sqrt() * ligand_pos_c + (1.0 - a_pos).sqrt() * noise
                lig_v_pert, _ = model.q_v_sample(log_v0, time_step, batch_ligand)
                err_apo = ligand_x0_error(model, net_cond, batch, apo, lig_pert, lig_v_pert,
                                          hbap_protein0, hbap_ligand0, time_step, ligand_pos_c, batch_ligand)
                err_holo = ligand_x0_error(model, net_cond, batch, holo, lig_pert, lig_v_pert,
                                           hbap_protein0, hbap_ligand0, time_step, ligand_pos_c, batch_ligand)
                records.append({'case': int(data_id), 't': t_val,
                                'err_apo': err_apo, 'err_holo': err_holo})
        print(f'[{ci+1}/{len(ids)}] case {data_id} done', flush=True)

    # aggregate
    def agg(rs):
        ea = [r['err_apo'] for r in rs]
        eh = [r['err_holo'] for r in rs]
        deltas = [r['err_apo'] - r['err_holo'] for r in rs]  # >0 means holo better
        return {
            'n': len(rs),
            'err_apo_mean': st.mean(ea), 'err_holo_mean': st.mean(eh),
            'delta_mean': st.mean(deltas), 'delta_median': st.median(deltas),
            'holo_better_frac': sum(1 for d in deltas if d > 0) / len(deltas),
            'rel_improve_pct': 100.0 * st.mean(deltas) / st.mean(ea),
        }

    overall = agg(records)
    by_t = {t: agg([r for r in records if r['t'] == t]) for t in timesteps}
    out = {'overall': overall, 'by_timestep': by_t, 'n_records': len(records),
           'n_cases': len(ids), 'records': records}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=2)

    print('\n==== STAGE-2 POCKET-SENSITIVITY ====', flush=True)
    print(json.dumps({'overall': overall, 'by_timestep': by_t}, indent=2), flush=True)


if __name__ == '__main__':
    main()
