import math

import torch
import torch.nn as nn


class LigandResidueCrossAttnGate(nn.Module):
    """Ligand-as-query cross-attention into pocket residues.

    For each graph in the batch we compute attn = softmax(Q K^T / sqrt(d))
    with Q from ligand atoms and K from residues, then pool over the ligand
    axis with max to obtain a per-residue gate logit. The gate w_i is
    sigmoid(logit), aligned 1:1 with `protein_translations_batch` order.

    The gate is trained with an induced-fit supervision loss (see
    ScorePosNet3D._build_induced_fit_target / get_diffusion_loss): the logits
    are supervised against a continuous soft target derived from the
    ground-truth apo->holo residue displacement. The same w_i is also
    multiplied into residue_tr / chi / rot so train- and sampling-time
    behaviour stay consistent.
    """

    def __init__(self, ligand_dim, residue_dim, inner_dim=64):
        super().__init__()
        self.scale = 1.0 / math.sqrt(inner_dim)
        self.to_q = nn.Linear(ligand_dim, inner_dim, bias=False)
        self.to_k = nn.Linear(residue_dim, inner_dim, bias=False)
        self.bias = nn.Parameter(torch.zeros(1))

    def forward(self, ligand_h, residue_h, batch_ligand, prot_update_batch):
        device = residue_h.device
        num_residues = residue_h.shape[0]
        logits = residue_h.new_zeros(num_residues)

        if num_residues == 0:
            return logits

        q_all = self.to_q(ligand_h)
        k_all = self.to_k(residue_h)
        num_graphs = int(prot_update_batch.max().item()) + 1
        for g in range(num_graphs):
            res_idx = (prot_update_batch == g).nonzero(as_tuple=True)[0]
            lig_idx = (batch_ligand == g).nonzero(as_tuple=True)[0]
            if res_idx.numel() == 0 or lig_idx.numel() == 0:
                continue
            q = q_all[lig_idx]
            k = k_all[res_idx]
            sim = q @ k.transpose(0, 1) * self.scale
            # for each residue, take the strongest ligand atom's attention logit
            per_res = sim.max(dim=0).values + self.bias
            # under bf16 autocast the fp32 bias promotes per_res to fp32; cast
            # back to the logits buffer dtype so index_copy dtypes match.
            logits = logits.index_copy(0, res_idx, per_res.to(logits.dtype))
        return logits
