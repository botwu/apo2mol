import math

import torch
import torch.nn as nn


class SinusoidalTimeEmbedding(nn.Module):
    def __init__(self, dim, max_period=10000):
        super().__init__()
        assert dim % 2 == 0
        self.dim = dim
        half = dim // 2
        freqs = torch.exp(
            -math.log(max_period) * torch.arange(half, dtype=torch.float32) / half
        )
        self.register_buffer('freqs', freqs)

    def forward(self, t):
        t = t.float().unsqueeze(-1)
        args = t * self.freqs.unsqueeze(0)
        return torch.cat([torch.sin(args), torch.cos(args)], dim=-1)


class DistanceGaussianEmbedding(nn.Module):
    def __init__(self, start=0.0, stop=20.0, num_gaussians=16):
        super().__init__()
        assert num_gaussians >= 2
        offset = torch.linspace(start, stop, num_gaussians)
        coeff = -0.5 / (offset[1] - offset[0]).item() ** 2
        self.register_buffer('offset', offset)
        self.coeff = coeff

    def forward(self, dist):
        dist = dist.view(-1, 1) - self.offset.view(1, -1)
        return torch.exp(self.coeff * torch.pow(dist, 2))


class LigandConditionedActiveSet(nn.Module):
    """Per-residue gate logits for active-set pocket adaptation.

    Inputs are gathered in `protein_translations_batch` row order so the
    output aligns 1:1 with the residue tensors used by the score model.
    """

    def __init__(
        self,
        residue_h_dim,
        num_dist_gaussians=16,
        time_emb_dim=16,
        hidden_dim=128,
        dist_max=20.0,
    ):
        super().__init__()
        self.dist_smear = DistanceGaussianEmbedding(0.0, dist_max, num_dist_gaussians)
        self.time_emb = SinusoidalTimeEmbedding(time_emb_dim)
        in_dim = residue_h_dim + num_dist_gaussians + 1 + 1 + time_emb_dim
        self.mlp = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(
        self,
        residue_h,
        current_min_dist,
        pred_tr_norm,
        apo_drift_norm,
        time_step_per_residue,
    ):
        dist_feat = self.dist_smear(current_min_dist)
        time_feat = self.time_emb(time_step_per_residue)
        x = torch.cat(
            [
                residue_h,
                dist_feat,
                pred_tr_norm.unsqueeze(-1),
                apo_drift_norm.unsqueeze(-1),
                time_feat,
            ],
            dim=-1,
        )
        return self.mlp(x).squeeze(-1)
