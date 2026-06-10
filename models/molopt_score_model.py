import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_scatter import scatter_sum, scatter_mean
from tqdm.auto import tqdm
import kornia as K

from Bio.PDB import PDBParser, Superimposer
from Bio.SVDSuperimposer import SVDSuperimposer

from models.common import compose_context, ShiftedSoftplus
from models.egnn import EGNN
from models.uni_transformer import UniTransformerO2TwoUpdateGeneral
from models.attn import CrossAttention, RetAugmentationLinearAttention
from models.ligand_residue_cross_attn import LigandResidueCrossAttnGate
from utils.data import apply_transforms_tensor_batch

import time


def get_refine_net(refine_net_type, config):
    if refine_net_type == 'uni_o2':
        refine_net = UniTransformerO2TwoUpdateGeneral(
            num_blocks=config.num_blocks,
            num_layers=config.num_layers,
            hidden_dim=config.hidden_dim,
            n_heads=config.n_heads,
            k=config.knn,
            edge_feat_dim=config.edge_feat_dim,
            num_r_gaussian=config.num_r_gaussian,
            num_node_types=config.num_node_types,
            act_fn=config.act_fn,
            norm=config.norm,
            cutoff_mode=config.cutoff_mode,
            ew_net_type=config.ew_net_type,
            num_x2h=config.num_x2h,
            num_h2x=config.num_h2x,
            r_max=config.r_max,
            x2h_out_fc=config.x2h_out_fc,
            sync_twoup=config.sync_twoup
        )
    elif refine_net_type == 'egnn':
        refine_net = EGNN(
            num_layers=config.num_layers,
            hidden_dim=config.hidden_dim,
            edge_feat_dim=config.edge_feat_dim,
            num_r_gaussian=1,
            k=config.knn,
            cutoff_mode=config.cutoff_mode
        )
    else:
        raise ValueError(refine_net_type)
    return refine_net


def get_beta_schedule(beta_schedule, *, beta_start, beta_end, num_diffusion_timesteps):
    def sigmoid(x):
        return 1 / (np.exp(-x) + 1)

    if beta_schedule == "quad":
        betas = (
                np.linspace(
                    beta_start ** 0.5,
                    beta_end ** 0.5,
                    num_diffusion_timesteps,
                    dtype=np.float64,
                )
                ** 2
        )
    elif beta_schedule == "linear":
        betas = np.linspace(
            beta_start, beta_end, num_diffusion_timesteps, dtype=np.float64
        )
    elif beta_schedule == "const":
        betas = beta_end * np.ones(num_diffusion_timesteps, dtype=np.float64)
    elif beta_schedule == "jsd":  # 1/T, 1/(T-1), 1/(T-2), ..., 1
        betas = 1.0 / np.linspace(
            num_diffusion_timesteps, 1, num_diffusion_timesteps, dtype=np.float64
        )
    elif beta_schedule == "sigmoid":
        betas = np.linspace(-6, 6, num_diffusion_timesteps)
        betas = sigmoid(betas) * (beta_end - beta_start) + beta_start
    else:
        raise NotImplementedError(beta_schedule)
    assert betas.shape == (num_diffusion_timesteps,)
    return betas


def cosine_beta_schedule(timesteps, s=0.008):
    """
    cosine schedule
    as proposed in https://openreview.net/forum?id=-NEXDKk8gZ
    """
    steps = timesteps + 1
    x = np.linspace(0, steps, steps)
    alphas_cumprod = np.cos(((x / steps) + s) / (1 + s) * np.pi * 0.5) ** 2
    alphas_cumprod = alphas_cumprod / alphas_cumprod[0]
    alphas = (alphas_cumprod[1:] / alphas_cumprod[:-1])

    alphas = np.clip(alphas, a_min=0.001, a_max=1.)

    # Use sqrt of this, so the alpha in our paper is the alpha_sqrt from the
    # Gaussian diffusion in Ho et al.
    alphas = np.sqrt(alphas)
    return alphas


def get_distance(pos, edge_index):
    return (pos[edge_index[0]] - pos[edge_index[1]]).norm(dim=-1)


def calculate_tm_score(predicted_pos, reference_pos):
    device = predicted_pos.device
    # Convert to numpy arrays
    predicted_pos = predicted_pos.cpu().numpy()
    reference_pos = reference_pos.cpu().numpy()

    # Initialize the superimposer
    sup = SVDSuperimposer()
    sup.set(reference_pos, predicted_pos)
    sup.run()

    # Get the RMSD and rotation/translation matrices
    rmsd = sup.get_rms()
    rot, tran = sup.get_rotran()

    # Apply the rotation and translation to the predicted positions
    predicted_pos_aligned = np.dot(predicted_pos, rot) + tran

    L = predicted_pos_aligned.shape[0]  # Number of atoms
    d_0 = 1.24 * (L ** (1/3)) - 1.8  # TM-score normalization factor

    # Compute pairwise Euclidean distances
    distances = np.linalg.norm(predicted_pos_aligned - reference_pos, axis=1)

    # Compute TM-score
    tm_score = np.sum(1 / (1 + (distances / d_0) ** 2)) / L
    # move tm_score to torch cuda
    tm_score = torch.tensor(tm_score).to(device)

    return tm_score


def to_torch_const(x):
    x = torch.from_numpy(x).float()
    x = nn.Parameter(x, requires_grad=False)
    return x


def center_pos(protein_pos, protein_pos_holo, ligand_pos, batch_protein, batch_ligand, mode='protein'):
    if mode == 'none':
        offset = 0.
        pass
    elif mode == 'protein':
        offset = scatter_mean(protein_pos, batch_protein, dim=0)
        protein_pos = protein_pos - offset[batch_protein]
        protein_pos_holo = protein_pos_holo - offset[batch_protein]
        ligand_pos = ligand_pos - offset[batch_ligand]
    else:
        raise NotImplementedError
    return protein_pos, protein_pos_holo, ligand_pos, offset


def index_to_log_onehot(x, num_classes):
    assert x.max().item() < num_classes, f'Error: {x.max().item()} >= {num_classes}'
    x_onehot = F.one_hot(x, num_classes)
    log_x = torch.log(x_onehot.float().clamp(min=1e-30))
    return log_x


def log_onehot_to_index(log_x):
    return log_x.argmax(1)


def categorical_kl(log_prob1, log_prob2):
    kl = (log_prob1.exp() * (log_prob1 - log_prob2)).sum(dim=1)
    return kl


def log_categorical(log_x_start, log_prob):
    return (log_x_start.exp() * log_prob).sum(dim=1)


def normal_kl(mean1, logvar1, mean2, logvar2):
    """
    KL divergence between normal distributions parameterized by mean and log-variance.
    """
    kl = 0.5 * (-1.0 + logvar2 - logvar1 + torch.exp(logvar1 - logvar2) + (mean1 - mean2) ** 2 * torch.exp(-logvar2))
    return kl.sum(-1)


def log_normal(values, means, log_scales):
    var = torch.exp(log_scales * 2)
    log_prob = -((values - means) ** 2) / (2 * var) - log_scales - np.log(np.sqrt(2 * np.pi))
    return log_prob.sum(-1)


def log_sample_categorical(logits):
    uniform = torch.rand_like(logits)
    gumbel_noise = -torch.log(-torch.log(uniform + 1e-30) + 1e-30)
    sample_index = (gumbel_noise + logits).argmax(dim=-1)
    return sample_index


def log_1_min_a(a):
    return np.log(1 - np.exp(a) + 1e-40)


def log_add_exp(a, b):
    maximum = torch.max(a, b)
    return maximum + torch.log(torch.exp(a - maximum) + torch.exp(b - maximum))


# Time embedding
class SinusoidalPosEmb(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.dim = dim

    def forward(self, x):
        device = x.device
        half_dim = self.dim // 2
        emb = np.log(10000) / (half_dim - 1)
        emb = torch.exp(torch.arange(half_dim, device=device) * -emb)
        emb = x[:, None] * emb[None, :]
        emb = torch.cat((emb.sin(), emb.cos()), dim=-1)
        return emb


# Model
class ScorePosNet3D(nn.Module):

    def __init__(self, config, protein_atom_feature_dim, ligand_atom_feature_dim):
        super().__init__()
        self.config = config

        # variance schedule
        self.model_mean_type = config.model_mean_type  # ['noise', 'C0']
        self.loss_v_weight = config.loss_v_weight

        self.sample_time_method = config.sample_time_method  # ['importance', 'symmetric']

        if config.beta_schedule == 'cosine':
            alphas = cosine_beta_schedule(config.num_diffusion_timesteps, config.pos_beta_s) ** 2
            betas = 1. - alphas
        else:
            betas = get_beta_schedule(
                beta_schedule=config.beta_schedule,
                beta_start=config.beta_start,
                beta_end=config.beta_end,
                num_diffusion_timesteps=config.num_diffusion_timesteps,
            )
            alphas = 1. - betas
        if config.lambda_schedule == 'sigmoid':
            lambdas = np.linspace(-6, 6, config.num_diffusion_timesteps)
            lambdas = 1 - 1 / (1 + np.exp(-lambdas))
        alphas_cumprod = np.cumprod(alphas, axis=0)
        alphas_cumprod_prev = np.append(1., alphas_cumprod[:-1])

        self.betas = to_torch_const(betas)
        self.num_timesteps = self.betas.size(0)
        self.alphas_cumprod = to_torch_const(alphas_cumprod)
        self.alphas_cumprod_prev = to_torch_const(alphas_cumprod_prev)
        self.lambdas = to_torch_const(lambdas)

        self.sqrt_alphas_cumprod = to_torch_const(np.sqrt(alphas_cumprod))
        self.sqrt_one_minus_alphas_cumprod = to_torch_const(np.sqrt(1. - alphas_cumprod))
        self.sqrt_recip_alphas_cumprod = to_torch_const(np.sqrt(1. / alphas_cumprod))
        self.sqrt_recipm1_alphas_cumprod = to_torch_const(np.sqrt(1. / alphas_cumprod - 1))

        posterior_variance = betas * (1. - alphas_cumprod_prev) / (1. - alphas_cumprod)
        self.posterior_mean_c0_coef = to_torch_const(betas * np.sqrt(alphas_cumprod_prev) / (1. - alphas_cumprod))
        self.posterior_mean_ct_coef = to_torch_const(
            (1. - alphas_cumprod_prev) * np.sqrt(alphas) / (1. - alphas_cumprod))
        self.posterior_var = to_torch_const(posterior_variance)
        self.posterior_logvar = to_torch_const(np.log(np.append(self.posterior_var[1], self.posterior_var[1:])))

        if config.v_beta_schedule == 'cosine':
            alphas_v = cosine_beta_schedule(self.num_timesteps, config.v_beta_s)
        else:
            raise NotImplementedError
        log_alphas_v = np.log(alphas_v)
        log_alphas_cumprod_v = np.cumsum(log_alphas_v)
        self.log_alphas_v = to_torch_const(log_alphas_v)
        self.log_one_minus_alphas_v = to_torch_const(log_1_min_a(log_alphas_v))
        self.log_alphas_cumprod_v = to_torch_const(log_alphas_cumprod_v)
        self.log_one_minus_alphas_cumprod_v = to_torch_const(log_1_min_a(log_alphas_cumprod_v))

        self.custom_noise = to_torch_const(np.array([1.19691158e-04, 3.37258288e-01, 3.08534372e-01, 
                                        6.15134064e-02, 5.82621236e-02, 1.98585290e-01,
                                        1.93860432e-03, 9.45952575e-03, 8.06443701e-03,
                                        0.00000000e+00, 7.97025380e-03, 2.89809573e-03,
                                        5.39591284e-03]))

        self.register_buffer('Lt_history', torch.zeros(self.num_timesteps))
        self.register_buffer('Lt_count', torch.zeros(self.num_timesteps))

        self.hidden_dim = config.hidden_dim
        self.num_classes = ligand_atom_feature_dim
        if self.config.node_indicator:
            emb_dim = self.hidden_dim - 1
        else:
            emb_dim = self.hidden_dim

        self.protein_atom_emb = nn.Linear(protein_atom_feature_dim, emb_dim)

        self.center_pos_mode = config.center_pos_mode  # ['none', 'protein']

        self.time_emb_dim = config.time_emb_dim
        self.time_emb_mode = config.time_emb_mode  # ['simple', 'sin']
        if self.time_emb_dim > 0:
            if self.time_emb_mode == 'simple':
                self.ligand_atom_emb = nn.Linear(ligand_atom_feature_dim + 1, emb_dim)
            elif self.time_emb_mode == 'sin':
                self.time_emb = nn.Sequential(
                    SinusoidalPosEmb(self.time_emb_dim),
                    nn.Linear(self.time_emb_dim, self.time_emb_dim * 4),
                    nn.GELU(),
                    nn.Linear(self.time_emb_dim * 4, self.time_emb_dim)
                )
                self.ligand_atom_emb = nn.Linear(ligand_atom_feature_dim + self.time_emb_dim, emb_dim)
            else:
                raise NotImplementedError
        else:
            self.ligand_atom_emb = nn.Linear(ligand_atom_feature_dim, emb_dim)

        self.refine_net_type = config.model_type
        self.refine_net = get_refine_net(self.refine_net_type, config)
        self.v_inference = nn.Sequential(
            nn.Linear(self.hidden_dim, self.hidden_dim),
            ShiftedSoftplus(),
            nn.Linear(self.hidden_dim, ligand_atom_feature_dim),
        )
        self.res_inference = nn.Sequential(
            # nn.Linear(self.hidden_dim//2 + 4*3, self.hidden_dim),
            nn.Linear(self.hidden_dim + 3, self.hidden_dim),
            nn.LeakyReLU(),
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.LeakyReLU(),
            nn.Linear(self.hidden_dim, 3 + 4 + 5)
        )

        self.cond_dim = config.cond_dim
        self.topk_prompt = config.topk_prompt
        self.emb_mlp_aug = nn.Linear(emb_dim + self.cond_dim, self.cond_dim)
        self.emb_mlp = nn.Linear(emb_dim + self.cond_dim * 2, emb_dim)
        self.prompt_protein_mlp = nn.Linear(self.cond_dim * (self.topk_prompt + 1), self.cond_dim)
        self.prompt_ligand_attn = RetAugmentationLinearAttention(in_dim=self.cond_dim, d=self.cond_dim, context_dim=self.cond_dim)

        self.protein_update_schedule = getattr(config, 'protein_update_schedule', 'static5')
        self.protein_update_interval = int(getattr(config, 'protein_update_interval', 50))
        self.protein_update_min_t = int(getattr(config, 'protein_update_min_t', 10))
        self.protein_update_residual_threshold = float(
            getattr(config, 'protein_update_residual_threshold', 0.5)
        )
        self.pocket_router_mode = getattr(config, 'pocket_router_mode', 'none')
        self.pocket_router_topk = int(getattr(config, 'pocket_router_topk', 0))
        self.pocket_router_shell_radius = float(getattr(config, 'pocket_router_shell_radius', 0.0))
        self.pocket_router_shell_weight = float(getattr(config, 'pocket_router_shell_weight', 0.0))
        self.pocket_router_background_weight = float(getattr(config, 'pocket_router_background_weight', 0.0))
        self.protein_update_steps = self._build_protein_update_steps()

        if self.pocket_router_mode == 'cross_attn_gate':
            gate_cfg = getattr(config, 'cross_attn_gate', None)
            inner_dim = int(getattr(gate_cfg, 'inner_dim', 64)) if gate_cfg is not None else 64
            self.cross_attn_gate = LigandResidueCrossAttnGate(
                ligand_dim=self.hidden_dim,
                residue_dim=self.hidden_dim + 3,
                inner_dim=inner_dim,
            )
            # Induced-fit supervision: the gate is trained to predict which
            # residues are in the ligand-induced active set, supervised by the
            # ground-truth apo->holo residue displacement (a continuous, soft,
            # per-residue target). This gives the gate a correct gradient even
            # with a frozen backbone (the protein L1 loss alone would otherwise
            # collapse the gate into a scalar shrinkage of the frozen head).
            self.gate_loss_weight = float(getattr(gate_cfg, 'gate_loss_weight', 1.0)) if gate_cfg is not None else 1.0
            self.gate_target_center_A = float(getattr(gate_cfg, 'gate_target_center_A', 0.75)) if gate_cfg is not None else 0.75
            self.gate_target_sharpness = float(getattr(gate_cfg, 'gate_target_sharpness', 4.0)) if gate_cfg is not None else 4.0
            self.gate_rot_radius_A = float(getattr(gate_cfg, 'gate_rot_radius_A', 3.0)) if gate_cfg is not None else 3.0
            self.gate_chi_radius_A = float(getattr(gate_cfg, 'gate_chi_radius_A', 1.5)) if gate_cfg is not None else 1.5
        else:
            self.cross_attn_gate = None

    def _build_protein_update_steps(self):
        _num = 5
        _T = self.num_timesteps - 1
        if self.protein_update_schedule == 'static5':
            return {int(_T * (1 - i / _num)) for i in range(1, _num)} | {10}
        if self.protein_update_schedule == 'none':
            return set()
        if self.protein_update_schedule == 'uniform10':
            return {int(_T * (1 - i / 10)) for i in range(1, 10)} | {10}
        if self.protein_update_schedule in {'late_dense', 'residual_adaptive'}:
            interval = max(1, self.protein_update_interval)
            return set(range(self.protein_update_min_t, (_T // 2) + 1, interval)) | {10}
        raise ValueError(f'Unknown protein_update_schedule: {self.protein_update_schedule}')

    def _should_update_protein(self, t, pred_ligand_pos=None, ligand_pos=None):
        step = int(t[0])
        if self.protein_update_schedule != 'residual_adaptive':
            return step in self.protein_update_steps
        if step not in self.protein_update_steps:
            return False
        if pred_ligand_pos is None or ligand_pos is None:
            return False
        residual = torch.sqrt(((pred_ligand_pos - ligand_pos) ** 2).sum(-1)).mean()
        return residual.item() >= self.protein_update_residual_threshold

    def _residue_centers(self, pos, atom_to_residue, num_residues):
        return scatter_mean(pos, atom_to_residue, dim=0, dim_size=num_residues)

    def _residue_motion_scores(self, pos_a, pos_b, atom_to_residue, num_residues):
        atom_motion = torch.sqrt(((pos_b - pos_a) ** 2).sum(-1))
        return scatter_mean(atom_motion, atom_to_residue, dim=0, dim_size=num_residues)

    def _build_pocket_router_weights(
            self, data, protein_pos, protein_pos_apo, protein_pos_holo,
            ligand_pos, batch_protein, batch_ligand, pred_residue_tr=None,
            residue_h=None, ligand_h=None, time_step=None, return_logits=False):
        prot_update_batch = data.protein_translations_batch.view(-1)
        num_residues = prot_update_batch.numel()
        if self.pocket_router_mode == 'cross_attn_gate':
            assert self.cross_attn_gate is not None
            assert residue_h is not None and ligand_h is not None, (
                'cross_attn_gate requires residue_h and ligand_h')
            logits = self.cross_attn_gate(
                ligand_h=ligand_h,
                residue_h=residue_h,
                batch_ligand=batch_ligand,
                prot_update_batch=prot_update_batch,
            )
            if return_logits:
                return logits
            return torch.sigmoid(logits)
        if self.pocket_router_mode == 'none' or self.pocket_router_topk <= 0:
            return torch.ones(num_residues, dtype=protein_pos.dtype, device=protein_pos.device)

        weights = torch.full(
            (num_residues,),
            self.pocket_router_background_weight,
            dtype=protein_pos.dtype,
            device=protein_pos.device,
        )
        num_graphs = int(batch_protein.max().item()) + 1
        for graph_idx in range(num_graphs):
            residue_global = (prot_update_batch == graph_idx).nonzero(as_tuple=True)[0]
            if residue_global.numel() == 0:
                continue
            atom_mask = batch_protein == graph_idx
            ligand_mask = batch_ligand == graph_idx
            if ligand_mask.sum() == 0:
                continue

            atom_to_residue = data.protein_atom_to_aa_group[atom_mask].long()
            num_local_residues = residue_global.numel()
            current_centers = self._residue_centers(
                protein_pos[atom_mask], atom_to_residue, num_local_residues)
            apo_centers = self._residue_centers(
                protein_pos_apo[atom_mask], atom_to_residue, num_local_residues)
            holo_centers = self._residue_centers(
                protein_pos_holo[atom_mask], atom_to_residue, num_local_residues)
            ligand_pos_i = ligand_pos[ligand_mask]

            current_dist = torch.cdist(current_centers, ligand_pos_i).min(dim=1).values
            apo_dist = torch.cdist(apo_centers, ligand_pos_i).min(dim=1).values
            holo_dist = torch.cdist(holo_centers, ligand_pos_i).min(dim=1).values

            if self.pocket_router_mode == 'distance':
                scores = -current_dist
            elif self.pocket_router_mode == 'random':
                scores = torch.rand(num_local_residues, device=protein_pos.device)
            elif self.pocket_router_mode == 'motion_oracle':
                scores = self._residue_motion_scores(
                    protein_pos_apo[atom_mask], protein_pos_holo[atom_mask],
                    atom_to_residue, num_local_residues)
            elif self.pocket_router_mode == 'contact_oracle':
                scores = -holo_dist
            elif self.pocket_router_mode == 'contact_change_oracle':
                scores = torch.clamp(apo_dist - holo_dist, min=0.0)
            elif self.pocket_router_mode == 'predicted_motion':
                if pred_residue_tr is None:
                    scores = -current_dist
                else:
                    scores = pred_residue_tr[residue_global].norm(dim=-1)
            else:
                raise ValueError(f'Unknown pocket_router_mode: {self.pocket_router_mode}')

            k = min(max(1, self.pocket_router_topk), num_local_residues)
            selected_local = torch.topk(scores, k=k, largest=True).indices
            weights[residue_global[selected_local]] = 1.0

            if self.pocket_router_shell_radius > 0 and self.pocket_router_shell_weight > 0:
                selected_centers = current_centers[selected_local]
                shell_dist = torch.cdist(current_centers, selected_centers).min(dim=1).values
                shell_local = shell_dist <= self.pocket_router_shell_radius
                shell_global = residue_global[shell_local]
                shell_weight = torch.full_like(weights[shell_global], self.pocket_router_shell_weight)
                weights[shell_global] = torch.maximum(weights[shell_global], shell_weight)
                weights[residue_global[selected_local]] = 1.0
        return weights.clamp(0.0, 1.0)

    def forward(self, protein_pos, protein_v, batch_protein, init_ligand_pos, init_ligand_v, batch_ligand, protein_atom_to_aa_group,
                time_step=None, return_all=False, fix_x=False, hbap_protein=None, hbap_ligand=None
                ):

        batch_size = batch_protein.max().item() + 1
        init_ligand_v = F.one_hot(init_ligand_v, self.num_classes).float()
        if self.time_emb_dim > 0:
            if self.time_emb_mode == 'simple':
                input_ligand_feat = torch.cat([
                    init_ligand_v,
                    (time_step / self.num_timesteps)[batch_ligand].unsqueeze(-1)
                ], -1)
            elif self.time_emb_mode == 'sin':
                time_feat = self.time_emb(time_step)
                input_ligand_feat = torch.cat([init_ligand_v, time_feat], -1)
            else:
                raise NotImplementedError
        else:
            input_ligand_feat = init_ligand_v

        h_protein = self.protein_atom_emb(protein_v)
        init_ligand_h = self.ligand_atom_emb(input_ligand_feat)

        if hbap_protein is None:
            hbap_protein = torch.zeros([h_protein.shape[0], self.cond_dim]).to(h_protein.device)
        if hbap_ligand is None:
            hbap_ligand = torch.zeros([init_ligand_h.shape[0], self.cond_dim]).to(init_ligand_h.device)

        hbap_protein_aug = self.emb_mlp_aug(torch.cat([h_protein, hbap_protein.detach()], dim=1))
        hbap_ligand_aug = self.emb_mlp_aug(torch.cat([init_ligand_h, hbap_ligand.detach()], dim=1))

        protein_prompt_list = [hbap_protein_aug]

        hbap_protein_aug = self.prompt_protein_mlp(torch.cat(protein_prompt_list, dim=1))

        max_lig_len = 150
        hbap_ligand_aug_batch = torch.zeros([batch_size, max_lig_len, hbap_ligand_aug.shape[1]]).to(hbap_ligand_aug.device)

        valid_num_atom_list = []
        for bi in range(batch_size):
            hbap_ligand_aug_batch_i = hbap_ligand_aug[batch_ligand == bi]
            num_atom = hbap_ligand_aug_batch_i.shape[0]
            assert num_atom <= max_lig_len

            valid_num_atom_list.append(num_atom)
            hbap_ligand_aug_batch[bi, :num_atom] = hbap_ligand_aug_batch_i

        prompt_hbap_ligand_batch_all_list = []

        if prompt_hbap_ligand_batch_all_list != []:
            prompt_hbap_ligand_batch_all = torch.cat(prompt_hbap_ligand_batch_all_list, dim=1)
            hbap_ligand_aug_batch = self.prompt_ligand_attn(h=hbap_ligand_aug_batch, h_retrieved=prompt_hbap_ligand_batch_all)
        else:
            prompt_hbap_ligand_batch_all = None
            hbap_ligand_aug_batch = self.prompt_ligand_attn(h=hbap_ligand_aug_batch, h_retrieved=hbap_ligand_aug_batch)

        hbap_ligand_aug_list = []
        for bi in range(batch_size):
            hbap_ligand_aug_i = hbap_ligand_aug_batch[bi][:valid_num_atom_list[bi]]
            hbap_ligand_aug_list.append(hbap_ligand_aug_i)
        hbap_ligand_aug = torch.cat(hbap_ligand_aug_list, dim=0).to(hbap_ligand_aug.device)

        h_protein = self.emb_mlp(torch.cat([h_protein, hbap_protein, hbap_protein_aug], dim=1))
        init_ligand_h = self.emb_mlp(torch.cat([init_ligand_h, hbap_ligand, hbap_ligand_aug], dim=1))

        if self.config.node_indicator:
            h_protein = torch.cat([h_protein, torch.zeros(len(h_protein), 1).to(h_protein)], -1)
            init_ligand_h = torch.cat([init_ligand_h, torch.ones(len(init_ligand_h), 1).to(h_protein)], -1)

        outputs = self.refine_net(
            h_protein, init_ligand_h,
            protein_pos, init_ligand_pos,
            batch_protein, batch_ligand,
            protein_atom_to_aa_group,
        )

        final_ligand_pos = outputs['ligand_pos']
        final_ligand_h = outputs['ligand_h']
        residue_h = outputs['residue_h']
        final_ligand_v = self.v_inference(final_ligand_h)
        final_res_out = self.res_inference(residue_h)
        final_res_tr = final_res_out[:, :3]  # translation
        final_res_rot = final_res_out[:, 3:7]  # rotation
        final_res_chi = final_res_out[:, 7:]  # chi angles

        preds = {
            'pred_ligand_pos': final_ligand_pos,
            'pred_ligand_v': final_ligand_v,
            'pred_residue_tr': final_res_tr,
            'pred_residue_rot': final_res_rot,
            'pred_residue_chi': final_res_chi,
            'final_ligand_h': final_ligand_h,
            'residue_h': residue_h,
        }
        return preds

    def q_v_pred_one_timestep(self, log_vt_1, t, batch):
        log_alpha_t = extract(self.log_alphas_v, t, batch)
        log_1_min_alpha_t = extract(self.log_one_minus_alphas_v, t, batch)

        log_probs = log_add_exp(
            log_vt_1 + log_alpha_t,
            log_1_min_alpha_t - np.log(self.num_classes)
        )
        return log_probs

    def q_v_pred(self, log_v0, t, batch):
        log_cumprod_alpha_t = extract(self.log_alphas_cumprod_v, t, batch)
        log_1_min_cumprod_alpha = extract(self.log_one_minus_alphas_cumprod_v, t, batch)

        log_probs = log_add_exp(
            log_v0 + log_cumprod_alpha_t,
            log_1_min_cumprod_alpha - np.log(self.num_classes)
        )
        return log_probs

    def q_v_sample(self, log_v0, t, batch):
        log_qvt_v0 = self.q_v_pred(log_v0, t, batch)
        sample_index = log_sample_categorical(log_qvt_v0)
        log_sample = index_to_log_onehot(sample_index, self.num_classes)
        return sample_index, log_sample

    def q_v_posterior(self, log_v0, log_vt, t, batch):
        t_minus_1 = t - 1
        t_minus_1 = torch.where(t_minus_1 < 0, torch.zeros_like(t_minus_1), t_minus_1)
        log_qvt1_v0 = self.q_v_pred(log_v0, t_minus_1, batch)
        unnormed_logprobs = log_qvt1_v0 + self.q_v_pred_one_timestep(log_vt, t, batch)
        log_vt1_given_vt_v0 = unnormed_logprobs - torch.logsumexp(unnormed_logprobs, dim=-1, keepdim=True)
        return log_vt1_given_vt_v0

    def kl_v_prior(self, log_x_start, batch):
        num_graphs = batch.max().item() + 1
        log_qxT_prob = self.q_v_pred(log_x_start, t=[self.num_timesteps - 1] * num_graphs, batch=batch)
        log_half_prob = -torch.log(self.num_classes * torch.ones_like(log_qxT_prob))
        kl_prior = categorical_kl(log_qxT_prob, log_half_prob)
        kl_prior = scatter_mean(kl_prior, batch, dim=0)
        return kl_prior

    def _predict_x0_from_eps(self, xt, eps, t, batch):
        pos0_from_e = extract(self.sqrt_recip_alphas_cumprod, t, batch) * xt - \
                      extract(self.sqrt_recipm1_alphas_cumprod, t, batch) * eps
        return pos0_from_e

    def q_pos_posterior(self, x0, xt, t, batch):
        pos_model_mean = extract(self.posterior_mean_c0_coef, t, batch) * x0 + \
                         extract(self.posterior_mean_ct_coef, t, batch) * xt
        return pos_model_mean

    def q_pos_linear(self, x0, xt, t, batch):
        lambda_t = extract(self.lambdas, t, batch)
        lambda_t_1 = extract(self.lambdas, t - 1, batch)
        pos_model_mean = (1 - lambda_t_1) / (1 - lambda_t) * (xt - lambda_t * x0) + \
                         (lambda_t_1 * x0)

        return pos_model_mean

    def kl_pos_prior(self, pos0, batch):
        num_graphs = batch.max().item() + 1
        a_pos = extract(self.alphas_cumprod, [self.num_timesteps - 1] * num_graphs, batch)  # (num_ligand_atoms, 1)
        pos_model_mean = a_pos.sqrt() * pos0
        pos_log_variance = torch.log((1.0 - a_pos).sqrt())
        kl_prior = normal_kl(torch.zeros_like(pos_model_mean), torch.zeros_like(pos_log_variance),
                             pos_model_mean, pos_log_variance)
        kl_prior = scatter_mean(kl_prior, batch, dim=0)
        return kl_prior

    def sample_time(self, num_graphs, device, method):
        if method == 'importance':
            if not (self.Lt_count > 10).all():
                return self.sample_time(num_graphs, device, method='symmetric')

            Lt_sqrt = torch.sqrt(self.Lt_history + 1e-10) + 0.0001
            Lt_sqrt[0] = Lt_sqrt[1]
            pt_all = Lt_sqrt / Lt_sqrt.sum()

            time_step = torch.multinomial(pt_all, num_samples=num_graphs, replacement=True)
            pt = pt_all.gather(dim=0, index=time_step)
            return time_step, pt

        elif method == 'symmetric':
            time_step = torch.randint(
                0, self.num_timesteps, size=(num_graphs // 2 + 1,), device=device)
            time_step = torch.cat(
                [time_step, self.num_timesteps - time_step - 1], dim=0)[:num_graphs]
            pt = torch.ones_like(time_step).float() / self.num_timesteps
            return time_step, pt

        else:
            raise ValueError

    def compute_pos_Lt(self, pos_model_mean, x0, xt, t, batch):
        pos_log_variance = extract(self.posterior_logvar, t, batch)
        pos_true_mean = self.q_pos_posterior(x0=x0, xt=xt, t=t, batch=batch)
        kl_pos = normal_kl(pos_true_mean, pos_log_variance, pos_model_mean, pos_log_variance)
        kl_pos = kl_pos / np.log(2.)

        decoder_nll_pos = -log_normal(x0, means=pos_model_mean, log_scales=0.5 * pos_log_variance)
        assert kl_pos.shape == decoder_nll_pos.shape
        mask = (t == 0).float()[batch]
        loss_pos = scatter_mean(mask * decoder_nll_pos + (1. - mask) * kl_pos, batch, dim=0)
        return loss_pos

    def compute_v_Lt(self, log_v_model_prob, log_v0, log_v_true_prob, t, batch):
        kl_v = categorical_kl(log_v_true_prob, log_v_model_prob)  # [num_atoms, ]
        decoder_nll_v = -log_categorical(log_v0, log_v_model_prob)  # L0
        assert kl_v.shape == decoder_nll_v.shape
        mask = (t == 0).float()[batch]
        loss_v = scatter_mean(mask * decoder_nll_v + (1. - mask) * kl_v, batch, dim=0)
        return loss_v

    def _build_induced_fit_target(self, prot_tr, prot_rot, prot_chi_update, prot_chi_mask):
        """Continuous soft target for the active-set gate.

        The ground-truth apo->holo residue motion is mapped to an *effective
        displacement* in Angstrom (all three components converted to A), then to
        a soft per-residue label in [0,1]:

            eff_motion_i = ||tr_i||                          # backbone translation (A)
                         + rot_radius * angle(rot_i)         # frame rotation -> A
                         + chi_radius * ||wrap(chi_update_i)|| # side-chain torsion -> A
            target_i     = sigmoid(sharpness * (eff_motion_i - center_A))

        With the defaults (center 0.75 A, sharpness 4) a rigid residue (eff
        motion ~0) maps to ~0.05, a shell-like residue (~0.5 A) to ~0.27, and a
        core residue (>=1.5 A) to ~1 -- i.e. the core/shell/background gradient
        emerges from the real induced-fit field instead of hand-tuned thresholds.

        The chi delta is wrapped to (-pi, pi] before taking its magnitude:
        chi is stored on [0, 2pi) (utils/data.py), so a raw subtraction would
        turn e.g. 6.25 -> 0.03 rad into a spurious 6.22 rad and mislabel a rigid
        residue as core. (The diffusion chi loss is unaffected: it uses a cosine
        loss which is already periodic-invariant.)
        """
        tr_mag = prot_tr.norm(dim=-1)
        rot_w = prot_rot[:, 0].abs().clamp(max=1.0)
        rot_angle = 2.0 * torch.acos(rot_w)
        chi_delta = torch.atan2(torch.sin(prot_chi_update), torch.cos(prot_chi_update))
        chi_mag = (chi_delta * prot_chi_mask).norm(dim=-1)
        eff_motion = (
            tr_mag
            + self.gate_rot_radius_A * rot_angle
            + self.gate_chi_radius_A * chi_mag
        )
        target = torch.sigmoid(self.gate_target_sharpness * (eff_motion - self.gate_target_center_A))
        return target

    def slerp_identity_to_q(self, q, lambdas):
        """
        q:       [B, 4]  (quaternion, w,x,y,z)
        lambdas: [B, 1]  (1-t)
        return:    [B, 4]
        """
        q0 = torch.zeros_like(q)
        q0[:,0] = 1.0  # identity, w=1
        q1 = q / q.norm(dim=-1, keepdim=True).clamp(min=1e-8)  # normalize

        dot = (q0 * q1).sum(-1, keepdim=True)
        q1 = torch.where(dot < 0, -q1, q1)
        dot = dot.abs().clamp(-1, 1)

        theta = torch.acos(dot)
        sin_theta = torch.sin(theta)

        mask = sin_theta > 1e-6
        s0 = torch.where(mask, torch.sin(lambdas * theta) / sin_theta, lambdas)
        s1 = torch.where(mask, torch.sin((1 - lambdas) * theta) / sin_theta, 1 - lambdas)

        out = s0 * q0 + s1 * q1
        out = out / out.norm(dim=-1, keepdim=True).clamp(min=1e-8)

        return out

    def add_noise_to_quaternion(self, q: torch.Tensor, noise_scale: torch.Tensor) -> torch.Tensor:
        """
        Add small rotational noise to unit quaternions.
        Args:
            q: (..., 4) unit quaternions
            noise_scale: (..., 1) scalar controlling the noise level (e.g., beta.sqrt())
        Returns:
            (..., 4) noisy unit quaternions
        """
        # Generate small random axis-angle noise
        axis = torch.randn_like(q[..., 1:])  # (..., 3)
        axis = axis / axis.norm(dim=-1, keepdim=True).clamp(min=1e-6)
        angle = torch.randn(q.shape[0], device=q.device) * noise_scale.squeeze(-1)  # (...,)
        axis_angle = axis * angle.unsqueeze(-1)  # (..., 3)

        # Convert to quaternion
        delta_q = K.geometry.conversions.axis_angle_to_quaternion(axis_angle)

        # Compose original quaternion with the noise
        q_noisy = quaternion_product(q, delta_q)
        q_noisy = q_noisy / q_noisy.norm(dim=-1, keepdim=True).clamp(min=1e-8)  # re-normalize

        return q_noisy

    def get_diffusion_loss(
            self, net_cond, data,
            protein_pos_apo, protein_pos_holo, protein_v, batch_protein,
            ligand_pos, ligand_v, batch_ligand,
            time_step=None
    ):
        num_graphs = batch_protein.max().item() + 1
        protein_pos_apo, protein_pos_holo, ligand_pos, offset = center_pos(
            protein_pos_apo, protein_pos_holo, ligand_pos, batch_protein, batch_ligand, mode=self.center_pos_mode)

        if time_step is None:
            time_step, pt = self.sample_time(num_graphs, protein_pos_holo.device, self.sample_time_method)
        else:
            pt = torch.ones_like(time_step).float() / self.num_timesteps

        # Add noise to ligand position
        a = self.alphas_cumprod.index_select(0, time_step)  # (num_graphs, )

        a_pos = a[batch_ligand].unsqueeze(-1)  # (num_ligand_atoms, 1)
        pos_noise = torch.zeros_like(ligand_pos)
        pos_noise.normal_()
        ligand_pos_perturbed = a_pos.sqrt() * ligand_pos + (1.0 - a_pos).sqrt() * pos_noise  # pos_noise * std
        # Add noise to ligand atom types
        log_ligand_v0 = index_to_log_onehot(ligand_v, self.num_classes) # torch.log(one_hot(ligand_v))
        ligand_v_perturbed, log_ligand_vt = self.q_v_sample(log_ligand_v0, time_step, batch_ligand)

        # Add noise to pocket position
        prot_update_batch = data.protein_translations_batch
        beta_update = self.betas.index_select(0, time_step)  # (num_graphs, )
        beta_update = beta_update[prot_update_batch].unsqueeze(-1)  # (num_res, 1)
        lambdas = self.lambdas.index_select(0, time_step)
        lambdas_update = lambdas[prot_update_batch].unsqueeze(-1) # (num_res, 1)

        prot_tr = data.protein_translations
        prot_rot = data.protein_rotations
        prot_chi_apo = data.protein_chi_apo
        prot_chi_holo = data.protein_chi_holo
        prot_chi_update = (prot_chi_apo - prot_chi_holo)
        prot_chi_mask = data.protein_chi_mask
        # prot_rot_mat = K.geometry.conversions.axis_angle_to_rotation_matrix(prot_rot)
        prot_rot_mat = K.geometry.conversions.quaternion_to_rotation_matrix(prot_rot)
        prot_tr = torch.matmul(prot_rot_mat, offset[prot_update_batch].unsqueeze(-1)).squeeze(-1) - offset[prot_update_batch] + prot_tr
        prot_tr_t = (1 - lambdas_update) * prot_tr
        # prot_rot_t = (1 - lambdas_update) * prot_rot
        prot_rot_t = self.slerp_identity_to_q(prot_rot, lambdas_update)
        prot_chi_t = (1 - lambdas_update) * prot_chi_update
        prot_tr_t += torch.randn_like(prot_tr_t) * beta_update.sqrt() * 3
        # prot_rot_t += torch.randn_like(prot_rot_t) * beta_update.sqrt() * 2
        prot_rot_t = self.add_noise_to_quaternion(prot_rot_t, beta_update.sqrt() * 2)
        prot_chi_t += torch.randn_like(prot_chi_t) * beta_update.sqrt() * 1
        prot_chi_t = prot_chi_t * prot_chi_mask

        protein_pos_perturbed = apply_transforms_tensor_batch(
            protein_pos=protein_pos_holo,
            protein_atom_name=data.protein_atom_name,
            protein_atom_to_aa_name=data.protein_atom_to_aa_name,
            protein_atom_to_aa_group=data.protein_atom_to_aa_group,
            protein_element_batch=data.protein_element_batch,
            rotations=prot_rot_t,
            translations=prot_tr_t,
            chi_update=prot_chi_t,
            chi_mask=prot_chi_mask,
            protein_translations_batch=prot_update_batch,
        )

        ligand_hbap_flag = (time_step[batch_ligand] == self.num_timesteps)
        protein_hbap_flag = time_step[batch_protein] == self.num_timesteps
        ligand_hbap_mask = torch.ones_like(batch_ligand).float()
        ligand_hbap_mask[ligand_hbap_flag] = 0.
        ligand_hbap_mask = ligand_hbap_mask.unsqueeze(-1)
        protein_hbap_mask = torch.ones_like(batch_protein).float()
        protein_hbap_mask[protein_hbap_flag] = 0.
        protein_hbap_mask = protein_hbap_mask.unsqueeze(-1)

        hbap_ligand = None
        hbap_protein = None

        gt_protein_v = protein_v
        gt_protein_pos = protein_pos_holo
        gt_protein_a_h = torch.argmax(gt_protein_v[:, :6], dim=1)
        gt_protein_r_h = torch.argmax(gt_protein_v[:, 6:26], dim=1)

        if self.model_mean_type == 'noise':
            pass
        elif self.model_mean_type == 'C0':
            gt_ligand_v = ligand_v
            gt_ligand_pos = ligand_pos
            gt_lig_a_h = gt_ligand_v

            hbap_ligand, hbap_protein = net_cond.extract_features(gt_ligand_pos, gt_protein_pos, gt_lig_a_h, gt_protein_a_h, gt_protein_r_h, batch_ligand, batch_protein)
            hbap_ligand = hbap_ligand * ligand_hbap_mask
            hbap_protein = hbap_protein * protein_hbap_mask
        else:
            raise ValueError

        # 3. forward-pass NN, feed perturbed pos and v, output noise
        preds = self.forward(
            protein_pos=protein_pos_perturbed,
            protein_v=protein_v,
            batch_protein=batch_protein,

            init_ligand_pos=ligand_pos_perturbed,
            init_ligand_v=ligand_v_perturbed,
            batch_ligand=batch_ligand,
            protein_atom_to_aa_group=data.protein_atom_to_aa_group,
            time_step=time_step,

            hbap_protein=hbap_protein,
            hbap_ligand=hbap_ligand,
        )
        # Get ligand predictions

        pred_ligand_pos, pred_ligand_v = preds['pred_ligand_pos'], preds['pred_ligand_v']
        pred_ligand_pos_noise = pred_ligand_pos - ligand_pos_perturbed
        if self.model_mean_type == 'noise':
            pos0_from_e = self._predict_x0_from_eps(
                xt=ligand_pos_perturbed, eps=pred_ligand_pos_noise, t=time_step, batch=batch_ligand)
            pos_model_mean = self.q_pos_posterior(
                x0=pos0_from_e, xt=ligand_pos_perturbed, t=time_step, batch=batch_ligand)
        elif self.model_mean_type == 'C0':
            pos_model_mean = self.q_pos_posterior(
                x0=pred_ligand_pos, xt=ligand_pos_perturbed, t=time_step, batch=batch_ligand)
        else:
            raise ValueError
        # Get protein predictions
        pred_res_tr = preds['pred_residue_tr']
        pred_res_rot = preds['pred_residue_rot']
        pred_res_chi = preds['pred_residue_chi']

        # Ligand-conditioned cross-attention gate. The gate predicts a per-residue
        # active-set weight w_i in [0,1] from (ligand state, residue state) and is
        # supervised by the ground-truth apo->holo residue displacement (computed
        # below as loss_gate). The same w_i is multiplied into the residue updates
        # so train-time and sampling-time behaviour stay consistent.
        router_w_mean = None
        router_w_min = None
        router_w_max = None
        router_selected_per_graph = None
        loss_gate = pred_res_tr.new_zeros(())
        if self.pocket_router_mode == 'cross_attn_gate':
            gate_logits = self._build_pocket_router_weights(
                data=data,
                protein_pos=protein_pos_perturbed,
                protein_pos_apo=protein_pos_apo,
                protein_pos_holo=protein_pos_holo,
                ligand_pos=ligand_pos_perturbed,
                batch_protein=batch_protein,
                batch_ligand=batch_ligand,
                residue_h=preds['residue_h'],
                ligand_h=preds['final_ligand_h'],
                return_logits=True,
            )
            router_weights = torch.sigmoid(gate_logits)
            # Induced-fit soft target from ground-truth apo->holo residue motion.
            gate_target = self._build_induced_fit_target(
                prot_tr=prot_tr, prot_rot=prot_rot,
                prot_chi_update=prot_chi_update, prot_chi_mask=prot_chi_mask,
            ).detach()
            bce_per_res = F.binary_cross_entropy_with_logits(
                gate_logits, gate_target, reduction='none')
            bce_per_graph = scatter_mean(bce_per_res, data.protein_translations_batch.view(-1), dim=0)
            loss_gate = bce_per_graph.mean()
            with torch.no_grad():
                router_w_mean = router_weights.mean()
                router_w_min = router_weights.min()
                router_w_max = router_weights.max()
                sel = (router_weights > 0.5).float()
                sel_per_graph = scatter_sum(sel, data.protein_translations_batch.view(-1), dim=0)
                router_selected_per_graph = sel_per_graph.mean()
            w_f = router_weights.unsqueeze(-1)
            pred_res_tr = pred_res_tr * w_f
            pred_res_chi = pred_res_chi * w_f
            pred_res_rot = self.slerp_identity_to_q(pred_res_rot, 1.0 - w_f)

        if self.model_mean_type == 'C0':
            target_ligand, pred_ligand = ligand_pos, pred_ligand_pos
        elif self.model_mean_type == 'noise':
            target_ligand, pred_ligand = pos_noise, pred_ligand_pos_noise
        else:
            raise ValueError
        # Loss for ligand
        loss_ligand_pos = scatter_mean(((pred_ligand - target_ligand) ** 2).sum(-1), batch_ligand, dim=0)
        loss_ligand_pos = torch.mean(loss_ligand_pos)
        log_ligand_v_recon = F.log_softmax(pred_ligand_v, dim=-1)
        log_v_model_prob = self.q_v_posterior(log_ligand_v_recon, log_ligand_vt, time_step, batch_ligand)
        log_v_true_prob = self.q_v_posterior(log_ligand_v0, log_ligand_vt, time_step, batch_ligand)
        kl_v = self.compute_v_Lt(log_v_model_prob=log_v_model_prob, log_v0=log_ligand_v0,
                                 log_v_true_prob=log_v_true_prob, t=time_step, batch=batch_ligand)
        loss_v = torch.mean(kl_v)
        # Loss for protein
        inverse_rot = torch.cat([prot_rot_t[:, :1], -prot_rot_t[:, 1:]], dim=-1)  # Conjugate: [w, -x, -y, -z]
        prot_rot_mat_t = K.geometry.conversions.quaternion_to_rotation_matrix(prot_rot_t)
        prot_rot_mat_inv_t = prot_rot_mat_t.transpose(-2, -1)
        inverse_tr = -torch.matmul(prot_rot_mat_inv_t, prot_tr_t.unsqueeze(-1)).squeeze(-1)
        inverse_chi = -prot_chi_t
        loss_prot_tr = scatter_mean(nn.L1Loss(reduction='none')(pred_res_tr, inverse_tr).sum(-1), data.protein_translations_batch, dim=0)
        loss_prot_tr = torch.mean(loss_prot_tr)
        loss_prot_rot, pred_res_rot = self.calculate_quat_loss(pred_res_rot, inverse_rot, data.protein_translations_batch)
        chi_loss = 1 - (pred_res_chi - inverse_chi).cos()  # Cosine loss
        chi_loss = (chi_loss * data.protein_chi_mask).sum(dim=-1) / (data.protein_chi_mask.sum(dim=-1) + 1e-12)  # Apply mask and normalize
        loss_prot_chi = scatter_mean(chi_loss, data.protein_translations_batch, dim=0)
        loss_prot_chi = torch.mean(loss_prot_chi)

        loss = loss_ligand_pos + loss_v * self.loss_v_weight + loss_prot_tr + loss_prot_rot + 5*loss_prot_chi
        if self.pocket_router_mode == 'cross_attn_gate':
            loss = loss + self.gate_loss_weight * loss_gate

        return {
            'loss_ligang_pos': loss_ligand_pos,
            'loss_protein_tr': loss_prot_tr,
            'loss_protein_rot': loss_prot_rot,
            'loss_protein_chi': loss_prot_chi,
            'loss_v': loss_v,
            'loss': loss,
            'loss_gate': loss_gate,
            'router_w_mean': router_w_mean,
            'router_w_min': router_w_min,
            'router_w_max': router_w_max,
            'router_selected_per_graph': router_selected_per_graph,
            'x0': ligand_pos,
            'p0': protein_pos_holo,
            'pred_ligand_pos': pred_ligand_pos,
            'perturbed_protein_pos': protein_pos_perturbed,
            "protein_pos_apo": protein_pos_apo,
            'pred_res_tr': pred_res_tr,
            'pred_res_rot': pred_res_rot,
            'pred_res_chi': pred_res_chi,
            'pred_ligand_v': pred_ligand_v,
            'pred_pos_noise': pred_ligand_pos_noise,
            'ligand_v_recon': F.softmax(pred_ligand_v, dim=-1)
        }

    def calculate_quat_loss(self, pred_rot, target_rot, batch):
        """
        Calculate quaternion loss between predicted and target rotations.
        Args:
            pred_rot: Predicted rotations (B, 4)
            target_rot: Target rotations (B, 4)
        """
        scale_loss = abs(1.0 - pred_rot.norm(dim=-1, keepdim=True))  # Scale loss for normalization
        scale_loss = scatter_mean(scale_loss, batch, dim=0)
        scale_loss = torch.mean(scale_loss)
        # Calculate quaternion dot product
        pred_rot = pred_rot / (pred_rot.norm(dim=-1, keepdim=True).clamp_min(1e-8))  # Normalize predicted rotations
        quat_loss = ((pred_rot - target_rot) ** 2).sum(-1)  # L2 loss for quaternion
        quat_loss = scatter_mean(quat_loss, batch, dim=0)
        quat_loss = torch.mean(quat_loss) * 10

        return scale_loss + quat_loss, pred_rot

    @torch.no_grad()
    def sample_diffusion(self, data,
                         protein_pos_apo, protein_pos_holo, protein_v, batch_protein,
                         init_ligand_pos, init_ligand_v, batch_ligand,
                         num_steps=None, center_pos_mode=None, pos_only=False, net_cond=None, cond_dim=128):

        if num_steps is None:
            num_steps = self.num_timesteps
        num_graphs = batch_protein.max().item() + 1

        protein_pos_apo, protein_pos_holo, init_ligand_pos, offset = center_pos(
            protein_pos_apo, protein_pos_holo, init_ligand_pos, batch_protein, batch_ligand, mode=center_pos_mode)

        protein_pos = protein_pos_apo

        protein_pos_traj, ligand_pos_traj, v_traj = [], [], []
        v0_pred_traj, vt_pred_traj = [], []
        router_selected_counts = []
        ligand_pos, ligand_v = init_ligand_pos, init_ligand_v

        gt_protein_v = protein_v.detach()
        gt_protein_a_h = torch.argmax(gt_protein_v[:, :6], dim=1)
        gt_protein_r_h = torch.argmax(gt_protein_v[:, 6:26], dim=1)

        hbap_protein = None
        hbap_ligand = None

        time_seq = list(reversed(range(self.num_timesteps - num_steps, self.num_timesteps)))
        for i in tqdm(time_seq, desc='sampling', total=len(time_seq)):

            t = torch.full(size=(num_graphs,), fill_value=i, dtype=torch.long, device=protein_pos_apo.device)
            preds = self.forward(
                protein_pos=protein_pos,
                protein_v=protein_v,
                batch_protein=batch_protein,
                init_ligand_pos=ligand_pos,
                init_ligand_v=ligand_v,
                batch_ligand=batch_ligand,
                protein_atom_to_aa_group=data.protein_atom_to_aa_group,
                time_step=t,

                hbap_protein=hbap_protein,
                hbap_ligand=hbap_ligand,
            )
            if self.model_mean_type == 'noise':
                pred_pos_noise = preds['pred_ligand_pos'] - ligand_pos
                pos0_from_e = self._predict_x0_from_eps(xt=ligand_pos, eps=pred_pos_noise, t=t, batch=batch_ligand)
                v0_from_e = preds['pred_ligand_v']
            elif self.model_mean_type == 'C0':
                ligand_pos0_from_e = preds['pred_ligand_pos']
                v0_from_e = preds['pred_ligand_v']
            else:
                raise ValueError

            pred_ligand_pos = ligand_pos0_from_e.detach()
            pred_lig_a_h = torch.argmax(v0_from_e.detach(), dim=1)
            pred_residue_tr = preds['pred_residue_tr'].detach()
            pred_residue_rot = preds['pred_residue_rot'].detach()
            pred_residue_rot = pred_residue_rot / pred_residue_rot.norm(dim=-1, keepdim=True).clamp_min(1e-8)  # Normalize quaternion
            pred_residue_chi = preds['pred_residue_chi'].detach()

            # Update ligand pos and atom types
            ligand_pos_prev = ligand_pos
            ligand_pos_model_mean = self.q_pos_posterior(x0=ligand_pos0_from_e, xt=ligand_pos, t=t, batch=batch_ligand)
            ligand_pos_log_variance = extract(self.posterior_logvar, t, batch_ligand)
            nonzero_mask = (1 - (t == 0).float())[batch_ligand].unsqueeze(-1)
            ligand_pos_next = ligand_pos_model_mean + nonzero_mask * (0.5 * ligand_pos_log_variance).exp() * torch.randn_like(ligand_pos)
            ligand_pos = ligand_pos_next
            if self._should_update_protein(t, pred_ligand_pos=pred_ligand_pos, ligand_pos=ligand_pos_prev):
                # Update protein pos
                prot_update_batch = data.protein_translations_batch
                beta_update = self.betas.index_select(0, t)
                beta_update = beta_update[prot_update_batch].unsqueeze(-1)
                residue_tr_t = (1 / t[0] + 1) * pred_residue_tr + torch.randn_like(pred_residue_tr) * beta_update.sqrt() * 3
                residue_chi_t = (1 / t[0] + 1) * pred_residue_chi + torch.randn_like(pred_residue_chi) * beta_update.sqrt() * 1
                residue_rot_t = self.slerp_identity_to_q(pred_residue_rot, (1 / t + 1)[prot_update_batch].unsqueeze(-1))
                residue_rot_t = self.add_noise_to_quaternion(residue_rot_t, beta_update.sqrt() * 2)
                router_weights = self._build_pocket_router_weights(
                    data=data,
                    protein_pos=protein_pos,
                    protein_pos_apo=protein_pos_apo,
                    protein_pos_holo=protein_pos_holo,
                    ligand_pos=ligand_pos_prev,
                    batch_protein=batch_protein,
                    batch_ligand=batch_ligand,
                    pred_residue_tr=pred_residue_tr,
                    residue_h=preds['residue_h'].detach(),
                    ligand_h=preds['final_ligand_h'].detach(),
                    time_step=t,
                )
                router_active = self.pocket_router_mode != 'none' and (
                    self.pocket_router_topk > 0
                    or self.pocket_router_mode == 'cross_attn_gate'
                )
                if router_active:
                    router_weights_f = router_weights.unsqueeze(-1)
                    residue_tr_t = residue_tr_t * router_weights_f
                    residue_chi_t = residue_chi_t * router_weights_f
                    residue_rot_t = self.slerp_identity_to_q(
                        residue_rot_t,
                        1.0 - router_weights_f,
                    )
                if self.pocket_router_mode == 'cross_attn_gate':
                    # per-graph mean count so batch size does not inflate the number
                    sel = (router_weights > 0.5).float()
                    sel_per_graph = scatter_sum(sel, prot_update_batch.view(-1), dim=0)
                    router_selected_counts.append(float(sel_per_graph.mean().item()))
                else:
                    router_selected_counts.append(int((router_weights > 0).sum().item()))
                pred_protein_pos = apply_transforms_tensor_batch(
                    protein_pos=protein_pos,
                    protein_atom_name=data.protein_atom_name,
                    protein_atom_to_aa_name=data.protein_atom_to_aa_name,
                    protein_atom_to_aa_group=data.protein_atom_to_aa_group,
                    protein_element_batch=data.protein_element_batch,
                    rotations=residue_rot_t,
                    translations=residue_tr_t,
                    chi_update=residue_chi_t,
                    chi_mask=data.protein_chi_mask,
                    protein_translations_batch=prot_update_batch,
                )
                pred_protein_pos_offset = scatter_mean(pred_protein_pos, batch_protein, dim=0)
                pred_protein_pos -= pred_protein_pos_offset[batch_protein]
                protein_pos = pred_protein_pos
            else:
                pred_protein_pos = protein_pos

            if not pos_only:
                log_ligand_v_recon = F.log_softmax(v0_from_e, dim=-1)
                log_ligand_v = index_to_log_onehot(ligand_v, self.num_classes)
                log_model_prob = self.q_v_posterior(log_ligand_v_recon, log_ligand_v, t, batch_ligand)
                ligand_v_next = log_sample_categorical(log_model_prob)

                v0_pred_traj.append(log_ligand_v_recon.clone().cpu())
                vt_pred_traj.append(log_model_prob.clone().cpu())
                ligand_v = ligand_v_next

            hbap_ligand, hbap_protein = net_cond.extract_features(pred_ligand_pos, pred_protein_pos, pred_lig_a_h, gt_protein_a_h, gt_protein_r_h, batch_ligand, batch_protein)
            hbap_ligand, hbap_protein = hbap_ligand.detach(), hbap_protein.detach()

            ori_ligand_pos = ligand_pos + offset[batch_ligand]
            ori_protein_pos = protein_pos + offset[batch_protein]
            ligand_pos_traj.append(ori_ligand_pos.clone().cpu())
            protein_pos_traj.append(ori_protein_pos.clone().cpu())
            v_traj.append(ligand_v.clone().cpu())

        ligand_pos = ligand_pos + offset[batch_ligand]
        protein_pos = protein_pos + offset[batch_protein]
        protein_pos_holo = protein_pos_holo + offset[batch_protein]

        # Calculate the RMSD between predicted protein_pos with protein_pos_holo for each graph
        protein_pos_rmsds = []
        for i in range(max(batch_protein) + 1):
            mask = batch_protein == i
            protein_pos_i = protein_pos[mask]
            protein_pos_holo_i = protein_pos_holo[mask]
            pos_rmsd = torch.sqrt(((protein_pos_i - protein_pos_holo_i) ** 2).sum(-1)).mean()
            protein_pos_rmsds.append(pos_rmsd)
        # Calculate the TM-score between predicted protein_pos with protein_pos_holo
        protein_pos_tmscores = []
        for i in range(max(batch_protein) + 1):
            mask = batch_protein == i
            protein_pos_i = protein_pos[mask]
            protein_pos_holo_i = protein_pos_holo[mask]
            pos_tmscore = calculate_tm_score(protein_pos_i, protein_pos_holo_i)
            protein_pos_tmscores.append(pos_tmscore)

        return {
            'ligand_pos': ligand_pos,
            'protein_pos': protein_pos,
            'protein_pos_rmsd': protein_pos_rmsds,
            'protein_pos_tmscore': protein_pos_tmscores,
            'v': ligand_v,
            'ligand_pos_traj': ligand_pos_traj,
            'protein_pos_traj': protein_pos_traj,
            'v_traj': v_traj,
            'v0_traj': v0_pred_traj,
            'vt_traj': vt_pred_traj,
            'router_selected_counts': router_selected_counts,
        }


def extract(coef, t, batch):
    out = coef[t][batch]
    return out.unsqueeze(-1)


def quaternion_product(q1: torch.Tensor, q2: torch.Tensor) -> torch.Tensor:
    """
    Hamilton product of two quaternions, [w, x, y, z] layout.
    Vectorised, avoids Python-side unbind/stack.
    """
    w1, w2 = q1[..., :1], q2[..., :1]          # (..., 1)
    v1, v2 = q1[..., 1:], q2[..., 1:]          # (..., 3)
    w = w1 * w2 - (v1 * v2).sum(dim=-1, keepdim=True)
    v = w1 * v2 + w2 * v1 + torch.cross(v1, v2, dim=-1)

    return torch.cat([w, v], dim=-1)
