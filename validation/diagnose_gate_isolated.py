"""Run the gate forward/backward in isolation with several alternative
implementations to localise which op produces the nan grad."""

import os, sys
sys.path.insert(0, ".")
import warnings; warnings.filterwarnings("ignore")
import math
import torch
import torch.nn as nn
import torch.nn.functional as F


def make_inputs(num_residues=77, num_lig=20, num_graphs=2, ligand_dim=128,
                residue_dim=131, seed=0):
    torch.manual_seed(seed)
    device = torch.device("cuda")
    ligand_h = torch.randn(num_lig, ligand_dim, device=device, requires_grad=False)
    residue_h = torch.randn(num_residues, residue_dim, device=device, requires_grad=False)
    # split residues + ligand atoms across graphs
    res_per = num_residues // num_graphs
    lig_per = num_lig // num_graphs
    prot_update_batch = torch.cat([
        torch.full((res_per,), 0, device=device),
        torch.full((num_residues - res_per,), 1, device=device),
    ])
    batch_ligand = torch.cat([
        torch.full((lig_per,), 0, device=device),
        torch.full((num_lig - lig_per,), 1, device=device),
    ])
    return ligand_h, residue_h, batch_ligand, prot_update_batch


def gate_v1_index_copy(ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias):
    """Original implementation."""
    inner = to_q.out_features
    scale = 1.0 / math.sqrt(inner)
    num_residues = residue_h.shape[0]
    logits = residue_h.new_zeros(num_residues)
    q_all = to_q(ligand_h)
    k_all = to_k(residue_h)
    num_graphs = int(prot_update_batch.max().item()) + 1
    for g in range(num_graphs):
        res_idx = (prot_update_batch == g).nonzero(as_tuple=True)[0]
        lig_idx = (batch_ligand == g).nonzero(as_tuple=True)[0]
        if res_idx.numel() == 0 or lig_idx.numel() == 0:
            continue
        q = q_all[lig_idx]
        k = k_all[res_idx]
        sim = q @ k.transpose(0, 1) * scale
        per_res = sim.max(dim=0).values + bias
        logits = logits.index_copy(0, res_idx, per_res.to(logits.dtype))
    return logits


def gate_v2_scatter(ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias):
    """index_copy -> scatter via index assignment."""
    inner = to_q.out_features
    scale = 1.0 / math.sqrt(inner)
    num_residues = residue_h.shape[0]
    logits = residue_h.new_zeros(num_residues)
    q_all = to_q(ligand_h)
    k_all = to_k(residue_h)
    num_graphs = int(prot_update_batch.max().item()) + 1
    for g in range(num_graphs):
        res_idx = (prot_update_batch == g).nonzero(as_tuple=True)[0]
        lig_idx = (batch_ligand == g).nonzero(as_tuple=True)[0]
        if res_idx.numel() == 0 or lig_idx.numel() == 0:
            continue
        q = q_all[lig_idx]
        k = k_all[res_idx]
        sim = q @ k.transpose(0, 1) * scale
        per_res = sim.max(dim=0).values + bias
        # use scatter_add on a fresh zero
        logits = logits + torch.zeros_like(logits).scatter(0, res_idx, per_res.to(logits.dtype))
    return logits


def gate_v3_global(ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias):
    """Compute Q*K^T globally then mask same-graph pairs and reduce per residue
    without per-graph python loop."""
    inner = to_q.out_features
    scale = 1.0 / math.sqrt(inner)
    q_all = to_q(ligand_h)
    k_all = to_k(residue_h)
    sim = q_all @ k_all.transpose(0, 1) * scale  # [L, R]
    same_graph = (batch_ligand[:, None] == prot_update_batch[None, :])
    sim = sim.masked_fill(~same_graph, float("-inf"))
    per_res = sim.max(dim=0).values + bias  # [R]
    return per_res


def run(name, fn, ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias, target):
    for p in (to_q.weight, to_k.weight, bias):
        if p.grad is not None:
            p.grad = None
    logits = fn(ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias)
    print(f"  {name}: logits shape={tuple(logits.shape)} min={logits.min().item():.4g} "
          f"max={logits.max().item():.4g} req_grad={logits.requires_grad}")
    bce = F.binary_cross_entropy_with_logits(logits, target, reduction="mean")
    print(f"    bce={bce.item():.4g}")
    bce.backward()
    print(f"    grad q.w nan={torch.isnan(to_q.weight.grad).any().item()} "
          f"norm={to_q.weight.grad.norm().item():.4g}")
    print(f"    grad k.w nan={torch.isnan(to_k.weight.grad).any().item()} "
          f"norm={to_k.weight.grad.norm().item():.4g}")
    print(f"    grad bias nan={torch.isnan(bias.grad).any().item()} "
          f"norm={bias.grad.norm().item():.4g}")


def main():
    torch.set_float32_matmul_precision("high")
    device = torch.device("cuda")

    ligand_h, residue_h, batch_ligand, prot_update_batch = make_inputs()
    inner = 64
    to_q = nn.Linear(128, inner, bias=False).to(device)
    to_k = nn.Linear(131, inner, bias=False).to(device)
    bias = nn.Parameter(torch.zeros(1, device=device))

    # binary target like induced_fit
    torch.manual_seed(1)
    target = torch.bernoulli(torch.full((residue_h.shape[0],), 0.3, device=device)).float()

    print("=== fp32 (no autocast) ===")
    run("v1_index_copy", gate_v1_index_copy, ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias, target)
    run("v2_scatter",    gate_v2_scatter,    ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias, target)
    run("v3_global",     gate_v3_global,     ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias, target)

    print()
    print("=== bf16 autocast ===")
    with torch.amp.autocast(device_type="cuda", dtype=torch.bfloat16):
        run("v1_index_copy", gate_v1_index_copy, ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias, target)
        run("v2_scatter",    gate_v2_scatter,    ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias, target)
        run("v3_global",     gate_v3_global,     ligand_h, residue_h, batch_ligand, prot_update_batch, to_q, to_k, bias, target)


if __name__ == "__main__":
    main()
