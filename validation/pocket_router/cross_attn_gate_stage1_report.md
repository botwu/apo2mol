# Ligand-Conditioned Induced-Fit Gate Stage-1 设计与进度汇报

更新时间：2026-06-10
项目位置：`/Users/wujian1/Downloads/a800_molecular/Apo2Mol`
相关设计文档：`validation/pocket_router/cross_attn_gate_design.md`

## 1. 当前一句话状态

当前实现的是 **stage-1 supervised induced-fit gate**：在预训练 Apo2Mol 生成器上增加一个轻量的 ligand→residue cross-attention 门控模块，预测每个 pocket residue 的可动权重 `w_i in [0, 1]`；门控由真实 apo→holo residue motion 构造的连续软标签监督，并在训练和采样时一致地乘到 residue translation / rotation / chi 更新上。

代码层面的关键问题已经修完，当前状态是：**可以上 A800 跑 50-step 训练 smoke，再跑 hard8 采样 smoke；暂不进入 stage-2 两遍 forward。**

## 2. 要解决的问题

原始 Apo2Mol 的 pocket 更新策略依赖 fixed rule，例如 top-k、shell radius、手调 shell weight。这类规则有两个局限：

1. active residue 的选择不是 ligand-conditioned 学出来的，而是人工规则。
2. fixed shell/top-k 的性能上限被规则本身限制，难以继续突破。

我们希望替换成一个可学习机制：给定当前 ligand 表示和 pocket residue 表示，由模型自己判断哪些 residue 应该释放、哪些 residue 应该锚定。

早期 v1 设想是“无辅助监督、完全端到端”，但 review 后发现它在当前单步去噪训练框架里不成立：门在 forward 后才乘到 protein update 上，ligand loss 已经由 gate 前的 `pred_ligand_pos` 计算完，因此 gate 拿不到 ligand loss 梯度。冻结 backbone 时，gate 主要从 protein loss 拿梯度，容易退化成“蛋白更新衰减器”。

stage-1 的修复是：给 gate 增加真实 induced-fit 位移软监督 `L_gate`，先让 gate 本身学稳。

## 3. Cross-Attention Gate 机制

实现文件：`models/ligand_residue_cross_attn.py`
核心类：`LigandResidueCrossAttnGate`

对 batch 内每个 complex 独立计算：

```text
Q = W_q(ligand_h)         # ligand atom as query
K = W_k(residue_h)        # pocket residue as key
sim = Q K^T / sqrt(d)

logit_i = max_over_ligand_atoms(sim[:, i]) + bias
w_i = sigmoid(logit_i)
```

设计含义：

- ligand atom 作为 query，residue 作为 key。
- 对每个 residue，取所有 ligand atom 中最强的 attention logit。
- `w_i` 是该 residue 的释放强度：越接近 1 越允许动，越接近 0 越锚定。
- 模块很轻：两个无 bias linear projection 加一个 scalar bias。
- 当前维度：`ligand_h=128`，`residue_h=131`，`inner_dim=64`。

这不是完整 Transformer cross-attention 融合层，而是 **cross-attention scoring gate**：用 ligand-conditioned attention 分数产生 active-set 权重。

## 4. Induced-Fit 软监督

训练数据里本来就有真实 apo→holo residue 运动：

- `prot_tr`：residue translation。
- `prot_rot`：residue frame rotation quaternion。
- `prot_chi_update`：side-chain chi torsion delta。

我们把三类运动统一换算成 Å 等效位移：

```text
chi_delta_i = wrap(prot_chi_update_i)  # wrap to (-pi, pi]

eff_motion_i = ||prot_tr_i||
             + gate_rot_radius_A * angle(prot_rot_i)
             + gate_chi_radius_A * ||chi_delta_i * chi_mask_i||

target_i = sigmoid(gate_target_sharpness * (eff_motion_i - gate_target_center_A))
L_gate   = BCEWithLogits(logit_i, target_i)
```

当前默认：

```yaml
gate_target_center_A: 0.75
gate_target_sharpness: 4.0
gate_rot_radius_A: 3.0
gate_chi_radius_A: 1.5
gate_loss_weight: 1.0
```

这组参数形成一个连续的 core / shell / background 软结构：

| effective motion | target | 含义 |
|---:|---:|---|
| `0.0 Å` | `0.047` | 背景 residue，基本锚定 |
| `0.5 Å` | `0.269` | shell residue，弱松弛，接近旧实验最优 shell weight 0.25 |
| `1.5 Å` | `0.953` | core residue，强释放 |

关键点：这不是手写 core/shell/background 三档权重，而是由真实 apo→holo induced-fit 位移场连续生成的软标签。旧 fixed-rule 的经验先验被内生化为可学习目标。

## 5. 门控如何作用到 protein update

在 `get_diffusion_loss` 和 `sample_diffusion` 中，得到 `w_i` 后应用到 residue update：

```text
pred_res_tr  = pred_res_tr  * w_i
pred_res_chi = pred_res_chi * w_i
pred_res_rot = slerp(identity, pred_res_rot, lambda = 1 - w_i)
```

语义：

- `w_i=1`：完整应用该 residue 的 translation / chi / rotation。
- `w_i=0`：translation 和 chi 归零，rotation 插回 identity，即该 residue 基本不动。
- 训练和采样共用同一套门控逻辑，避免 train / inference 行为不一致。

总损失：

```text
L_total = L_ligand_pos
        + 100 * L_v
        + L_prot_tr
        + L_prot_rot
        + 5 * L_prot_chi
        + gate_loss_weight * L_gate
```

stage-1 中，`L_gate` 是 gate 的主驱动；protein losses 作为一致性约束。ligand loss 暂时不会反传到 gate，这留到 stage-2。

## 6. 冻结和 checkpoint 策略

训练入口：`train_pl.py`

推荐 stage-1 配置：

```yaml
model:
  pocket_router_mode: cross_attn_gate
  cross_attn_gate:
    inner_dim: 64
    gate_loss_weight: 1.0
    gate_target_center_A: 0.75
    gate_target_sharpness: 4.0
    gate_rot_radius_A: 3.0
    gate_chi_radius_A: 1.5

train:
  freeze_backbone: true
  freeze_residue_head: true
  pretrained_ckpt: <apo2mol_pretrained_ckpt>
```

冻结策略：

- `cross_attn_gate.*` 保持 trainable。
- 默认 `res_inference.*` 冻结。
- 其余 Apo2Mol backbone 全冻结。
- optimizer 只接收 `requires_grad=True` 的参数。

fail-fast 保护：

- 如果 `freeze_backbone=true` 但没有设置 `train.pretrained_ckpt`，`train_pl.py` 会直接报错。
- 这样避免“冻结随机初始化 backbone，只训练 gate”的无效训练。

采样 checkpoint 保护：

- `sample_split.py` 在 `cross_attn_gate` 模式下会检查 checkpoint 是否包含 `cross_attn_gate.*` 权重。
- 如果使用原始 Apo2Mol checkpoint 去采样，会直接报错，避免随机 gate 被用于评估。

## 7. 已完成修复清单

| # | 问题 | 当前修复 |
|---|---|---|
| 1 | v1 gate 拿不到 ligand loss 梯度，冻结时退化成 protein shrinkage | 加 `L_gate`，用真实 induced-fit motion 做软监督 |
| 2 | `freeze_backbone=true` 但无 checkpoint 也能跑 | `train_pl.py` fail-fast |
| 3 | sampling selected count 被 batch size 放大 | cross-attn gate 分支改为 per-graph mean selected count |
| 4 | quaternion normalize 可能除零 NaN | `slerp_identity_to_q`、sampling normalize、noise normalize 加 clamp |
| 5 | 50-step smoke 看不到 router 诊断 | `loss_gate/router_w_min/max/mean/router_selected_per_graph` 每步 log |
| 6 | target 中 chi 角差跨 `0/2π` 周期边界会误标 core | gate target 内使用 `atan2(sin(delta), cos(delta))` wrap 到 `(-π, π]` |
| 7 | 文档和代码表述不一致 | docstring、设计文档公式、checkpoint 加载说明已同步 |

## 8. 当前验证状态

本地已做：

- `py_compile` 通过。
- 合成 BCE backward 通过，`cross_attn_gate.to_q / to_k / bias` 均有非零梯度。
- target 数值验证通过：
  - `0.0 Å -> 0.047`
  - `0.5 Å -> 0.269`
  - `1.5 Å -> 0.953`
- chi wrap 验证通过：
  - raw delta `6.22 rad`
  - wrapped delta `-0.063 rad`

本地没有完整 GPU 训练；下一步必须在 A800 上做 smoke。

## 9. A800 下一步验证计划

### 9.1 训练 smoke

目标：冻结 backbone，只训练 gate，跑约 50 steps。

观察指标：

- `train/loss_gate` 应该下降。
- `train/router_w_mean` 不应贴近 0 或 1，健康区间先看 `[0.1, 0.5]`。
- `train/router_w_min` 和 `train/router_w_max` 不应全部饱和。
- `train/router_selected_per_graph` 应落在个位数。
- 无 NaN。

### 9.2 hard8 采样 smoke

目标：用训练出的 gate checkpoint 跑 hard8 中少量 case。

观察指标：

- `router_selected_counts` 是 per-graph 平均，不再被 batch size 放大。
- 建议健康范围先看 `[3, 10]`。
- 生成流程 complete，不出现 NaN、不出现明显 clash 恶化。

### 9.3 hard8 / hard20 A/B

先 hard8，再 hard20。对比：

- original static5
- distance_top4_hard
- active_set_shell4_w0_25
- pocket_router_cross_attn_gate

hard20 对标 fixed shell4 w0.25：

- protein RMSD ≤ `1.5935`
- selected_per_graph 在 `[4, 8]`
- complete ≥ `0.8500`
- clash ≤ `1.90`
- boundary jump < static5 的 `0.8569`

## 10. 风险和预案

| 风险 | 现象 | 处理 |
|---|---|---|
| gate 过稀疏 | `router_w_mean` 贴 0，selected 很少 | 降低 `gate_target_center_A` 或增大 `gate_loss_weight` |
| gate 过密 | `router_w_mean` 贴 1，selected 过多 | 增大 `gate_target_center_A` |
| 选择过散 | selected 很多但质量不稳定 | 增大 `gate_target_sharpness` |
| protein loss 与 gate loss 拉扯 | `loss_gate` 降但 protein RMSD 不改善 | stage-1 后可考虑解冻 `res_inference`，不要立刻上全 backbone |
| ligand 质量未提升 | protein 指标改善但 ligand 不变 | 这是 stage-1 预期边界，stage-2 两遍 forward 再解决 |

## 11. 为什么现在不上 stage-2

stage-2 设想：两遍 forward。

```text
first forward -> gate -> update pocket
second forward on gated pocket -> ligand prediction
ligand loss -> gate
```

这会让 ligand loss 真正反传到 gate，实现完整端到端 induced-fit。但它会同时引入：

- 约 2 倍 forward 计算和显存压力。
- 训练分布变化。
- gate、protein head、ligand head 的耦合不稳定。
- 排错复杂度显著上升。

当前最稳妥路线是先让 stage-1 supervised gate 在 hard8/hard20 上证明有效，再决定是否进入 stage-2。

## 12. 对外汇报口径

我们目前完成了一个 ligand-conditioned induced-fit active-set gate。它不是简单复刻 fixed shell rule，而是用 ligand→residue cross-attention 预测 residue 释放强度，并用真实 apo→holo 位移构造连续软监督。这个软监督天然形成 background / shell / core 梯度：不动 residue 被锚定，shell residue 弱松弛，core residue 强释放。与原来的 fixed-rule active set 相比，新方法把手写规则转化为从真实 induced-fit field 学到的连续可学习门控。

当前代码已完成 stage-1：gate 训练、冻结策略、checkpoint 保护、采样 per-graph 统计、NaN 防护、chi 周期 wrap 和日志诊断都已补齐。本地通过语法检查和合成梯度验证。下一步是在 A800 上跑 50-step smoke 和 hard8 采样，确认 `loss_gate` 下降、门控不塌缩、selected_per_graph 合理；通过后再进入 hard20 A/B。stage-2 的两遍 forward 暂缓，避免一次引入过多不确定性。

## 13. 关键代码位置

| 文件 | 作用 |
|---|---|
| `models/ligand_residue_cross_attn.py` | cross-attention gate 模块 |
| `models/molopt_score_model.py` | gate 初始化、target 构造、训练 loss、采样 gate 应用 |
| `models/pl_model.py` | 冻结策略、optimizer 参数选择、每步 router 诊断日志 |
| `train_pl.py` | `pretrained_ckpt` 加载和 freeze fail-fast |
| `configs/training.yaml` | gate 参数和 freeze 参数 |
| `sample_split.py` | 采样 checkpoint gate 权重校验 |
| `validation/run_new_method_ab.py` | A/B arm：`pocket_router_cross_attn_gate` |
| `validation/pocket_router/cross_attn_gate_design.md` | 详细技术设计文档 |
