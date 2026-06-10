# Ligand-Conditioned Induced-Fit Gate（pocket router 可学习方案 v2）

更新时间：2026-06-10
位置：`/Users/wujian1/Downloads/a800_molecular/Apo2Mol`

## 0. 一句话结论

> 用一个**ligand→residue cross-attention** 预测每个 pocket residue 的 active-set 门控 `w_i ∈ [0,1]`，并用**训练里已有的真实 apo→holo 残基位移**作连续软监督，教会门"给定当前 ligand，哪些 residue 会被诱导移动"。门同时乘到 residue 的 translation / chi / rotation 更新上，使训练和采样行为一致。

这是对 v1（纯端到端、无监督）的修订。v1 在代码 review 中被发现有一个架构级硬伤（门拿不到 ligand 梯度，会退化成"蛋白更新衰减器"），本版从机制上修掉它，并对齐近两年 induced-fit 文献的主流做法。

## 1. 为什么 v1（纯端到端）不 work

v1 的设想是"门乘到 residue 更新上，梯度从下游 ligand/protein 损失端到端流回"。但实际代码里：

- `forward` 一次性算出 `pred_ligand_pos` 和 `pred_residue_*`（`molopt_score_model.py:849`）。
- 门在 `molopt_score_model.py:885` 之后才乘到 `pred_res_tr/chi/rot`。
- 而 **ligand loss 用的是门之前的 `pred_ligand_pos`**（`molopt_score_model.py:910`）。

后果：训练是**单步去噪**，没有"门改 pocket → 重新 forward → 影响 ligand"的闭环。冻结 backbone 时，门**唯一的梯度来自 protein L1 loss**，等价于学一个标量收缩 `w·frozen_pred ≈ inverse_target`。大多数 residue 目标运动很小，会把 `w` 往低值推 → 门塌缩成"蛋白更新衰减器"，而不是 ligand-conditioned active-set selector。

所以 v1 文档里"通过 backbone 间接影响 ligand 损失"是不成立的。

## 2. 借鉴的近年工作

- **DynamicBind**(Nat. Commun. 2024)：用等变扩散网络学 ligand 特异的蛋白构象变化，**直接用真实 holo 构象作监督**（morph-like 变换造 decoy），从 apo 恢复 holo。→ 启示：induced-fit 应该用真实 apo→holo 位移监督，而不是靠间接梯度。
- **FABFlex**(ICLR 2025)：明确"**只聚焦 pocket core 区域的构象变化，而非整个蛋白**"，因为 pocket 只占蛋白一小部分，能大幅降低计算量。→ 启示：sparse active-set（少量 residue）是已发表的、被认可的设计哲学。

我们的定位：在一个**已预训练**的 Apo2Mol 生成器之上，加一个**轻量、可学习、ligand-conditioned 的稀疏门**，预测 induced-fit active set——结合 DynamicBind 的"真实位移监督"和 FABFlex 的"pocket 聚焦"，但不重训整个网络。

## 3. 方案设计

### 3.1 门控模块（cross-attention）

文件：`models/ligand_residue_cross_attn.py`，类 `LigandResidueCrossAttnGate`。逐 graph 计算：

```
Q = W_q · ligand_h        # ligand 作 query
K = W_k · residue_h        # residue 作 key
sim = Q · Kᵀ / sqrt(d)
logit_i = max_over_ligand( sim[:, i] ) + b      # 每个 residue 取最强 ligand 原子
w_i = sigmoid(logit_i)
```

- aggregation 用 **max over ligand atoms**：residue 只要被任一 ligand 原子强烈关注就该被释放。
- 输出 **logit**，sigmoid 在调用点做（BCE-with-logits 数值稳定）。
- 极轻量：`W_q`、`W_k` 两个无 bias 线性层 + 标量 bias。
- 维度：`ligand_h = hidden_dim = 128`，`residue_h = hidden_dim + 3 = 131`，`inner_dim = 64`。

### 3.2 induced-fit 软监督（核心修复）

训练里 `get_diffusion_loss` 已有 ground-truth 的 apo→holo 残基运动：`prot_tr`（平移）、`prot_rot`（旋转四元数）、`prot_chi_update`（侧链）。由此构造每个 residue 的运动幅度并映射成连续软目标：

把三种运动**统一折算成 Å 等效位移**（平移本就是 Å，旋转角×`rot_radius`、侧链扭转×`chi_radius`），再映射成软目标：

```
eff_motion_i = ||prot_tr_i|| + rot_radius · angle(prot_rot_i) + chi_radius · ||wrap(prot_chi_update_i) · mask||
target_i     = sigmoid( sharpness · (eff_motion_i − center_A) )   # 软标签 ∈ (0,1)，detach
loss_gate    = BCEWithLogits( logit_i, target_i )                  # 逐 graph 平均后再平均
```

- 这是**真实的 induced-fit 位移场**，不是 v1 LASC 那套"4 个硬编 Å 阈值去模仿 shell4"。
- 默认 `center_A=0.75`、`sharpness=4.0` 下，软目标天然形成 core/shell/background 梯度：
  - 刚性 residue（eff_motion≈0）→ **~0.05**（背景，锚定）
  - shell 级 residue（≈0.5Å）→ **~0.27**（弱松弛，正好对应实证最优 shell4 w0.25）
  - core residue（≥1.5Å）→ **~1**（强释放）
- 单位严格统一为 Å 等效位移（`rot_radius`/`chi_radius` 把弧度折算成 Å），所以 `center_A` 是一个真正的 Å 尺度，论文里可写成 "effective induced-fit displacement"。
- 门由此拿到**正确、强、不塌缩**的梯度，彻底解决第 1 节的退化问题。

### 3.3 门作用到 pocket 更新

```
pred_res_tr  = pred_res_tr  · w_i
pred_res_chi = pred_res_chi · w_i
pred_res_rot = slerp(identity, pred_res_rot, λ = 1 − w_i)
```

训练时同样乘进去（再算 protein loss），保证**训练/采样行为一致**；这条 protein-loss 通路与 induced-fit 监督方向一致（两者都想让真动的 residue `w→1`、不动的 `w→0`），起耦合/一致性作用，监督 loss 是主驱动。

### 3.4 总损失

```
L_total = L_ligand_pos + 100·L_v + L_prot_tr + L_prot_rot + 5·L_prot_chi + gate_loss_weight·L_gate
```

冻结 backbone + residue head 时，前面的生成损失里只有 `L_prot_*` 随门变化，`L_gate` 是主驱动。

## 4. 第一阶段训练

冻结策略（`pl_model.py` `_apply_freeze_policy`）：只放开 `cross_attn_gate.*`（可选 `res_inference.*`），其余全冻。优化器只收 `requires_grad=True` 的参数。

加载：训练入口 `train_pl.py` 通过 `train.pretrained_ckpt` + `load_state_dict(strict=False)`（`load_pretrained_weights`）加载 Apo2Mol backbone，`cross_attn_gate.*` 从零初始化。

**强制约束（已加 fail-fast）**：`freeze_backbone=true` 但未给 `pretrained_ckpt` 时 `train_pl.py` 直接报错，避免冻结随机初始化 backbone 只训门这种"看似能跑实则全废"的情况。

## 5. 代码 review 修复清单

针对 review 提出的 5 个问题：

| # | 问题 | 修复 |
|---|---|---|
| 1 | 门拿不到 ligand 梯度，退化成蛋白衰减器 | 加 induced-fit 软监督 `L_gate`（3.2），门有了正确梯度来源 |
| 2 | `freeze_backbone=true` 但无 ckpt 不报错 | `train_pl.py` 加 fail-fast |
| 3 | 采样 selected count 是 batch 总数 | `sample_diffusion` 改 per-graph 平均（`scatter_sum` over `protein_translations_batch`） |
| 4 | slerp `q/q.norm()` 可能 NaN | `slerp_identity_to_q` 两处 norm 加 `.clamp(min=1e-8)` |
| 5 | 50 iter 看不到 router 日志 | router 诊断移出 `train_report_iter` 块，**每步**都 log |

## 6. 配置

```yaml
model:
  pocket_router_mode: cross_attn_gate
  cross_attn_gate:
    inner_dim: 64
    gate_loss_weight: 1.0       # induced-fit 监督权重
    gate_target_center_A: 0.75  # eff 位移软目标过 0.5 的位置 (A)
    gate_target_sharpness: 4.0  # 位移->软目标的软硬（刚性 residue->~0.05）
    gate_rot_radius_A: 3.0      # 旋转角(rad)->Å 等效位移
    gate_chi_radius_A: 1.5      # 侧链扭转(rad)->Å 等效位移
train:
  freeze_backbone: true
  freeze_residue_head: true
  pretrained_ckpt: <apo2mol_ckpt>   # 冻结时必填，否则报错
```

## 7. 怎么跑（GPU 机）

本地无 torch，只做 AST/lint。实际在 A800：

1. **训练烟雾**（冻结 backbone，~50 iters）：观察每步 log 的 `train/loss_gate` 下降；`train/router_w_mean` 不塌缩到 0 或 1（健康 `[0.1, 0.5]`）；`train/router_selected_per_graph` 落在个位数；无 NaN。
2. **采样烟雾**（1 个 hard8 case）：`router_selected_counts`（per-graph 平均，门 >0.5）均值落 `[3, 10]`。
3. **A/B**：先 hard8，过了再 hard20，对比 `original static5` / `distance_top4_hard` / `active_set_shell4_w0_25` / `pocket_router_cross_attn_gate`。

## 8. 第一阶段成功标准

对标 fixed shell4 w0.25（hard20）：

- protein RMSD ≤ `1.5935`（有望突破，因为不再被规则封顶）
- selected_per_graph（门 >0.5）∈ `[4, 8]`
- complete 不低于 `0.8500`
- clash 不高于 `1.90`
- boundary jump < static5 的 `0.8569`

健康信号异常时的调参（不写死、可配置）：
- `router_w_mean` 贴 0 → 调小 `gate_target_center_A` 或调大 `gate_loss_weight`。
- `router_w_mean` 贴 1 / 选太多 → 调大 `gate_target_center_A`。
- 选择过散 → 调大 `gate_target_sharpness`。

## 9. 论文叙事

- DiffSBDD / TargetDiff：pocket 完全冻结，不选。
- DynamicBind：全 residue 预测构象变化，无 gate。
- PocketGen / FABFlex：固定半径/core 区域，pocket 聚焦但不学选择。
- 本方案：**用 ligand→residue cross-attention 学习 ligand 诱导的稀疏 active-set，由真实 induced-fit 位移监督**，在预训练生成器上即插即用。

> 创新点不是"加了一个 router"，而是"在预训练 SBDD 生成器上，用 ligand-conditioned cross-attention 学习 ligand 诱导的稀疏 pocket 自由度释放（induced-fit active set）"。

## 10. 后续阶段（先不做）

- stage-2：解冻 `res_inference` 联合微调；并可改成**两遍 forward**（门改 pocket → 再 forward 预测 ligand），让 ligand loss 也真正流回门，实现完整端到端 induced-fit。
- stage-3：迁移到 ligand-conditioned active-set flow matching（见 `active_set_flow_framework_full.md`）。

第一阶段未在 hard20 过线前，不进 stage-2 / flow matching。
