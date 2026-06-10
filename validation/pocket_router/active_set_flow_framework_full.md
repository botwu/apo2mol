# Ligand-Conditioned Active-Set Flow Matching：策略、实验和最终框架

更新时间：2026-06-03 12:15 CST

这份文档用于系统整理当前项目的完整思路，包括：

- 目前我们到底在验证什么；
- 为什么这个方向不是简单的“只更新部分 residues”；
- 当前 hard8 和 hard20 实验怎么设计；
- 目前实验结果说明了什么；
- 最后真正要做成什么样的 flow matching 框架；
- 如果要冲 AAAI 风格顶会，后面还缺哪些关键实验和工程。

核心结论先放在前面：

> 我们的 idea 不是“少更新几个残基”，而是把 apo-to-holo pocket adaptation 建模为 ligand-conditioned active-set optimization。
>
> 大白话说：ligand 生成过程中，模型要判断哪些 pocket residues 需要认真动，哪些邻居需要轻轻跟着动，哪些背景区域必须保持稳定。最终方法应该用 flow matching 来统一建模 ligand 生成和 pocket adaptation，而不是停留在手写规则的 sparse update。

## 1. 当前研究问题

### 1.1 我们真正想解决什么问题

结构生成式药物设计里面，很多方法把 protein pocket 当成一个静态条件：

```text
给定 pocket -> 生成 ligand
```

但在 apo-to-holo 场景下，这个假设是不够真实的。

真实情况往往是：

- apo pocket 不是天然处于最适合 ligand 结合的状态；
- ligand 结合会诱导 pocket 局部构象变化；
- 不是整个 protein 都应该动；
- 不是只有离 ligand 最近的几个 residues 会有影响；
- selected residues 周围的邻居可能也需要小幅松弛；
- 背景结构必须稳定，否则 protein 会被生成过程带偏。

所以我们最后要解决的任务不是：

```text
apo pocket -> ligand
```

而应该是：

```text
apo pocket -> ligand + adapted holo-like pocket
```

也就是说，模型不只是生成一个 ligand，还要同时给出一个 ligand 诱导后的 pocket adaptation。

大白话：

> ligand 不是在一个死的 pocket 里生成。ligand 生成过程中，pocket 也应该发生合理的局部适应。

### 1.2 为什么不能简单全 pocket 更新

一种直接想法是：

```text
每一步都让 ligand 和所有 pocket residues 做 dense interaction
```

问题是，这会带来很多噪声。

因为当前 ligand fragment 可能只和 pocket 的一小部分区域真正相关。如果把整个 pocket 都放开：

- 远处 residue 可能被不必要地扰动；
- 背景 protein 会漂移；
- protein RMSD 可能变差；
- ligand 生成会受到无关 pocket motion 的干扰；
- 模型更难学习到“哪些自由度真的该动”。

大白话：

> ligand 只碰到了 pocket 的一小块区域，却让整个 pocket 都跟着动，这不符合物理直觉，也容易把结构搞乱。

### 1.3 为什么不能只更新 top-k residues

另一种直接想法是：

```text
只更新离当前 ligand 最近的 top-k residues
```

这个想法比 dense update 更合理，但仍然不够。

因为只动 selected residues，会造成几个问题：

- selected residues 动了，周围 unselected residues 不动，局部结构可能不连续；
- core residue 和 background residue 之间可能出现硬切边界；
- ligand-protein clash 可能变多；
- protein RMSD 可能改善，但 ligand geometry 或 binding contact 可能变差；
- 只看 ligand 附近 residues，可能忽略二级影响。

大白话：

> 只动最近的几个残基，就像只拽动衣服上的几个点，周围布料完全不跟着动，局部会很别扭。

### 1.4 为什么 random sparse 是危险负例

我们之前实验发现，random sparse 有时 protein RMSD 看起来并不差。

但 hard8 的几何诊断显示：

```text
random_top4 protein RMSD 看起来接近 distance_top4
但 ligand-protein clash 非常高
```

这说明：

> 只看 protein RMSD 会被 random sparse 骗。

大白话：

> protein 看起来被调得更像 holo 了，但 ligand 可能已经被 pocket 挤坏了。

这也是为什么后续论文不能只汇报 protein RMSD，必须同时汇报：

- ligand stability；
- ligand-protein clash；
- contact recall / precision；
- boundary jump；
- QED / SA / Lipinski；
- docking 或 binding surrogate 指标。

## 2. 当前核心 idea

### 2.1 一句话版本

当前最核心的创新点应该表述为：

> Apo-to-holo SBDD should be formulated as ligand-conditioned active-set pocket adaptation.

中文：

> apo-to-holo 分子生成不应该把 pocket 当成静态条件，也不应该让整个 pocket 自由乱动，而应该让当前 ligand 状态动态决定 pocket 哪些自由度释放、哪些邻域弱松弛、哪些背景保持稳定。

### 2.2 Active-set 的三个区域

我们把 pocket residues 分成三类：

| 区域 | 含义 | 更新方式 |
|---|---|---|
| core | 被当前 ligand 强烈诱导的核心 residues | 强更新 |
| shell | core 周围的邻居 residues | 弱松弛 |
| background | 当前不相关或远处 residues | 锚定或近似不动 |

公式上可以写成：

```text
pocket_update_i =
  core_i * strong_update_i
  + shell_i * weak_update_i
  + background_i * anchor_update_i
```

大白话：

```text
核心残基：认真动
邻居残基：轻轻动
背景残基：别乱动
```

### 2.3 为什么这个 idea 比 sparse update 更强

“只更新 top-k residues”只是一个 hard mask。

我们的 active-set 思想多了两个关键点：

1. shell relaxation

   core 周围的 residues 不能完全僵住，否则会有边界不连续。

2. background anchoring

   背景不能乱动，否则 whole-pocket noise 会破坏结构稳定。

所以我们不是在说：

```text
只动少数 residues 就行
```

而是在说：

```text
ligand-conditioned core release
+ local shell relaxation
+ background anchoring
```

这才是完整的 pocket adaptation 机制。

## 3. 当前实验策略

### 3.1 当前实验不是最终方法

目前正在跑的实验是 fixed-rule validation，不是最终投稿方法。

它的作用是回答：

> active-set 这个第一性原理思想是否值得继续做成 learnable model？

更具体地说，当前实验在验证：

1. sparse pocket update 是否比 dense pocket update 更合理；
2. ligand-near sparse 是否比 random sparse 更安全；
3. top4 hard update 是否会有局部不连续或 clash 风险；
4. 加 shell weak relaxation 是否能缓解这些风险；
5. protein RMSD、ligand stability、clash、contact 是否会给出不同结论；
6. `core + shell + background` 是否比单纯 top-k 更像合理机制。

### 3.2 当前实验和最终方法的关系

当前实验：

```text
手写规则 active-set
distance top4 core
固定 shell 半径
固定 shell 权重
```

最终方法：

```text
learnable active-set router
ligand-conditioned active weights
flow matching joint ligand-pocket generation
core/shell/background 动态预测
```

所以当前实验只能证明机制是否有信号，不能作为最终顶会方法。

大白话：

> 现在是在证明“这条路值得走”，不是证明“最终模型已经完成”。

## 4. hard20 实验细节

### 4.1 实验目录

当前 hard20 实验目录：

```text
validation/ab_runs/hard20_active_set_candidate_repeat_steps1000_n3
```

主日志：

```text
validation/ab_runs/hard20_active_set_candidate_repeat_steps1000_n3/hard20_tmux_cpu.log
```

tmux 会话：

```text
apo2mol_hard20
```

### 4.2 启动命令

当前实验使用的命令是：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --active-set-candidate-repeat \
  --run-dir validation/ab_runs/hard20_active_set_candidate_repeat_steps1000_n3 \
  --num-cases 20 \
  --num-samples 3 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none \
  --device cpu
```

这里固定 `--device cpu` 的原因是：

- `--device auto` 在 tmux 里可能选到 Apple MPS；
- 但当前项目里的 `torch_cluster.knn` 路径要求 CPU tensor；
- 如果用 MPS，会出现 `x must be CPU tensor` 错误；
- 所以为了稳定运行，当前实验固定 CPU。

### 4.3 为什么这个实验这么慢

完整 hard20 是：

```text
6 arms x 20 cases x 3 samples x 1000 diffusion steps
= 360,000 denoising steps
```

每一步不是简单算一个指标，而是要跑：

- ligand-protein 图构建；
- KNN graph；
- equivariant message passing；
- PMINet 条件特征；
- protein-ligand joint update；
- diffusion denoising。

大白话：

> 这不是在算 120 个 RMSD，而是在生成 120 组 ligand-pocket 轨迹，每组还有 3 个 samples，每个 sample 要走 1000 步。

### 4.4 hard20 的 6 个 arms

| arm | 目的 |
|---|---|
| `pocket_router_random_top4` | 随机选 4 个 residues，负对照，用来证明 random sparse 是否危险 |
| `pocket_router_distance_top4_hard` | 只强更新 ligand 最近的 top4 residues，强手工 baseline |
| `active_set_distance_top4_shell4_w025` | 主候选，top4 core 强更新，4A shell 以 0.25 权重弱松弛 |
| `active_set_distance_top4_shell3_w010` | 更保守的 shell，3A shell，0.10 权重 |
| `active_set_distance_top4_shell3_w025` | 3A shell，0.25 权重 |
| `active_set_distance_top4_shell4_w010` | 4A shell，0.10 权重 |

### 4.5 当前 hard20 进度

截至 2026-06-03 12:15 CST：

```text
已完成 result files: 78 / 120
总体进度: 65.0%
```

分 arm 进度：

| arm | 当前进度 |
|---|---:|
| `pocket_router_random_top4` | 20 / 20 |
| `pocket_router_distance_top4_hard` | 20 / 20 |
| `active_set_distance_top4_shell4_w025` | 20 / 20 |
| `active_set_distance_top4_shell3_w010` | 18 / 20 |
| `active_set_distance_top4_shell3_w025` | 0 / 20 |
| `active_set_distance_top4_shell4_w010` | 0 / 20 |

因为全量还没跑完，所以最终汇总文件还没有生成：

```text
results.json                 未生成
geometry_diagnostics.json    未生成
ligand_quality.json          未生成
```

## 5. 当前 hard20 阶段性结果

目前 hard20 已经完整 eval 的 arm 有三个。

| arm | samples | mol stable | atom stable | recon success | complete |
|---|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 60 | 0.2833 | 0.9191 | 0.9000 | 0.7833 |
| `pocket_router_distance_top4_hard` | 60 | 0.4667 | 0.9188 | 0.8833 | 0.7667 |
| `active_set_distance_top4_shell4_w025` | 60 | 0.4500 | 0.9199 | 0.9167 | 0.8500 |

### 5.1 hard20 当前能说明什么

从目前已经完成的三个 arm 看：

1. `random_top4` 的 molecule stability 明显较低：

   ```text
   mol_stable = 0.2833
   ```

2. `distance_top4_hard` 明显比 random sparse 稳定：

   ```text
   mol_stable = 0.4667
   ```

3. `active_set_shell4_w025` 的 molecule stability 接近 hard distance top4：

   ```text
   mol_stable = 0.4500
   ```

4. `active_set_shell4_w025` 的 complete 最好：

   ```text
   complete = 0.8500
   ```

5. `active_set_shell4_w025` 的 reconstruction success 也最好：

   ```text
   recon_success = 0.9167
   ```

当前阶段性解释：

> active-set shell4 w025 没有破坏 ligand 生成质量，反而在 complete 和 reconstruction 上更稳。它的 mol_stable 略低于 hard distance top4，但明显优于 random top4。

还不能下最终结论的原因：

- hard20 还没完成全部 6 个 arms；
- geometry diagnostics 还没跑；
- ligand quality diagnostics 还没跑；
- protein RMSD / clash / contact / boundary 还没有最终表格；
- docking 当前也没跑。

## 6. 已完成 hard8 repeat 结果

hard8 repeat 已经完整完成，并且已有：

```text
results.json
geometry_diagnostics.json
ligand_quality.json
```

所以 hard8 是目前最完整的证据来源。

实验目录：

```text
validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3
```

### 6.1 hard8 主生成指标

| arm | result files | protein RMSD | TM-score | selected residues/update | mol stable | atom stable | recon | complete |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 8 | 1.8153 | 0.9253 | 4.00 | 0.2500 | 0.9094 | 0.8750 | 0.7083 |
| `pocket_router_distance_top4_hard` | 8 | 1.8148 | 0.9245 | 4.00 | 0.5000 | 0.9108 | 0.8333 | 0.6667 |
| `active_set_distance_top4_shell4_w025` | 8 | 1.8207 | 0.9243 | 6.08 | 0.5000 | 0.9173 | 0.8750 | 0.7500 |
| `active_set_distance_top4_shell3_w010` | 8 | 1.8171 | 0.9244 | 4.25 | 0.5000 | 0.9147 | 0.8333 | 0.6667 |
| `active_set_distance_top4_shell3_w025` | 8 | 1.8145 | 0.9245 | 4.26 | 0.4583 | 0.9081 | 0.8333 | 0.6250 |
| `active_set_distance_top4_shell4_w010` | 8 | 1.8168 | 0.9244 | 6.08 | 0.4583 | 0.9094 | 0.8333 | 0.7500 |

这里最重要的点是：

> distance-based arms 的 protein RMSD 非常接近，差距小到不能作为强论文结论。

例如：

```text
distance_top4_hard RMSD = 1.8148
shell4_w025 RMSD        = 1.8207
差距                    = 0.0059 A
```

这个差距太小，不能说谁在 protein RMSD 上碾压谁。

所以真正需要看的，是 ligand 侧和 geometry 侧指标。

### 6.2 hard8 几何诊断

| arm | clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 9.7500 | 1.3000 | 0.6592 | 0.6821 | 0.5024 | 0.6222 |
| `pocket_router_distance_top4_hard` | 1.7500 | 1.7124 | 0.6880 | 0.7131 | 0.5349 | 0.6143 |
| `active_set_distance_top4_shell4_w025` | 1.3750 | 1.7239 | 0.6886 | 0.7557 | 0.5540 | 0.5167 |
| `active_set_distance_top4_shell3_w010` | 1.7500 | 1.7170 | 0.6964 | 0.7155 | 0.5414 | 0.6124 |
| `active_set_distance_top4_shell3_w025` | 1.7500 | 1.7073 | 0.6964 | 0.7210 | 0.5437 | 0.6097 |
| `active_set_distance_top4_shell4_w010` | 1.8750 | 1.6569 | 0.6940 | 0.7346 | 0.5487 | 0.6023 |

这个表是目前最有价值的证据。

关键观察：

1. `random_top4` 的 protein RMSD 看起来不差，但 clash 高达 `9.75`。

   说明 random sparse 会骗过 protein RMSD。

2. `distance_top4_hard` 是很强的 protein baseline，但 clash 和 contact precision 不是最优。

3. `shell4_w025` 的综合几何安全性最好：

   ```text
   clashes            = 1.3750   # 最低
   contact precision  = 0.7557   # 最高
   contact Jaccard    = 0.5540   # 最高
   boundary jump      = 0.5167   # 最低
   ```

大白话解释：

> hard top4 像是只把最关键的几个 residues 强行拉过去。shell4_w025 则像加了一层缓冲垫，让核心周围的 residues 轻轻跟着动，所以 clash 更少、contact 更干净、边界也更平滑。

### 6.3 hard8 ligand quality

| arm | complete | QED | SA | LogP | Lipinski |
|---|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 0.7083 | 0.4441 | 0.5147 | 4.7667 | 4.2941 |
| `pocket_router_distance_top4_hard` | 0.6667 | 0.4856 | 0.5025 | 4.8934 | 4.2500 |
| `active_set_distance_top4_shell4_w025` | 0.7500 | 0.5017 | 0.5106 | 4.8415 | 4.3333 |
| `active_set_distance_top4_shell3_w010` | 0.6667 | 0.4802 | 0.5025 | 5.0589 | 4.1875 |
| `active_set_distance_top4_shell3_w025` | 0.6250 | 0.4973 | 0.5080 | 4.4610 | 4.3333 |
| `active_set_distance_top4_shell4_w010` | 0.7500 | 0.4889 | 0.4956 | 4.7371 | 4.3333 |

这里 `shell4_w025` 仍然比较稳：

- QED 最高；
- complete 并列最高；
- Lipinski 稳定；
- geometry 指标也优于 hard top4。

所以目前最适合作为 fixed active-set baseline 的是：

```text
active_set_distance_top4_shell4_w025
```

## 7. 当前实验结论

### 7.1 已经支持的结论

目前已经比较有把握的结论：

1. sparse pocket update 是有价值的。

2. random sparse 不是有效方法。

   它可能让 protein RMSD 看起来不错，但 clash 很差。

3. hard distance top4 是一个很强的手工 baseline。

   后续 learned router 必须超过或至少明显更安全。

4. shell weak relaxation 是有必要的。

   它能缓解 hard top4 的边界和 clash 风险。

5. `shell4_w025` 是目前最均衡的 fixed active-set candidate。

6. 不能只看 protein RMSD。

   后续论文主表必须同时包括 protein、ligand、ligand-pocket、geometry safety 指标。

### 7.2 还没有证明的东西

目前还不能说已经证明：

1. learned active-set router 一定能超过 distance top4；
2. flow matching 一定比当前 diffusion 生成质量更好；
3. docking 或 binding 指标一定提升；
4. 方法能在更大数据集上稳定泛化；
5. 当前 fixed-rule 方法已经足够投稿。

更准确的判断是：

```text
active-set 机制：有正信号，值得继续
fixed-rule shell：可作为 baseline，不是最终贡献
learnable active-set flow matching：最终需要实现的主方法
```

一句话：

> 这个 idea 作为机制是 work 的，但还不是最终 paper method。

## 8. 最终框架设计

最终方法应该是：

```text
Ligand-Conditioned Active-Set Flow Matching
```

可作为论文题目方向：

```text
Active-Set Flow Matching for Ligand-Induced Apo-to-Holo Pocket Adaptation
```

中文可以叫：

```text
配体条件化主动集流匹配，用于配体诱导的 apo-to-holo 口袋适应
```

### 8.1 最终任务定义

输入：

```text
P_apo: apo pocket
```

输出：

```text
L: generated ligand
P_adapted: adapted holo-like pocket
```

训练目标：

```text
L_holo: 真实 holo ligand
P_holo: 真实 holo pocket
```

也就是：

```text
P_apo -> (L_holo-like, P_holo-like)
```

### 8.2 Flow matching 的基本状态

定义 joint state：

```text
z_t = {L_t, P_t}
```

其中：

```text
L_t = ligand atom coordinates / features at time t
P_t = pocket residue or atom coordinates at time t
```

flow matching 学的是：

```text
v_theta(z_t, t, P_apo)
```

意思是：

> 当前 ligand-pocket 状态应该往哪个方向走，才能变成目标 ligand + holo-like pocket。

### 8.3 Source 和 target

一个最简单的训练设置：

```text
z_0 = {noise ligand, P_apo}
z_1 = {holo ligand, P_holo}
```

插值：

```text
z_t = (1 - t) z_0 + t z_1
```

目标 velocity：

```text
u_t = z_1 - z_0
```

模型预测：

```text
v_theta(z_t, t)
```

基础 flow matching loss：

```text
L_flow = ||v_theta(z_t, t) - u_t||^2
```

但这只是基础。我们的创新不在“用了 flow matching”，而在：

> pocket 的 vector field 由 ligand-conditioned active-set 控制。

### 8.4 Active-set vector field

对每个 pocket residue `i`，模型预测：

```text
a_i_core(t)
a_i_shell(t)
a_i_background(t)
```

这些是连续权重，而不是固定 hard mask。

最终 pocket velocity：

```text
v_P_i =
  a_i_core * v_i_core
  + a_i_shell * v_i_shell
  + a_i_background * v_i_anchor
```

其中：

```text
v_i_core       = 强 pocket adaptation flow
v_i_shell      = 弱局部松弛 flow
v_i_anchor     = 背景锚定 flow，可以是 0 或拉回 apo 的力
```

大白话：

> flow matching 负责告诉系统“怎么走到目标状态”，active set 负责告诉 pocket“哪些地方可以认真走，哪些地方只能小步走，哪些地方不能乱走”。

## 9. 模型模块设计

最终模型建议拆成 5 个模块：

```text
1. ApoPocketEncoder
2. LigandPocketJointEncoder
3. ActiveSetRouter
4. EquivariantFlowDecoder
5. OptionalConfidenceHead
```

### 9.1 ApoPocketEncoder

输入：

```text
residue type
atom coordinates
residue frame
side-chain geometry
local neighborhood
apo pocket shape
```

输出：

```text
h_P_apo
```

作用：

> 给模型一个稳定的 apo pocket memory，告诉模型当前 pocket 的初始结构是什么。

### 9.2 LigandPocketJointEncoder

输入：

```text
L_t
P_t
t
P_apo context
ligand-pocket distance
current contact pattern
```

输出：

```text
h_L(t), h_P(t)
```

作用：

> 建模当前 ligand 状态和 pocket 状态之间的相互作用。

### 9.3 ActiveSetRouter

router 是最终方法的核心。

输入可以包括：

```text
当前 noisy ligand 坐标
当前 ligand partial graph / atom features
residue-ligand distance
residue context
apo pocket geometry
time t
local contact features
```

输出：

```text
a_core_i
a_shell_i
a_background_i
```

注意：

> router 必须是 ligand-conditioned，不能只依赖静态 apo-to-holo displacement。

原因是当前 denoising step 最需要更新的 residues，未必等于最终 holo 中位移最大的 residues。

### 9.4 EquivariantFlowDecoder

输出：

```text
v_L: ligand coordinate velocity
v_P: pocket coordinate / residue-frame velocity
v_atom_type: ligand atom type flow or categorical prediction
```

这个 decoder 最好是 E(3)-equivariant 或至少保持坐标一致性。

### 9.5 OptionalConfidenceHead

可以额外预测：

```text
contact confidence
clash risk
active-set uncertainty
```

这些对论文可解释性很有用：

- 可以画 active residues；
- 可以画 shell relaxation；
- 可以解释为什么某些 case 成功或失败；
- 可以证明模型不是黑箱。

## 10. Loss 设计

总 loss 可以写成：

```text
L =
  L_ligand_flow
  + lambda_pocket * L_pocket_flow
  + lambda_active * L_active_regularization
  + lambda_anchor * L_background_anchor
  + lambda_shell * L_shell_smoothness
  + lambda_contact * L_contact_recovery
  + lambda_clash * L_clash_penalty
```

### 10.1 Ligand flow loss

让 ligand 生成方向正确：

```text
L_ligand_flow = ||v_L - u_L||^2
```

### 10.2 Pocket flow loss

让 pocket adaptation 方向正确：

```text
L_pocket_flow = ||v_P - u_P||^2
```

但这个 loss 不能简单平均所有 residues，否则背景 residues 太多，会稀释真正 active residues。

所以可以用 active-set weight 进行加权：

```text
L_pocket_flow = sum_i w_i ||v_P_i - u_P_i||^2
```

其中：

```text
w_i 可以来自 contact、displacement、active-set probability
```

### 10.3 Active sparsity loss

避免模型把整个 pocket 都激活：

```text
L_active_regularization = mean_i a_i_core
```

或者使用 entropy / top-k style regularization。

目标：

> 让模型只释放必要的 pocket 自由度。

### 10.4 Background anchor loss

保持背景稳定：

```text
L_background_anchor =
  sum_i a_i_background ||P_t_i - P_apo_i||^2
```

大白话：

> 如果一个 residue 被认为是背景，它就不应该被 ligand 生成过程带着乱动。

### 10.5 Shell smoothness loss

避免 core 和 background 之间出现硬切：

```text
L_shell_smoothness =
  sum_(i,j neighbors) shell_edge_weight_ij ||v_P_i - v_P_j||^2
```

目标：

> 让 core 周围有自然过渡，而不是 selected residues 动、旁边 residues 完全不动。

### 10.6 Contact recovery loss

保留 holo-like binding contact：

```text
L_contact_recovery =
  BCE(predicted_contact, holo_contact)
```

目标：

> 生成出来的 ligand-pocket contact 要像真实 holo contact。

### 10.7 Clash penalty

避免 ligand 被 protein 挤坏：

```text
L_clash =
  sum max(0, d_min - d_ligand_protein)^2
```

目标：

> ligand 和 protein 原子不能过近。

## 11. 最终论文 claim 应该怎么写

不要写：

```text
We use flow matching for SBDD.
```

这个太弱，因为 flow matching 本身已经不是新概念。

应该写：

```text
We formulate apo-to-holo SBDD as ligand-conditioned active-set flow matching,
where ligand generation dynamically releases pocket degrees of freedom through
core activation, shell relaxation, and background anchoring.
```

中文：

> 我们把 apo-to-holo SBDD 建模为配体条件化的 active-set flow matching。在生成过程中，当前 ligand 状态动态释放核心 pocket 自由度，弱松弛邻域 residues，并锚定背景结构。

更像论文摘要的说法：

> Existing SBDD methods often condition on static pockets or allow uncontrolled pocket motion. We propose an active-set flow formulation for ligand-induced pocket adaptation. During generation, the current ligand state selects a residue core for strong adaptation, relaxes a local shell for geometric continuity, and anchors the background for structural stability.

中文解释：

> 现有方法要么把 pocket 当静态条件，要么没有明确控制 pocket 运动。我们提出 active-set flow，把 ligand-induced pocket adaptation 拆成核心释放、邻域松弛和背景锚定三个部分，并在 flow matching 框架中联合生成 ligand 和 adapted pocket。

## 12. 最终实验设计

### 12.1 主对比实验

最终 learnable active-set flow matching 至少要对比：

| baseline | 作用 |
|---|---|
| static apo-pocket generation | 证明静态 apo pocket 不够 |
| holo-pocket oracle generation | 上界参考 |
| dense pocket update | 证明全 pocket 更新噪声大 |
| random top4 | 负对照，证明 random sparse 不安全 |
| hard distance top4 | 强手工 sparse baseline |
| fixed active-set shell4 w025 | 当前最强 fixed active-set baseline |
| ligand-only flow matching | 证明 pocket adaptation 有必要 |
| learnable active-set flow matching | 我们最终方法 |

### 12.2 Ablation

| ablation | 要回答的问题 |
|---|---|
| no active-set router | learned routing 是否必要 |
| no shell relaxation | shell 是否缓解边界不连续 |
| no background anchor | 背景锚定是否防止漂移 |
| no contact loss | contact supervision 是否重要 |
| no clash penalty | clash loss 是否保护 ligand |
| hard top-k active set | soft active weights 是否更好 |
| static router | dynamic ligand-conditioned router 是否必要 |

### 12.3 评价指标

Protein / pocket 指标：

```text
pocket RMSD
TM-score
active residue displacement
background drift
boundary jump
```

Ligand 指标：

```text
mol_stable
atom_stable
recon_success
complete
QED
SA
LogP
Lipinski
```

Ligand-pocket 指标：

```text
ligand-protein clash
minimum ligand-protein distance
contact recall
contact precision
contact Jaccard
docking score
```

效率指标：

```text
sampling steps
wall-clock time
number of function evaluations
```

flow matching 的一个潜在优势是：

> 有机会把 1000-step diffusion 采样减少到几十步 ODE / flow steps。

## 13. 当前风险判断

### 13.1 有希望的地方

目前比较积极的信号：

- hard8 已证明 `random_top4` 会产生严重 clash；
- hard8 中 `shell4_w025` 的 clash、contact precision、contact Jaccard、boundary jump 最好；
- hard20 中 `shell4_w025` 的 complete 当前最高；
- active-set 不是单纯追 protein RMSD，而是更符合 ligand-pocket joint quality；
- 当前结果已经能支撑“为什么不能只看 protein RMSD”的论文动机。

### 13.2 风险

主要风险：

1. learned router 可能短期内不容易超过 hard distance top4；
2. flow matching 的坐标处理和 atom type 处理需要设计清楚；
3. apo/holo pocket alignment 必须可靠；
4. docking 工具当前本机还不可用；
5. CPU 实验太慢，不能无限做 sweep；
6. 如果只做 fixed-rule active-set，顶会贡献会不够。

### 13.3 降低风险的策略

建议：

1. hard20 继续跑完，得到完整 fixed-rule evidence；
2. 不再无限增加手工 shell sweep；
3. 尽快实现 minimal learnable active-set router；
4. 先做一个小规模 flow matching prototype；
5. 用 hard8/hard20 作为调试和快速验证集；
6. 只有 prototype 有信号后，再扩大训练和评估。

## 14. 后续工程计划

### Phase 1：完成 hard20 验证

需要做：

```text
等待 hard20 完成 120/120 result files
生成 results.json
跑 geometry diagnostics
跑 ligand quality diagnostics
写 hard20 summary
```

预期输出：

```text
validation/pocket_router/hard20_candidate_repeat_n3.md
```

### Phase 2：实现 minimal learnable active-set router

目标：

```text
把 fixed distance/shell rules 换成可训练 residue router
```

router 输出：

```text
a_core_i
a_shell_i
a_background_i
```

初期可以用 weak labels：

```text
distance-to-ligand
apo-to-holo displacement
contact change
clash/contact diagnostics
```

但最终不能只依赖静态标签，必须是 ligand-conditioned dynamic routing。

### Phase 3：实现 flow matching prototype

第一版可以简化：

```text
continuous ligand coordinates
continuous pocket coordinates
fixed atom count or prior atom count
simple atom type prediction
```

后续再升级：

```text
variable atom count
discrete atom/bond flow
side-chain-aware pocket DOFs
more complete active-set losses
```

### Phase 4：完整评估

需要跑：

```text
main comparison
ablation
geometry diagnostics
ligand quality
docking or docking surrogate
efficiency comparison
case study visualization
active-set visualization
```

### Phase 5：论文写作

论文故事线：

1. apo-to-holo SBDD 需要 pocket adaptation；
2. static pocket 不够真实；
3. dense pocket motion 有噪声；
4. hard sparse update 有边界和 ligand clash 风险；
5. ligand-conditioned active-set 是更合理的抽象；
6. flow matching 给出了连续生成路径；
7. 我们的方法联合生成 ligand 和 adapted pocket；
8. 实验证明它在 ligand-pocket safety、contact、pocket adaptation 和效率上更好。

## 15. 最终一句话框架

英文版：

> We propose ligand-conditioned active-set flow matching for apo-to-holo SBDD, where ligand generation dynamically releases a core set of pocket degrees of freedom, weakly relaxes neighboring residues for geometric continuity, and anchors the background to preserve structural stability.

中文版：

> 我们提出配体条件化 active-set flow matching，用于 apo-to-holo 分子生成。模型在 ligand 生成过程中动态释放核心 pocket 自由度，弱松弛邻近 residues 以保持几何连续性，并锚定背景结构以维持整体稳定。

大白话：

> ligand 不是在死的 pocket 里生成的。ligand 生成过程中会诱导 pocket 变化，但 pocket 不能全局乱动。模型要学会判断哪些 residues 该认真动，哪些邻居该轻轻动，哪些背景必须稳住。flow matching 负责生成路径，active set 负责控制 pocket 自由度释放。

## 16. 当前最终判断

目前最准确的判断是：

```text
active-set 机制：值得继续，已有实验证据支持
fixed-rule shell：可以作为强 baseline，但不是最终贡献
learnable active-set flow matching：必须成为最终投稿方法
```

所以现在不应该继续把主要精力放在手写规则上。

下一步真正关键的里程碑是：

```text
实现并训练第一个 learnable active-set flow matching prototype
```

如果这个 prototype 能超过 hard distance top4 和 fixed shell4 w025，或者在 clash/contact/complete/efficiency 上给出明确优势，这个方向就有进入 AAAI 风格论文主线的可能。
