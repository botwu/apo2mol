# Hard8 Active-Set Shell Suite

日期：2026-05-29

运行目录：

`validation/ab_runs/hard8_active_set_shell_only_steps1000_n1`

结构化结果：

- `validation/ab_runs/hard8_active_set_shell_only_steps1000_n1/results.json`
- `validation/ab_runs/hard8_active_set_shell_only_steps1000_n1/geometry_diagnostics.json`

这轮实验是在 hard2 active-set pilot 之后做的 hard8 扩展。它专门回答一个问题：

> 只更新 ligand 附近少数 residues 很强，但会不会因为 unselected residues 不动而产生局部不连续、clash，或者伤害 ligand contact？

大白话：

> 之前我们发现“只修 ligand 旁边最相关的一小块 pocket”很有效。现在要检查的是：只修这一小块，会不会把它旁边没有动的区域挤坏？如果给旁边一圈 residues 一点点松弛自由度，会不会更自然？

---

## 1. 本轮核心想法

我们不再把 pocket update 看成一个简单的二选一：

```text
selected residues: update
unselected residues: freeze
```

而是把它看成 ligand-conditioned active-set optimization：

```text
current ligand state
  -> 选择最相关的 core residues
  -> core residues 强更新
  -> core 周围 shell residues 弱松弛
  -> background residues 保持稳定
```

大白话：

> ligand 当前靠近哪里，哪里就是核心工作区；核心旁边的一圈可以轻轻跟着让位，避免硬切造成不连续；更远处的背景不要乱动，否则整个 pocket 会被带偏。

这比“只更新部分 residues”更像一个论文级主张，因为它不是简单省计算，而是在定义 induced-fit 的自由度释放规则。

---

## 2. 本轮具体改了什么

### 2.1 模型侧已有 active-set weighting

上一轮已经在 `models/molopt_score_model.py` 里加入了 active-set 权重机制：

| pocket 区域 | update weight | 含义 |
|---|---:|---|
| core selected residues | 1.0 | 正常执行模型预测的 translation / rotation / chi update |
| shell residues | 0.25 或 0.50 | 只吃一部分更新，相当于弱松弛 |
| background residues | 0.0 | 基本锚住，不让它每步乱动 |

采样时不是重新训练模型，而是在 protein update step 对模型预测出的 residue 更新做 gating：

- translation 乘 active-set weight；
- side-chain chi update 乘 active-set weight；
- residue rotation 用 slerp 缩回 identity，weight 越小，旋转越接近不动；
- `router_selected_counts` 记录本轮真正 weight > 0 的 residues 数量。

大白话：

> 模型还是会给所有 residues 提修改建议，但我们只完整采纳 core 的建议，shell 的建议只采纳一点点，背景的建议基本不采纳。

### 2.2 实验脚本新增 shell-only 模式

这轮为了只扩展关键 active-set arms，在 `validation/run_new_method_ab.py` 里加了一个更窄的入口：

```text
--active-set-shell-only
```

它只跑 4 个 active-set arms：

| arm | 目的 |
|---|---|
| `active_set_distance_top4_shell4_w025` | distance top4 core，加 4A 弱 shell，测试小 shell 是否安全 |
| `active_set_distance_top4_shell6_w025` | distance top4 core，加 6A 弱 shell，测试更大 shell 是否改善连续性 |
| `active_set_distance_top4_shell6_w050` | shell6 但权重提高到 0.50，测试更强 shell 松弛是否有用 |
| `active_set_random_top4_shell6_w025` | random core + shell，对照“是不是随便多动一点就行” |

大白话：

> 这次不再重复跑所有 dense / oracle arms，而是专门比较“小 shell、大 shell、强一点的大 shell、随机 shell”。

---

## 3. 实验设置

运行命令：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --active-set-shell-only \
  --run-dir validation/ab_runs/hard8_active_set_shell_only_steps1000_n1 \
  --num-cases 8 \
  --num-samples 1 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none
```

固定设置：

| 设置 | 值 |
|---|---|
| cases | 8 hardest apo-to-holo test cases |
| samples per case | 1 |
| denoising steps | 1000 |
| initialization | realistic apo |
| docking | disabled |
| protein update schedule | late-stage update schedule |

hard8 cases：

| test position | original index | metadata apo-holo RMSD |
|---:|---:|---:|
| 477 | 24500 | 4.0188 |
| 327 | 24350 | 3.9616 |
| 390 | 24413 | 3.3143 |
| 310 | 24333 | 2.9210 |
| 377 | 24400 | 2.7470 |
| 342 | 24365 | 2.6896 |
| 365 | 24388 | 2.6712 |
| 347 | 24370 | 2.4535 |

总采样量：

```text
4 arms x 8 cases x 1 sample = 32 result_*.pt
```

---

## 4. 主结果

| arm | mean protein RMSD | mean TM-score | updated residues/update | mol stable | atom stable | recon | complete | evaluated mols |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell4_w025` | 1.8589 | 0.9213 | 6.11 | 0.875 | 0.9675 | 0.875 | 0.750 | 0 |
| `active_set_distance_top4_shell6_w025` | 1.8666 | 0.9217 | 16.38 | 0.375 | 0.9747 | 1.000 | 1.000 | 0 |
| `active_set_distance_top4_shell6_w050` | 1.9301 | 0.9192 | 15.51 | 0.625 | 0.9892 | 1.000 | 0.875 | 0 |
| `active_set_random_top4_shell6_w025` | 1.8524 | 0.9249 | 15.75 | 0.250 | 0.8556 | 0.875 | 0.875 | 0 |

对照 hard8 full router suite 里的关键 baseline：

| arm | mean protein RMSD | mean TM-score | updated residues/update | mol stable | complete |
|---|---:|---:|---:|---:|---:|
| `baseline_realistic_static5` | 2.4365 | 0.8977 | 52.25 | 0.250 | 0.875 |
| `control_realistic_late_dense` | 2.4142 | 0.9005 | 52.25 | 0.250 | 0.875 |
| `pocket_router_distance_top4` | 1.7924 | 0.9249 | 4.00 | 0.375 | 0.750 |
| `pocket_router_random_top4` | 1.8305 | 0.9248 | 4.00 | 0.375 | 0.625 |
| `active_set_distance_top4_shell4_w025` | 1.8589 | 0.9213 | 6.11 | 0.875 | 0.750 |
| `active_set_distance_top4_shell6_w025` | 1.8666 | 0.9217 | 16.38 | 0.375 | 1.000 |
| `active_set_random_top4_shell6_w025` | 1.8524 | 0.9249 | 15.75 | 0.250 | 0.875 |

大白话读法：

> active-set shell 仍然远好于 dense pocket update，但没有超过 hard `distance_top4`。小 shell 最稳，大 shell 不一定更好，随机 shell 在 protein RMSD 上看着不差，但 ligand 几何和 clash 很糟。

---

## 5. 几何诊断

为了专门回答“局部不连续 / clash / contact”问题，我补了诊断脚本：

`validation/analyze_geometry_diagnostics.py`

它计算：

| 指标 | 意思 | 大白话 |
|---|---|---|
| `ligand_protein_clashes` | ligand 和 protein 原子之间过近的 clash 数 | ligand 有没有被 pocket 挤到 |
| `ligand_protein_min_dist` | ligand-protein 最近原子距离 | 最近有没有贴得太近 |
| `contact_recall` | holo 参考 contact 有多少被保留 | 真正该接触的位置有没有找回来 |
| `contact_precision` | 预测 contact 中有多少是参考 contact | 生成出来的 contact 有没有乱碰 |
| `contact_jaccard` | contact set overlap | contact 整体重合度 |
| `neighbor_displacement_jump_mean` | active/background 边界位移差 | 动的区域和不动区域之间是否有断层 |
| `active_displacement_mean` | active residues 平均移动幅度 | 被释放的区域动了多少 |
| `background_displacement_mean` | background 平均移动幅度 | 背景有没有被带着乱动 |

hard8 active-set shell 诊断：

| arm | RMSD | ligand-protein clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell4_w025` | 1.8589 | 0.25 | 2.3097 | 0.7228 | 0.7185 | 0.5501 | 0.7971 |
| `active_set_distance_top4_shell6_w025` | 1.8666 | 0.25 | 2.3752 | 0.7430 | 0.7237 | 0.5777 | 0.4775 |
| `active_set_distance_top4_shell6_w050` | 1.9301 | 0.375 | 2.3211 | 0.7495 | 0.7501 | 0.5957 | 0.5336 |
| `active_set_random_top4_shell6_w025` | 1.8524 | 9.50 | 1.3783 | 0.7121 | 0.7423 | 0.5742 | 0.3965 |

关键对照：

| arm | RMSD | ligand-protein clashes | min dist | contact recall | boundary jump |
|---|---:|---:|---:|---:|---:|
| `baseline_realistic_static5` | 2.4365 | 3.375 | 1.8340 | 0.6434 | 0.9078 |
| `control_realistic_late_dense` | 2.4142 | 2.000 | 1.7417 | 0.6579 | 0.7434 |
| `pocket_router_distance_top4` | 1.7924 | 1.875 | 1.8062 | 0.7215 | 1.0044 |
| `pocket_router_random_top4` | 1.8305 | 12.000 | 1.2323 | 0.7492 | 0.5505 |
| `active_set_distance_top4_shell4_w025` | 1.8589 | 0.250 | 2.3097 | 0.7228 | 0.7971 |
| `active_set_distance_top4_shell6_w025` | 1.8666 | 0.250 | 2.3752 | 0.7430 | 0.4775 |
| `active_set_random_top4_shell6_w025` | 1.8524 | 9.500 | 1.3783 | 0.7121 | 0.3965 |

---

## 6. 怎么解释这些结果

### 6.1 小 shell 是安全的，但不是 RMSD 最强

`active_set_distance_top4_shell4_w025`：

- protein RMSD 1.8589，比 dense baselines 好很多；
- mol stable 0.875，是本轮 active-set 里最好；
- ligand-protein clashes 0.25，明显低于 hard `distance_top4` 的 1.875；
- boundary jump 0.7971，比 hard `distance_top4` 的 1.0044 更平滑。

大白话：

> 小 shell 像一个安全缓冲垫。它让 ligand-protein clash 少很多，也让动/不动的边界没那么硬，但代价是 protein RMSD 比纯 hard top4 稍差一点。

### 6.2 shell6 让边界更连续，但会伤害 ligand 稳定性

`active_set_distance_top4_shell6_w025`：

- boundary jump 从 shell4 的 0.7971 降到 0.4775；
- contact recall 从 0.7228 提到 0.7430；
- complete 从 0.750 提到 1.000；
- 但 mol stable 从 0.875 掉到 0.375。

大白话：

> 6A shell 确实让结构过渡更平滑，contact 也稍好，但 ligand 自身稳定性坏了。也就是说，“更连续”不等于“更好的生成结果”。

### 6.3 shell weight 0.50 不是更好

`active_set_distance_top4_shell6_w050`：

- contact 指标略升；
- RMSD 变差到 1.9301；
- clash 从 0.25 升到 0.375；
- complete 低于 shell6 w0.25。

大白话：

> shell 不是越强越好。邻居只是应该轻轻松一下，不应该跟 core 一样大幅参与更新。

### 6.4 random shell 是危险对照

`active_set_random_top4_shell6_w025` 的 protein RMSD 是 1.8524，看起来甚至略好于 distance shell4。但几何诊断揭示了问题：

- ligand-protein clashes = 9.50；
- min dist = 1.3783 A；
- mol stable = 0.250；
- atom stable = 0.8556。

大白话：

> 单看 protein RMSD 会被误导。random shell 可能把 protein 调得像 holo，但 ligand 被挤坏了。真正的 active set 不能随机选，必须 ligand-conditioned。

这正好回答了用户提出的担心：

> 局部 pocket update 可能改善 protein RMSD，但损害 ligand geometry 或 binding contact。

是的，random shell 就是这个负例。

### 6.5 只动部分 residues 会不会造成其他 residues 冲突？

从本轮看，答案更细：

- hard `distance_top4` 有较高 boundary jump 和一定 clash；
- distance shell4 / shell6 显著降低 ligand-protein clash；
- shell6 显著降低 boundary jump；
- random shell 虽然 boundary jump 更低，但 clash 爆炸。

大白话：

> 只动 core 确实可能有边界问题；加 shell 能缓解。但 shell 必须围绕正确的 core 加，不能随机加。

### 6.6 只看 ligand 附近 residues 会不会忽略二级影响？

本轮没有证明二级影响不重要。它证明的是：

> 简单扩大 shell 半径不是解决二级影响的正确方式。

更合理的做法应该是：

```text
core: 当前 fragment 直接相关 residues
inner shell: 小半径弱松弛
long-range residues: learned router 判断是否需要释放
background: 大多数时候锚住，只做低频弱 relaxation
```

大白话：

> 远处影响可能重要，但不能把一大片 pocket 都放出来一起动。应该让模型学会哪条传导链真的重要。

---

## 7. 当前最稳结论

这轮 hard8 结果支持：

> ligand-conditioned active-set optimization 是合理方向，但默认应该是小 shell、弱松弛、强背景锚定。

更具体：

- dense update 仍然不是好选择；
- hard `distance_top4` 仍是当前 protein RMSD 最强 hand-crafted baseline；
- distance shell4 是更安全、更物理的折中：clash 少、稳定性好、边界较平滑；
- distance shell6 说明“更连续”会带来 ligand 稳定性代价；
- random shell 是强负例：protein RMSD 不差，但 ligand-protein clash 很差；
- 单看 protein RMSD 不够，必须同时看 ligand stability、clash、contact 和 boundary。

一句大白话：

> 论文不能只讲“让更少 residues 动”。真正要讲的是：当前 ligand 决定哪些自由度释放，旁边哪些自由度轻轻松动，哪些背景必须稳住；并且这个选择错了会直接把 ligand 挤坏。

---

## 8. 对顶会方向的影响

这轮把 idea 又收紧了一步。

不应该把主贡献写成：

> sparse residue update

更合理的主贡献是：

> fragment-conditioned active-set pocket adaptation

核心 claim 可以是：

> During ligand generation, each partial ligand state induces an active set of pocket degrees of freedom: a contact-relevant residue core is updated strongly, its local shell is weakly relaxed to preserve geometric continuity, and the background is anchored to prevent denoising noise.

中文大意：

> ligand 生成到每一步时，不是让整个 pocket 一起动，而是释放当前 ligand 真正需要的一组 pocket 自由度。核心区域强更新，邻近区域弱松弛，背景区域稳定住。

这个表述比 hard mask 更完整，因为它同时解释：

- 为什么 dense update 有噪声；
- 为什么 top4 很强；
- 为什么只 hard freeze 可能有边界问题；
- 为什么 shell 不能随机加；
- 为什么只看 protein RMSD 不够。

---

## 9. 下一步实验

### 9.1 训练前还应补的算法实验

1. Multi-sample repeat

```text
hard8, n=3 or n=5 samples per case
```

目的：

> 当前 n=1，单个 sample 偶然性还比较大。需要看均值、方差和胜率。

2. Shell hyperparameter sweep

建议只扫保守范围：

```text
radius: 3A, 4A, 5A
weight: 0.10, 0.25
```

不建议继续大力扫 6A / 0.50，因为当前已经显示风险。

3. Core top-k + shell 联合 sweep

```text
core topk: 3, 4, 6
shell radius: 3A, 4A
shell weight: 0.10, 0.25
```

目的：

> 现在只验证了 top4 core。训练前要知道 core budget 和 shell budget 的合理组合。

### 9.2 ligand 侧必须补

必须补：

- reconstruction success；
- mol stable / atom stable；
- ligand-protein clash；
- contact recall / precision / Jaccard；
- docking score；
- QED / SA；
- PoseBusters-style geometry checks。

本轮 `evaluated_mols = 0`，因为 `--docking-mode none`，所以 QED / SA / docking 不能作为结论。

大白话：

> pocket 修得像 holo 不代表药物分子真的好。random shell 已经证明，protein RMSD 可以好看，但 ligand 被挤坏。

### 9.3 训练版建议

训练版不建议直接把 shell 当固定规则写死。更适合：

```text
core router:
  predicts ligand-contact relevance and update necessity

shell router:
  predicts weak relaxation weights for neighboring residues

background stabilizer:
  penalizes unnecessary far-field motion

losses:
  ligand diffusion/reconstruction loss
  pocket displacement/rotamer loss
  contact recovery loss
  clash penalty
  sparsity budget loss
  background stability loss
```

大白话：

> 让模型学会“哪里必须动、哪里轻轻动、哪里不能动”，而不是只学“离 ligand 最近”。

---

## 10. 一句话结论

> Hard8 active-set shell suite 说明：小范围弱 shell 可以缓解 hard top4 的局部连续性和 clash 风险，但大 shell 或 random shell 会损害 ligand 稳定性；顶会方向应从 hard sparse update 升级为 fragment-conditioned active-set optimization，并且后续必须用 ligand/clash/contact 指标证明这个 active set 选得对。

---

## 11. 后续 conservative shell sweep 已完成

完整记录：

`validation/pocket_router/hard8_conservative_shell_sweep.md`

运行目录：

`validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1`

该 follow-up 跑完了 4 个更保守的 distance top4 active-set arms：

| arm | RMSD | mol stable | complete | clashes | contact recall |
|---|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell3_w010` | 1.8571 | 0.750 | 0.875 | 0.250 | 0.7342 |
| `active_set_distance_top4_shell3_w025` | 1.8584 | 0.750 | 0.875 | 0.250 | 0.7342 |
| `active_set_distance_top4_shell4_w010` | 1.8512 | 0.750 | 0.750 | 0.250 | 0.7366 |
| `active_set_distance_top4_shell5_w025` | 1.8563 | 0.625 | 0.750 | 0.375 | 0.7537 |

更新后的理解：

> `shell3_w010/w025` 是更均衡的保守方案；`shell5_w025` contact recall 更高，但 ligand 稳定性和 clash 变差。固定 shell 参数没有全赢，下一步应该训练 learned active-set policy，而不是继续手写更大的半径。

---

## 12. Hard8 n=3 candidate repeat 已完成

完整记录：

`validation/pocket_router/hard8_candidate_repeat_n3.md`

运行目录：

`validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3`

本轮把 6 个关键候选 arm 每个 hard case 跑 3 个 sample：

```text
6 arms x 8 cases x 3 samples = 144 generated samples
48 / 48 result_*.pt completed
```

主结果：

| arm | protein RMSD | TM-score | mol stable | complete | clashes | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 1.8153 | 0.9253 | 0.2500 | 0.7083 | 9.7500 | 0.6821 | 0.5024 | 0.6222 |
| `pocket_router_distance_top4_hard` | 1.8148 | 0.9245 | 0.5000 | 0.6667 | 1.7500 | 0.7131 | 0.5349 | 0.6143 |
| `active_set_distance_top4_shell4_w025` | 1.8207 | 0.9243 | 0.5000 | 0.7500 | 1.3750 | 0.7557 | 0.5540 | 0.5167 |
| `active_set_distance_top4_shell3_w010` | 1.8171 | 0.9244 | 0.5000 | 0.6667 | 1.7500 | 0.7155 | 0.5414 | 0.6124 |
| `active_set_distance_top4_shell3_w025` | 1.8145 | 0.9245 | 0.4583 | 0.6250 | 1.7500 | 0.7210 | 0.5437 | 0.6097 |
| `active_set_distance_top4_shell4_w010` | 1.8168 | 0.9244 | 0.4583 | 0.7500 | 1.8750 | 0.7346 | 0.5487 | 0.6023 |

更新后的理解：

> n=3 repeat 后，几个 distance-based candidates 的 protein RMSD 差距非常小，不能只凭 RMSD 选方案。`random_top4` 再次证明 sparse 本身会误导：RMSD 接近 distance top4，但 clash 高到 9.75。`shell4_w025` 不是 RMSD 最低，但 clash 最低、contact precision/Jaccard 最高、boundary jump 最低，是目前最像 active-set optimization 的安全候选。
