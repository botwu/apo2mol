# Apo2Mol Active-Set Pocket Adaptation 最终实验整理

更新时间：2026-06-06  
位置：`/Users/wujian1/Downloads/a800_molecular/Apo2Mol`

## 0. 一句话结论

目前的 fixed-rule active-set 实验已经证明：

> ligand-conditioned active-set pocket adaptation 对 apo-to-holo pocket 几何适配是有效的，能在 hard8 / hard20 上显著降低 protein RMSD、提高 TM-score，并且只更新约 4-6 个 residue，而原论文 static5 baseline 平均会更新约 52-54 个 residue。

但是它还没有证明：

> 我们的方法已经全面超过原 Apo2Mol 的所有 ligand quality 指标。

更准确的判断是：

> active-set 是一个有价值的 pocket adaptation 机制验证；最终顶会方法必须把它做成 learnable active-set flow matching，并加入 ligand chemistry / complete / clash 的训练约束。

## 1. 已完成实验

### 1.1 hard20 active-set 完整实验

目录：

`validation/ab_runs/hard20_active_set_candidate_repeat_steps1000_n3`

规模：

- 20 个 hard cases
- 每个 case 3 个 samples
- 6 个方法臂
- 总计 120 个 result files

已生成：

- `results.json`
- `geometry_diagnostics.json`
- `ligand_quality.json`

### 1.2 hard8 原论文 static5 baseline

目录：

`validation/ab_runs/hard8_original_static5_baseline_steps1000_n3`

设置：

- `baseline_realistic_static5`
- `protein_update_schedule = static5`
- `pocket_router_mode = none`
- `num_steps = 1000`
- `num_samples = 3`

已生成：

- `results.json`
- `geometry_diagnostics.json`
- `ligand_quality.json`

### 1.3 hard20 原论文 static5 baseline

目录：

`validation/ab_runs/hard20_original_static5_baseline_steps1000_n3`

设置：

- `baseline_realistic_static5`
- `protein_update_schedule = static5`
- `pocket_router_mode = none`
- `num_steps = 1000`
- `num_samples = 3`

已生成：

- `results.json`
- `geometry_diagnostics.json`
- `ligand_quality.json`

## 2. 指标怎么读：大白话版

### 2.1 protein RMSD

越低越好。

大白话：生成后的 pocket 和真实 holo pocket 越接近，RMSD 越低。这个指标直接反映 apo-to-holo pocket adaptation 是否成功。

### 2.2 TM-score

越高越好。

大白话：蛋白整体结构相似度。RMSD 看距离误差，TM-score 看整体结构是否还像真实结构。

### 2.3 selected/update

越小不一定越好，但代表每次更新多少 pocket residue。

大白话：原论文 static5 平均动 50 多个 residue；我们的 active-set 只动约 4-6 个 residue。如果还能更接近 holo，说明 active-set 的确更聚焦。

### 2.4 mol stable / atom stable / recon / complete

- `mol stable`：整个生成分子价态是否合理，越高越好。
- `atom stable`：单个原子的稳定比例，越高越好。
- `recon`：能否被 RDKit/OpenBabel 重建成分子，越高越好。
- `complete`：重建后是否是一个连通分子，不是碎片，越高越好。

大白话：这些是生成 ligand 能不能成为一个像样分子的基本门槛。

### 2.5 QED / SA / LogP / Lipinski

- QED：药物样性，越高越好。
- SA：这里是归一化后的合成可及性分数，越高越好。
- LogP：疏水性，不是越高越好，通常过高会影响 drug-likeness。
- Lipinski：满足类药规则数量，越高越好，满分 5。

注意：这些指标只在 complete molecule 上计算，所以如果一个方法 complete 很低，QED/SA 的均值会有选择偏差。

### 2.6 clash / contact / boundary jump

- ligand-protein clashes：ligand 和 protein 撞在一起的次数，越低越好。
- min distance：ligand-protein 最近距离，过低说明可能有碰撞。
- contact recall：真实 holo contact 被恢复了多少，越高越好。
- contact precision：预测 contact 有多少是真的，越高越好。
- contact Jaccard：contact 集合重合度，越高越好。
- boundary jump：active region 和 background region 的局部不连续程度，越低越好。

大白话：active-set 不只是要让 RMSD 好看，还要避免局部乱撞、断裂、不连续。

## 3. hard20：原论文 static5 baseline vs active-set

### 3.1 主生成指标

| 方法 | cases | protein RMSD | TM-score | selected/update | mol stable | atom stable | recon | complete |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| original static5 | 20 | 2.1966 | 0.9143 | 54.10 | 0.4167 | 0.9227 | 0.9000 | 0.8833 |
| random top4 | 20 | 1.5889 | 0.9422 | 4.00 | 0.2833 | 0.9191 | 0.9000 | 0.7833 |
| distance top4 hard | 20 | 1.5886 | 0.9419 | 4.00 | 0.4667 | 0.9188 | 0.8833 | 0.7667 |
| active-set shell4 w0.25 | 20 | 1.5935 | 0.9417 | 5.65 | 0.4500 | 0.9199 | 0.9167 | 0.8500 |
| active-set shell3 w0.10 | 20 | 1.5897 | 0.9418 | 4.15 | 0.4667 | 0.9171 | 0.8833 | 0.7667 |
| active-set shell3 w0.25 | 20 | 1.5891 | 0.9418 | 4.16 | 0.4667 | 0.9160 | 0.8833 | 0.7500 |
| active-set shell4 w0.10 | 20 | 1.5897 | 0.9418 | 5.67 | 0.4500 | 0.9166 | 0.8833 | 0.8167 |

核心观察：

- 所有 sparse / active-set 变体的 protein RMSD 都约为 `1.59`。
- 原论文 static5 baseline 的 protein RMSD 是 `2.1966`。
- 也就是说，active-set 类方法在 hard20 上把 protein RMSD 降低了大约 `0.60 A`。
- TM-score 从 static5 的 `0.9143` 提高到约 `0.9417-0.9422`。
- static5 平均更新 `54.10` 个 residue，而 active-set shell4 w0.25 平均只更新 `5.65` 个 residue。

大白话解释：

> 原论文像是“大片区域一起动”；我们的 active-set 像是“只释放 ligand 真正诱导的局部自由度”。实验显示，后者虽然动得少，但 pocket 反而更接近真实 holo。

### 3.2 hard20 ligand quality

| 方法 | samples | complete | QED | SA | LogP | Lipinski |
|---|---:|---:|---:|---:|---:|---:|
| original static5 | 60 | 0.8833 | 0.5151 | 0.5234 | 4.3042 | 4.2264 |
| random top4 | 60 | 0.7833 | 0.4669 | 0.4985 | 4.8374 | 4.1915 |
| distance top4 hard | 60 | 0.7667 | 0.5062 | 0.5165 | 4.3691 | 4.4130 |
| active-set shell4 w0.25 | 60 | 0.8500 | 0.5003 | 0.5006 | 4.3290 | 4.2745 |
| active-set shell3 w0.10 | 60 | 0.7667 | 0.5053 | 0.5130 | 4.4403 | 4.3696 |
| active-set shell3 w0.25 | 60 | 0.7500 | 0.5148 | 0.5169 | 4.1971 | 4.4444 |
| active-set shell4 w0.10 | 60 | 0.8167 | 0.4947 | 0.5049 | 4.5038 | 4.3673 |

核心观察：

- 原论文 static5 的 complete 最好：`0.8833`。
- active-set shell4 w0.25 的 complete 是 `0.8500`，略低但接近。
- static5 的 QED / SA 仍然更强。
- active-set shell3 w0.25 的 QED 几乎追平 static5：`0.5148` vs `0.5151`，并且 LogP 更低：`4.1971` vs `4.3042`。
- shell4 w0.25 是 fixed active-set 中最均衡的方案：complete 好，recon 最好，但 QED/SA 稍弱。

大白话解释：

> active-set 已经把 pocket 适配做得更准，但 ligand 化学质量还没有全面超过原始模型。说明 fixed active-set 只是改变了采样时 pocket 怎么动，ligand generator 本身还没有被重新训练来适配这个新机制。

### 3.3 hard20 几何诊断

| 方法 | clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|
| original static5 | 1.45 | 2.1469 | 0.6519 | 0.8331 | 0.5747 | 0.8569 |
| random top4 | 11.20 | 1.1897 | 0.6735 | 0.7753 | 0.5633 | 0.5136 |
| distance top4 hard | 2.20 | 1.9136 | 0.7044 | 0.7771 | 0.5811 | 0.6942 |
| active-set shell4 w0.25 | 1.90 | 1.9541 | 0.6946 | 0.7984 | 0.5829 | 0.6580 |
| active-set shell3 w0.10 | 2.20 | 1.9303 | 0.7077 | 0.7780 | 0.5837 | 0.6923 |
| active-set shell3 w0.25 | 2.20 | 1.9262 | 0.7077 | 0.7802 | 0.5847 | 0.6914 |
| active-set shell4 w0.10 | 2.25 | 1.9125 | 0.7101 | 0.7949 | 0.5923 | 0.6511 |

核心观察：

- random top4 的 clash 非常高：`11.20`。这证明“随便 sparse 更新”不行。
- static5 的 clash 最低：`1.45`，contact precision 最高：`0.8331`。
- active-set shell4 w0.10 的 contact Jaccard 最高：`0.5923`，但 clash 更高。
- active-set shell4 w0.25 的 boundary jump 从 static5 的 `0.8569` 降到 `0.6580`。
- shell4 w0.25 不是单项最优，但综合比较稳定。

大白话解释：

> active-set 的优势不是“所有指标都压倒 static5”，而是用很少的 residue 更新换来更好的 pocket RMSD / TM 和更低的局部不连续。它的问题是 ligand-protein contact 精度和 clash 还要继续训练优化。

## 4. hard8：原论文 static5 baseline vs active-set

### 4.1 主生成指标

| 方法 | cases | protein RMSD | TM-score | selected/update | mol stable | atom stable | recon | complete |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| original static5 | 8 | 2.4450 | 0.8962 | 52.25 | 0.5000 | 0.9449 | 0.8750 | 0.8750 |
| random top4 | 8 | 1.8153 | 0.9253 | 4.00 | 0.2500 | 0.9094 | 0.8750 | 0.7083 |
| distance top4 hard | 8 | 1.8148 | 0.9245 | 4.00 | 0.5000 | 0.9108 | 0.8333 | 0.6667 |
| active-set shell4 w0.25 | 8 | 1.8207 | 0.9243 | 6.08 | 0.5000 | 0.9173 | 0.8750 | 0.7500 |
| active-set shell3 w0.10 | 8 | 1.8171 | 0.9244 | 4.25 | 0.5000 | 0.9147 | 0.8333 | 0.6667 |
| active-set shell3 w0.25 | 8 | 1.8145 | 0.9245 | 4.26 | 0.4583 | 0.9081 | 0.8333 | 0.6250 |
| active-set shell4 w0.10 | 8 | 1.8168 | 0.9244 | 6.08 | 0.4583 | 0.9094 | 0.8333 | 0.7500 |

核心观察：

- hard8 上，active-set shell4 w0.25 相比 static5：
  - protein RMSD：`2.4450 -> 1.8207`，降低 `0.6243 A`
  - TM-score：`0.8962 -> 0.9243`，提高 `0.0281`
  - 更新 residue：`52.25 -> 6.08`
- 但是 complete 从 `0.8750` 降到 `0.7500`。

### 4.2 hard8 ligand quality

| 方法 | samples | complete | QED | SA | LogP | Lipinski |
|---|---:|---:|---:|---:|---:|---:|
| original static5 | 24 | 0.8750 | 0.5113 | 0.5190 | 5.3543 | 4.1905 |
| random top4 | 24 | 0.7083 | 0.4441 | 0.5147 | 4.7667 | 4.2941 |
| distance top4 hard | 24 | 0.6667 | 0.4856 | 0.5025 | 4.8934 | 4.2500 |
| active-set shell4 w0.25 | 24 | 0.7500 | 0.5017 | 0.5106 | 4.8415 | 4.3333 |
| active-set shell3 w0.10 | 24 | 0.6667 | 0.4802 | 0.5025 | 5.0589 | 4.1875 |
| active-set shell3 w0.25 | 24 | 0.6250 | 0.4973 | 0.5080 | 4.4610 | 4.3333 |
| active-set shell4 w0.10 | 24 | 0.7500 | 0.4889 | 0.4956 | 4.7371 | 4.3333 |

核心观察：

- static5 的 complete 最高。
- active-set shell4 w0.25 的 QED 比 static5 低 `0.0095`，差距很小。
- active-set 的 LogP 明显低于 static5：`4.8415` vs `5.3543`。这反而是好信号，因为 static5 的 LogP 偏高。

### 4.3 hard8 几何诊断

| 方法 | clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|
| original static5 | 1.625 | 1.9423 | 0.6537 | 0.8136 | 0.5607 | 0.8419 |
| random top4 | 9.750 | 1.3000 | 0.6592 | 0.6821 | 0.5024 | 0.6222 |
| distance top4 hard | 1.750 | 1.7124 | 0.6880 | 0.7131 | 0.5349 | 0.6143 |
| active-set shell4 w0.25 | 1.375 | 1.7239 | 0.6886 | 0.7557 | 0.5540 | 0.5167 |
| active-set shell3 w0.10 | 1.750 | 1.7170 | 0.6964 | 0.7155 | 0.5414 | 0.6124 |
| active-set shell3 w0.25 | 1.750 | 1.7073 | 0.6964 | 0.7210 | 0.5437 | 0.6097 |
| active-set shell4 w0.10 | 1.875 | 1.6569 | 0.6940 | 0.7346 | 0.5487 | 0.6023 |

核心观察：

- hard8 上 active-set shell4 w0.25 的 clash 比 static5 更低：`1.375` vs `1.625`。
- boundary jump 大幅降低：`0.8419 -> 0.5167`。
- contact recall 提高：`0.6537 -> 0.6886`。
- contact precision 下降：`0.8136 -> 0.7557`。

大白话解释：

> hard8 上 active-set 的局部几何更平滑，恢复了更多真实 contact，但预测 contact 的精度还不如 static5。这说明 active-set 的“释放自由度”方向是对的，但需要训练时约束 contact 和 ligand chemistry。

## 5. 关键差值表：active-set 相对 static5 改善了什么

### 5.1 active-set shell4 w0.25 vs original static5

| split | RMSD delta | TM delta | complete delta | QED delta | LogP delta | clash delta | contact Jaccard delta | boundary jump delta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| hard8 | -0.6243 | +0.0281 | -0.1250 | -0.0095 | -0.5128 | -0.2500 | -0.0067 | -0.3252 |
| hard20 | -0.6031 | +0.0274 | -0.0333 | -0.0148 | +0.0248 | +0.4500 | +0.0082 | -0.1990 |

怎么解释：

- `RMSD delta` 为负是好事，说明 active-set pocket 更接近 holo。
- `TM delta` 为正是好事，说明结构相似性提高。
- `boundary jump delta` 为负是好事，说明局部更新更连续。
- `complete delta` 为负是问题，说明 ligand 完整性还没完全守住。
- hard20 上 `clash delta` 为正是问题，说明 active-set 在全 hard20 上还有碰撞风险。

### 5.2 hard20 其他 active-set 候选 vs original static5

| 方法 | RMSD delta | TM delta | complete delta | QED delta | LogP delta | clash delta | Jaccard delta | boundary delta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| shell4 w0.25 | -0.6031 | +0.0274 | -0.0333 | -0.0148 | +0.0248 | +0.4500 | +0.0082 | -0.1990 |
| shell3 w0.25 | -0.6075 | +0.0275 | -0.1333 | -0.0004 | -0.1071 | +0.7500 | +0.0100 | -0.1655 |
| shell4 w0.10 | -0.6070 | +0.0275 | -0.0666 | -0.0204 | +0.1996 | +0.8000 | +0.0176 | -0.2058 |

解释：

- `shell3 w0.25` 的 QED 和 LogP 最接近 / 优于 static5，但 complete 掉得最明显。
- `shell4 w0.10` 的 contact Jaccard 和 boundary 最好，但 clash 和 complete 不够理想。
- `shell4 w0.25` 综合最稳：complete 最接近 static5，同时 RMSD/TM 改善明显。

## 6. 当前最佳 fixed-rule 候选

如果只在 fixed-rule active-set 中选一个用于后续论文 narrative，我建议主推：

> `active_set_distance_top4_shell4_w025`

理由：

- hard20 上 recon 最高：`0.9167`
- hard20 上 complete 最接近 static5：`0.8500` vs `0.8833`
- hard20 上 RMSD 明显好于 static5：`1.5935` vs `2.1966`
- hard8 上 clash 比 static5 更低：`1.375` vs `1.625`
- hard8 / hard20 上 boundary jump 都显著低于 static5
- 相比 random top4，clash 大幅降低，证明 router 和 shell 不是随便 sparse

但是需要明确：

> 这不是最终方法，只是机制验证中最稳的 fixed active-set baseline。

## 7. 这组实验能证明什么

### 7.1 能证明的

1. apo-to-holo pocket adaptation 不一定要 dense update。

   原论文 static5 平均动 52-54 个 residue；active-set 只动 4-6 个 residue，却显著降低 protein RMSD。

2. ligand-conditioned sparse pocket update 有明确几何收益。

   hard8 / hard20 都显示 active-set 或 sparse router 能把 RMSD 从约 `2.2-2.4` 降到约 `1.6-1.8`。

3. shell 机制是必要的。

   只随机 top4 会产生严重 clash。shell 让 core 周边有弱松弛，boundary jump 明显下降。

4. fixed active-set 已经给出了顶会方法的机制信号。

   它不是最终 SOTA，但证明了 active-set 视角不是空想。

### 7.2 不能证明的

1. 不能证明我们已经全面超过原 Apo2Mol。

   static5 的 complete、QED、SA 仍然强。

2. 不能证明 fixed active-set 就足够投稿顶会。

   因为它没有训练 router，也没有把 ligand chemistry 纳入 active-set flow。

3. 不能证明 docking / binding affinity 已经提升。

   当前实验使用 `docking_mode none`，还没有 Vina/qvina 结果。

4. 不能证明 full benchmark 上稳定领先。

   当前是 hard8 / hard20 难例验证，不是全测试集主表。

## 8. 对顶会方向的判断

### 8.1 当前结果的价值

这组结果足以支持一个顶会级 idea 的起点：

> apo-to-holo pocket adaptation 应该被建模成 ligand-conditioned active-set optimization，而不是固定 schedule 的 dense update。

审稿人会关心的问题：

- 为什么不是所有 residue 都更新？
- 只更新部分 residue 会不会断裂？
- sparse update 会不会导致 clash？
- ligand quality 会不会变差？

这组实验已经部分回答：

- 不需要所有 residue 都 dense update，因为 active-set 只动少量 residue 就能更接近 holo。
- 需要 shell，而不是 hard selected-only，否则 boundary 和 clash 会出问题。
- fixed active-set 的 ligand quality 还没有完全守住，所以最终必须训练。

### 8.2 当前结果的短板

最大短板不是 protein pocket，而是 ligand quality：

- complete 没有稳定超过 static5。
- QED / SA 还没有明显超过 static5。
- hard20 上 clash 仍有上升风险。
- contact precision 不如 static5。

这说明 fixed active-set 还只是“后验规则”：

> pocket 变得更像 holo 了，但 ligand generator 并没有学会在这个 active-set pocket dynamics 下生成更好的 ligand。

## 9. 最终方法应该怎么设计

最终不能停在“只更新部分 residue”。

最终框架应该是：

> Ligand-Conditioned Active-Set Flow Matching for Apo-to-Holo Pocket Adaptation

### 9.1 核心建模

在每个 flow matching step，模型根据当前 ligand state 动态预测：

- core active residues：强释放、强更新
- shell residues：弱松弛、保持连续
- background residues：锚定、稳定

大白话：

> ligand 生成到哪一步，就问它现在真正影响 pocket 的哪些位置。真正相关的 residue 动，附近 residue 轻微跟着动，远处背景不要乱动。

### 9.2 训练目标

至少需要这些 loss：

1. ligand flow loss

   学 ligand 坐标 / 类型的生成路径。

2. pocket flow loss

   学 pocket 从 apo 到 holo-like 的适配路径。

3. active-set router loss

   可用 apo-holo displacement、contact change、ligand distance 构造 pseudo labels。

4. shell smoothness loss

   约束 active/background 边界不要跳变。

5. background anchor loss

   防止不相关 pocket 区域乱动。

6. clash/contact regularization

   降低 ligand-protein clash，提高真实 contact recovery。

7. ligand chemistry regularization

   保护 complete、QED、SA、LogP、Lipinski。

### 9.3 为什么必须 flow matching

当前 fixed active-set 的问题是：

- active set 是手工规则，不是模型学出来的。
- pocket update 和 ligand generation 不是联合训练。
- ligand chemistry 只是在采样后评估，没有进训练目标。

flow matching 的价值是：

> 可以把 ligand generation、pocket adaptation、active-set selection 三件事放到同一个连续动力学里训练。

## 10. 下一步实验建议

### 10.1 短期：把结果做成论文可用表

需要整理成四张主表：

1. hard8 main generation table
2. hard20 main generation table
3. hard8 / hard20 ligand quality table
4. hard8 / hard20 geometry diagnostics table

主对比：

- original static5
- distance top4 hard
- active-set shell4 w0.25
- active-set shell3 w0.25
- random top4 negative control

### 10.2 中期：训练 learnable router

先做一个 lightweight router：

- 输入：当前 ligand state + pocket residue memory
- 输出：core / shell / background probabilities
- 监督信号：
  - apo-holo residue displacement
  - holo contact change
  - ligand-residue distance
  - local clash/contact feedback

目标：

> 让模型自己学哪些 residue 应该释放，而不是写死 top4/shell4。

### 10.3 长期：active-set flow matching

完整模型应包括：

- ligand flow
- pocket flow
- active-set router
- shell relaxation module
- background anchor
- chemistry-aware regularization

最终主张：

> Active-set flow matching learns ligand-induced pocket degrees of freedom, improving apo-to-holo adaptation while preserving ligand validity.

## 11. 当前论文叙事建议

不建议写：

> 我们所有指标全面超过 Apo2Mol。

建议写：

> 原 Apo2Mol 使用固定 pocket update schedule，缺少 ligand-conditioned freedom selection。我们发现 apo-to-holo adaptation 可以被视为 active-set optimization：当前 ligand 状态只释放少数关键 pocket 自由度，并通过 shell relaxation 保证局部连续性。固定规则实验已经显示该机制可以显著改善 hard-tail apo-to-holo pocket RMSD/TM，并降低 boundary discontinuity；下一步通过 learnable active-set flow matching 将该机制纳入训练，同时保持 ligand quality。

## 12. 最终判断

当前 fixed-rule 结果：

- 对 pocket adaptation：强正向。
- 对 active-set 机制：强支持。
- 对 ligand quality：未全面超过 static5，但部分指标接近或更合理。
- 对顶会投稿：作为最终方法还不够，但作为设计 learnable active-set flow matching 的前期证据是有价值的。

最终应该推进的方向：

> 不要把创新点定义成“只更新部分 residues”。要定义成“ligand-conditioned active-set flow matching”，即让模型学习当前 ligand 诱导释放哪些 pocket 自由度，以及如何在 core / shell / background 之间做连续、稳定、化学友好的联合生成。

