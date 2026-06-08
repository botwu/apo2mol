# Hard8 Active-Set Candidate Repeat n=3

日期：2026-05-31

运行目录：

`validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3`

结构化结果：

- `validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3/results.json`
- `validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3/geometry_diagnostics.json`

这轮实验接在 hard8 conservative shell sweep 后面做，目的不是继续扩大搜索空间，而是把最有价值的几个 candidate 用 n=3 重复采样重新比一遍，降低单次采样偶然性。

大白话：

> 前面 n=1 的结果已经告诉我们，小 shell 比大 shell 更安全。但 n=1 可能有运气成分，所以这轮把几个候选方案每个 case 跑 3 个 sample，看看结论还稳不稳。

---

## 1. 这轮要回答什么

用户提出的核心担心是：

- selected residues 动了，unselected residues 不动，会不会造成局部结构不连续；
- 只动 ligand 附近 residues，会不会忽略二级影响；
- 局部 pocket 更新可能改善 protein RMSD，但损害 ligand geometry 或 binding contact；
- 只看 protein RMSD 可能把 ligand 挤坏。

这轮实验专门围绕这些问题做 repeat。

我们比较 6 个 arms：

| arm | 作用 |
|---|---|
| `pocket_router_random_top4` | 随机选 4 个 residues 更新，对照“稀疏本身是否就有效” |
| `pocket_router_distance_top4_hard` | 只更新距离当前 ligand 最近的 top4 residues |
| `active_set_distance_top4_shell4_w025` | top4 core 强更新，4A shell 以 0.25 弱松弛 |
| `active_set_distance_top4_shell3_w010` | top4 core 强更新，3A shell 以 0.10 弱松弛 |
| `active_set_distance_top4_shell3_w025` | top4 core 强更新，3A shell 以 0.25 弱松弛 |
| `active_set_distance_top4_shell4_w010` | top4 core 强更新，4A shell 以 0.10 弱松弛 |

大白话：

> 我们不再问“要不要 sparse”，而是问：只动核心够不够？核心旁边要不要轻轻让位？shell 是 3A 还是 4A？权重要 0.10 还是 0.25？

---

## 2. 代码改动

本轮只改实验编排脚本：

`validation/run_new_method_ab.py`

新增函数：

```python
make_active_set_candidate_repeat_arms()
```

新增命令行参数：

```text
--active-set-candidate-repeat
```

这个入口固定生成上述 6 个 candidate arms，并保持下面这些设置一致：

| 设置 | 值 |
|---|---|
| `protein_update_schedule` | `late_dense` |
| `protein_update_interval` | 50 |
| `protein_update_min_t` | 10 |
| `protein_update_residual_threshold` | 0.5 |
| `pocket_router_topk` | 4 |
| `pocket_router_background_weight` | 0.0 |
| `sample_num_atoms` | `prior` |
| `init_center_mode` | `apo` |
| `num_steps` | 1000 |

大白话：

> 所有方法都用同样的采样步数、同样的 protein update 时机、同样的 top4 core。唯一变化是：core 怎么选、core 外面要不要加弱 shell、shell 多大、shell 多强。

---

## 3. 运行命令

准备 / preflight：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --active-set-candidate-repeat \
  --run-dir validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3 \
  --num-cases 8 \
  --num-samples 3 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none
```

正式运行：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --active-set-candidate-repeat \
  --run-dir validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3 \
  --num-cases 8 \
  --num-samples 3 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none
```

总采样量：

```text
6 arms x 8 hard cases x 3 samples = 144 generated samples
```

结果文件数量：

```text
6 arms x 8 cases = 48 result_*.pt
```

最终确认：

```text
48 / 48 result files completed
```

注意：

> 每个 `result_*.pt` 对应一个 hard case，里面包含 3 个 generated samples。所以 `result_files=8`，但 sampling 统计里的 `n=24`。

---

## 4. Hard8 cases

| test position | original index | metadata apo-holo RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |
| 390 | 24413 | 3.3143 | `8u6b__1__1.A__1.C/1.C.sdf` |
| 310 | 24333 | 2.9210 | `8ow3__1__1.A__1.B/1.B.sdf` |
| 377 | 24400 | 2.7470 | `8sfu__2__1.B__1.E/1.E.sdf` |
| 342 | 24365 | 2.6896 | `8qn5__3__1.C__1.L/1.L.sdf` |
| 365 | 24388 | 2.6712 | `8pqh__1__1.A__1.B/1.B.sdf` |
| 347 | 24370 | 2.4535 | `8sbv__2__1.B__1.G/1.G.sdf` |

---

## 5. 主指标结果

| arm | result files | protein RMSD mean | RMSD median | RMSD min | RMSD max | TM-score mean | updated residues/update | mol stable | atom stable | recon | complete | evaluated mols |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 8 | 1.8153 | 1.8123 | 0.9891 | 2.9326 | 0.9253 | 4.00 | 0.2500 | 0.9094 | 0.8750 | 0.7083 | 0 |
| `pocket_router_distance_top4_hard` | 8 | 1.8148 | 1.7952 | 0.9509 | 3.0098 | 0.9245 | 4.00 | 0.5000 | 0.9108 | 0.8333 | 0.6667 | 0 |
| `active_set_distance_top4_shell4_w025` | 8 | 1.8207 | 1.7991 | 0.9529 | 2.9943 | 0.9243 | 6.08 | 0.5000 | 0.9173 | 0.8750 | 0.7500 | 0 |
| `active_set_distance_top4_shell3_w010` | 8 | 1.8171 | 1.8001 | 0.9509 | 3.0122 | 0.9244 | 4.25 | 0.5000 | 0.9147 | 0.8333 | 0.6667 | 0 |
| `active_set_distance_top4_shell3_w025` | 8 | 1.8145 | 1.8034 | 0.9509 | 2.9964 | 0.9245 | 4.26 | 0.4583 | 0.9081 | 0.8333 | 0.6250 | 0 |
| `active_set_distance_top4_shell4_w010` | 8 | 1.8168 | 1.7945 | 0.9518 | 2.9891 | 0.9244 | 6.08 | 0.4583 | 0.9094 | 0.8333 | 0.7500 | 0 |

大白话读法：

> 单看 protein RMSD，几个 distance-based arms 几乎贴在一起，差距只有 0.006 A 左右。这个级别不能作为强论文结论。真正有区分度的是 ligand 稳定性、clash、contact 和边界连续性。

---

## 6. 几何诊断结果

诊断命令：

```bash
.venv310/bin/python validation/analyze_geometry_diagnostics.py \
  --run-dir validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3
```

结果：

| arm | geometry RMSD | ligand-protein clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump | active disp | background disp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell3_w010` | 1.7752 | 1.7500 | 1.7170 | 0.6964 | 0.7155 | 0.5414 | 0.6124 | 0.5902 | 0.2871 |
| `active_set_distance_top4_shell3_w025` | 1.7760 | 1.7500 | 1.7073 | 0.6964 | 0.7210 | 0.5437 | 0.6097 | 0.5924 | 0.2887 |
| `active_set_distance_top4_shell4_w010` | 1.7756 | 1.8750 | 1.6569 | 0.6940 | 0.7346 | 0.5487 | 0.6023 | 0.5352 | 0.2776 |
| `active_set_distance_top4_shell4_w025` | 1.7777 | 1.3750 | 1.7239 | 0.6886 | 0.7557 | 0.5540 | 0.5167 | 0.4149 | 0.3091 |
| `pocket_router_distance_top4_hard` | 1.7749 | 1.7500 | 1.7124 | 0.6880 | 0.7131 | 0.5349 | 0.6143 | 0.5950 | 0.2868 |
| `pocket_router_random_top4` | 1.8668 | 9.7500 | 1.3000 | 0.6592 | 0.6821 | 0.5024 | 0.6222 | 0.5608 | 0.4201 |

指标解释：

| 指标 | 大白话 |
|---|---|
| `ligand-protein clashes` | ligand 和 protein 有没有挤到一起，越低越好 |
| `min dist` | ligand-protein 最近原子距离，太小通常危险 |
| `contact recall` | holo 参考里该有的接触找回了多少 |
| `contact precision` | 预测 contact 里有多少不是乱碰 |
| `contact Jaccard` | contact 整体重合度 |
| `boundary jump` | active 区域和背景之间有没有明显断层，越低越平滑 |
| `active disp` | 被释放区域平均动了多少 |
| `background disp` | 背景区域平均动了多少 |

---

## 7. Ligand-quality 补充评估

由于本机没有 `qvina2`，Python 环境里也没有 `vina` 包，本轮暂时不能跑真正 docking。为避免 `eval_split.py --docking-mode none` 下 QED/SA 为空，我补了一个不依赖 docking 的轻量评估脚本：

`validation/analyze_ligand_quality.py`

运行命令：

```bash
.venv310/bin/python validation/analyze_ligand_quality.py \
  --run-dir validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3
```

输出：

`validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3/ligand_quality.json`

结果：

| arm | generated | mol stable | recon | complete | chem success | QED mean | SA mean | LogP mean | Lipinski mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 24 | 0.2500 | 0.8750 | 0.7083 | 0.7083 | 0.4441 | 0.5147 | 4.7667 | 4.2941 |
| `pocket_router_distance_top4_hard` | 24 | 0.5000 | 0.8333 | 0.6667 | 0.6667 | 0.4856 | 0.5025 | 4.8934 | 4.2500 |
| `active_set_distance_top4_shell4_w025` | 24 | 0.5000 | 0.8750 | 0.7500 | 0.7500 | 0.5017 | 0.5106 | 4.8415 | 4.3333 |
| `active_set_distance_top4_shell3_w010` | 24 | 0.5000 | 0.8333 | 0.6667 | 0.6667 | 0.4802 | 0.5025 | 5.0589 | 4.1875 |
| `active_set_distance_top4_shell3_w025` | 24 | 0.4583 | 0.8333 | 0.6250 | 0.6250 | 0.4973 | 0.5080 | 4.4610 | 4.3333 |
| `active_set_distance_top4_shell4_w010` | 24 | 0.4583 | 0.8333 | 0.7500 | 0.7500 | 0.4889 | 0.4956 | 4.7371 | 4.3333 |

大白话：

> 在不跑 docking 的化学质量指标上，`shell4_w025` 仍然比较稳：QED 最高，complete/chem success 也是并列最高。`shell4_w010` 的 SA 最低，但前面的几何诊断显示它 clash 更差。`random_top4` 的 QED 最低、mol stable 最低、clash 最高，进一步说明随机 sparse 不是可发表的方法。

注意：

> 这还不是 docking 结论。它只是补上 QED/SA/LogP/Lipinski 这些 RDKit 化学指标。真正 binding-quality 还需要安装/配置 vina 或 qvina 后再跑。

---

## 8. 关键结论

### 8.1 Random sparse 不能作为方法

`pocket_router_random_top4` 的主 RMSD 看起来不差：

```text
protein RMSD mean = 1.8153
```

但几何诊断很差：

```text
ligand-protein clashes = 9.75
min dist = 1.30 A
contact recall = 0.6592
background disp = 0.4201
```

大白话：

> 随机 sparse 可以把 protein RMSD 做得看起来还行，但 ligand 被 pocket 挤坏了。这说明只汇报 protein RMSD 会误导；顶会论文必须汇报 ligand clash/contact/validity。

### 8.2 Hard distance top4 仍是强 baseline，但不是最安全

`pocket_router_distance_top4_hard` 的 RMSD 很好：

```text
protein RMSD mean = 1.8148
geometry RMSD = 1.7749
```

但它的 ligand 几何不是最安全：

```text
ligand-protein clashes = 1.75
contact precision = 0.7131
contact Jaccard = 0.5349
boundary jump = 0.6143
complete = 0.6667
```

大白话：

> 只动 top4 core 能很好地修 protein，但它像硬切：核心动了，旁边不一定跟得上，所以 ligand 和局部 contact 还有风险。

### 8.3 Shell4 w0.25 是这轮最均衡的安全 candidate

`active_set_distance_top4_shell4_w025` 不是 protein RMSD 最低：

```text
protein RMSD mean = 1.8207
```

但它在几何指标上最均衡：

```text
ligand-protein clashes = 1.3750  # 本轮最低
min dist = 1.7239
contact precision = 0.7557       # 本轮最高
contact Jaccard = 0.5540         # 本轮最高
boundary jump = 0.5167           # 本轮最低
complete = 0.7500
```

大白话：

> shell4_w025 像一个缓冲层。它牺牲一点点 protein RMSD，但换来更少 clash、更平滑边界、更靠谱的 contact precision。这更接近我们想要的 active-set optimization。

### 8.4 Shell3 variants 更像保守 baseline，不是最终答案

`shell3_w010/w025` 的 contact recall 比 hard top4 好：

```text
hard top4 contact recall = 0.6880
shell3 contact recall = 0.6964
```

但 clash 没有降下来：

```text
hard top4 clashes = 1.75
shell3 clashes = 1.75
```

大白话：

> 3A shell 太保守，能稍微补一点 contact，但对缓解 ligand-protein clash 不够明显。

### 8.5 Shell4 w0.10 说明“更平滑”不等于“更安全”

`shell4_w010` 的 boundary jump 较低：

```text
boundary jump = 0.6023
```

但 ligand-protein clash 最高于 distance active-set candidates：

```text
clashes = 1.8750
min dist = 1.6569
```

大白话：

> 让边界平滑本身不是最终目标。真正目标是：core/shell/background 的自由度释放要同时服务 ligand geometry 和 binding contact。

---

## 9. 对第一性原理问题的回答

### selected residues 动了，unselected residues 不动，会不会冲突？

会，有风险。

证据：

- hard top4 的 clash 是 1.75；
- random top4 的 clash 是 9.75；
- shell4_w025 把 clash 降到 1.375。

结论：

> 只做 hard mask 不够。core 外面需要可控的弱松弛，尤其要用 clash/contact 指标约束，而不是只看 protein RMSD。

### 会不会造成局部结构不连续？

会有边界风险，但 shell 能缓解。

证据：

- hard top4 boundary jump = 0.6143；
- shell4_w025 boundary jump = 0.5167；
- shell3_w025 boundary jump = 0.6097。

结论：

> shell 半径和权重要足够覆盖核心邻域。3A 太保守，4A/0.25 更像有效缓冲层。

### 会不会改善 protein RMSD 但损害 ligand geometry？

会，random top4 是最清楚的反例。

证据：

- random top4 protein RMSD = 1.8153，和 distance top4 很接近；
- 但 random top4 clash = 9.75，min dist = 1.30 A。

结论：

> protein RMSD 不能单独作为主指标。后续论文实验必须把 ligand geometry、clash、contact、docking 同时放进主表。

### 只看 ligand 附近 residues，会不会忽略二级影响？

会，所以 active-set 不能只是 top4 core。

证据：

- hard top4 RMSD 很强，但 contact precision/Jaccard 不如 shell4_w025；
- shell4_w025 的 boundary jump 和 clash 更好。

结论：

> 第一性原理上更合理的形式不是“只更新部分 residues”，而是“core 强更新 + neighbor weak relaxation + background anchoring”。

---

## 10. 当前最合理的 paper idea

这轮实验后，最准确的主张应该是：

> Apo-to-holo pocket adaptation should be formulated as ligand-conditioned active-set optimization: the current ligand state releases a sparse core set of pocket degrees of freedom, softly relaxes neighboring residues, and anchors the far-field background.

大白话：

> ligand 长到哪里，哪里附近的 pocket 核心就强动；核心旁边轻轻让位；远处不要乱动。

这比“只更新 top-k residues”更好，因为它明确处理了三个问题：

- 哪些 residues 必须释放自由度；
- 哪些邻居需要弱松弛来避免硬切；
- 哪些背景必须保持稳定，避免全局漂移。

---

## 11. 下一步实验建议

### 10.1 训练前最应该补的实验

1. 对少数 finalist arms 开 docking/eval：

```text
distance_top4_hard
active_set_shell4_w025
active_set_shell3_w025
active_set_shell4_w010
```

目的：

> 现在 `evaluated_mols = 0`，QED/SA/docking 都不能下结论。必须补 ligand-side evaluation。

2. 做 paired case-level 统计：

```text
hard8, n=3 已完成
下一步：hard8 n=5 或 hard20 n=3
```

目的：

> 现在各 distance-based arms 的 protein RMSD 太接近，需要看 case-level 胜率和方差。

3. 小范围 core/shell 联合 sweep：

```text
core topk: 3, 4, 6
shell radius: 3A, 4A
shell weight: 0.10, 0.25
```

目的：

> 确认 shell4_w025 是稳定规律，还是只在 top4 core 下成立。

### 10.2 训练版应该怎么做

训练版不建议固定写死 shell4_w025，而应该让模型学习：

```text
core release score:
  哪些 residues 当前必须强更新

shell relaxation score:
  哪些邻居需要弱松弛，松弛多强

background anchoring score:
  哪些区域必须保持稳定
```

推荐 loss：

| loss | 作用 |
|---|---|
| ligand diffusion / reconstruction loss | 保证生成分子基本正确 |
| pocket displacement / rotamer loss | 学会 apo-to-holo pocket adaptation |
| contact recovery loss | 保留 holo-like binding contact |
| clash penalty | 避免 ligand 被 pocket 挤坏 |
| sparsity budget loss | 避免全 pocket 都被释放 |
| background stability loss | 防止远处 residues 被带着乱动 |
| shell smoothness loss | 避免 core/background 边界硬切 |

大白话：

> 不要手写“永远 4A shell、0.25 权重”。顶会方法应该让模型根据当前 ligand 状态判断：哪里强动，哪里轻动，哪里别动。

---

## 12. 一句话结论

> hard8 n=3 repeat 说明：protein RMSD 上 hard distance top4 和几个 active-set shell candidates 差距很小；真正有价值的信号来自几何诊断和 ligand-quality 补充评估。random top4 证明 sparse 本身会误导，shell4_w025 则给出目前最均衡的 active-set 证据：更少 clash、更高 contact precision/Jaccard、更低边界跳变，并且 QED/complete 也较稳。下一步应从固定 shell 规则推进到 learned ligand-conditioned active-set policy，并补真正 docking / PoseBusters-style ligand-side evaluation。
