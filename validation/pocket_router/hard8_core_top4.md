# Hard8 Core PocketRouter Top4 实验记录

日期：2026-05-28

这份文档记录这次 hard8 core 实验。它的目的不是再证明“代码能跑”，而是回答一个更接近顶会标准的问题：

> 前面 hard2 看到的 `distance_top4` 信号，放到 8 个 hardest cases 上还稳不稳？

大白话说：

> 如果只在一两个 case 上好，可能是巧合；如果在 8 个最难 case 上仍然好，说明这个 idea 至少有一个值得继续打磨的核心现象。

---

## 1. 这次实验想验证什么

前面实验已经看到一个现象：

> 在 Apo2Mol 采样过程中，不让整个 pocket 每次都动，而是只更新当前 ligand 附近很少几个 residues，protein RMSD 明显下降。

但这件事有几个风险：

1. 可能只是某一两个 case 的偶然结果；
2. 可能只是因为 top-k sparse mask 有正则化效果，和“选得准”关系不大；
3. 可能 top-k 不是越小越好，也不是越大越好，需要找到预算；
4. 可能 dense late update 已经足够，router 其实没必要。

所以这次 hard8 core 实验专门做一个小而关键的对照：

| arm | 大白话目的 |
|---|---|
| `baseline_realistic_static5` | 原始风格：少数 timestep 更新整个 pocket |
| `control_realistic_late_dense` | 同样更多次更新，但每次仍然更新整个 pocket |
| `pocket_router_random_top4` | 只更新 4 个 residue，但随机选，看“稀疏本身”有没有用 |
| `pocket_router_distance_top4` | 只更新当前 ligand 最近的 4 个 residue，这是主打手工 router |
| `pocket_router_distance_top8` | 只更新最近 8 个 residue，看 top-k 变大是否更好 |

这组实验最重要的判断标准是：

> `distance_top4` 要同时超过 dense baseline、late-dense control、random top4，并且 top4/top8 的趋势要合理。

---

## 2. 我们具体改动了什么

这次没有改模型主干，也没有重新训练。

主要是在实验编排脚本里增加一个 hard-tail core 模式：

- 文件：`validation/run_new_method_ab.py`
- 新增 arm 组：`HARDTAIL_CORE_ROUTER_ARMS`
- 新增命令参数：`--router-hardtail-core`
- 新增汇总字段：
  - `hardtail_core_comparison`
  - `distance_topk_comparison`

这次的核心 arms 是：

```text
baseline_realistic_static5
control_realistic_late_dense
pocket_router_random_top4
pocket_router_distance_top4
pocket_router_distance_top8
```

大白话解释：

> 我们把“只动最近 4 个 residues”放到一个更公平的擂台里，不只和原始 baseline 比，也和“随机动 4 个 residues”“动最近 8 个 residues”“后期更频繁地动整个 pocket”比。

---

## 3. 选择了哪些数据

实验选择 test set 里 apo/holo pocket RMSD 最大的 8 个 hard cases。

| test position | original index | metadata apo/holo RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |
| 390 | 24413 | 3.3143 | `8u6b__1__1.A__1.C/1.C.sdf` |
| 310 | 24333 | 2.9210 | `8ow3__1__1.A__1.B/1.B.sdf` |
| 377 | 24400 | 2.7470 | `8sfu__2__1.B__1.E/1.E.sdf` |
| 342 | 24365 | 2.6896 | `8qn5__3__1.C__1.L/1.L.sdf` |
| 365 | 24388 | 2.6712 | `8pqh__1__1.A__1.B/1.B.sdf` |
| 347 | 24370 | 2.4535 | `8sbv__2__1.B__1.G/1.G.sdf` |

为什么选这些？

> 因为它们是 pocket 变化最大的难题。这里能改善，比在简单 case 上改善更有说服力。

---

## 4. 采样协议

所有 arm 都使用同样的现实设置：

| 设置 | 值 | 大白话 |
|---|---|---|
| `num_steps` | 1000 | 跑完整 reverse diffusion，不是 smoke test |
| `num_samples` | 1 | 每个 case 生成 1 个样本 |
| `batch_size` | 1 | CPU 本地稳定跑 |
| `sample_num_atoms` | `prior` | 分子原子数按 prior 采样 |
| `init_center_mode` | `apo` | 从 apo pocket 初始化，不偷看 holo center |
| `docking_mode` | `none` | 本轮只看生成和 protein adaptation，不跑 docking |

实际命令：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --router-hardtail-core \
  --run-dir validation/ab_runs/hard8_core_top4_steps1000 \
  --num-cases 8 \
  --num-samples 1 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none
```

输出目录：

```text
validation/ab_runs/hard8_core_top4_steps1000
```

核心结果文件：

```text
validation/ab_runs/hard8_core_top4_steps1000/results.json
validation/ab_runs/hard8_core_top4_steps1000/preflight_report.md
```

---

## 5. 每个 arm 到底怎么更新 pocket

### 5.1 baseline_realistic_static5

配置：

- update schedule：`static5`
- router：`none`
- 更新 timestep：`799, 599, 399, 199, 10`

含义：

> 原始风格。只在 5 个关键 timestep 更新 protein，但每次更新时整个 pocket 都可以动。

这里的 `router_selected_counts` 平均是 52.25，不代表 router 选了 52 个，而是没有 sparse router 时，整个 pocket 的 residues 都参与更新。

### 5.2 control_realistic_late_dense

配置：

- update schedule：`late_dense`
- router：`none`
- 更新 timestep：`460, 410, 360, 310, 260, 210, 160, 110, 60, 10`

含义：

> 后半程 ligand 更成形以后，让 protein 更频繁地整体更新。

这个 arm 用来排除一个很关键的问题：

> 如果 `distance_top4` 变好，只是因为它用了 10 次 update，而不是 5 次 update 呢？

所以 late-dense 也用 10 次 update，但不做 top-k sparse。

### 5.3 pocket_router_random_top4

配置：

- update schedule：`late_dense`
- router：`random`
- top-k：4

含义：

> 每次只允许 4 个 residues 更新，但这 4 个 residues 是随机选的。

这个 arm 是负控制。

如果 random top4 和 distance top4 一样好，说明：

> 可能不是“选 ligand 附近 residues”有用，而只是“少动几个 residues”有用。

### 5.4 pocket_router_distance_top4

配置：

- update schedule：`late_dense`
- router：`distance`
- top-k：4

含义：

> 每次看当前 ligand 在哪里，只更新离当前 ligand 最近的 4 个 residues。

这是当前最强的手工 router。

大白话：

> ligand 长到哪里，就只让它附近最相关的一小块 pocket 动起来。

### 5.5 pocket_router_distance_top8

配置：

- update schedule：`late_dense`
- router：`distance`
- top-k：8

含义：

> 和 top4 一样，但预算扩大到 8 个 residues。

这个 arm 用来判断：

> 是不是多给一些 residues 更好？

如果 top8 比 top4 差，说明这个现象更像“非常局部的精修”，不是“越多越好”。

---

## 6. 采样完成情况

这次最终每个 arm 都有 8 个 result 文件：

| arm | result files |
|---|---:|
| `baseline_realistic_static5` | 8 |
| `control_realistic_late_dense` | 8 |
| `pocket_router_random_top4` | 8 |
| `pocket_router_distance_top4` | 8 |
| `pocket_router_distance_top8` | 8 |

总计：

```text
40 result_*.pt files
```

运行过程中复用了前面相同配置已经完成的部分样本，并补齐缺失样本；主控脚本会检查已有 `result_*.pt`，已齐全的 arm 会跳过采样，只做评估和汇总。

---

## 7. Protein adaptation 结果

核心结果如下：

| arm | mean protein RMSD | mean TM-score | selected residues/update | delta vs static5 | delta vs late-dense |
|---|---:|---:|---:|---:|---:|
| `baseline_realistic_static5` | 2.4365 | 0.8977 | 52.25 | 0.0000 | +0.0223 |
| `control_realistic_late_dense` | 2.4142 | 0.9005 | 52.25 | -0.0223 | 0.0000 |
| `pocket_router_random_top4` | 1.8305 | 0.9248 | 4.00 | -0.6060 | -0.5837 |
| `pocket_router_distance_top4` | 1.7924 | 0.9249 | 4.00 | -0.6441 | -0.6218 |
| `pocket_router_distance_top8` | 1.8371 | 0.9229 | 8.00 | -0.5994 | -0.5772 |

最优 arm：

```text
pocket_router_distance_top4
```

对应结果：

```text
mean protein RMSD = 1.7924 A
mean TM-score     = 0.9249
```

相对 baseline：

```text
2.4365 A -> 1.7924 A
improvement = 0.6441 A
```

相对 late-dense：

```text
2.4142 A -> 1.7924 A
improvement = 0.6218 A
```

大白话：

> 在 8 个最难 case 上，只更新当前 ligand 附近 4 个 residues，比让整个 pocket 后期密集更新更好。

---

## 8. Top4 和 Top8 的结论

distance router 的 top-k 对比：

| distance router | mean protein RMSD | mean TM-score |
|---|---:|---:|
| `top4` | 1.7924 | 0.9249 |
| `top8` | 1.8371 | 0.9229 |

top4 更好。

这和 hard2 top-k sweep 的趋势一致：

| hard2 distance router | mean protein RMSD |
|---|---:|
| `top4` | 2.0624 |
| `top8` | 2.1050 |
| `top12` | 2.1731 |
| `top16` | 2.2407 |
| `top24` | 2.3321 |

大白话：

> 不是更新越多 residues 越好。当前结果更像是：更新太多会把无关区域也带着动，反而增加噪声。

这对论文故事很有用，因为它支持一个明确机制：

> fragment-conditioned very-local pocket update

也就是：

> 当前 fragment 只需要一个很小的局部 pocket 工作区，不需要每一步都唤醒整个 pocket。

---

## 9. Random top4 的重要信号

random top4 的结果：

```text
mean protein RMSD = 1.8305 A
```

它也明显好于 dense baseline：

```text
baseline static5 = 2.4365 A
late-dense       = 2.4142 A
random top4      = 1.8305 A
```

这说明：

> sparse mask 本身就有很强的 regularization 效果。

大白话：

> 哪怕随机少动几个 residues，也比整个 pocket 都动要稳很多。

但 distance top4 仍然比 random top4 更好：

```text
random top4   = 1.8305 A
distance top4 = 1.7924 A
difference    = 0.0381 A
```

这个差距目前不大。

所以现在不能过度宣称：

> 我们已经证明了 contact router 本身非常强。

更准确的说法是：

> 我们证明了 sparse local update 是强信号；distance/contact routing 是当前最好的手工版本，但它相对 random 的优势还需要更多 case、多 seed、以及 learned router 来放大和稳定。

---

## 10. 分子评估结果

本轮 `docking_mode=none`，所以没有 docking score。

并且当前 eval 里：

```text
evaluated_mols = 0
QED = null
SA = null
```

这意味着：

> 这轮结果主要能支持 protein adaptation 结论，不能支持 ligand quality / docking 结论。

分子 reconstruction / completeness 结果如下：

| arm | mol stable | atom stable | recon success | eval success | complete |
|---|---:|---:|---:|---:|---:|
| `baseline_realistic_static5` | 0.2500 | 0.8750 | 0.8750 | 0.8750 | 0.8750 |
| `control_realistic_late_dense` | 0.2500 | 0.8583 | 0.8750 | 0.8750 | 0.8750 |
| `pocket_router_random_top4` | 0.3750 | 0.8881 | 0.8750 | 0.6250 | 0.6250 |
| `pocket_router_distance_top4` | 0.3750 | 0.9292 | 0.8750 | 0.7500 | 0.7500 |
| `pocket_router_distance_top8` | 0.2500 | 0.8417 | 0.8750 | 0.7500 | 0.7500 |

这里要非常谨慎：

- `distance_top4` 的 atom stability 最好；
- 但 complete 比 dense baseline 低；
- QED/SA/docking 都没有有效数值；
- 每个 case 只有 1 个 sample，统计还很粗。

大白话：

> pocket 结构变好了，但药物分子质量还没有被完整证明。

---

## 11. 这次实验能说明什么

可以比较有把握地说：

1. sparse pocket update 在 hard-tail cases 上是强信号；
2. dense late update 不是答案，单纯多更新整个 pocket 只带来很小改善；
3. very-local top4 比 top8 更好，支持“局部精修”而不是“扩大更新范围”；
4. distance/contact top4 是当前最好的手工 router；
5. random top4 很强，说明 top-k 稀疏本身是一个强正则化 baseline。

最重要的结论：

> 现在这个方向值得继续做，但顶会主张必须从“手工 distance top4 有效”推进到“learned fragment-conditioned router 稳定超过 distance top4 和 random top4”。

---

## 12. 这次实验不能说明什么

还不能说：

1. 这已经是顶会级别方法；
2. learned router 一定会超过 distance top4；
3. ligand 质量已经变好；
4. docking score 已经变好；
5. 这个结果在 full test set、多 seed、多 sample 下仍然稳定。

原因很直接：

- 当前只是 8 个 hard cases；
- 每个 case 只有 1 个 sample；
- 没有 docking；
- QED/SA 没有有效统计；
- random top4 非常接近 distance top4；
- router 还不是训练出来的。

---

## 13. 对顶会方向的影响

这次结果把路线变得更清楚了。

不要把论文核心写成：

> 我们发现 distance top4 很好。

这太像工程 trick。

更有潜力的写法是：

> 生成过程中的 pocket adaptation 应该是 state-conditioned sparse memory routing：当前 fragment 只激活少量 contact-relevant pocket memory slots，避免 dense pocket update 的噪声传播。

大白话：

> 不是所有 pocket residues 都应该在每一步参与决策。当前分子片段只需要叫醒它附近最相关的一小撮 residues。

这次实验给这个故事提供了两个关键证据：

1. sparse top4 大幅超过 dense update；
2. top4 超过 top8，说明“少而准”比“多而全”更好。

但顶会门槛还差一块：

> learned router 必须显著超过 random top4 和 distance top4。

如果 learned router 只能超过 dense baseline，不能超过 distance top4，那么论文贡献会弱很多，因为最强 baseline 只是一个简单几何规则。

---

## 14. 下一步实验建议

### 14.1 先做多 seed / 多 sample 稳定性

当前 random top4 和 distance top4 差距只有 0.0381 A。

下一步要问：

> 这个差距是真信号，还是采样随机性？

建议：

- hard8；
- 每个 case 3-5 samples；
- 至少 3 个 random seeds；
- 统计 mean / median / std / best-of-N。

### 14.2 做 learned router v1

目标不是复杂，而是先打败强基线。

最小 learned router 可以输入：

- residue center；
- current ligand atom center；
- residue-ligand distance features；
- residue type embedding；
- denoising timestep；
- 当前 ligand residual / clash proxy。

输出：

- 每个 residue 一个 relevance score；
- top-k 或 soft top-k gate。

训练信号可以先用弱监督：

- current/holo contact labels；
- apo-to-holo displacement；
- contact-change labels；
- ligand-near labels。

评价标准：

> learned_top4 必须超过 distance_top4 和 random_top4。

### 14.3 补 ligand quality / docking

需要补：

- valid molecule count；
- QED；
- SA；
- vina score；
- clash / contact / PoseBusters-style checks；
- ligand RMSD 或 pocket-ligand interaction similarity。

否则当前结果只能讲 pocket adaptation，不能完整讲 SBDD。

### 14.4 扩到 full hard-tail / full test

hard8 是强信号，但还不够。

建议分三层：

1. hard8：快速迭代机制；
2. hard32 / hard64：确认 hard-tail 稳定性；
3. full test：最终论文表格。

---

## 15. 当前一句话结论

这次 hard8 core 实验支持继续推进 PocketRouter：

> 在 8 个 hardest apo-to-holo cases 上，`distance_top4` 把 mean protein RMSD 从 `2.4365 A` 降到 `1.7924 A`，比 late-dense 整体更新也低 `0.6218 A`，并且 top4 优于 top8。

但它还不是最终顶会证据：

> random top4 也很强，ligand quality/docking 还没补齐，所以真正的下一步是训练 learned fragment-conditioned router，并要求它稳定超过 `distance_top4` 这个强手工基线。
