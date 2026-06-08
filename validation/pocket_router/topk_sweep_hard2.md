# PocketRouter Top-k Sweep: Hard2

日期：2026-05-28

这份报告记录 `distance` router 的 top-k sweep。它回答一个很关键的问题：

> PocketRouter 的收益到底是不是来自“合理稀疏预算”？还是随便设一个 top12 就碰巧有效？

大白话：

> ligand 附近到底让几个 residues 动最合适？4 个、8 个、12 个，还是越多越好？

---

## 1. 实验设置

运行目录：

`validation/ab_runs/hard2_distance_topk_sweep_steps1000`

数据：

| test position | original index | metadata apo-holo RMSD | ligand path |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |

共同设置：

- sampling steps：1000
- samples per case：1
- batch size：1
- init center：apo
- sample atom number：prior
- update schedule：late_dense
- router mode：distance
- docking：none

比较的 top-k：

- top4
- top8
- top12
- top16
- top24

其中 `top12` 复用了上一轮 hard2 中完全相同设置下的结果，避免重复计算。

---

## 2. 结果

结果文件：

`validation/ab_runs/hard2_distance_topk_sweep_steps1000/results.json`

| distance router | mean protein RMSD | min-max RMSD | mean TM-score | mean selected residues |
|---|---:|---:|---:|---:|
| `top4` | 2.0624 | 2.0402-2.0846 | 0.8796 | 4 |
| `top8` | 2.1050 | 2.1046-2.1053 | 0.8784 | 8 |
| `top12` | 2.1731 | 2.1205-2.2258 | 0.8764 | 12 |
| `top16` | 2.2407 | 2.1333-2.3482 | 0.8745 | 16 |
| `top24` | 2.3321 | 2.1631-2.5010 | 0.8695 | 24 |

最优：

> `distance_top4`

平均 protein RMSD：

> 2.0624 A

---

## 3. 怎么读

### 3.1 top4 最好

这次 sweep 里，top4 比 top8、top12、top16、top24 都好。

大白话：

> 不是“多选一点更保险”，而是只让最贴近 ligand 的极小局部动，反而最稳。

### 3.2 top-k 越大，平均 RMSD 越差

这次结果几乎是单调变差：

```text
top4  -> 2.0624
top8  -> 2.1050
top12 -> 2.1731
top16 -> 2.2407
top24 -> 2.3321
```

大白话：

> 更新范围越大，更多不该动的 residues 被拉进来，模型反而更容易把 pocket 搅乱。

### 3.3 这加强了论文主线

之前 hard2 已经说明：

- dense update 不够好；
- random sparse 不如 distance sparse；
- motion oracle 和 contact-change oracle 不稳。

现在 top-k sweep 进一步说明：

> 关键不是“能 sparse 就行”，而是“在当前 ligand contact 附近用很小预算做精准更新”。

这比单独说 top12 有效更强，因为它展示了一个可解释的预算规律。

---

## 4. 对顶会叙事的影响

之前可以说：

> distance top12 是强 baseline。

现在应该改成：

> best distance top-k 是强 baseline，当前 hard2 上 top4 最强。

这会提高 learned router 的门槛：

> learned router 不能只打过 random 或 top12；它必须打过 best hand-crafted distance top-k，当前就是 top4。

大白话：

> 如果手写规则只选最近 4 个 residues 已经很好，学习出来的 router 至少要比这个更聪明，否则不能算顶会主方法。

---

## 5. 下一步建议

### 5.1 hard8 上优先跑 top4

下一轮 hard8 不建议先全量 sweep。更合理的是：

- static5
- late_dense
- random_top4
- distance_top4
- distance_top8 或 distance_top12 作为预算对照

目的：

> 验证 top4 在更多 hard cases 上是否稳定。

### 5.2 训练 learned router 时用 top4 做强 baseline

learned router 的预算可以先固定为 4 或 8。

训练目标应该围绕：

- 当前 ligand-fragment 距离/contact；
- residue chemistry compatibility；
- denoising timestep；
- local clash / residual proxy；
- contact-change 和 displacement 作为辅助信号，而不是唯一标签。

### 5.3 增加 per-case 和 multi-seed

当前每个 case 只有一个 sample，所以还不能过度解释 top4 的绝对优势。

需要补：

- hard8；
- 每个 case 3 seeds；
- 统计均值、方差、胜率；
- ligand 侧指标和 clash/docking 指标。

---

## 6. 一句话结论

top-k sweep 把 PocketRouter 的主线进一步收紧：

> 当前最有证据的机制是 very sparse current-contact pocket update。对 hard2 来说，只更新当前 ligand 附近 top4 residues 比 top8/top12/top16/top24 都更好。

大白话：

> 不是 pocket 动得越多越好，而是要非常克制地只动 ligand 当前真正碰到的那一小块。

