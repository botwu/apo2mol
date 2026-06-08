# PocketRouter Model Validation: Hard2 Follow-up

日期：2026-05-27

这份报告记录第二轮模型层验证。它不是重新训练模型，而是在已经下载好的 Apo2Mol 数据和已有 checkpoint 上，继续检验 PocketRouter 这个方向是否值得往顶会投稿推进。

大白话总结：

> 这次实验是在问：只让 ligand 附近的一小圈 pocket residues 动，是不是比“整个 pocket 一起动”更靠谱？如果随机挑 residues 也一样好，那这个 idea 就不成立；如果 contact/distance 挑出来的 residues 明显更好，说明 router 真的抓到了有用局部。

---

## 1. 这次实验想回答什么

hard1 pilot 已经看到一个正信号：

- `distance_top12` 在一个 hard case 上把 protein RMSD 从 3.1191 降到 2.2258；
- `motion_oracle_top12` 和 `random_top12` 没有明显跟上；
- 说明“离当前 ligand 近的 residues”可能比“apo-to-holo 位移最大的 residues”更适合作为 sparse pocket update 区域。

但是 hard1 只有一个 case，不能支撑顶会级结论。

所以 hard2 follow-up 主要回答三个问题：

1. `distance_top12` 的优势能不能在第二个 hard case 上继续成立？
2. `random_top12` 会不会也差不多好？如果是，那收益可能只是 sparse regularization，不是 router。
3. `contact_change_oracle_top12` 会不会比 distance 更强？如果它更强，后续 learned router 就应该优先学习 contact-change；如果它不强，主线就要收紧到 current-fragment distance/contact relevance。

---

## 2. 使用的数据

运行目录：

`validation/ab_runs/hard2_steps1000_router_v2`

选择的是 test split 里 apo-holo RMSD 最大的两个 hard cases：

| test position | original index | metadata apo-holo RMSD | ligand path |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |

大白话：

> 我们不是挑简单样本，而是挑了 pocket 变化最大的两个样本。这样更能看出方法有没有处理 induced-fit 的能力。

---

## 3. 共同采样设置

所有 arms 使用相同的基础采样条件：

- checkpoint：`./apo2mol_dataset/apo2mol_checkpoint.ckpt`
- 采样步数：1000 reverse denoising steps
- 每个 case 采样数：1
- batch size：1
- ligand atom number：`prior`
- 初始化中心：`apo`
- docking：未开启
- 设备：本机自动落到 CPU

大白话：

> 这不是训练新模型，而是拿同一个模型，用不同的 pocket update 策略去采样。这样差异主要来自 update / router，而不是训练随机性。

---

## 4. 比较了哪些 arms

### 4.1 `baseline_realistic_static5`

原始风格的稀疏时间更新：

- 只在少数 timestep 更新 protein；
- 一旦更新，就是对整个 pocket 做全局更新。

大白话：

> pocket 偶尔整体动一下。

### 4.2 `control_realistic_late_dense`

更频繁的 late-stage dense update：

- 在后期 denoising 中更密集地更新；
- 每次仍然是整个 pocket 更新。

大白话：

> 不是只偶尔动，而是后期多动几次；但每次还是整块 pocket 一起动。

这个 arm 的作用：

> 排除“distance router 只是因为 update 次数更多才变好”的可能。

### 4.3 `pocket_router_distance_top12`

同样使用 late-stage update schedule，但每次只允许 top-12 residues 做高精度更新。

选择规则：

> 选当前 ligand 附近最近的 12 个 residues。

大白话：

> ligand 现在靠近谁，就让谁动；离 ligand 远的 residue 先别乱动。

### 4.4 `pocket_router_motion_oracle_top12`

选择 apo-to-holo 位移最大的 12 个 residues。

大白话：

> 事后看 ground truth，谁最终动得最大，就选谁。

这个 arm 看起来像 oracle，但它不是我们最终想要的方法。它回答的是：

> 只选“最会动”的 residues 是否足够？

### 4.5 `pocket_router_contact_change_oracle_top12`

选择 apo/holo contact 变化最大的 12 个 residues。

大白话：

> 哪些 residues 的接触关系在 apo 到 holo 之间变化最大，就选哪些。

这个 arm 是本轮新补的关键对照。它回答的是：

> 如果我们知道最终 contact 怎么变，只更新这些 residues 会不会最好？

### 4.6 `pocket_router_random_top12`

随机选择 12 个 residues 更新。

大白话：

> 随便挑 12 个 residues。

这个 arm 是负对照。它回答的是：

> 是不是只要 sparse 一下就会变好？如果 random 也一样好，那 router 没有说服力。

---

## 5. 代码层到底改了什么

本轮之前已经把 PocketRouter 接到了采样链路里，本轮补充的是 `contact_change_oracle` arm。

核心改动在：

- `models/molopt_score_model.py`
- `sample_split.py`
- `configs/training.yaml`
- `validation/run_new_method_ab.py`

### 5.1 在模型里加 router mask

模型每次准备更新 protein 时，会构造一个 residue-level mask。

如果是 dense update：

> 所有 pocket residues 都允许更新。

如果是 router update：

> 只有被 router 选中的 top-k residues 更新，其它 residues 保持原状。

大白话：

> 每次 protein 要动之前，先问 router：“这一步到底让谁动？”

### 5.2 非选中 residues 怎么处理

对于没有被选中的 residues：

- rigid transform 不更新；
- translation 设为 0；
- side-chain chi update 设为 0；
- 等价于保持上一时刻状态。

大白话：

> 没选中的 residues 不是删掉，而是先按兵不动。

### 5.3 在采样脚本里记录 router 使用情况

`sample_split.py` 会把每次更新选了多少 residues 保存到结果里。

所以结果里会看到：

- dense arms：平均约 48 个 residues 更新；
- router arms：每次固定 12 个 residues 更新。

大白话：

> 我们不仅看结果，还看它到底是不是按预算只动了 12 个 residues。

### 5.4 在 A/B 脚本里补 contact-change oracle

`validation/run_new_method_ab.py` 新增：

`pocket_router_contact_change_oracle_top12`

配置上它和其它 router arms 一样：

- `protein_update_schedule: late_dense`
- `protein_update_interval: 50`
- `protein_update_min_t: 10`
- `protein_update_residual_threshold: 0.5`
- `pocket_router_topk: 12`

唯一不同是：

`pocket_router_mode: contact_change_oracle`

---

## 6. 运行命令

语法检查：

```bash
.venv310/bin/python -m py_compile validation/run_new_method_ab.py models/molopt_score_model.py sample_split.py validation/validate_pocket_router.py
```

preflight：

```bash
.venv310/bin/python validation/run_new_method_ab.py --router-validation --run-dir validation/ab_runs/hard2_steps1000_router_v2 --num-cases 2 --num-samples 1 --batch-size 1 --num-steps 1000 --docking-mode none
```

正式运行：

```bash
.venv310/bin/python validation/run_new_method_ab.py --run --router-validation --run-dir validation/ab_runs/hard2_steps1000_router_v2 --num-cases 2 --num-samples 1 --batch-size 1 --num-steps 1000 --docking-mode none
```

---

## 7. 核心结果

结果文件：

`validation/ab_runs/hard2_steps1000_router_v2/results.json`

| arm | result files | mean protein RMSD | min-max RMSD | mean TM-score | delta vs static5 | delta vs late-dense |
|---|---:|---:|---:|---:|---:|---:|
| `baseline_realistic_static5` | 2 | 2.8381 | 2.5501-3.1261 | 0.8483 | 0.0000 | +0.0514 |
| `control_realistic_late_dense` | 2 | 2.7867 | 2.4542-3.1191 | 0.8496 | -0.0514 | 0.0000 |
| `pocket_router_distance_top12` | 2 | 2.1731 | 2.1205-2.2258 | 0.8764 | -0.6650 | -0.6135 |
| `pocket_router_motion_oracle_top12` | 2 | 2.6050 | 2.1548-3.0552 | 0.8596 | -0.2331 | -0.1817 |
| `pocket_router_contact_change_oracle_top12` | 2 | 2.8072 | 2.2100-3.4044 | 0.8493 | -0.0309 | +0.0206 |
| `pocket_router_random_top12` | 2 | 2.5913 | 2.1264-3.0563 | 0.8574 | -0.2468 | -0.1953 |

最强 arm：

> `pocket_router_distance_top12`

平均 protein RMSD：

> 2.1731 A

相对原始 static5 改善：

> -0.6650 A

相对 late-dense control 改善：

> -0.6135 A

---

## 8. 怎么读这个结果

### 8.1 distance router 复现了 hard1 的正信号

hard1：

- `distance_top12`: 2.2258
- `late_dense`: 3.1191
- 改善：0.8933 A

hard2：

- `distance_top12`: 2.1731
- `late_dense`: 2.7867
- 改善：0.6135 A

大白话：

> 第二个 hard case 加进来以后，distance router 还是最强。这个方向不是只在一个样本上碰巧赢。

### 8.2 late-dense 不是答案

`late_dense` 相比 static5 只从 2.8381 到 2.7867，改善 0.0514 A。

大白话：

> 只是让整个 pocket 后期多动几次，基本没有解决问题。

这说明 distance router 的收益不是简单来自“update 次数变多”。

### 8.3 random 有改善，但明显不如 distance

`random_top12` 平均 RMSD 是 2.5913，比 static5 好 0.2468 A，但比 distance 差 0.4182 A。

大白话：

> sparse 本身确实有一点 regularization 效果，因为少动一些 residues 会减少噪声；但随便 sparse 不如选 ligand 附近的 residues。

这对论文很重要：

> 我们不能只说 sparse 有用，而要证明“选哪些 residues”有用。

### 8.4 motion oracle 不是主线

`motion_oracle_top12` 平均 RMSD 是 2.6050，和 random 接近，明显不如 distance。

大白话：

> 最终动得最大的 residues，不一定是当前生成过程中最该更新的 residues。

原因可能是：

- 大位移 residues 可能是远端柔性区域；
- 它们最终会动，但当前 fragment 还没有接触它们；
- 硬选 mobile residues 会让模型在不合适的时间更新不合适的位置。

### 8.5 contact-change oracle 这轮不稳

`contact_change_oracle_top12` 平均 RMSD 是 2.8072，接近 baseline，并且比 late-dense 略差。

这不是说 contact-change 永远没用，而是说明：

> 静态 ground-truth contact-change 标签，不一定等于当前 denoising step 最该更新的 residue 集合。

大白话：

> 最终接触关系会变的地方，不一定在 ligand 当前长到一半的时候就应该动。

可能原因：

- contact-change 是最终 holo 相对 apo 的事后标签，和当前 partial ligand 状态不完全对齐；
- top-12 contact-change residues 可能覆盖了最终构象变化，但没有覆盖当前 fragment 的局部接触；
- contact-change label 对 threshold / ligand pose / residue granularity 很敏感；
- 当前模型没有训练过这种 hard mask，直接换 mask 可能导致分布不匹配。

---

## 9. 当前最稳的结论

这轮 hard2 后，最稳妥的结论是：

> PocketRouter 不是“选最会动的 residues”，也不是“随便 sparse 一下”。目前最有信号的是 current-ligand distance/contact conditioned sparse pocket update。

更像顶会主线的表述是：

> fragment-conditioned contact-sparse pocket adaptation

大白话：

> ligand 当前长到哪里，就让附近真正相关的 pocket 局部动；不要让全 pocket 每次一起动。

---

## 10. 对顶会投稿意味着什么

### 10.1 现在还不能直接说能投顶会

原因：

- hard2 仍然只有两个 cases；
- 还没有 learned router；
- ligand 侧评估不完整，`QED/SA/PoseBusters/docking` 还没跑通；
- 目前最强的是 hand-crafted distance router，不是一个完整可学习方法；
- 还没有和其它 SBDD / apo-to-holo baselines 做系统比较。

### 10.2 但它已经给了一个值得推进的主创新点

现在最有价值的创新点不是“加 sparse attention”，而是：

> state-conditioned sparse pocket routing for fragment-conditioned apo-to-holo adaptation

大白话：

> 把 pocket 看成一个很大的记忆库。每一步 ligand 生成时，只去调用当前 fragment 真正相关的一小块 pocket 记忆。

### 10.3 顶会级要求

如果要往 NeurIPS / ICML / ICLR / ACL / ECCV 这种风格推进，后续必须让方法满足：

1. 不是工程 patch，而是一个清楚的计算原则；
2. 有 learned router，而不是只靠 distance heuristic；
3. learned router 必须打过 `random_top12`、`late_dense`，以及 best hand-crafted distance top-k；
4. 在 hard-tail 大样本上稳定提升；
5. ligand validity、clash、docking、pocket adaptation 都不能掉；
6. 有清楚可解释的 router behavior，例如 router 是否真的选中 current-fragment contact residues。

---

## 11. 下一步实验建议

### 11.1 先做 top-k sweep

对 distance router 跑：

- top4
- top8
- top12
- top16
- top24

目的：

> 找到“动太少不够、动太多变吵”的 sweet spot。

已完成的 hard2 top-k sweep 显示，当前两个 hardest cases 上 `top4` 最好：

| distance router | mean protein RMSD | mean TM-score |
|---|---:|---:|
| `top4` | 2.0624 | 0.8796 |
| `top8` | 2.1050 | 0.8784 |
| `top12` | 2.1731 | 0.8764 |
| `top16` | 2.2407 | 0.8745 |
| `top24` | 2.3321 | 0.8695 |

这个结果说明：

> 在当前 hard2 上，不是选越多 residues 越好。只更新 ligand 最近的极小局部，即 top4，反而最稳。

所以 learned router 后续不能只打过 `distance_top12`，而要打过当前 best hand-crafted baseline：`distance_top4`。

### 11.2 跑 hard8 / hard-tail

至少扩到：

- hard8：沿用数据层验证的 8 个 hardest cases；
- hard20：如果 CPU 时间允许；
- 每个 case 至少 3 个 samples，减少单次采样随机性。

### 11.3 做 learned router

学习目标不要只用 motion，也不要只用 contact-change。

更合理的是混合目标：

- 当前 fragment 到 residue 的距离/contact；
- residue 与 ligand 的化学兼容性；
- apo-to-holo contact-change；
- apo-to-holo displacement；
- denoising timestep；
- 当前预测残差或 clash proxy。

大白话：

> 不是问“哪个 residue 最终动得最大”，而是问“当前这一步，哪个 residue 最值得花精力更新”。

### 11.4 补 ligand 侧指标

现在 results 里的 protein RMSD / TM-score 是可用的，但 ligand 侧评估还不完整：

- QED / SA 是 null；
- evaluated molecules 是 0；
- docking 没开；
- PoseBusters / clash 没跑。

顶会投稿不能只看 protein RMSD。

需要补：

- reconstruction success；
- atom stability；
- molecule stability；
- QED / SA；
- docking score；
- clash/contact validity；
- generated ligand diversity。

---

## 12. 一句话更新结论

hard2 follow-up 把结论从：

> PocketRouter 可能有用。

推进到：

> contact/distance-conditioned sparse pocket update 有可复现实验信号，而且 hard2 top-k sweep 显示 very sparse top4 当前最强；最终顶会方法必须把这个 hand-crafted distance router 升级成 learned state-conditioned router，并且稳定打过 best distance-only baseline。
