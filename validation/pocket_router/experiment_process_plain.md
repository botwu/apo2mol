# PocketRouter 实验全过程说明

日期：2026-05-27

最近更新：2026-05-29，加入 hard8 conservative active-set shell sweep。

这份文档用尽量大白话的方式解释我们刚刚做的 PocketRouter 实验：

- 我们到底想验证什么；
- 改了哪些代码；
- 数据层验证怎么做；
- 模型层 A/B 怎么跑；
- 每个结果代表什么；
- 为什么现在的结论不是“已经能投顶会”，而是“方向有信号，但还要继续扩展”。

---

## 0. 一句话总结

我们验证的是：

> ligand 生成过程中，protein pocket 不一定每次都要整体更新。也许只更新当前 ligand 附近最相关的一小部分 residues，反而更稳定、更准。

实验结果显示：

> hard1 上，只更新离当前 ligand 最近的 top-12 residues，protein RMSD 从 3.1191 降到 2.2258，改善约 0.8933 A。

扩展到两个 hardest cases 后：

> `distance_top12` 平均 protein RMSD 为 2.1731，优于 static5 的 2.8381 和 late-dense 的 2.7867。

继续扩展到 8 个 hardest cases 后：

> `distance_top4` 平均 protein RMSD 为 1.7924，优于 static5 的 2.4365、late-dense 的 2.4142、random top4 的 1.8305、distance top8 的 1.8371。

但也要注意：

> 这不是证明“选最会动的 residues 就行”。实验反而显示，真正有用的是当前 ligand/contact 附近 residues，而不是单纯 apo-to-holo 位移最大的 residues，也不是静态 contact-change oracle。并且 random top4 也很强，所以 learned router 后续必须超过 random top4 和 distance top4，不能只超过 dense baseline。

所以目前最准确的 idea 应该叫：

> Fragment-Guided Contact-Sparse Pocket Memory / PocketRouter

大白话：

> 药物分子长到哪里，就重点让那个局部 pocket 动起来；不要让整个 pocket 每一步都乱动。

---

## 1. 原始问题是什么

Apo2Mol 做的是 apo-to-holo / SBDD 方向的问题。

简单说：

- apo protein：没有 ligand 结合时的蛋白构象；
- holo protein：ligand 结合后的蛋白构象；
- ligand 结合后，pocket 往往会发生 induced-fit 变化；
- 模型希望从 apo pocket 出发，生成 ligand，同时让 pocket 往合理的 holo 状态靠近。

原始 Apo2Mol 的采样过程中，会在少数几个 denoising timestep 更新 protein pocket。

可以理解为：

> ligand 一边生成，protein pocket 偶尔跟着动一下。

但是这里有一个潜在问题：

> 每次更新 protein 时，如果让整个 pocket 都参与高精度更新，可能会引入很多无关扰动。

因为实际 ligand 当前只接触 pocket 的一小部分区域。很多 residues 虽然在 pocket 里，但当前 fragment 根本还没碰到它们。

所以我们提出一个问题：

> 能不能根据当前 ligand / fragment 状态，只挑最相关的一小部分 residues 更新？

这就是 PocketRouter 的核心。

---

## 2. PocketRouter 的直觉

### 2.1 大白话解释

想象 protein pocket 是一个很大的工作台，ligand 是正在拼装的零件。

原来的做法类似：

> 每装一个零件，都把整张工作台重新调整一遍。

PocketRouter 的想法是：

> 当前零件碰到工作台哪里，就只精修那一小块区域。其他地方先别动，避免越调越乱。

这里的“当前零件”就是当前 ligand / fragment 状态。

这里的“一小块区域”就是 top-k selected residues。

这里的“精修”就是 protein residue 的 rigid transform / side-chain chi update。

### 2.2 更研究一点的说法

我们把 pocket 看作一个 memory bank：

- 全部 residues = memory bank；
- 当前 ligand 状态 = query；
- router 根据 query 选出 top-k residues；
- 只对这些 residues 做高精度 pocket update；
- 其他 residues 保持不动，或者以后再用更轻量的全局上下文更新。

注意：这不是简单“加 attention”。

真正的创新点应该是：

> state-conditioned sparse pocket update

也就是：

> 根据当前生成状态，动态决定 pocket 哪些局部区域该被激活。

---

## 3. 这轮实验的总体设计

我们没有直接训练一个新模型。

这轮做的是 pilot validation，分两层：

1. 数据层验证；
2. 模型层采样验证。

### 3.1 数据层验证在问什么

问题：

> 如果只允许一小部分 residues 从 apo 变到 holo，能不能解释大部分 pocket 变化？

这里不跑生成模型，只看数据本身。

做法：

- 找 hard test cases；
- 对每个 case 计算 apo/holo residue movement；
- 用不同策略选 top-k residues；
- 只把这些 residues 替换成 holo，其他保持 apo；
- 看 replay RMSD 能降多少。

大白话：

> 假设我们提前知道该动哪些位置，只动这些位置，能不能把 apo pocket 修得像 holo？

如果可以，说明“稀疏更新”有物理可行性。

### 3.2 模型层验证在问什么

问题：

> 真正在 Apo2Mol 采样过程中，只更新 top-k residues，会不会比更新整个 pocket 更好？

这里跑完整 1000-step 采样。

大白话：

> 不是只看静态数据，而是真的让模型生成一次，看最终 protein RMSD 是否变好。

---

## 4. 我们改了哪些代码

### 4.1 配置文件：增加 update schedule 和 router 参数

文件：

- `configs/training.yaml`

新增配置：

```yaml
protein_update_schedule: static5
protein_update_interval: 50
protein_update_min_t: 10
protein_update_residual_threshold: 0.5
pocket_router_mode: none
pocket_router_topk: 0
```

这些参数控制两件事：

1. protein 什么时候更新；
2. 如果开启 router，每次更新时选几个 residues。

大白话：

> 我们给采样过程加了两个旋钮：什么时候让 pocket 动，以及每次让哪些 residues 动。

---

### 4.2 模型采样逻辑：支持不同 protein update schedule

文件：

- `models/molopt_score_model.py`

新增的 schedule：

| schedule | 含义 | 大白话 |
|---|---|---|
| `static5` | 原始风格，少数固定 timestep 更新 protein | 偶尔整体动一下 |
| `late_dense` | 后半程更密集更新 protein | ligand 逐渐成形后多动几次 |
| `uniform10` | 均匀 10 次 protein update | 均匀安排 10 次调整 |
| `residual_adaptive` | 根据 ligand residual 决定是否更新 | ligand 预测不稳时才动 |
| `none` | 不更新 protein | pocket 完全不动 |

为什么需要这些 schedule？

因为我们不能只比较 router 和原始 baseline。

如果 router 好了，可能不是因为“选 residue”好，而只是因为 update 次数变多了。

所以要有 `late_dense` 这种对照：

> 同样是 10 次 update，但不做 top-k，整个 pocket 都更新。

这样才能回答：

> 改善到底来自更多 update，还是来自 sparse router？

---

### 4.3 核心改动：构建 PocketRouter mask

文件：

- `models/molopt_score_model.py`

新增核心函数：

```python
_build_pocket_router_mask(...)
```

这个函数做的事情：

1. 把 protein atoms 聚合成 residue center；
2. 根据当前 ligand 位置和 router mode 给每个 residue 打分；
3. 每个样本选 top-k residues；
4. 返回一个 boolean mask：
   - `True`：这个 residue 允许更新；
   - `False`：这个 residue 本轮不动。

大白话：

> 先给每个 residue 排个队，然后只让前 k 个参与这次 pocket 更新。

---

### 4.4 Router mode 有哪些

这轮实现了几种 router：

| router mode | 怎么选 residues | 它验证什么 |
|---|---|---|
| `none` | 不选 top-k，所有 residues 都可更新 | dense baseline |
| `distance` | 选离当前 ligand 最近的 residues | contact/locality 是否有用 |
| `random` | 随机选 top-k | 负控制，看是不是随便选也行 |
| `motion_oracle` | 选真实 apo-to-holo 位移最大的 residues | 最会动的 residues 是否最有用 |
| `contact_oracle` | 选 holo 中最靠近 ligand 的 residues | 真实 holo contact 是否有用 |
| `contact_change_oracle` | 选 apo 到 holo contact 变化最大的 residues | induced-fit contact change 是否有用 |
| `predicted_motion` | 选模型预测 translation 最大的 residues | 未来 learned/predicted router 的雏形 |

这轮模型层主要跑了：

- `distance`
- `motion_oracle`
- `random`

为什么需要 random？

> 如果 random 也一样好，那说明 top-k 机制本身可能只是正则化，不说明 router 真的选对了。

为什么需要 motion_oracle？

> 如果真实最会动的 residues 都没明显更好，那就不能把论文故事讲成“找到 mobile residues”。

---

### 4.5 真正 sparse update 是怎么做的

在每个 protein update timestep，模型本来会预测每个 residue 的：

- translation；
- rotation；
- side-chain chi。

我们没有改网络结构，也没有让模型只预测 top-k。

我们做的是：

1. 模型照常预测所有 residues；
2. router 选出 top-k residues；
3. 对没选中的 residues，把预测更新清零：
   - translation = 0；
   - chi = 0；
   - rotation = identity；
4. 再调用原来的 protein transform。

大白话：

> 模型虽然给所有 residues 都提了修改建议，但我们只采纳 top-k 的建议，其他建议先忽略。

这样做的好处：

- 改动小；
- 不需要重新训练；
- 可以快速验证 sparse update 有没有潜力。

缺点：

- router 不是训练出来的；
- hard mask 比较粗糙；
- 后续如果要成论文方法，需要 learned router 和更平滑的 gating。

---

### 4.6 采样脚本改动

文件：

- `sample_split.py`

主要改动：

1. 支持 `device=auto`；
2. 支持 CPU / MPS / CUDA；
3. 支持只跑指定 test ids；
4. 支持 `init_center_mode=apo`；
5. 保存 `router_selected_counts`。

为什么要支持指定 test ids？

> 因为完整 test set 很大，本地 CPU 跑 1000-step 非常慢。我们先挑 hard case 做 pilot。

为什么要支持 `init_center_mode=apo`？

> 如果用 holo center 初始化，会给模型偷看答案的味道。真实 apo-only SBDD 更应该从 apo pocket 出发。

为什么保存 `router_selected_counts`？

> 为了确认实验真的只更新了 top-12，而不是代码没生效。

结果里 `router_selected_counts = 10 x 12`，说明：

- 一共 10 次 protein update；
- 每次选 12 个 residues；
- router 确实生效。

---

### 4.7 实验编排脚本

文件：

- `validation/run_new_method_ab.py`

这个脚本负责：

1. 选择 hard cases；
2. 为每个 arm 生成 config；
3. 跑 sampling；
4. 跑 eval；
5. 汇总 `results.json`；
6. 写 preflight/report。

这轮 router validation arms：

| arm | protein update | router |
|---|---|---|
| `baseline_realistic_static5` | static5 | none |
| `control_realistic_late_dense` | late_dense | none |
| `pocket_router_distance_top12` | late_dense | distance top12 |
| `pocket_router_motion_oracle_top12` | late_dense | motion_oracle top12 |
| `pocket_router_random_top12` | late_dense | random top12 |

大白话：

> 我们不是只跑一个新方法，而是跑了一组对照，看看 improvement 到底来自哪里。

---

## 5. 数据层验证：具体怎么做

文件：

- `validation/validate_pocket_router.py`

输出：

- `validation/pocket_router/hard8_data/pocket_router_validation.md`
- `validation/pocket_router/hard8_data/pocket_router_validation.json`

命令大致是：

```bash
.venv310/bin/python validation/validate_pocket_router.py \
  --run-dir validation/pocket_router/hard8_data \
  --num-cases 8 \
  --topk 4,8,12,16 \
  --random-trials 64
```

### 5.1 选了哪些 cases

选的是 test set 里 apo/holo RMSD 最大的一批 hard cases。

前 8 个包括：

| test position | original index | metadata RMSD |
|---:|---:|---:|
| 477 | 24500 | 4.0188 |
| 327 | 24350 | 3.9616 |
| 390 | 24413 | 3.3143 |
| 310 | 24333 | 2.9210 |
| 377 | 24400 | 2.7470 |
| 342 | 24365 | 2.6896 |
| 365 | 24388 | 2.6712 |
| 347 | 24370 | 2.4535 |

大白话：

> 不挑简单题，先挑 pocket 变化最大的难题。

### 5.2 数据层指标是什么意思

#### motion coverage

选中的 residues 覆盖了多少真实 apo-to-holo movement。

大白话：

> 真正该动的幅度里，有多少落在我们选中的 residues 上？

#### mobile recall

能不能抓住最 mobile 的那批 residues。

大白话：

> 最会动的 residues，我们有没有选中？

#### holo-contact recall

能不能抓住 holo 状态下真正接触 ligand 的 residues。

大白话：

> ligand 最后真正碰到的 residues，我们有没有选中？

#### replay RMSD

只把选中的 residues 替换成 holo，其他保持 apo，看最终距离 holo 多远。

大白话：

> 如果只修这些 residues，pocket 能变好多少？

### 5.3 数据层关键结果

hard8 mean apo atom RMSD：

```text
1.7027 A
```

top12 结果：

| router | motion coverage | mobile recall | holo-contact recall | replay RMSD |
|---|---:|---:|---:|---:|
| distance | 0.212 | 0.255 | 0.897 | 1.3415 |
| motion_oracle | 0.644 | 1.000 | 0.139 | 0.6426 |
| contact_oracle | 0.180 | 0.198 | 1.000 | 1.3569 |
| contact_change_oracle | 0.473 | 0.488 | 0.302 | 0.9225 |
| random | 0.228 | 0.225 | 0.227 | 1.3115 |

### 5.4 数据层结论

数据层说明两件事：

1. 稀疏 residue set 的确有可能解释 apo-to-holo pocket adaptation；
2. 但是单纯静态 replay 不足以证明方法有效，因为 random replay 也不差。

特别注意：

- `motion_oracle` replay RMSD 很好，说明少数 mobile residues 的确解释了很多结构变化；
- `distance` holo-contact recall 很高，说明 ligand 附近 residues 选得准；
- 但 `distance` replay RMSD 和 `random` 接近，所以还必须跑真实模型采样。

大白话：

> 静态看，确实存在“只动少数 residues 就能修好 pocket”的可能。但到底生成时有没有用，必须上模型跑。

---

## 6. 模型层实验：具体怎么跑

输出目录：

- `validation/ab_runs/hard1_steps1000_ab_gated`

结果文件：

- `validation/ab_runs/hard1_steps1000_ab_gated/results.json`
- `validation/pocket_router/model_validation_hard1.md`

### 6.1 为什么只先跑 hard1

完整 1000-step CPU 采样很慢。

一个 arm 约 6-7 分钟。

这轮有 5 个 arms，因此跑一次完整 hard1 已经要比较久。

我们先用最难 case 做 pilot：

- test position：477
- original index：24500
- metadata apo/holo RMSD：4.0188
- ligand：`8pyx__2__1.B__1.J/1.J.sdf`

大白话：

> 先拿最难的一题验方向，不一上来就烧完整 benchmark。

### 6.2 为什么必须跑 1000-step

当前 sampler 是 reverse diffusion chain。

如果 `num_steps < 1000`，它不是完整从噪声走到最终结构，而是只走高噪声尾段。

所以短步数只能用来 smoke test，看代码通不通。

不能用短步数下科研结论。

大白话：

> 短步数相当于只看考试开头几分钟，不能代表最终答案。

### 6.3 这轮采样设置

| 设置 | 值 |
|---|---|
| num_steps | 1000 |
| num_samples | 1 |
| batch_size | 1 |
| init center | apo |
| atom number | prior |
| docking mode | none |
| selected case | test position 477 |

late-dense / router arms 的 protein update timestep：

```text
460, 410, 360, 310, 260, 210, 160, 110, 60, 10
```

router top-k：

```text
12 residues per update
```

---

## 7. 模型层每个 arm 在验证什么

### 7.1 `baseline_realistic_static5`

含义：

> 原始风格，固定少数几次 protein update，没有 router。

作用：

> 作为最基本 baseline。

### 7.2 `control_realistic_late_dense`

含义：

> 10 次 late-stage protein update，但每次整个 pocket 都能更新。

作用：

> 控制 update 次数。

如果这个 arm 比 baseline 好很多，说明“多更新几次”本身就有用。

如果 router 比这个还好，才说明“选 residues”有额外价值。

### 7.3 `pocket_router_distance_top12`

含义：

> 10 次 late-stage protein update，但每次只更新离当前 ligand 最近的 top-12 residues。

作用：

> 验证 ligand/contact locality 是否是有效 sparse route。

大白话：

> ligand 附近的 residues 优先动。

### 7.4 `pocket_router_motion_oracle_top12`

含义：

> 10 次 late-stage protein update，但每次只更新真实 apo-to-holo 位移最大的 top-12 residues。

注意：

这个 arm 用了 holo 信息，所以不是现实可用方法，而是 oracle / upper-bound control。

作用：

> 验证“最会动的 residues”是不是最应该被更新。

大白话：

> 如果提前知道哪些 residues 最会动，只动它们会不会最好？

### 7.5 `pocket_router_random_top12`

含义：

> 随机选 top-12 residues。

作用：

> 负控制。

如果 random 也很好，说明 sparse mask 可能只是普通正则化，不一定说明 router 真会选。

---

## 8. 模型层结果

最终结果：

| arm | protein RMSD | TM-score | validity/recon | router usage |
|---|---:|---:|---|---:|
| baseline_realistic_static5 | 3.1261 | 0.8148 | recon 1.0, complete 1.0 | n/a |
| control_realistic_late_dense | 3.1191 | 0.8130 | recon 1.0, complete 1.0 | n/a |
| pocket_router_distance_top12 | 2.2258 | 0.8517 | recon 1.0, complete 1.0 | 10 x 12 |
| pocket_router_motion_oracle_top12 | 3.0552 | 0.8199 | recon 1.0, complete 1.0 | 10 x 12 |
| pocket_router_random_top12 | 3.0563 | 0.8160 | recon 1.0, complete 1.0 | 10 x 12 |

Best arm：

```text
pocket_router_distance_top12
```

关键差值：

| comparison | protein RMSD delta |
|---|---:|
| distance_top12 vs static5 | -0.9003 A |
| distance_top12 vs late_dense | -0.8933 A |
| motion_oracle_top12 vs late_dense | -0.0639 A |
| random_top12 vs late_dense | -0.0628 A |

解释：

- `distance_top12` 明显优于 dense update；
- `motion_oracle_top12` 只小幅优于 dense；
- `random_top12` 和 `motion_oracle_top12` 几乎一样；
- 所以当前结果不支持“只选最 mobile residues”；
- 当前结果支持“ligand/contact-near sparse update”。

大白话：

> 真正大幅变好的是“动 ligand 附近的 residues”，不是“动真实最会动的 residues”。

---

## 9. 为什么 distance router 反而比 motion oracle 好

这个现象很重要。

直觉上，motion oracle 用了真实 apo-to-holo 信息，好像应该最好。

但模型层不是静态 replay。

模型采样时，protein update 是为了帮助当前 ligand denoising 稳定下来。

当前 ligand 最需要的是：

> 它附近的 pocket 不要乱、不冲突、能配合当前接触。

而不是：

> 整个 apo-to-holo 位移最大的 residues 一定要马上动。

有些最 mobile residues 可能：

- 离当前 ligand 还远；
- 是全局构象变化的一部分；
- 当前 fragment 阶段还不该动；
- 被 hard mask 更新后反而扰乱当前局部几何。

所以这轮实验提示：

> generation-time routing 应该更关注当前 ligand/contact context，而不是静态 apo-holo displacement。

大白话：

> 最会动的地方，不一定是当前最该动的地方。当前最该动的，是 ligand 正在碰的地方。

---

## 10. 这轮实验支持什么

支持：

1. Sparse pocket update 有真实信号；
2. 只更新 top-12 residues 可以比 dense update 更好；
3. ligand/contact-conditioned routing 是一个值得继续推的方向；
4. 这个方向比简单 update schedule 更像一个主创新点。

可以形成的 paper idea：

> Fragment-Guided Contact-Sparse Pocket Memory for Apo-to-Holo Molecular Generation

核心 claim 可以是：

> 当前 ligand fragment 不需要激活整个 pocket，而是动态路由到 contact-relevant pocket memory，并只在该局部执行高精度 induced-fit update。

---

## 11. 这轮实验不支持什么

不支持：

1. 不支持“已经可以投顶会”；
2. 不支持“mobile residue routing 就是核心”；
3. 不支持“oracle sparse 一定比 learned/distance 强”；
4. 不支持“只要 sparse 就一定好”，因为 random 也有轻微改善；
5. 不支持 ligand property / docking 指标上的强结论，因为这轮样本太少，docking 也关了。

大白话：

> 有苗头，但还不是论文证据。现在只能说方向值得继续，不是已经赢了。

---

## 12. 当前证据等级

### 12.1 强证据

- 完整 1000-step，不是短步 smoke test；
- realistic apo initialization；
- 有 dense update control；
- 有 random top-k control；
- 有 motion oracle control；
- hard case 上 distance top12 明显改善。

### 12.2 弱点

- 模型层只有 n=1；
- 没有 learned router；
- 没有 hard-tail 平均和显著性；
- ligand docking / PoseBusters / clash 指标还没系统评估；
- distance router 是手工规则，不够顶会方法化；
- 当前 sparse update 是 hard mask，不够优雅。

---

## 13. 如果要往顶会推进，下一步怎么做

### 13.1 先扩实验

必须把 hard1 扩到 hard8：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --router-validation \
  --run-dir validation/ab_runs/hard8_steps1000_router \
  --num-cases 8 \
  --num-samples 1 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none
```

目标：

> 看 distance_top12 是否在多个 hard cases 上稳定优于 late_dense 和 random。

### 13.2 加 `contact_change_oracle_top12`

当前模型层还没跑 contact-change oracle。

应该新增 arm：

```text
pocket_router_contact_change_oracle_top12
```

它验证：

> 真实 apo-to-holo contact 变化最大的 residues 是否比 distance-only 更好。

如果 contact-change oracle 明显更好，说明 learned router 应该学习 contact-change，而不是只学距离。

### 13.3 训练 learned router

最终不能靠 hand-crafted distance rule。

顶会方法需要 learned router。

可学习目标可以包括：

- 当前 ligand-residue distance/contact；
- apo-to-holo contact change；
- residue displacement；
- residue type；
- ligand atom/fragment feature；
- denoising timestep；
- clash/residual proxy。

大白话：

> 让模型自己学“当前该动哪几个 residues”，而不是我们手写最近距离规则。

### 13.4 必须超过哪些 baseline

learned router 至少要超过：

- original static5；
- late_dense；
- random_top12；
- distance_top12；
- motion_oracle_top12；
- uniform10；
- residual_adaptive。

尤其是：

> learned router 必须超过 distance_top12，否则主创新不够强。

如果 learned router 不能超过 distance_top12，那论文最多只能说：

> 一个强 hand-crafted sparse induced-fit baseline。

这不太像顶会主贡献。

---

## 14. 推荐的论文叙事

不要这样讲：

> We add sparse attention to Apo2Mol.

这个太泛，也太像工程 patch。

应该这样讲：

> We introduce fragment-guided contact-sparse pocket memory for realistic apo-only SBDD. During generation, each partial ligand state routes to a small set of contact-relevant pocket memories, enabling high-fidelity induced-fit updates only where the current fragment demands them.

中文大意：

> 我们提出 fragment 引导的 contact-sparse pocket memory。生成过程中，每个 ligand 中间状态只激活与当前接触相关的小部分 pocket memory，从而只在真正需要的局部区域进行高精度 induced-fit 更新。

主创新点：

> state-conditioned sparse pocket routing for fragment-conditioned apo-to-holo adaptation

大白话：

> 药物分子长到哪，pocket 就重点修哪里。

---

## 15. 当前文件索引

### 15.1 想法与结论

- `validation/top_conf_single_idea.md`
- `validation/pocket_router/model_validation_hard1.md`
- `validation/pocket_router/model_validation_hard2.md`
- `validation/pocket_router/topk_sweep_hard2.md`
- `validation/pocket_router/experiment_process_plain.md`

### 15.2 数据层验证

- `validation/validate_pocket_router.py`
- `validation/pocket_router/hard8_data/pocket_router_validation.md`
- `validation/pocket_router/hard8_data/pocket_router_validation.json`

### 15.3 模型层 A/B

- `validation/run_new_method_ab.py`
- `validation/ab_runs/hard1_steps1000_ab_gated/results.json`
- `validation/ab_runs/hard1_steps1000_ab_gated/preflight_report.md`
- `validation/ab_runs/hard2_steps1000_router_v2/results.json`
- `validation/ab_runs/hard2_steps1000_router_v2/preflight_report.md`
- `validation/ab_runs/hard2_distance_topk_sweep_steps1000/results.json`
- `validation/ab_runs/hard2_distance_topk_sweep_steps1000/preflight_report.md`

### 15.4 核心代码改动

- `configs/training.yaml`
- `models/molopt_score_model.py`
- `sample_split.py`

---

## 16. Hard2 扩展实验

为了继续往顶会方向推进，我们又跑了第二轮模型层 A/B。

这轮不是只看一个 hard case，而是看两个 test split 里 apo-holo RMSD 最大的 hard cases：

| test position | original index | metadata apo-holo RMSD | ligand path |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |

运行目录：

`validation/ab_runs/hard2_steps1000_router_v2`

新增关键 arm：

`pocket_router_contact_change_oracle_top12`

它的作用是测试：

> 如果我们知道 apo 到 holo 哪些 residue contact 变化最大，只更新这些 residues 会不会最好？

结果如下：

| arm | mean protein RMSD | mean TM-score | delta vs static5 | delta vs late-dense |
|---|---:|---:|---:|---:|
| `baseline_realistic_static5` | 2.8381 | 0.8483 | 0.0000 | +0.0514 |
| `control_realistic_late_dense` | 2.7867 | 0.8496 | -0.0514 | 0.0000 |
| `pocket_router_distance_top12` | 2.1731 | 0.8764 | -0.6650 | -0.6135 |
| `pocket_router_motion_oracle_top12` | 2.6050 | 0.8596 | -0.2331 | -0.1817 |
| `pocket_router_contact_change_oracle_top12` | 2.8072 | 0.8493 | -0.0309 | +0.0206 |
| `pocket_router_random_top12` | 2.5913 | 0.8574 | -0.2468 | -0.1953 |

大白话读法：

- `distance_top12` 还是最强，说明“更新 ligand 附近 residues”这个信号复现了；
- `late_dense` 几乎没带来收益，说明不是“后期多更新整个 pocket”就行；
- `random_top12` 有一定提升，说明 sparse 本身确实能减少一些噪声；
- 但 random 明显不如 distance，说明“选对 residues”仍然重要；
- `motion_oracle_top12` 不够强，说明“最终动得最大”的 residues 不是当前生成过程中最该更新的 residues；
- `contact_change_oracle_top12` 这轮不稳，说明静态 ground-truth contact-change 标签不一定等于当前 denoising step 最该更新的局部区域。

所以 hard2 后，结论更精确了：

> 当前最有信号的不是 mobile-residue routing，也不是 contact-change oracle routing，而是 current-fragment distance/contact-conditioned sparse pocket update。

一句更直白的话：

> ligand 现在靠近哪里，就重点修哪里；不要提前去修那些“最终可能会变”的远端区域。

## 17. 更新后的最终结论

这轮实验给出的最稳妥结论是：

> PocketRouter 方向有正信号，但应该从“mobile residue sparse memory”进一步收紧成“current-fragment contact-conditioned sparse pocket adaptation”。

具体证据是：

- hard1 中，`distance_top12` 把 protein RMSD 从 3.1191 降到 2.2258；
- hard2 中，`distance_top12` 在两个 hard cases 上平均 protein RMSD 为 2.1731，优于 static5 的 2.8381 和 late-dense 的 2.7867；
- 改善不是来自更多 dense update，因为 `late_dense` 只比 static5 好 0.0514 A；
- 改善也不只是随便 sparse，因为 `random_top12` 虽然有提升，但仍明显差于 `distance_top12`；
- 改善不像来自“选最 mobile residues”，因为 `motion_oracle_top12` 接近 random；
- `contact_change_oracle_top12` 当前不稳，说明最终 apo-holo contact-change 标签不能直接当作每一步 denoising 的最佳 update 区域。

一句大白话：

> 不要让整个 pocket 每次都跟着 ligand 乱动；让 ligand 当前旁边真正相关的 residues 动，效果反而更好。

这可以作为后续顶会投稿方向的核心种子，但还需要：

1. hard8 / hard-tail 多样本验证；
2. top-k sweep，比如 top4、top8、top12、top16、top24；
3. learned router，并且必须打过 best hand-crafted distance top-k，目前 hard2 上是 `distance_top4`；
4. docking / clash / validity / QED / SA 系统指标；
5. 与现有 SBDD / apo-to-holo baselines 的完整比较。

最关键的一句话：

> 如果 learned router 不能稳定打过 best distance-only router，这个工作还不能作为顶会主方法；但 distance-only router 已经是一个很强的发现，可以作为设计 learned router 的核心老师和强 baseline。

## 18. Hard2 Top-k Sweep

继续往顶会方向推进时，我们又补了一个关键实验：

> distance router 到底应该选几个 residues？

运行目录：

`validation/ab_runs/hard2_distance_topk_sweep_steps1000`

比较：

- `distance_top4`
- `distance_top8`
- `distance_top12`
- `distance_top16`
- `distance_top24`

结果：

| distance router | mean protein RMSD | mean TM-score | mean selected residues |
|---|---:|---:|---:|
| `top4` | 2.0624 | 0.8796 | 4 |
| `top8` | 2.1050 | 0.8784 | 8 |
| `top12` | 2.1731 | 0.8764 | 12 |
| `top16` | 2.2407 | 0.8745 | 16 |
| `top24` | 2.3321 | 0.8695 | 24 |

最重要的现象：

> top4 最好，并且 top-k 越大，平均 protein RMSD 越差。

大白话：

> 不是让更多 pocket residues 一起动更好，而是只让 ligand 当前最近的极小局部动更好。

这对论文很有价值，因为它说明：

- PocketRouter 不是随便调 top12；
- 稀疏预算本身有规律；
- 当前最像“机制”的说法是 very sparse current-contact pocket update；
- learned router 的门槛变高：必须打过 `distance_top4`，不能只打过 random 或 top12。

更新后的最强结论：

> 药物分子长到哪，pocket 就只修最近的一小块；修太多反而会把无关残基拉进来，增加噪声。

## 19. Hard8 Core Top4 验证

在 hard2 top-k sweep 之后，我们继续跑了一个更接近顶会证据链的 hard8 core 实验。

详细记录见：

`validation/pocket_router/hard8_core_top4.md`

运行目录：

`validation/ab_runs/hard8_core_top4_steps1000`

这次实验选择 8 个 hardest apo-to-holo cases，每个 arm 都跑完整 1000-step 采样。

比较的 arms：

| arm | 作用 |
|---|---|
| `baseline_realistic_static5` | 原始风格，5 次整体 pocket update |
| `control_realistic_late_dense` | 后半程 10 次整体 pocket update |
| `pocket_router_random_top4` | 10 次 update，每次随机只动 4 个 residues |
| `pocket_router_distance_top4` | 10 次 update，每次只动当前 ligand 最近 4 个 residues |
| `pocket_router_distance_top8` | 10 次 update，每次只动当前 ligand 最近 8 个 residues |

核心结果：

| arm | mean protein RMSD | mean TM-score | selected residues/update |
|---|---:|---:|---:|
| `baseline_realistic_static5` | 2.4365 | 0.8977 | 52.25 |
| `control_realistic_late_dense` | 2.4142 | 0.9005 | 52.25 |
| `pocket_router_random_top4` | 1.8305 | 0.9248 | 4 |
| `pocket_router_distance_top4` | 1.7924 | 0.9249 | 4 |
| `pocket_router_distance_top8` | 1.8371 | 0.9229 | 8 |

最强结果：

> `distance_top4` 把 mean protein RMSD 从 static5 的 2.4365 A 降到 1.7924 A，相对 late-dense 也低 0.6218 A。

大白话：

> 在 8 个最难的 pocket 变化案例上，还是“只修 ligand 当前最近的一小块”最好；整体多修 pocket 没有明显帮助。

这次实验也暴露了一个很重要的问题：

> `random_top4` 也很强，RMSD = 1.8305 A，只比 `distance_top4` 差 0.0381 A。

这说明 sparse mask 本身就是强正则化，不能只拿“超过 dense baseline”当作顶会证据。

更新后的投稿判断：

- 可以继续推进 PocketRouter；
- 主创新点应该收紧为 current-fragment contact-sparse pocket memory；
- learned router 必须打过 `distance_top4` 和 `random_top4`；
- ligand quality、QED、SA、docking 还没有补齐，本轮只支持 protein adaptation 结论。

一句最准确的当前结论：

> PocketRouter 现在有一个很强的物理和模型层信号：dense pocket update 会带来噪声，very sparse current-contact update 更稳；但要到顶会标准，还必须证明 learned router 比强手工 distance top4 更好，并补齐 ligand/docking 指标。

## 20. Hard8 Full Router Suite 完整实验

完整记录见：

`validation/pocket_router/hard8_full_router_suite.md`

运行目录：

`validation/ab_runs/hard8_full_router_steps1000_n1`

这次我们把 hard8 上所有关键 router/control 都补齐了。

总共 14 个 arms，每个 arm 都是：

- 8 个 hardest cases；
- 每个 case 1 个 sample；
- 1000-step sampling；
- realistic apo initialization；
- docking disabled。

最终检查：

> 14 arms x 8 cases = 112 个 `result_*.pt`，全部完成。

完整结果：

| arm | mean protein RMSD | mean TM-score | complete |
|---|---:|---:|---:|
| `baseline_realistic_static5` | 2.4365 | 0.8977 | 0.875 |
| `control_realistic_late_dense` | 2.4142 | 0.9005 | 0.875 |
| `control_realistic_uniform10` | 3.1475 | 0.8514 | 0.875 |
| `adaptive_realistic_residual` | 2.0904 | 0.9141 | 0.875 |
| `pocket_router_random_top4` | 1.8305 | 0.9248 | 0.625 |
| `pocket_router_random_top12` | 1.9555 | 0.9196 | 1.000 |
| `pocket_router_motion_oracle_top12` | 2.0511 | 0.9136 | 0.750 |
| `pocket_router_contact_oracle_top12` | 2.1195 | 0.9065 | 0.750 |
| `pocket_router_contact_change_oracle_top12` | 2.2410 | 0.9037 | 0.875 |
| `pocket_router_distance_top4` | 1.7924 | 0.9249 | 0.750 |
| `pocket_router_distance_top8` | 1.8371 | 0.9229 | 0.750 |
| `pocket_router_distance_top12` | 1.8992 | 0.9206 | 0.750 |
| `pocket_router_distance_top16` | 1.9534 | 0.9180 | 0.875 |
| `pocket_router_distance_top24` | 2.0531 | 0.9140 | 0.750 |

最强结果：

> `distance_top4` 最好，mean protein RMSD = 1.7924 A。

它比：

- `static5` 低 0.6441 A；
- `late_dense` 低 0.6218 A；
- `random_top4` 低 0.0381 A；
- `distance_top8` 低 0.0446 A。

大白话：

> 在 hard8 上，还是“只修 ligand 当前最近的一小块 pocket”最好。修更多 residues 不是更好，反而越来越差。

distance top-k sweep 的趋势非常清楚：

| distance top-k | mean protein RMSD |
|---:|---:|
| 4 | 1.7924 |
| 8 | 1.8371 |
| 12 | 1.8992 |
| 16 | 1.9534 |
| 24 | 2.0531 |

这说明：

> PocketRouter 的关键不是“多更新 pocket”，而是“极小局部更新”。

但这次实验也把顶会门槛抬高了：

> `random_top4` 也很强，只比 `distance_top4` 差 0.0381 A。

所以不能只说：

> 我们比 dense baseline 好。

真正要说服顶会，需要后续 learned router 同时打过：

- dense baselines；
- `random_top4`；
- `distance_top4`。

另外，ligand 侧指标还没完成：

- 本轮 `docking-mode none`；
- `evaluated_mols = 0`；
- QED / SA 是 nan；
- 所以这轮只能支持 protein adaptation 结论，不能直接声称完整 SBDD 性能提升。

更新后的最准确结论：

> PocketRouter 的机制信号已经很清楚：dense pocket update 有噪声，very sparse current-contact update 更稳。下一步要做 learned fragment-conditioned router，并且必须超过 `distance_top4` 和 `random_top4`，同时补齐 ligand validity、docking、QED/SA。

---

## 21. Active-set shell pilot：从“只动 top-k”升级到“核心强动、邻居弱动”

完整记录：

`validation/pocket_router/active_set_shell_pilot.md`

运行目录：

`validation/ab_runs/hard2_active_set_shell_pilot_steps1000_n1`

这轮实验回应的是一个很关键的担心：

> 只更新 selected residues，其他 residues 完全不动，会不会造成局部结构不连续，甚至让 ligand contact 变差？

所以我们把 hard mask 改成 soft active set：

```text
selected core residues: weight = 1.0
nearby shell residues: weight = 0.25 or 0.50
background residues: weight = 0.0
```

大白话：

> 不再是“选中就动、没选中就死死固定”，而是“最相关的地方正常动，旁边一小圈轻轻跟着松，远处背景先稳住”。

本轮 hard2 结果：

| arm | mean protein RMSD | mean TM-score | updated residues/update | mol stable |
|---|---:|---:|---:|---:|
| `control_realistic_late_dense` | 2.7867 | 0.8496 | 48.00 | 0.5000 |
| `pocket_router_random_top4` | 2.1252 | 0.8795 | 4.00 | 0.5000 |
| `pocket_router_distance_top4_hard` | 2.0624 | 0.8796 | 4.00 | 1.0000 |
| `active_set_distance_top4_shell4_w025` | 2.0641 | 0.8796 | 5.95 | 1.0000 |
| `active_set_distance_top4_shell6_w025` | 2.0806 | 0.8794 | 17.80 | 0.5000 |
| `active_set_distance_top4_shell6_w050` | 2.0950 | 0.8802 | 17.50 | 0.5000 |
| `active_set_random_top4_shell6_w025` | 2.1625 | 0.8786 | 15.50 | 0.5000 |

怎么理解：

- `shell4_w025` 几乎不损伤 hard `distance_top4`，RMSD 只差 0.0017 A；
- `shell6` 反而变差一点，说明邻居不是越多越好；
- `random shell6` 比 `distance shell6` 更差，说明不是“多放松几个 residues 就行”，还是要动对地方；
- 这轮 `evaluated_mols = 0`，所以不能下 ligand 质量结论，只能支持 protein adaptation 侧判断。

更新后的设计判断：

> 顶会方向不应该停留在“只更新部分 residues”，而应该表述为 ligand-conditioned active-set pocket optimization：当前 ligand 释放核心 pocket 自由度，小范围邻居弱松弛，背景保持稳定。

当时下一步最该补：

- hard8 active-set shell suite；
- 每个 case 多采样 n=3 / n=5；
- selected-unselected boundary clash / contact / strain 指标；
- ligand validity、docking、QED/SA、contact recovery。

---

## 22. Hard8 active-set shell suite：回答“局部不连续”和“ligand 被挤坏”问题

完整记录：

`validation/pocket_router/hard8_active_set_shell_suite.md`

运行目录：

`validation/ab_runs/hard8_active_set_shell_only_steps1000_n1`

这轮是在上一节 hard2 pilot 后继续做的 hard8 扩展。

它专门回答三个问题：

1. 只让 selected residues 动，会不会导致旁边 unselected residues 不连续？
2. 加一圈 weak shell 能不能缓解这个问题？
3. protein RMSD 变好时，会不会反而伤害 ligand geometry / binding contact？

本轮只跑 4 个关键 active-set arms：

| arm | 解释 |
|---|---|
| `active_set_distance_top4_shell4_w025` | ligand-near top4 core 强更新，4A 邻居 shell 以 0.25 权重弱松弛 |
| `active_set_distance_top4_shell6_w025` | 同上，但 shell 扩到 6A |
| `active_set_distance_top4_shell6_w050` | 6A shell，但 shell 权重提高到 0.50 |
| `active_set_random_top4_shell6_w025` | 随机 core + 6A shell，作为负控制 |

主结果：

| arm | protein RMSD | TM-score | updated residues/update | mol stable | atom stable | recon | complete |
|---|---:|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell4_w025` | 1.8589 | 0.9213 | 6.11 | 0.875 | 0.9675 | 0.875 | 0.750 |
| `active_set_distance_top4_shell6_w025` | 1.8666 | 0.9217 | 16.38 | 0.375 | 0.9747 | 1.000 | 1.000 |
| `active_set_distance_top4_shell6_w050` | 1.9301 | 0.9192 | 15.51 | 0.625 | 0.9892 | 1.000 | 0.875 |
| `active_set_random_top4_shell6_w025` | 1.8524 | 0.9249 | 15.75 | 0.250 | 0.8556 | 0.875 | 0.875 |

局部几何诊断：

| arm | ligand-protein clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell4_w025` | 0.25 | 2.3097 | 0.7228 | 0.7185 | 0.5501 | 0.7971 |
| `active_set_distance_top4_shell6_w025` | 0.25 | 2.3752 | 0.7430 | 0.7237 | 0.5777 | 0.4775 |
| `active_set_distance_top4_shell6_w050` | 0.375 | 2.3211 | 0.7495 | 0.7501 | 0.5957 | 0.5336 |
| `active_set_random_top4_shell6_w025` | 9.50 | 1.3783 | 0.7121 | 0.7423 | 0.5742 | 0.3965 |

怎么读：

- `shell4_w025` 是最稳的 active-set 版本，mol stable 最高，clash 很低；
- `shell6_w025` 边界最平滑，contact recall 也略高，但 mol stable 掉到 0.375；
- `shell6_w050` 说明 shell 权重不能太大；
- `random_shell6` 是最重要的负例：protein RMSD 看起来很好，但 ligand-protein clash 平均 9.50，min distance 只有 1.3783 A。

大白话：

> 只看 protein RMSD 会被 random shell 骗。它能把 protein 调得像 holo，但 ligand 被挤坏了。真正的顶会方法必须证明 active set 选得对，而不是只证明 pocket RMSD 低。

这轮给出的实验结论：

> active-set shell 方向成立，但默认应该是小 shell、弱权重、强背景锚定。更大的 shell 可以缓解边界不连续，却可能伤 ligand 稳定性；随机 shell 会造成严重 clash。

更新后的顶会主张：

> apo-to-holo pocket adaptation 应该被建模为 ligand-conditioned active-set optimization：当前 ligand 激活核心 pocket 自由度，邻域弱松弛，背景稳定；并用 clash/contact/ligand validity 约束 active set 是否真的合理。

还不能 claim 的内容：

- 本轮 `docking-mode none`；
- `evaluated_mols = 0`；
- QED / SA / docking score 仍不可用；
- 所以现在只能支持 protein adaptation + local geometry 结论，不能说完整 SBDD 性能已经提升。

下一步实验：

- hard8 multi-sample：n=3 或 n=5；
- 小范围 shell sweep：3A/4A/5A，weight 0.10/0.25；
- core top-k + shell 联合 sweep：top3/top4/top6；
- 开启 ligand-side evaluation：docking、QED/SA、PoseBusters-style checks；
- 训练 learned router，并要求它打过 hard `distance_top4`、`random_top4` 和 active-set shell4。

---

## 23. Hard8 conservative shell sweep：把 shell 调得更保守

完整记录：

`validation/pocket_router/hard8_conservative_shell_sweep.md`

运行目录：

`validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1`

上一轮 shell suite 说明：

> 小 shell 安全，大 shell 可以改善边界和 contact，但可能伤 ligand 稳定性。

所以这一轮没有继续加大 shell，而是反过来做保守 sweep：

```text
shell radius: 3A, 4A, 5A
shell weight: 0.10, 0.25
core: distance top4
background: anchored, weight 0.0
```

新增命令：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --active-set-conservative-shell-sweep \
  --run-dir validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1 \
  --num-cases 8 \
  --num-samples 1 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none
```

本轮完成：

```text
4 arms x 8 cases x 1 sample = 32 result_*.pt
```

主结果：

| arm | protein RMSD | TM-score | updated residues/update | mol stable | atom stable | recon | complete |
|---|---:|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell3_w010` | 1.8571 | 0.9215 | 4.33 | 0.750 | 0.9928 | 1.000 | 0.875 |
| `active_set_distance_top4_shell3_w025` | 1.8584 | 0.9214 | 4.31 | 0.750 | 0.9928 | 1.000 | 0.875 |
| `active_set_distance_top4_shell4_w010` | 1.8512 | 0.9218 | 6.10 | 0.750 | 0.9495 | 0.875 | 0.750 |
| `active_set_distance_top4_shell5_w025` | 1.8563 | 0.9219 | 9.98 | 0.625 | 0.9386 | 0.875 | 0.750 |

几何诊断：

| arm | ligand-protein clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell3_w010` | 0.250 | 2.2840 | 0.7342 | 0.7142 | 0.5588 | 0.8375 |
| `active_set_distance_top4_shell3_w025` | 0.250 | 2.2759 | 0.7342 | 0.7189 | 0.5611 | 0.8018 |
| `active_set_distance_top4_shell4_w010` | 0.250 | 2.2827 | 0.7366 | 0.7082 | 0.5534 | 0.7814 |
| `active_set_distance_top4_shell5_w025` | 0.375 | 2.2842 | 0.7537 | 0.7270 | 0.5811 | 0.7830 |

怎么读：

- `shell3_w010` 和 `shell3_w025` 是这轮最均衡的 conservative variants；
- `shell4_w010` protein RMSD 最低，但 complete 下降到 0.750；
- `shell5_w025` contact recall 最高，但 mol stable 下降到 0.625，clash 升到 0.375；
- 更大 shell 能保留更多 contact，但不是免费收益；
- 固定 shell 参数没有全赢，说明后续应该训练 learned active-set policy。

大白话：

> 这轮证明了“轻一点、小一点”的 shell 更像安全缓冲垫。它不能保证每个指标都最好，但能避免大 shell 把 ligand 稳定性弄坏。真正的论文方法应该让模型自己判断：哪里强动、哪里轻动、哪里不动。

当前最准确结论：

> conservative shell sweep 支持 ligand-conditioned active-set optimization，但固定半径/固定权重只能作为 baseline。顶会方向要继续推进到 learned router：预测 core release、shell relaxation 和 background anchoring，并且用 ligand stability、clash、contact、docking 一起验证。

---

## 24. Hard8 n=3 candidate repeat：降低采样偶然性

完整记录：

`validation/pocket_router/hard8_candidate_repeat_n3.md`

运行目录：

`validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3`

上一轮 conservative shell sweep 是 n=1，也就是每个 case 只生成一次。n=1 可以看方向，但容易被随机性影响。

所以这轮把最关键的 6 个候选方案拿出来，每个 hard case 生成 3 次：

```text
6 arms x 8 cases x 3 samples = 144 generated samples
48 / 48 result files completed
```

比较对象：

| arm | 大白话 |
|---|---|
| `pocket_router_random_top4` | 随便选 4 个 residues，看 sparse 本身有多强 |
| `pocket_router_distance_top4_hard` | 只强更新 ligand 最近的 top4 residues |
| `active_set_distance_top4_shell4_w025` | top4 强更新，4A 邻居轻轻动 25% |
| `active_set_distance_top4_shell3_w010` | top4 强更新，3A 邻居轻轻动 10% |
| `active_set_distance_top4_shell3_w025` | top4 强更新，3A 邻居轻轻动 25% |
| `active_set_distance_top4_shell4_w010` | top4 强更新，4A 邻居轻轻动 10% |

主指标：

| arm | protein RMSD | TM-score | mol stable | complete |
|---|---:|---:|---:|---:|
| `pocket_router_random_top4` | 1.8153 | 0.9253 | 0.2500 | 0.7083 |
| `pocket_router_distance_top4_hard` | 1.8148 | 0.9245 | 0.5000 | 0.6667 |
| `active_set_distance_top4_shell4_w025` | 1.8207 | 0.9243 | 0.5000 | 0.7500 |
| `active_set_distance_top4_shell3_w010` | 1.8171 | 0.9244 | 0.5000 | 0.6667 |
| `active_set_distance_top4_shell3_w025` | 1.8145 | 0.9245 | 0.4583 | 0.6250 |
| `active_set_distance_top4_shell4_w010` | 1.8168 | 0.9244 | 0.4583 | 0.7500 |

几何诊断：

| arm | clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump |
|---|---:|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 9.7500 | 1.3000 | 0.6592 | 0.6821 | 0.5024 | 0.6222 |
| `pocket_router_distance_top4_hard` | 1.7500 | 1.7124 | 0.6880 | 0.7131 | 0.5349 | 0.6143 |
| `active_set_distance_top4_shell4_w025` | 1.3750 | 1.7239 | 0.6886 | 0.7557 | 0.5540 | 0.5167 |
| `active_set_distance_top4_shell3_w010` | 1.7500 | 1.7170 | 0.6964 | 0.7155 | 0.5414 | 0.6124 |
| `active_set_distance_top4_shell3_w025` | 1.7500 | 1.7073 | 0.6964 | 0.7210 | 0.5437 | 0.6097 |
| `active_set_distance_top4_shell4_w010` | 1.8750 | 1.6569 | 0.6940 | 0.7346 | 0.5487 | 0.6023 |

怎么读：

- 只看 protein RMSD，所有 distance-based 方法都很接近，差距小到不能作为强结论；
- `random_top4` 的 RMSD 看起来也不错，但 clash 高到 9.75，说明 sparse 本身会骗过 protein RMSD；
- `distance_top4_hard` 是强 baseline，但 complete 只有 0.6667，clash 也不是最低；
- `shell4_w025` 的 RMSD 稍差一点，但 clash 最低、contact precision/Jaccard 最高、boundary jump 最低；
- `shell3` 更保守，contact recall 稍高，但没有明显降低 clash；
- `shell4_w010` 边界稍平滑，但 clash 变差，说明 shell 不是越平滑越好。

大白话：

> 这轮最重要的结论不是“谁 RMSD 第一”，而是“protein RMSD 已经不够用了”。真正合理的 active-set 方法必须同时看：protein 有没有修好、ligand 有没有被挤坏、contact 有没有保留、core 和 background 之间有没有硬断层。

更新后的顶会方向：

> 固定 active-set shell 可以作为强 baseline，但不应该作为最终方法。最终方法应该训练 learned active-set policy，让模型根据当前 ligand 状态预测 core release、neighbor relaxation 和 background anchoring。

仍然没解决的限制：

- 本轮 `docking-mode none`；
- `evaluated_mols = 0`；
- QED / SA / docking score 仍不可用；
- 所以下一步必须补 ligand-side evaluation 或直接进入 learned router 训练并把这些指标纳入验证。

补充更新：

我又加了一个不依赖 docking 的轻量 ligand-quality 评估脚本：

`validation/analyze_ligand_quality.py`

它已经跑完 hard8 n=3：

`validation/ab_runs/hard8_active_set_candidate_repeat_steps1000_n3/ligand_quality.json`

结果：

| arm | complete | QED | SA | LogP | Lipinski |
|---|---:|---:|---:|---:|---:|
| `pocket_router_random_top4` | 0.7083 | 0.4441 | 0.5147 | 4.7667 | 4.2941 |
| `pocket_router_distance_top4_hard` | 0.6667 | 0.4856 | 0.5025 | 4.8934 | 4.2500 |
| `active_set_distance_top4_shell4_w025` | 0.7500 | 0.5017 | 0.5106 | 4.8415 | 4.3333 |
| `active_set_distance_top4_shell3_w010` | 0.6667 | 0.4802 | 0.5025 | 5.0589 | 4.1875 |
| `active_set_distance_top4_shell3_w025` | 0.6250 | 0.4973 | 0.5080 | 4.4610 | 4.3333 |
| `active_set_distance_top4_shell4_w010` | 0.7500 | 0.4889 | 0.4956 | 4.7371 | 4.3333 |

大白话：

> `shell4_w025` 不只是 clash/contact 更均衡，QED 也是这组最高，complete 也是并列最高。它目前最适合作为 fixed active-set baseline。真正 docking 还没跑，因为本机还缺 `vina/qvina`。
