# Hard8 Conservative Active-Set Shell Sweep

日期：2026-05-29

运行目录：

`validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1`

结构化结果：

- `validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1/results.json`
- `validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1/geometry_diagnostics.json`

这轮实验接着上一轮 hard8 active-set shell suite 做。它专门回答一个更细的问题：

> 如果只更新 core residues 可能造成局部不连续，而大 shell 又可能伤 ligand，那么更保守的小半径、低权重 shell 是不是更合理？

大白话：

> 之前我们试了“核心强更新，旁边一圈弱松弛”。结果发现小 shell 安全，大 shell 不一定好。这轮就是把 shell 调得更保守，看看能不能在 protein RMSD、ligand 稳定性、clash、contact 之间找到更好的平衡。

---

## 1. 为什么要做这轮

用户提出的关键担心是对的：

- selected residues 动了，unselected residues 不动，会不会局部结构不连续；
- 只看 ligand 附近 residues，可能漏掉二级影响；
- 局部 pocket 更新可能改善 protein RMSD，但损害 ligand geometry 或 binding contact；
- 只优化 protein RMSD 可能把 ligand 挤坏。

上一轮 hard8 shell-only 已经给出两个信号：

| 现象 | 说明 |
|---|---|
| `shell4_w025` mol stable 最高 | 小 shell 确实更安全 |
| `shell6_w025` boundary jump 更低，但 mol stable 大幅下降 | 更连续不等于更好 |
| `random_shell6_w025` protein RMSD 看着不错，但 clash 很高 | 只看 protein RMSD 会被误导 |

所以这轮不继续盲目放大 shell，而是转向更保守的 active-set：

```text
core residues: 强更新
near shell: 更小半径、更低权重弱松弛
background: 锚住
```

大白话：

> 核心区域该动就动；旁边一圈只轻轻让位；远处不要乱动。我们不是追求“动更多”，而是追求“该动的动，别的稳住”。

---

## 2. 代码改动

本轮只改实验编排脚本：

`validation/run_new_method_ab.py`

新增函数：

```python
make_active_set_conservative_shell_sweep_arms()
```

新增命令行参数：

```text
--active-set-conservative-shell-sweep
```

它生成 4 个更保守的 active-set arms：

| arm | core router | core top-k | shell radius | shell weight | background weight |
|---|---|---:|---:|---:|---:|
| `active_set_distance_top4_shell3_w010` | distance | 4 | 3A | 0.10 | 0.0 |
| `active_set_distance_top4_shell3_w025` | distance | 4 | 3A | 0.25 | 0.0 |
| `active_set_distance_top4_shell4_w010` | distance | 4 | 4A | 0.10 | 0.0 |
| `active_set_distance_top4_shell5_w025` | distance | 4 | 5A | 0.25 | 0.0 |

固定设置：

| 设置 | 值 |
|---|---|
| `protein_update_schedule` | `late_dense` |
| `protein_update_interval` | 50 |
| `protein_update_min_t` | 10 |
| `protein_update_residual_threshold` | 0.5 |
| `sample_num_atoms` | `prior` |
| `init_center_mode` | `apo` |
| `num_steps` | 1000 |

大白话：

> 每次模型提出 pocket update 时，我们不再让所有 pocket residues 同等接受更新。离当前 ligand 最近的 top4 residue 完整接受更新，旁边一小圈只接受 10% 或 25% 的更新，背景不接受更新。

---

## 3. 运行命令

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

总采样量：

```text
4 arms x 8 hard cases x 1 sample = 32 result_*.pt
```

最终确认：

```text
32 result_*.pt
```

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

---

## 4. 主指标结果

| arm | result files | mean protein RMSD | mean TM-score | updated residues/update | mol stable | atom stable | recon | complete | evaluated mols |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell3_w010` | 8 | 1.8571 | 0.9215 | 4.33 | 0.750 | 0.9928 | 1.000 | 0.875 | 0 |
| `active_set_distance_top4_shell3_w025` | 8 | 1.8584 | 0.9214 | 4.31 | 0.750 | 0.9928 | 1.000 | 0.875 | 0 |
| `active_set_distance_top4_shell4_w010` | 8 | 1.8512 | 0.9218 | 6.10 | 0.750 | 0.9495 | 0.875 | 0.750 | 0 |
| `active_set_distance_top4_shell5_w025` | 8 | 1.8563 | 0.9219 | 9.98 | 0.625 | 0.9386 | 0.875 | 0.750 | 0 |

大白话读法：

> `shell3_w010` 和 `shell3_w025` 是这轮最均衡的：protein RMSD 还在 1.86 A 左右，分子稳定性 0.75，reconstruction 1.0，complete 0.875。`shell5_w025` 放开更多邻居后，contact 更强一点，但分子稳定性掉了。

按 case 看 protein RMSD 最优的 arm：

| original index | best arm | best RMSD |
|---:|---|---:|
| 24333 | `active_set_distance_top4_shell3_w010` | 2.0989 |
| 24350 | `active_set_distance_top4_shell4_w010` | 2.0842 |
| 24365 | `active_set_distance_top4_shell3_w010` | 1.5050 |
| 24370 | `active_set_distance_top4_shell4_w010` | 0.9490 |
| 24388 | `active_set_distance_top4_shell5_w025` | 1.6921 |
| 24400 | `active_set_distance_top4_shell5_w025` | 1.3724 |
| 24413 | `active_set_distance_top4_shell4_w010` | 3.0318 |
| 24500 | `active_set_distance_top4_shell3_w010` | 2.0414 |

胜率：

| arm | RMSD wins / 8 cases |
|---|---:|
| `active_set_distance_top4_shell3_w010` | 3 |
| `active_set_distance_top4_shell3_w025` | 0 |
| `active_set_distance_top4_shell4_w010` | 3 |
| `active_set_distance_top4_shell5_w025` | 2 |

大白话：

> 从 paired case 看，没有哪个 conservative shell 稳赢所有 case。`shell3_w010` 和 `shell4_w010` 各赢 3 个 case，`shell5_w025` 赢 2 个 case，但它的 ligand 稳定性更差。这再次说明固定规则不是终点。

---

## 5. 几何诊断

诊断命令：

```bash
.venv310/bin/python validation/analyze_geometry_diagnostics.py \
  --run-dir validation/ab_runs/hard8_active_set_conservative_shell_sweep_steps1000_n1
```

诊断结果：

| arm | RMSD | ligand-protein clashes | min dist | contact recall | contact precision | contact Jaccard | boundary jump | active disp | background disp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell3_w010` | 1.8571 | 0.250 | 2.2840 | 0.7342 | 0.7142 | 0.5588 | 0.8375 | 1.0538 | 0.3460 |
| `active_set_distance_top4_shell3_w025` | 1.8584 | 0.250 | 2.2759 | 0.7342 | 0.7189 | 0.5611 | 0.8018 | 0.9771 | 0.3559 |
| `active_set_distance_top4_shell4_w010` | 1.8512 | 0.250 | 2.2827 | 0.7366 | 0.7082 | 0.5534 | 0.7814 | 0.8248 | 0.3344 |
| `active_set_distance_top4_shell5_w025` | 1.8563 | 0.375 | 2.2842 | 0.7537 | 0.7270 | 0.5811 | 0.7830 | 1.0130 | 0.3617 |

指标解释：

| 指标 | 大白话 |
|---|---|
| ligand-protein clashes | ligand 有没有被 pocket 挤到 |
| min dist | ligand 和 protein 最近距离，太小通常危险 |
| contact recall | 参考 holo 里该接触的地方找回了多少 |
| contact precision | 生成 contact 里有多少不是乱碰 |
| contact Jaccard | contact 整体重合度 |
| boundary jump | active 区域和 background 之间是否有突变 |
| active/background disp | 被释放区域和背景分别动了多少 |

关键现象：

- 所有保守 distance-shell arms 的 ligand-protein clash 都很低，前三个是 0.25；
- `shell5_w025` contact recall 最高，0.7537；
- 但 `shell5_w025` 的 clash 升到 0.375，mol stable 降到 0.625；
- `shell4_w010` protein RMSD 最低，1.8512，但 complete 只有 0.750；
- `shell3_w010/w025` complete 更好，0.875，recon 也都是 1.000。

大白话：

> 放大 shell 能多找回一点 binding contact，但也更容易影响 ligand 自身稳定性。小 shell 不一定 contact 最高，但更稳。

---

## 6. 和上一轮 hard8 shell-only 对照

上一轮运行目录：

`validation/ab_runs/hard8_active_set_shell_only_steps1000_n1`

关键对照：

| arm | RMSD | selected/update | mol stable | complete | clashes | contact recall | boundary jump |
|---|---:|---:|---:|---:|---:|---:|---:|
| old `shell4_w025` | 1.8589 | 6.11 | 0.875 | 0.750 | 0.250 | 0.7228 | 0.7971 |
| old `shell6_w025` | 1.8666 | 16.38 | 0.375 | 1.000 | 0.250 | 0.7430 | 0.4775 |
| old `shell6_w050` | 1.9301 | 15.51 | 0.625 | 0.875 | 0.375 | 0.7495 | 0.5336 |
| new `shell3_w010` | 1.8571 | 4.33 | 0.750 | 0.875 | 0.250 | 0.7342 | 0.8375 |
| new `shell3_w025` | 1.8584 | 4.31 | 0.750 | 0.875 | 0.250 | 0.7342 | 0.8018 |
| new `shell4_w010` | 1.8512 | 6.10 | 0.750 | 0.750 | 0.250 | 0.7366 | 0.7814 |
| new `shell5_w025` | 1.8563 | 9.98 | 0.625 | 0.750 | 0.375 | 0.7537 | 0.7830 |

怎么读：

1. 如果只看 mol stable，旧 `shell4_w025` 仍然最好，0.875。
2. 如果看 complete 和 recon，新 `shell3_w010/w025` 更均衡，complete 0.875，recon 1.000。
3. 如果看 contact recall，新 `shell5_w025` 最高，但它付出了更多 clash 和更低 mol stable。
4. 如果看 boundary jump，旧 `shell6_w025` 最平滑，但 ligand 稳定性最差。

大白话：

> 没有一个固定 shell 参数全赢。小 shell 更稳，大 shell 更会找 contact，但更容易伤 ligand。这个结果支持“训练一个会预测 active-set 权重的模型”，而不是手写一个固定半径。

---

## 7. 回答几个关键担心

### 7.1 只更新部分 residues 会不会造成其他 residues 冲突？

会有这个风险，但不是无解。

本轮结果显示，distance-conditioned 小 shell 的 ligand-protein clash 很低：

```text
shell3_w010: 0.25
shell3_w025: 0.25
shell4_w010: 0.25
```

这说明：

> 只要 core 选得对，并且 shell 不乱放大，局部更新不一定导致严重 clash。

但上一轮 random shell 也说明：

> 如果 active set 选错，protein RMSD 可以好看，但 ligand 会被挤坏。

### 7.2 selected residues 动了，unselected residues 不动，会不会不连续？

会，所以需要 shell。

但是 shell 的目标不是越大越好，而是：

```text
core: 大幅适配当前 ligand
shell: 小幅吸收 core/background 边界差
background: 保持结构稳定
```

本轮 `shell3/4` 给出了一个更稳的折中，上一轮 `shell6` 说明过度追求边界平滑会伤 ligand。

### 7.3 局部 pocket 更新会不会改善 protein RMSD，但损害 ligand geometry？

会。上一轮 random shell 是最明显负例。

本轮保守 distance shell 的意义是：

> 不只看 protein RMSD，而是同时看 mol stable、clash、contact、complete。

这也是后续论文必须强调的评估方式。

### 7.4 只看 ligand 附近 residues，会不会忽略二级影响？

可能会。

但实验说明，直接扩大 shell 不是最好答案。更合理的是：

```text
直接接触区域: 高置信强更新
一阶邻居: 低权重弱松弛
二级/远程区域: learned router 判断是否需要释放
背景区域: 默认锚住
```

大白话：

> 二级影响不是不要管，而是不能靠“半径越大越好”来管。应该让模型学会哪条传导链真的需要释放。

---

## 8. 当前结论

这轮实验支持：

> fragment-conditioned active-set optimization 是比 hard sparse update 更合理的顶会方向，但固定 shell 规则只能当 baseline，不能当最终方法。

更具体：

- `shell3_w010/w025` 是当前最均衡的保守方案；
- `shell4_w010` protein RMSD 最低，但 ligand complete 略差；
- `shell5_w025` contact recall 最高，但 mol stable 和 clash 变差；
- 旧 `shell4_w025` mol stable 单项最好，但 complete 较低；
- 大 shell 确实能改善 contact/boundary，但会带来 ligand 代价；
- 后续不能只用 protein RMSD 选模型，必须同时看 ligand 几何。

一句大白话：

> 我们已经看到“怎么动 pocket”比“动不动 pocket”更关键：核心要敢动，邻居要轻动，背景要稳住；shell 过大或选错区域，ligand 会付代价。

---

## 9. 对训练版的直接启发

训练版不应该把 `shell3` 或 `shell4` 写死成最终规则。更应该学三个东西：

1. Core release score

```text
当前 ligand 真的需要哪个 residue 强更新？
```

2. Shell relaxation weight

```text
哪些邻居只需要轻轻松一下？权重是 0.1、0.25，还是更小？
```

3. Background anchor

```text
哪些 residues 不能被当前 fragment 带着乱动？
```

建议的训练目标：

| loss / metric | 目的 |
|---|---|
| ligand diffusion / reconstruction loss | 保住生成主任务 |
| pocket displacement / rotamer loss | 学会 apo-to-holo adaptation |
| contact recovery loss | 保留真实 binding contact |
| ligand-protein clash penalty | 防止 pocket 把 ligand 挤坏 |
| active-set sparsity budget | 防止模型动太多 residues |
| background stability loss | 防止远处背景漂移 |
| boundary smoothness loss | 防止 core/background 断层 |

大白话：

> 不是让模型学“最近的几个 residues 动”，而是学“哪些必须动、哪些轻轻动、哪些绝不能动”。

---

## 10. 下一步实验建议

### 10.1 先做 n=3 或 n=5 repeat

候选 arms 不要太多，建议：

```text
pocket_router_distance_top4
active_set_distance_top4_shell4_w025
active_set_distance_top4_shell3_w010
active_set_distance_top4_shell3_w025
active_set_distance_top4_shell4_w010
```

目的：

> n=1 只能看方向，n=3/n=5 才能看均值、方差、胜率。

### 10.2 加 ligand-side evaluation

目前 `--docking-mode none`，所以：

```text
evaluated_mols = 0
QED / SA / docking = nan 或未评估
```

后续必须补：

- docking score；
- QED / SA；
- PoseBusters-style geometry；
- bond/angle validity；
- ligand-protein clash；
- contact precision/recall；
- per-case failure analysis。

### 10.3 再进入训练

建议训练从轻量版开始：

```text
distance/contact features + residue features + denoising time
  -> predict core probability
  -> predict shell relaxation weight
  -> regularize active-set size and background motion
```

如果 learned router 不能超过：

```text
distance_top4
random_top4
safe shell baseline
```

那就还不能作为顶会主方法。

---

## 11. 一句话结论

> Conservative shell sweep 说明：更保守的小 shell 能在低 clash、较好 ligand 稳定性和较强 protein adaptation 之间取得平衡；但固定 shell 参数没有全局最优，真正顶会方向应转向 learned fragment-conditioned active-set：核心强释放、邻居弱松弛、背景锚定，并用 ligand geometry/contact/clash 指标证明 active set 选得对。
