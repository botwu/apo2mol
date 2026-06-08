# Active-Set Shell Pilot

日期：2026-05-28

运行目录：

`validation/ab_runs/hard2_active_set_shell_pilot_steps1000_n1`

本轮实验的核心问题是：

> 如果把 apo-to-holo pocket adaptation 看成 ligand-conditioned active-set optimization，是否应该从“只更新被选中的 residues”升级成“核心 residues 强更新，邻居 residues 弱松弛，背景 residues 保持稳定”？

大白话：

> ligand 当前碰到 pocket 哪一块，哪一块就应该重点动；它旁边的一圈可以轻轻跟着松一下，避免结构断裂；更远的背景不要乱动，避免整个 pocket 被带偏。

---

## 1. 为什么要做这轮实验

前面 hard2 / hard8 实验证明了一件事：

> very sparse pocket update 很强，尤其是 `distance_top4`。

但它也留下一个合理担心：

> 如果只让 top-k residues 动，其他 residues 完全不动，会不会造成局部结构不连续？比如 selected residues 往 holo 方向动了，但旁边 unselected residues 被固定住，两者之间产生冲突。

这就是 active-set shell pilot 要回答的问题。

我们没有直接训练新模型，而是先做一个采样阶段的 algorithmic ablation。原因很简单：

> 如果这个机制在不训练时就完全不稳定，那不值得投入训练；如果它稳定，并且指标有信号，再把它做成 learned router。

---

## 2. 本轮改动是什么

### 2.1 从 hard mask 变成 soft active set

之前的 router 是 hard mask：

```text
selected residues: update
unselected residues: do not update
```

现在改成 active-set weights：

```text
core selected residues: weight = 1.0
neighbor shell residues: weight = 0.25 or 0.50
background residues: weight = 0.0
```

大白话：

> 以前是开关：选中就动，没选中就完全不动。
>
> 现在是旋钮：最相关的 residue 正常动，旁边的 residue 轻轻动，远处不动。

### 2.2 代码层改动

主要文件：

- `models/molopt_score_model.py`
- `configs/training.yaml`
- `validation/run_new_method_ab.py`

新增配置：

```yaml
pocket_router_shell_radius: 0.0
pocket_router_shell_weight: 0.0
pocket_router_background_weight: 0.0
```

含义：

| 参数 | 作用 | 大白话 |
|---|---|---|
| `pocket_router_shell_radius` | selected residues 周围多少 A 内算邻居 shell | 旁边多大一圈可以跟着轻轻动 |
| `pocket_router_shell_weight` | shell residues 的更新权重 | 邻居跟着动多少 |
| `pocket_router_background_weight` | 背景 residues 的更新权重 | 远处背景是否允许动 |

采样时的更新逻辑：

- residue translation 乘以 active-set weight；
- side-chain chi update 乘以 active-set weight；
- residue rotation 用 slerp 向 identity 缩回，weight 越小，旋转越接近不动；
- `router_selected_counts` 统计 weight > 0 的 residues 数量。

大白话：

> 平移、侧链角、旋转都不是硬切了，而是按权重缩放。核心 residue 完整执行模型预测，shell residue 只吃一小部分预测，背景基本锚住。

---

## 3. 实验设计

本轮还是 hard2 pilot：

- hard cases: 2 个；
- 每个 arm 每个 case 生成 1 个样本；
- 1000-step sampling；
- realistic apo initialization；
- docking disabled；
- docking / QED / SA 本轮不作为有效结论；
- 目标先看 protein adaptation 和基础稳定性。

启动命令：

```bash
.venv310/bin/python validation/run_new_method_ab.py \
  --run \
  --active-set-ablation \
  --run-dir validation/ab_runs/hard2_active_set_shell_pilot_steps1000_n1 \
  --num-cases 2 \
  --num-samples 1 \
  --batch-size 1 \
  --num-steps 1000 \
  --docking-mode none
```

本轮 7 个 arms：

| arm | 目的 |
|---|---|
| `control_realistic_late_dense` | dense update baseline，全 pocket 更新 |
| `pocket_router_random_top4` | 随机 hard top4，对照 sparse 本身是否有用 |
| `pocket_router_distance_top4_hard` | 当前最强 hard router baseline |
| `active_set_distance_top4_shell4_w025` | distance top4 + 4A shell，shell 权重 0.25 |
| `active_set_distance_top4_shell6_w025` | distance top4 + 6A shell，shell 权重 0.25 |
| `active_set_distance_top4_shell6_w050` | distance top4 + 6A shell，shell 权重 0.50 |
| `active_set_random_top4_shell6_w025` | random top4 + 6A shell，对照 shell 是否只是多动 residues |

---

## 4. 完整结果

共完成：

> 7 arms x 2 cases = 14 个 `result_*.pt`

结果表：

| arm | mean protein RMSD | mean TM-score | updated residues/update | complete | mol stable | atom stable | evaluated mols |
|---|---:|---:|---:|---:|---:|---:|---:|
| `control_realistic_late_dense` | 2.7867 | 0.8496 | 48.00 | 1.0000 | 0.5000 | 0.9844 | 0 |
| `pocket_router_random_top4` | 2.1252 | 0.8795 | 4.00 | 1.0000 | 0.5000 | 0.9844 | 0 |
| `pocket_router_distance_top4_hard` | 2.0624 | 0.8796 | 4.00 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `active_set_distance_top4_shell4_w025` | 2.0641 | 0.8796 | 5.95 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `active_set_distance_top4_shell6_w025` | 2.0806 | 0.8794 | 17.80 | 1.0000 | 0.5000 | 0.9844 | 0 |
| `active_set_distance_top4_shell6_w050` | 2.0950 | 0.8802 | 17.50 | 1.0000 | 0.5000 | 0.9844 | 0 |
| `active_set_random_top4_shell6_w025` | 2.1625 | 0.8786 | 15.50 | 1.0000 | 0.5000 | 0.9844 | 0 |

相对 control / hard baseline 的变化：

| arm | RMSD vs dense control | RMSD vs hard distance top4 |
|---|---:|---:|
| `pocket_router_random_top4` | -0.6615 | +0.0628 |
| `pocket_router_distance_top4_hard` | -0.7242 | 0.0000 |
| `active_set_distance_top4_shell4_w025` | -0.7225 | +0.0017 |
| `active_set_distance_top4_shell6_w025` | -0.7060 | +0.0182 |
| `active_set_distance_top4_shell6_w050` | -0.6916 | +0.0326 |
| `active_set_random_top4_shell6_w025` | -0.6241 | +0.1001 |

---

## 5. 怎么读这些结果

### 5.1 active-set soft shell 是能跑通的

所有 soft shell arms 都完整生成，没有 shape、rotation、broadcast 或数值崩溃问题。

这说明：

> active-set weighting 作为算法机制是可运行的。

大白话：

> 不是一改成“部分强动、邻居弱动”模型就直接炸了。这个方向可以继续做。

### 5.2 shell4 几乎不损伤 hard top4

对比：

| arm | RMSD | TM-score | mol stable |
|---|---:|---:|---:|
| `pocket_router_distance_top4_hard` | 2.0624 | 0.8796 | 1.0000 |
| `active_set_distance_top4_shell4_w025` | 2.0641 | 0.8796 | 1.0000 |

差别：

- RMSD 只差 0.0017 A；
- TM-score 几乎相同；
- mol stable / atom stable 没变差。

解释：

> 加一个很窄、很弱的 shell，基本不会破坏当前最强 hard top4 结果。

大白话：

> 如果我们担心 hard top4 太硬，给它旁边加一小圈缓冲垫是安全的，至少这 2 个 case 没看到明显副作用。

### 5.3 shell6 开始变差

对比：

| arm | updated residues/update | RMSD |
|---|---:|---:|
| hard top4 | 4.00 | 2.0624 |
| shell4 w0.25 | 5.95 | 2.0641 |
| shell6 w0.25 | 17.80 | 2.0806 |
| shell6 w0.50 | 17.50 | 2.0950 |

趋势：

> shell 半径越大，参与更新的 residues 越多，RMSD 反而稍微变差。

这和 hard8 full suite 的结论一致：

> 多更新 residues 不一定更好，很多额外 residues 可能只是带来 update noise。

大白话：

> pocket 不是越多人一起动越协调。拉进来的人太多，反而会把局部精修带乱。

### 5.4 random shell 更差，说明 route 还是有意义

对比：

| arm | RMSD | selected/update |
|---|---:|---:|
| `active_set_distance_top4_shell6_w025` | 2.0806 | 17.80 |
| `active_set_random_top4_shell6_w025` | 2.1625 | 15.50 |

两者 selected 数量差不多，但 distance shell 更好。

解释：

> 不是随便多放松几个 residues 就行；从 ligand 当前位置出发选 active set 仍然更合理。

大白话：

> 不是“多动一点就好”，而是“动对地方才好”。

### 5.5 本轮还不能支持 ligand 侧 claim

本轮 `evaluated_mols = 0`，因为 docking disabled 且 eval 没有产生可用的 QED / SA / docking 统计。

所以这轮只能说：

> active-set shell 对 protein adaptation 的影响。

不能说：

> active-set shell 已经改善完整 SBDD 质量。

大白话：

> 现在只能说 pocket 修得怎么样，还不能说药物分子本身整体更好了。

---

## 6. 对用户提出的几个担心的回答

### 6.1 只更新部分 residues 会不会造成其他 residues 冲突？

可能会。

所以本轮加入 shell 的目的，就是给 selected residues 周围一圈邻居一个弱松弛自由度。

实验看到：

- 4A shell 没有损伤 hard top4；
- 6A shell 稍微变差；
- 因此最合理的默认不是“完全 hard mask”，也不是“大范围 shell”，而是“小半径、低权重 shell”。

当前建议：

> 训练版可以保留 soft shell，但初始设置应该偏保守，比如 shell radius 4A、weight 0.1-0.25。

### 6.2 selected residues 动了，unselected residues 不动，会不会局部不连续？

可能会，但需要用更直接的 clash/contact/bonded-neighbor metric 来测。

本轮只通过 protein RMSD/TM 和基础稳定性间接看，没有做局部几何连续性检查。

后续应补：

- selected-unselected boundary distance change；
- side-chain clash count；
- residue-residue contact break；
- ligand-protein steric clash；
- pocket local strain proxy。

大白话：

> RMSD 看的是整体像不像 holo，但看不出局部有没有挤在一起。这个必须补专门的冲突指标。

### 6.3 局部 pocket 更新会不会改善 protein RMSD，但损害 ligand geometry 或 binding contact？

可能会。

本轮没有 docking / QED / SA / PoseBusters，所以不能排除这个风险。

后续必须补：

- ligand validity；
- reconstruction success；
- QED / SA；
- docking score；
- ligand-protein contact recovery；
- steric clash；
- PoseBusters-style checks。

大白话：

> pocket 看起来更像 holo，不代表 ligand 就真的放得更好。药物侧指标必须补。

### 6.4 只看 ligand 附近 residues，会不会忽略二级影响？

可能会。

本轮 shell6 的结果提醒我们：

> 二级影响不能简单靠“扩大半径”解决。

更合理的做法可能是：

- core active set：当前 fragment 直接相关 residues；
- shell：小范围弱松弛；
- background：用低频、低强度 global relaxation，而不是每步高精度 dense update；
- learned router：预测哪些远端 residues 虽然离 ligand 不近，但会通过 contact network 间接受影响。

大白话：

> 远处影响确实可能重要，但不能把一大片 pocket 都强行拉进来一起动。更好的方式是让模型学会“哪条传导链真的重要”。

---

## 7. 当前结论

本轮最准确结论：

> active-set soft shell 是可运行的；小 shell 基本安全；大 shell 会引入轻微噪声；distance-conditioned active set 仍优于 random active set。

更具体：

- hard top4 仍是 hard2 上最强结果；
- shell4 w0.25 几乎等同于 hard top4，提供了一个“更物理、更连续”的安全变体；
- shell6 w0.25 / w0.50 没有带来收益，反而稍微变差；
- random shell 更差，说明 ligand-conditioned / distance-conditioned routing 仍有价值；
- 本轮还不是 ligand 质量结论，因为 ligand 侧有效指标不足。

---

## 8. 对顶会方向的影响

这轮实验让 idea 更精确了。

不应该把贡献说成：

> 只更新一部分 residues。

更合理的说法是：

> ligand-conditioned active-set pocket optimization。

也就是：

```text
current ligand state
  -> select active pocket degrees of freedom
  -> strong update core residues
  -> weak relax local shell
  -> anchor background
  -> optionally perform low-frequency global relaxation
```

大白话：

> 不是把 pocket 粗暴切成“动/不动”，而是判断当前 ligand 让哪些自由度释放，哪些自由度轻轻跟着松，哪些自由度必须稳住。

这比 hard mask 更像一个顶会主张，因为它回答的是第一性原理问题：

> induced fit 不是全 pocket 同时发生，而是 ligand 当前状态激活一组局部自由度。

---

## 9. 下一步实验建议

### 9.1 立即补的实验

1. hard8 active-set shell suite

把本轮 hard2 的 4 个关键 active-set arms 扩到 hard8：

- `pocket_router_distance_top4_hard`
- `active_set_distance_top4_shell4_w025`
- `active_set_distance_top4_shell6_w025`
- `active_set_random_top4_shell6_w025`

目的：

> 看 shell4 是否在 hard8 仍然不损伤，random shell 是否仍然更差。

2. multi-sample statistics

每个 case 不只采 1 个样本，而是 n=3 或 n=5。

目的：

> 排除单次采样偶然性。

3. local conflict metrics

补 selected-unselected boundary 和 clash/contact 指标。

目的：

> 专门回答“局部不连续/冲突”问题。

4. ligand-side metrics

开启可用 docking / QED / SA / contact / clash 评估。

目的：

> 判断 protein RMSD 改善是否伤害 ligand。

### 9.2 训练前的算法定型

建议训练版采用：

```text
core: learned or distance-initialized top-k, weight 1.0
inner shell: 3-4A, weight 0.1-0.25
background: 0 or very small low-frequency relaxation
router loss: contact relevance + displacement uncertainty + sparsity budget
```

不建议一开始采用：

```text
large shell 6A+
high shell weight 0.50+
dense background update every protein update step
```

原因：

> 当前结果已经看到：扩大 active set 容易引入噪声。

---

## 10. 一句话结论

> 本轮支持把 PocketRouter 从 hard sparse update 升级为 ligand-conditioned active-set optimization，但默认应该是“小 shell、弱松弛、强背景锚定”，不是大范围放松；下一步要在 hard8/multi-sample/local clash/ligand metrics 上验证它是否真正达到顶会级证据。

---

## 11. Hard8 扩展已经完成

日期：2026-05-29

完整记录：

`validation/pocket_router/hard8_active_set_shell_suite.md`

运行目录：

`validation/ab_runs/hard8_active_set_shell_only_steps1000_n1`

这次把 hard2 pilot 里的关键 active-set arms 扩到 8 个 hardest cases，并补了 local geometry diagnostics。

主结果：

| arm | mean protein RMSD | mean TM-score | updated residues/update | mol stable | complete |
|---|---:|---:|---:|---:|---:|
| `active_set_distance_top4_shell4_w025` | 1.8589 | 0.9213 | 6.11 | 0.875 | 0.750 |
| `active_set_distance_top4_shell6_w025` | 1.8666 | 0.9217 | 16.38 | 0.375 | 1.000 |
| `active_set_distance_top4_shell6_w050` | 1.9301 | 0.9192 | 15.51 | 0.625 | 0.875 |
| `active_set_random_top4_shell6_w025` | 1.8524 | 0.9249 | 15.75 | 0.250 | 0.875 |

几何诊断：

| arm | ligand-protein clashes | min dist | contact recall | boundary jump |
|---|---:|---:|---:|---:|
| `active_set_distance_top4_shell4_w025` | 0.25 | 2.3097 | 0.7228 | 0.7971 |
| `active_set_distance_top4_shell6_w025` | 0.25 | 2.3752 | 0.7430 | 0.4775 |
| `active_set_distance_top4_shell6_w050` | 0.375 | 2.3211 | 0.7495 | 0.5336 |
| `active_set_random_top4_shell6_w025` | 9.50 | 1.3783 | 0.7121 | 0.3965 |

新的结论：

- shell4 w0.25 是当前最稳的 active-set 折中：RMSD 明显优于 dense baseline，同时 ligand-protein clash 很低，mol stability 最好；
- shell6 w0.25 让边界更平滑，也略提高 contact recall，但 mol stability 明显变差；
- shell6 w0.50 没有收益，说明 shell 不能太强；
- random shell 是关键负例：protein RMSD 不差，但 ligand-protein clash 爆炸，说明不能只看 protein RMSD；
- hard `distance_top4` 仍是 protein RMSD 最强 baseline，但 active-set shell 提供了更物理的连续性/安全性方向。

大白话：

> 加一小圈弱 shell 像给 hard top4 加缓冲垫；加太大一圈会把 ligand 稳定性拖坏；随机加 shell 会直接把 ligand 挤坏。所以顶会方向要强调“ligand-conditioned active set”，不是“多放松一点 residues”。
