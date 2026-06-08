# Apo2Mol Research Idea Memo

这份 memo 目标不是罗列灵感，而是把 Apo2Mol 放进一个可复制的 idea 生产线里：先从代码和论文证据定位真实缺口，再用三套公式生成可投稿方向。

## 0. Evidence Snapshot

本地代码和公开资料共同说明，Apo2Mol 的核心卖点是从 apo pocket 同时生成 ligand 和 holo-like pocket conformation。README 明确写到它区别于 holo-conditioned 方法，直接在 ligand-free apo 设置下建模 pocket flexibility。公开 arXiv/AAAI 摘要也强调：多数 SBDD 方法默认 rigid pocket，而 Apo2Mol 用约 24k experimental apo-holo pairs 训练 joint generation。

本地 `selected_index_apo_druglike.pkl` 有 24,501 条 apo-holo-ligand 记录；RMSD 分布显示 median 约 0.846 A，mean 约 1.001 A，但长尾很明显：RMSD >= 2 A 有 1,435 条，>= 3 A 有 362 条，>= 5 A 有 90 条。这意味着普通指标容易被“小构象变化”主导，真正有药物发现价值的 induced-fit / cryptic-pocket 场景反而是少数长尾。

代码层面有几个直接可用的切口：

- `configs/sampling.yaml` 默认 `sample_num_atoms: ref`，测试时可直接使用参考 ligand 原子数；真实 de novo 设计通常没有这个信息。
- `configs/training.yaml` 和模型里都有 `topk_prompt`，但默认是 0；`models/molopt_score_model.py` 里 prompt list 实际为空，检索增强能力没有真正用起来。
- 蛋白 pocket 在采样中只在少数几个 diffusion time steps 更新，而 ligand 每一步更新；这是一个天然的时空分解入口。
- 评估主要集中在 stability/reconstruction/Vina/QED/SA/bond distance 等，缺少 pocket 物理合理性、side-chain rotamer、clash、interaction fingerprint、PoseBusters 类检查和不确定性校准。

近年趋势也支持这些缺口：DynamicBind、FlexSBDD、DynamicFlow、FlowDock 都在把 rigid docking/generation 推向 flexible protein 或 flow matching；PoseBusters 则提醒只看 RMSD/Vina 不足以证明物理合理性。

## 1. Production-Line Formula A: Mature Hammer on Emerging Nail

### Idea A1: RA-Apo2Mol, Retrieval-Augmented Dynamic Pocket Prompting

一句话：把成熟的 RAG / case-based reasoning 用到 Apo2Mol 的 apo-to-holo transition 上，检索相似 pocket 的已知 apo-holo 变形作为 conformational prompt，指导当前 ligand-pocket joint diffusion。

为什么值钱：Apo2Mol 代码已经有 `topk_prompt` 和 prompt attention 的骨架，但没有启用。不是硬加一个新模块，而是补全作者留下的自然延长线。这个方向容易写成“memory-augmented dynamic SBDD”，故事清楚：实验 apo-holo 数据不只是训练样本，也可以成为生成时的结构记忆库。

最小实现：

1. 为每个 pocket 建索引：residue type histogram、pocket shape descriptor、apo geometry embedding、BAPNet/EGNN pooled embedding、apo-holo RMSD bucket。
2. 对测试 apo pocket 检索 top-k 相似 apo-holo pairs，取它们的 residue-level translation/rotation/chi delta 或 BAPNet hidden features。
3. 把检索结果接进当前空着的 prompt path：`protein_prompt_list` 和 `prompt_hbap_ligand_batch_all_list`。
4. 做 k=0/1/3/5/10 ablation，重点看 RMSD >= 2 A 和 >= 3 A 的 long-tail bucket。

核心实验：

- Main metrics: Vina min, QED, validity, novelty, high-affinity rate。
- Pocket metrics: generated-vs-holo RMSD/JSD、pocket volume JSD、side-chain clash、rotamer plausibility。
- Stress subset: apo-holo RMSD >= 2 A、>= 3 A、kinase/GPCR/ion-channel family split。
- Negative control: random retrieval、same-family retrieval、geometry-only retrieval、sequence-only retrieval。

投稿叙事：Apo2Mol 证明“apo-holo paired data 能训练 dynamic generator”；RA-Apo2Mol 进一步证明“apo-holo paired data 还能作为 non-parametric conformational memory，在长尾大构象变化上显著提升”。

风险：检索相似不等于变形相似。应加入 delta-consistency filter 和 uncertainty gate，避免错误 prompt 伤害简单样本。

优先级：最高。原因是创新点贴合本仓库，工程成本中等，容易做出 ablation。

### Idea A2: Foundation-Model Verifier / Preference-Guided Apo2Mol

一句话：把 AlphaFold3/Boltz/Chai 一类 biomolecular structure model 当成 verifier 或 preference model，对 Apo2Mol 生成的 protein-ligand complex 做 reranking、DPO 或 guidance。

为什么值钱：2024 以后结构生物学 foundation model 已经成为新锤子。Apo2Mol 仍主要依赖 Vina/QED/SA 等传统评价；可以把“生成器”和“结构/亲和力 verifier”分成两阶段，形成现实药物发现工作流。

最小实现：

1. Apo2Mol 生成 N 个 candidate complexes。
2. 用 Boltz-1/2、Chai-1 或可用的开源 verifier 打分结构可信度、clash、affinity/confidence。
3. 构造 preference pairs：高 verifier 分数 vs 低 verifier 分数。
4. 先做 reranking paper baseline，再做 DPO/contrastive fine-tuning。

核心实验：

- 生成 5/20/50 samples 后，比较 top-1 by Vina、top-1 by verifier、combined score。
- 检查 verifier 是否只偏好训练集相似 molecule，需要 novelty/diversity control。
- 对比 wet-lab 不可得时，用 PoseBusters、ProLIF interaction fingerprint、MMFF/UFF relaxation 后 energy 作为 secondary evidence。

风险：AF3 访问和许可证可能限制；优先选开源 verifier 或只做 post-hoc benchmark。

优先级：中高。适合做增强版工作，但对外部工具依赖较强。

### Idea A3: Physics-Guided Dynamic Pocket Diffusion

一句话：把经典物理/药化先验作为 guidance 加到 Apo2Mol 采样：steric clash、hydrogen bond、salt bridge、pi-stacking、rotamer、torsion strain、ligand internal strain。

为什么值钱：PoseBusters 已经指出 deep docking/generation 只看 RMSD 会产生不物理结构；Apo2Mol 的动态 pocket 更容易引入 side-chain clash 或 rotamer 不合理。把成熟 physics checks 放进 flexible generation，是“老锤子砸新钉子”的典型方向。

最小实现：

1. 先做无训练的 energy-guided reranking：Vina + clash + rotamer + interaction fingerprint。
2. 再做 differentiable 或 approximate guidance，只在最后 50-100 steps 施加，避免采样早期过约束。
3. 对比 Apo2Mol 原始采样、Vina-only rerank、physics-composite rerank。

风险：物理项若设计粗糙，容易只是后处理；论文贡献要落在“dynamic pocket validity”而不是普通 docking cleanup。

优先级：中高。适合与 A1 或 C1 合并，作为 robustness/value amplifier。

## 2. Production-Line Formula B: Temporal-Spatial Reframing

### Idea B1: Adaptive Asynchronous Co-Denoising

一句话：把 Apo2Mol 的 1000-step joint diffusion 改成 adaptive ligand-pocket co-denoising：ligand 高频更新，pocket 只在不确定性、clash、interaction mismatch 高的时候更新；不同 pocket residues 也按局部重要性异步更新。

为什么值钱：现有代码已经是隐式时空分解：ligand 每步动，protein 只在少数固定步更新。但固定 schedule 很粗糙，不知道哪些 target 需要大 pocket motion。把“固定 5 次更新”改成“事件触发/残基级局部更新”，既能提升速度，也能提升大构象变化样本。

最小实现：

1. 记录每步 ligand-pocket clash、distance-map error proxy、BAPNet confidence/entropy。
2. 训练一个 lightweight update controller，输出是否更新 pocket、更新哪些 residues。
3. 设计 1000->200->50 step 的 progressive distillation 或 consistency distillation。
4. 与固定更新 schedule 对比质量/速度。

核心实验：

- Wall-clock speed、GPU memory、samples/sec。
- 质量不降或在 hard subset 提升：Vina min、validity、pocket RMSD、clash。
- Ablation: fixed 5 updates vs fixed 10/20 updates vs adaptive global vs adaptive residue-local。

投稿叙事：Apo2Mol 解决“动态 pocket 能不能生成”；本工作解决“动态 pocket 如何高效、按需地生成”。这符合大规模筛选场景。

风险：速度论文需要硬指标；如果质量提升不明显，仍可定位成 efficient dynamic SBDD。

优先级：高，但工程量大于 A1。

### Idea B2: Bridge-Apo2Mol, From Linear Interpolation to Stochastic Conformational Bridges

一句话：Apo2Mol 用 apo-holo 之间的 residue-level 插值构造 pseudo trajectory；可以换成 Schrödinger bridge / optimal transport / flow matching，把 apo-to-holo 过渡建模为多路径、多模态分布。

为什么值钱：真实蛋白 conformational transition 通常不只有一条直线，尤其是 cryptic pocket 或 induced fit。DynamicFlow 走的是 MD trajectory；Apo2Mol 走的是 experimental apo-holo pair。Bridge-Apo2Mol 可以结合二者优点：不需要完整 MD，也不假设线性路径。

最小实现：

1. 先不改主网络，只改 forward process：从 deterministic interpolation 变成 stochastic bridge noise schedule。
2. 对有多个 holo ligands 的同一 apo/protein family，学习 conditional transition distribution。
3. 训练时预测 bridge velocity/score，采样时输出 pocket ensemble，而不是单一 holo-like pocket。

核心实验：

- Pocket ensemble coverage: generated pocket RMSD distribution 是否覆盖多个 holo states。
- Ligand diversity vs affinity tradeoff。
- 对 high-RMSD/多配体 protein family 单独分析。

风险：理论包装容易过重；必须控制成可实现的 bridge loss 和明确的 pocket ensemble metric。

优先级：中高。适合作为更主线的顶会方法，但时间成本较高。

### Idea B3: Two-Stage Pocket-First Then Ligand Generation

一句话：先生成一组 plausible holo-like pockets，再对每个 pocket 生成 ligand；不要强迫 ligand 和 pocket 在一个 1000-step loop 里完全同步。

为什么值钱：这是空间分治。动态 pocket 本身是一个结构预测问题，ligand generation 是另一个问题。先产生 pocket ensemble，再用成熟 SBDD generator/docking/guidance 生成 ligand，可以提升可解释性和可控性。

最小实现：

1. 从 Apo2Mol 抽出 pocket refinement branch，训练/采样 pocket ensemble。
2. 固定每个 generated pocket，用现有 ligand generator 或 Apo2Mol ligand branch 生成 ligand。
3. 做 ensemble selection：哪个 pocket state 更支持高亲和、低 clash ligand。

风险：拆开后可能损失 ligand-induced fit 的耦合优势；需要用 interaction-aware feedback loop 补回来。

优先级：中。适合走“实用系统/benchmark”路线。

## 3. Production-Line Formula C: Extreme-Pressure Method

### Idea C1: Apo2Mol-Hard, Long-Tail Dynamic Pocket Benchmark and Robust Training

一句话：围绕 Apo2Mol 数据里的大构象变化长尾构建 hard benchmark，并提出 robust training/evaluation。不要再平均所有测试样本，而是问：RMSD >= 2 A / 3 A / 5 A 时，模型还行不行？

为什么值钱：本地数据统计直接支持这个方向。24,501 条里，大部分 apo-holo RMSD 小于 2 A；平均指标容易掩盖真正难的 dynamic SBDD。Benchmark + robust baseline 是很稳的论文生产线。

最小实现：

1. 按 apo-holo RMSD 分桶：0-0.5、0.5-1、1-1.5、1.5-2、2-3、3-5、5-10 A。
2. 复现或运行 Apo2Mol checkpoint，对每桶分别报告 ligand quality、pocket quality、physical validity。
3. 训练 robust variant：long-tail reweighting、RMSD-conditioned diffusion timestep、hard-example curriculum、uncertainty-aware sampling。
4. 加入 stress perturbation：apo pocket 坐标噪声、side-chain rotamer noise、missing side chain、wrong protonation proxy、pocket center shift。

核心实验：

- Bucket-wise performance table，而不是单一平均。
- Severe shift subset case studies。
- 与 DynamicFlow/FlexSBDD/FlowDock 等可运行方法在同一 hard split 上比较。

投稿叙事：现有 dynamic SBDD 评价不区分小幅 breathing 和大幅 induced-fit；我们提出 Apo2Mol-Hard，证明长尾构象变化是当前方法的真实瓶颈，并给出 robust training baseline。

风险：如果无法跑全模型，可先做 benchmark design + partial evaluation；但最好至少跑 Apo2Mol checkpoint 和一两个 baseline。

优先级：最高。原因是数据证据强，idea 不依赖大改模型，最适合研究生快速形成可发表工作。

### Idea C2: Realistic Apo Stress Test, Remove the Reference Crutches

一句话：系统测试 Apo2Mol 在真实使用条件下的退化：没有参考 ligand 原子数、没有精确 pocket、只有 AlphaFold/低置信度 apo、side chain 缺失、pocket residue 选择错误。

为什么值钱：当前 sampling 默认 `sample_num_atoms: ref`。这对 benchmark 可以理解，但真实 de novo design 没有 reference ligand size。只要证明 ref atom count 对结果影响很大，就能形成一个强问题；再提出 joint atom-count/pocket-volume prior 就是方法贡献。

最小实现：

1. 对比 `sample_num_atoms=ref` vs `prior` vs learned atom-count predictor。
2. 学一个 pocket-volume + residue-composition conditioned atom-count distribution。
3. 加 pocket detection noise：随机删/加边界 residues，口袋中心偏移，side chain 去除后重建。
4. 指标分为 ligand validity、affinity、pocket plausibility 和 sample efficiency。

投稿叙事：Apo2Mol 解决 ligand-free holo absence，但仍隐含 benchmark-only 信息；本工作把它推进到 real-world apo-only SBDD。

风险：如果 prior 已经不差，创新要转向“uncertainty-aware size ensemble + reranking”。

优先级：高。容易和 C1 合并成一个完整 benchmark/method paper。

### Idea C3: Dynamic Pocket Physical Plausibility Audit

一句话：给 generated holo-like pocket 建一套专门的 validity tests：side-chain rotamer、backbone geometry、steric clash、ligand strain、interaction fingerprint、pocket volume continuity、protein-ligand contact realism。

为什么值钱：Apo2Mol 评估里已有 pocket RMSD/TM-score 记录，但这些不能证明 pocket 可物理存在。PoseBusters 对 docking 的批评可以直接迁移到 dynamic pocket generation。

最小实现：

1. 整合 RDKit/PoseBusters-style checks for ligand。
2. 加 protein side-chain rotamer 和 clash checks，可用 Dunbrack/Rosetta/PyRosetta/OpenMM 简化版。
3. 用 ProLIF 统计 H-bond、hydrophobic、salt bridge、pi interactions。
4. 报告：高 Vina 但物理不合格的比例，作为现有方法崩溃点。

投稿叙事：dynamic SBDD 需要“pocket validity”，不是只需要 ligand validity。

风险：如果只是评估工具，贡献偏 benchmark；最好再加 guidance/reranking 改善结果。

优先级：中高。适合与 A3 或 C1 组合。

## 4. Best Three Paper Directions

### Paper 1: Apo2Mol-Hard + Realistic Apo Stress Test

最稳。以 C1+C2 为主，少量方法改进。标题可以是：

> Beyond Average Apo-to-Holo Generation: A Stress Benchmark for Dynamic Structure-Based Molecular Design

卖点：

- 本地数据天然存在 long-tail pocket dynamics。
- 指标从 average Vina 扩展到 bucket-wise dynamic robustness。
- 揭示 `ref atom count`、pocket noise、large RMSD shift 对模型的影响。
- 提供 robust baseline：RMSD-conditioned reweighting + learned atom-count prior + physical reranking。

适合投：NeurIPS/ICLR workshop、Bioinformatics、JCIM、Journal of Cheminformatics；如果 baseline 强，也可冲主会。

### Paper 2: RA-Apo2Mol Retrieval-Augmented Dynamic Pocket Diffusion

最贴合代码。以 A1 为主，C1 hard subset 作为主要验证场。

标题可以是：

> Retrieval-Augmented Apo-to-Holo Pocket Prompting for Dynamic 3D Molecular Generation

卖点：

- 不是泛泛 RAG，而是检索 apo-holo conformational deltas。
- Apo2Mol 代码已有 prompt path，补全自然。
- 长尾大构象变化上最容易显示收益。

适合投：ICLR/NeurIPS/AAAI/ICML workshop 到主会梯度；药化方向可投 JCIM/Bioinformatics。

### Paper 3: Adaptive Asynchronous Dynamic Co-Denoising

最工程但也最有系统价值。以 B1 为主。

标题可以是：

> Adaptive Asynchronous Co-Denoising for Efficient Flexible-Pocket Molecular Generation

卖点：

- 解决 1000-step 采样贵的问题。
- 从固定 protein update schedule 变成 event-driven residue-local updates。
- 对大规模 virtual screening 有直接价值。

适合投：机器学习会议 workshop/main 或系统型药物发现期刊。

## 5. One-Month MVP Plan

第一周：做 Apo2Mol-Hard 数据切分和 stress protocol。

- 从 `selected_index_apo_druglike.pkl` 生成 RMSD bucket split。
- 输出每桶数量、年份、protein family/ligand size 统计。
- 明确 ref atom count vs prior atom count 的评估设置。

第二周：跑原始 Apo2Mol checkpoint 小规模评估。

- 每个 bucket 选 20-50 个样本。
- 每个样本 5-20 个 ligands。
- 记录 time、validity、Vina、pocket RMSD、clash proxy。

第三周：实现一个小方法。

- 快线：learned atom-count predictor + physical reranker。
- 中线：RMSD-conditioned sampling/reweighting。
- 高线：top-k retrieval prompt 接入已有 prompt path。

第四周：写出 paper story。

- 图 1：Apo2Mol average benchmark 掩盖 long-tail。
- 图 2：stress protocol。
- 表 1：bucket-wise breakdown。
- 表 2：ref vs prior vs learned atom count。
- 表 3：baseline vs robust variant。
- case study：RMSD >= 3 A 的 induced-fit examples。

## 6. Source Links

- Apo2Mol arXiv: https://arxiv.org/abs/2511.14559
- Apo2Mol dataset card: https://huggingface.co/datasets/AIDD-LiLab/Apo2Mol_Dataset
- Apo2Mol GitHub: https://github.com/AIDD-LiLab/Apo2Mol
- DynamicBind, Nature Communications 2024: https://www.nature.com/articles/s41467-024-45461-2
- FlexSBDD, NeurIPS 2024: https://proceedings.neurips.cc/paper_files/paper/2024/hash/60fb8cf8000f0386063fb24ead366330-Abstract-Conference.html
- DynamicFlow, ICLR 2025: https://arxiv.org/abs/2503.03989
- FlowDock: https://arxiv.org/abs/2412.10966
- DecompDiff: https://arxiv.org/abs/2403.07902
- PoseBusters: https://arxiv.org/abs/2308.05777
- AlphaFold 3, Nature 2024: https://www.nature.com/articles/s41586-024-07487-w
