# multi-k 迭代组装：unitig 图跨轮验证（借鉴 metaMDBG）

> 2026-08-14。动机：anchr = Assembler of **N-free** CHRomosomes，目标是
> 拼出**无 gap、无嵌合的完整染色体**。用户裁定（2026-08-14）：metaMDBG 的
> **multi-k 迭代 + 跨接验证**是本项目核心方向，"只有这一条路能做到"。
> 机制细节见 `references/metaMDBG.md` §4.1.1（2026-08-14 重读补充）；
> 本文把它映射到 anchr 的碱基 k-mer 空间，形成可落地设计。
>
> **状态：v4 已实现（2026-08-14）**：`anchr asm multik`（`libs/asm/multik.rs`
> + `cmd/asm/multik.rs`），4 单测 + 6 集成测试，361 全量测试绿。
> * v1（跨接验证 + 嵌合清理 + 压实）：Lambda 20k k=21,51,81 → 64 contigs /
>   最长 7924 / N50 3790 / 总长 46670；
> * v2（+ 渐进丰度过滤 + recompact）：同数据 → **49 contigs / 最长 19207 /
>   N50 16676 / 总长 49579**（≈ 参考 48502，无 OLC 的 2.4× 冗余）；
>   OLC A 实验对照：52 / 19035 / 3409 / 116 kb（含冗余）。实测偏差见 §4.5。
> * v3（+ k 序列自适应）：`--kmer auto` 按读长 N50 生成 k 序列；**合成长读
>   端到端验证无 N 可达**（见 §4.7）。
> * v4（+ 主路径保护）：渐进过滤只删分支/孤立节点（直链永远保留），
>   recompact 后用 `compute_links` 重算边——**Lambda 20k → 38 contigs /
>   最长 46467 / N50 46467**（≈ 参考 48502 的 95.8%，reads 100% 覆盖、
>   零缺口；OLC 对照最长 19035）；无重复长读保持单条 98566（见 §4.8）。
>   **决定性测试**：`command_asm_multik_full_coverage_single_contig`——
>   完整环状覆盖长读 → 单条 = 基因组 100%（k-mer 多重集一致）。
>   核心算法单元测试：`bridge_kmer` 三个方向（正链/rc 目标/左端出）、
>   `progressive_filter` 主路径保护 + 孤立删除、`recompact_graph` 链合并
>   + 假边消失——374 全量测试绿。
> * v5（真实数据 + 工程修复，2026-08-14）：**G37 真实数据端到端验证**
>   （首个真实基因组）：覆盖参考 96.9% k-mer、总长 571,370（98.5% 参考）、
>   最长 91,246 bp 单条（旧路线 bcalm Sum 561.2K）。工程修复 4 处：
>   pass 0 启用 supermer+DFA（真实数据 0.66s→秒级）、环状图 head-walk
>   死循环、渐进过滤时机（每轮→最终一次）与 cutoff 上限（max→中位 25%）、
>   remove_unsupported 容错（允许 <2% 内部 k-mer 缺失）。见 §4.10。
> * v6（被删分支回灌，2026-08-14）：渐进过滤每轮删的低丰度分支**作为
>   下一轮 unitigs 继续参与**（megahit bubble 回灌 + metaMDBG unitig
>   反馈）——菌株/多态序列参与组装而非只输出。G37：misassemblies 保持 0、
>   最长 43,790 → 52,807（分支连接）、Genome fraction 95.99% → 95.33%
>   （权衡）。见 §4.11。
> * v7（megahit 清洗，2026-08-14）：tip_remover + weak_link_remover 移到
>   最终阶段（每轮做会拆主路径）。见 §4.12。
> * v8（严格链唯一性，2026-08-14，借鉴 SKESA）：merge/recompact 从"先到先得
>   占用"改为**严格两端唯一**（汇点不吞进链，含对称 link 去重修正）。G37：
>   misassemblies 保持 0、mismatches 27.7/100kbp 历史最佳、N50 24.4K（-8%
>   宁断勿嵌合）、Genome fraction 95.86%。见 §4.13。

## 1. 为什么现有路线不能保证无 N

现有 `asm olc`（`asm-olc.md`）是**多 k 并行 + 启发式拼接**：

```
k1/k2/k3 各自独立生成 unitigs → 合并为伪 reads → 精确 overlap
  → greedy best-edge layout（互惠 + top2 重复标记）→ consensus 缝合
```

这条路线有三处硬伤，正好对应无 N 的三个必要条件（连续 / 正确 / 完整）：

1. **选边是启发式，不是证据**：layout 选"最长 overlap + 互惠 best edge"。
   在重复区/菌株分歧（bubble）处，两条候选 overlap 长度可能相近，启发式
   要么选错（嵌合）、要么因 top2 重复标记断开（gap）——**没有 reads 计数
   证据区分主路径**。
2. **k 之间无反馈**：各 k 独立出 unitig，合并时靠 overlap 找关系。低覆盖
   区只有小 k 能连通，但小 k 的特异性差；大 k 特异性好但连不通。并行合并
   把"哪个 k 的连接可信"交给 overlap 长度，而不是图结构本身。
3. **无嵌合清理**：拼错/上一轮低 k 产生的错误 unitig 会留在路径里，没有
   "更大 k 下内部 k-mer 不 solid → 删除"的收敛机制。

metaMDBG 的迭代路线正好逐一补上：**选边用当前更大 k 的 solid k-mer 计数
（≥2）验证，k 之间 unitig 图结构反馈，嵌合 unitig 在下一轮被剔除**——
每一步都指向"路径由 reads 证据唯一决定"，这是无 N 的核心。

## 2. 借鉴机制回顾（metaMDBG §4.1.1 摘要）

```
每轮 k+1（k′-min-mer 长度 +1）：
  count(k)        : reads + 上一轮 unitig 序列 → 当前 k 的 solid 计数（≥2 才进表）
  solveEdges      : 每条相邻 unitig 边构造跨接 k′-min-mer（doublet），
                    查当前计数表 → 存在保留（压实成 edge node）/ 不存在剪边
  removeUnsupportedUnitigs : unitig 内部 k′-min-mer 有缺失 → 整条删（嵌合清理）
  solveSmallUnitigs        : 长度恰为单个 k′-min-mer 的小 unitig，
                    前向/后向 triplet 验证，支持的前驱×后继两两建边后删节点
```

关键洞察：**bubble 的各候选分支跨接 k′-min-mer 互不相同，只有被 reads 以
丰度 ≥ 2 支撑的分支保留**——更大 k 的计数天然区分主路径与菌株分支，无需
启发式。

## 3. 碱基空间等价映射

| metaMDBG（minimizer 空间） | anchr（碱基空间） | 现状 |
|---|---|---|
| k′-min-mer（minimizer 序列） | 碱基 k-mer（FastK 字节键，`Kmer::MAX_K=128`） | ✅ pgr |
| 每轮 k′ +1（≈ +200 bp 跨度） | 用户给定递增 k 序列（如 21→51→81） | ✅ 参数化 |
| reads + unitig 序列计数 | 多文件计数：reads + 上一轮 unitig FASTA | ✅ `TadpoleTable::build_supermer` 多序列输入 |
| unitig 图（`unitigGraph_prev`） | `asm unitig` 输出 + `compute_links` 的 L: 边 | ✅ `libs/asm/assemble.rs` |
| solveEdges（doublet 查 solid） | 每条 L: 边构造跨接 k-mer（u 末尾 k−1 碱基 + v 延续碱基，见 §4.2），查当前 k 计数 ≥ 2 | ➕ 新逻辑 |
| removeUnsupportedUnitigs | unitig 内部所有当前 k 的 k-mer 必须 solid，缺失 → 删整条 | ➕ 新逻辑 |
| solveSmallUnitigs（triplet） | 短 unitig（len < k）：前向 triplet（pred 末尾 k−len + small）与后向 triplet（small + succ 开头 k−len）验证 | ➕ 新逻辑 |
| createDoubletNode（压实） | 沿保留边合并 unitig（overlap k−1 碱基）→ 长 unitig | ➕ 新逻辑 |
| 渐进丰度过滤（同轮） | unitig 按 `cov=` 渐进剔除（v1 已有布局前过滤思路） | 挂账，先不做 |

碱基空间比 minimizer 空间简单的地方：**不需要丰度传播**。metaMDBG 的
`getAbundance`（取相邻 2 个 prev k-min-mer 的最小丰度）是因为 k′ 每轮 +1
导致 k-min-mer 集变化；碱基空间每轮直接对当前 k 重计即可，计数表就是
当前 k 的 solid 表。

## 4. 算法设计：`anchr asm multik`

### 4.1 总流程

```
输入: reads（1+ 文件）, k 序列 [k1 < k2 < ... < km], min_count_seed, min_count_extend
输出: 长 unitigs（无 N 染色体候选），可选 GFA

pass 0（k = k1）:
    unitigs ← asm unitig(reads, k1)          # 现有语义，最大无分支路径
    links ← compute_links(unitigs, k1)       # 现有 L: 边（共享 (k1-1)-mer）
    G ← (unitigs, links)

for pass i = 1..=m-1（k = k_{i+1}）:
    # 1. 计数（reads + 上一轮 unitig 序列）
    C ← count(k, reads ∪ G.unitigs)          # solid: count ≥ min_count_extend(=2)

    # 2. 跨接验证：G 的每条边
    for edge (u → v) in G.links:
        jk ← u 末尾 (k−1) 碱基 + v 的延续碱基    # 覆盖连接点的当前 k 的 k-mer
        if C[jk] ≥ 2: 保留
        else:         断开（remove_successor 双向）

    # 3. 嵌合清理
    for u in G.unitigs:
        if 任一 u 内部 k-mer ∉ C: 删 u（含其边）   # 上一轮拼错的在更大 k 下现形

    # 4. 小 unitig 处理（len(u) < k）
    for u in G.unitigs where len(u) < k:
        supported_pred ← {p | triplet(p,u) ∈ C}    # p 末尾 (k−len) + u
        supported_succ ← {s | triplet(u,s) ∈ C}    # u + s 开头 (k−len)
        if supported_pred / supported_succ 非空:
            支持的前驱 × 后继两两建边
        删 u

    # 5. 压实：沿保留边合并（u + v[(k_i−1)..]，去掉共享重叠），重新编号
    G ← merge_chains(G)

输出 G.unitigs（长 unitigs，按长度降序，头带 len/cov）
```

### 4.2 关键语义与细节

**跨接 k-mer 的构造**：`compute_links` 的 L: 边语义 = u 末尾与 v 开头共享
`(k_i−1)` 个碱基（建图 k 的重叠）。拼接 u→v 时 v 的**延续碱基**是
`v[k_i−1]`（v 中第一个不参与共享的碱基）。跨接 k-mer =
`u 末尾 (k−1) 碱基 + v[k_i−1]`（长度恰为当前 k）——即覆盖连接点的当前 k
的窗口，等价于 metaMDBG `getDoublet2` 的
`pred 末尾 (k−1) 个 minimizer + succ 延续 minimizer`
（`CreateMdbg.hpp:3328`，注意 succ 开头 `(k_i−1)` 与 pred 末尾共享，延续
从 `v[k_i−1]` 起，**不是** `v[0]`）。该 k-mer 不完整出现在 u 或 v 的序列里
（跨接窗口含 u 特有前缀），因此 unitig 自身计数不自动支撑它——验证才有
区分力。

**跨接 k-mer 的方向**：`compute_links` 已给每条 L: 边标注
`from_rc`/`to_rc`（u 的哪端、v 的哪端）。跨接 k-mer 必须在 u 的出端
（`from_rc=false` → u 右端；`from_rc=true` → u 左端取 rc）与 v 的入端
构造，且 canonical 归一化后查表——与 `libs/olc/overlap.rs` 的 seed 方向
处理同构。

**合并（压实）**：沿保留边合并时 u + v[(k_i−1)..]（去掉共享重叠），
与 `asm unitig` 的 `L` 边 overlap `(k−1)M` 语义一致。

**为什么 unitig 序列参与计数**：两重作用——(a) 上一轮确认的 unitig 内部
k-mer 在下一轮保持 solid，不会被 `removeUnsupportedUnitigs` 误删（低丰度
宏基因组物种尤其需要）；(b) 跨接 k-mer 获得"unitig 自身延续"的证据，让
确认的连接单调增长，不反复切分。这与 metaMDBG 的 `parser3(unitig_data.txt,
IndexKminmerFunctor(_extractingContigs=true))` 同义。

**单调性 vs 重新 unitig 化**：迭代**不重新跑 unitig 压缩**，只做"验证边 +
删坏点 + 合并"，保证已确认的连接不被更大 k 的新图结构破坏；新 k 的图结构
只用于验证连接与剔除嵌合。这是与"每轮直接重跑 `asm unitig`"的本质区别
（后者会把已合并的 unitig 重新切碎）。

**确定性**：计数表（排序键）、L: 边（排序去重）、unitig 集合（长度排序）
全部确定性；验证/合并按 (u,v) 键序处理，输出与并行线程数无关。

**复杂度**：每轮一次计数（supermer，与 `asm unitig` 同量级）+ 一次
O(unitig 总长 × log|C|) 的验证扫描。m 轮约 m × 单轮成本；k 序列默认
3 轮（21/51/81）→ 成本可接受。

### 4.3 与现有命令的关系

* `asm multik` 是**新的迭代驱动器**，内部复用 `assemble_unitigs_core`（pass 0）
  与 `compute_links`；不经过 overlap/layout/cns。
* 输出长 unitigs 后，**可选**接现有 `asm ovlp → layout → cns`（或直接输出）。
  v1 以迭代输出本身为结果（无 N 的验证对象），OLC 收尾留真实数据评估后定。
* `asm olc` 保持现状（并行多 k），`asm multik` 是串行迭代路线——两者并存，
  文档注明迭代路线是长期方向。

### 4.4 参数

```
anchr asm multik <infiles>... -o unitigs.fa \
    --kmer 21,51,81        # 递增 k 序列（逗号分隔，默认 21,51,81）
    --min-count-seed 3     # pass 0 的 solid 阈值（透传 asm unitig）
    --min-count-extend 2   # 跨接验证/内部验证的 solid 阈值（默认 2）
    --parallel <int|auto>  # 计数并行（supermer/streamed）
    --gfa                  # 输出 unitig 图（S/L 行，可选）
    --keep-dir DIR         # 每轮中间产物（unitigs/links/计数表统计）
```

### 4.5 v1 实现偏差（2026-08-14 实测）

* **小 unitig 验证从简**：设计里的小 unitig triplet 在碱基空间退化为与
  doublet 相同的"覆盖连接点的当前 k 窗口"（metaMDBG minimizer 空间
  doublet==triplet）；实现统一用 `bridge_kmer`（实际序列匹配方向，不解释
  `to_rc` 符号）。**当 unitig 短于 `k_prev` 时无法构造窗口，其边一律断开**
  （保守，短 unitig 保留为独立输出）——合成重复序列（unitig 21 bp）在
  21→41 大步长下因此不合并；随机基因组（unitig 长）正常合并。
* **大步长限制**：碱基空间 unitig 最短 = `k_prev`（建图 k 的单 k-mer），
  metaMDBG minimizer 空间最短 = `k−1` minimizer。故大步长（如 21→81）时
  短 unitig 无法参与验证。缓解：k 步长收敛（21/51/81 比 21/81 好）、短
  unitig 自然隔离。v2 可考虑"链太短则跳过验证"或小步长自动插值。
* **`--gfa` / `--keep-dir` 未实现**（v1 只输出 FASTA）；`--parallel` 接受
  但计数走 supermer 内存路径（unitig 序列无质量，统一无质量门控）。
* **合并只沿"两端唯一"的边**：bubble 两侧不合并（保守），输出保留分支；
  metaMDBG 的渐进丰度过滤（同轮删低丰度分支）未并入——bubble 内低覆盖
  分支在跨接验证下因 51-mer 不 solid 自然断开（Lambda 实测有效），但
  高覆盖菌株差异（两条都 solid）需 v2 丰度过滤。
* **空/过短输入**：输出空文件（exit 0），与 `asm unitig` 一致。

### 4.6 v2：渐进丰度过滤 + recompact（2026-08-14）

**背景**：v1 只处理"低覆盖分支"（51-mer 不 solid 自然断开）；**高覆盖菌株
差异**（两条分支都 solid）v1 全部保留，bubble 处主路径无法合并。Lambda
实测 v1 最长 7924 bp 卡住。

**实现**（metaMDBG `removeAbundanceNoQueue` + `recompact` 映射）：

1. **渐进过滤**（`progressive_filter`）：cutoff `t` 从 1.1 起步、每轮
   `t += min(t*0.1, 10)`，删 `coverage < t` 的 unitig（进低丰度输出列表，
   不丢失）；直到 `t >= max_abundance` 或无节点可删。每轮迭代后执行。
2. **recompact**（`recompact_graph`）：删除后立即把唯一链合并成更长
   unitig（覆盖度取平均）——**关键**：合并后主路径 unitig 继承侧翼的高
   丰度，不会被后续更高 cutoff 误删（metaMDBG 同轮 `recompact` 语义）。
3. **方向匹配统一用建图 k0**（不是当前验证 k）：unitig 边共享
   `(k0−1)-mer`，`bridge_kmer`/`recompact_graph`/`merge_chains` 全部以 k0
   为共享段长度——v1 曾误用当前 k 导致第二/三轮验证全失败（真实数据
   触发，已修）。
4. **被删的 unitig 保留为独立输出**：对应 metaMDBG 的 cutoff 快照（低丰度
   物种/菌株不丢失，作为独立 contig 输出）。

**实测**：合成"高覆盖菌株差异 bubble"（30× 主路径 + 5× 分支，共享侧翼）：
v1 输出 4 条碎片（主路径 500 bp 无法合并）；v2 输出 **500 bp 主路径** +
240 bp 菌株分支 ✓。Lambda 20k：最长 19207 / N50 16676（OLC 对照 3409）。

**阈值语义**：v2 沿用 metaMDBG 的"渐进到 max_abundance"（主路径优先），
低丰度内容靠独立输出保留；`0.25×` 类保护阈值未引入（metaMDBG superbubble
的 `currentCutoff/0.25` 保护留给 v3，若宏基因组低丰度物种连通性被破坏再
加）。

### 4.7 v3：k 序列自适应 + 长读端到端验证（2026-08-14）

**k 序列自适应**（`auto_ks`，`--kmer auto` 默认）：按读长 N50 生成——
`k_max = min(0.8×N50, 128)`，起点 `clamp(N50/10, 21, 31)`，步长
`clamp(N50/100, 20, 30)`。短读（108 bp）→ 21/41/61/81；长读（≥10 kb）→
31/61/91/121（对应 metaMDBG `computeLastK` 的"最后 k-min-mer 跨度 ≈ 2×N50"）。
Lambda 短读 auto 与手动 21/51/81 结果相同（49 contigs / N50 16676）。

**合成长读端到端验证（无 N 可达性的核心证据）**：

| 场景 | 输入 | 输出 | reads 回贴 |
|---|---|---|---|
| 无错长读 | 100 kb 随机基因组，15 kb reads ×30× | **单条 98,499 bp**（98.5%） | 98499/98499 覆盖、零缺口 |
| 真实 HiFi 错误率 | 同上 + 0.1% 替换错误 + 覆盖两端 | **单条 98,566 bp**（98.6%），N50 98566 | —（0.1% 错误下完美匹配 map 不适用） |
| **完整覆盖（环状）** | 100 kb 基因组 + 0.1% 错误 + 30× 随机 + 50 条跨起点/从起点 reads | **单条 100,000 bp = 基因组 100%**，N50 100000 | k-mer 覆盖 100.000%（99970/99970），**多重集与基因组完全一致（环状旋转等价）** |
| 重复区 + 0.1% 错误 | 104 kb 基因组 + 2 kb 精确重复 ×2 | 9 contigs / 最长 69702 / 总长 99510 | **k-mer 覆盖 97.4%**（缺失仅 reads 边缘；两个重复区 980/980 全覆盖） |
| 重复区 + 0.5% 错误 | 同上（过噪） | 28 contigs / 最长 30368 | 0.5% 下跨接 k-mer 45% 带错，重复区断开（真实 HiFi ≈0.1% 无此问题） |

结论：**multi-k 迭代路线在长读数据下能拼出单条近全长（≈98.5%）零缺口
contig**，**完整覆盖（含环状跨起点 reads）下拼出 100,000 bp = 基因组 100%
（k-mer 多重集一致 = 环状旋转等价）**——"无 N 染色体"可达性得到决定性
验证（这是用户"只有这一条路能做到"判断的直接证据）。残余差距仅为数据
覆盖边缘（reads 未覆盖的两端）；**重复区在真实错误率（0.1%）下被完整
覆盖**（k-mer 层面 980/980）。

### 4.8 v4：主路径保护（2026-08-14）

**Bug 实锤**：渐进过滤"删 `cov < max_abundance`"会误删主路径——同一基因组
的覆盖有波动（重复区 76.5× vs 主路径 47.3× vs 其他 23.8×），`max_abundance`
常被重复区抬高，主路径低覆盖段被当成"低丰度分支"删掉（Lambda 实测最长卡在
19 kb；合成长读 round k=91 曾把全部 unitigs 删光）。

**修复**（两处）：
1. **直链保护**（`progressive_filter`）：先统计每个 unitig 的端级出入边数；
   **两端唯一（直链）的 unitig 永远保留**（主路径），只有分支点（某端
   >1 边）和孤立节点按丰度删。这同时避免宏基因组里低丰度物种的直链被误杀。
2. **recompact 后重算边**：`recompact_graph` 合并链后不再手工 relink（链间
   边方向易错、端点被合并吸收后产生假边），改用 `compute_links` 按新端点
   的 `(k_build−1)-mer` 重算——假边消失，真实边由下一轮跨接验证重新把关。

**实测**：
* Lambda 20k：**38 contigs / 最长 46467 / N50 46467 / 总长 49359**——最长
  ≈ 参考 48502 的 95.8%，reads 回贴 100% 覆盖、零缺口（OLC A 实验最长
  19035）；这是**短读数据下第一条近整条染色体级无 N contig**；
* 无重复合成长读：单条 98566 保持（修复无回归）；
* 0.1% 错误 + 重复区合成长读：round k=91 保留主路径（修复前删光），最终
  9 contigs / 最长 69702 / k-mer 覆盖 97.4%——69702 与 28498 之间无桥接
  reads（真 gap，multik 正确断开）。

**语义修正**：渐进过滤从"纯丰度删"变为"丰度删 + 主路径保护"——更贴近
metaMDBG 的图结构语义（它的主路径是单 unitig 节点，天然不会被拆），同时
避免了碱基空间"同一基因组覆盖波动"导致的误删。

### 4.9 PE mate 桥接分析（2026-08-14，结论：不做 scaffold）

Lambda 20k 是 **PE reads**（R1/R2，insert 中位数 405.5 bp）。`asm multik`
输出 unitig_1（46467，≈ 参考 48502 的 95.8%，零缺口）+ unitig_2（1930）+
unitig_3（108）——三者 k-mer 零重叠、在参考上覆盖互补区域（主体 + 开头），
**正好可拼出 ≈ 完整 Lambda 参考**。但：

* **mate 桥接证据充足**：106 对 R1/R2 跨 contig（unitig_1↔unitig_2 79 对、
  unitig_2↔unitig_3 27 对），方向一致（RF 构型）、位置明确（unitig_1 尾
  ↔ unitig_2 头）；
* **无序列重叠**：unitig_1 尾与 unitig_2 全序列（正反链）零共享 k-mer——
  断点处是真实序列缺口（估计 ≈277 bp，insert 405 − 端距）；
* **结论**：PE mate 只能 **scaffold（gap 填 N）**，违背"无 N"目标；gap 恢复
  需要局部组装（gap filling，108 bp reads 无法单条覆盖 277 bp gap）——
  **不做 scaffold**。短读下 95.8% + 零缺口即 multik 的数据极限；完整无 N
  染色体依赖长读（合成长读已验证单条 98.5% 零缺口）或 gap filling（v5
  候选，需真实数据评估）。

**gap filling 实验（2026-08-14，结论：数据极限）**：收集锚定到断点两侧
的 215 对 reads（unitig_1 尾 `pos>46200` / unitig_2 头 `pos<300`），用
`asm multik` 局部组装 → 最长 640 bp contig，但**不含 unitig_1 尾或
unitig_2 头序列**（断点未桥接）。原因：断点处 reads 覆盖剖面——unitig_1
尾 46400-46467 覆盖 29→**4×**（急剧下降），unitig_2 头 0-10 bp **0×**、
10-30 bp 13×——**gap 中间 reads 覆盖不足，局部组装无原料**。结论：
Lambda 短读的 95.8% 是数据覆盖极限（非算法缺陷）；multik 正确输出零缺口
主 contig，gap 处保持断开（无 N 优先于完整）。

### 4.10 v5：真实数据验证 + 工程修复（2026-08-14）

**G37（Mycoplasma genitalium）真实数据端到端**（`results/model.md` 的
手动模拟例子，纠错后 reads `Q25L60X40P000/pe.cor.fa`，150 bp × 155k 条，
40×；参考 580,076 bp）：

* multik auto：395 contigs / 最长 **91,246 bp** / N50 24,527 / 总长
  **571,370**（参考 **98.5%**）；k-mer 覆盖参考 **96.9%**；
* 对比旧路线（`results/model.md` 手动流程）：bcalm 6 k + contained →
  Sum 561.2K / N50 ~15-25K；merge anchors → N50 55K / 17 条。multik
  总长更接近参考（571K vs 561K），有 91 kb 单条长 contig；碎片（395 条）
  多于旧路线（17-50 条）——旧路线用 `anchr contained` 激进去冗余，multik
  保留全部内容（含低丰度 dropped）。

**真实数据暴露的工程问题（全部修复）**：

1. **pass 0 性能**：multik 未启用 `use_supermer`/`use_dfa`（asm unitig
   默认有），真实数据（155k reads）单轮卡 30s+ → 启用后 0.66s；
2. **环状图死循环**：`merge_chains`/`recompact_graph` 的 head-walk 无环
   检测（细菌染色体首尾相接 → left_of 环 → 无限循环）→ 加 `seen` 检测；
3. **渐进过滤时机与 cutoff**：
   * 每轮迭代都跑渐进过滤（累积删除、单菌株覆盖波动误删）→ **只在最终
     merge_chains 前做一次**；迭代轮只做跨接验证 + remove_unsupported +
     recompact（主路径每轮增长）；
   * cutoff 上限从 `max_abundance`（重复区可到 600×）改为 **cov 中位数的
     25%**——单菌株正常覆盖（40×）不被删，只删显著低丰度；
4. **remove_unsupported 容错**：单菌株覆盖波动（个别内部 k-mer <2×）会
   误删整条主 unitig（Lambda 46,467 曾因此消失）→ **允许 <2% 内部 k-mer
   缺失**（`max_missing = max(1, n_kmers/50)`），只删真正嵌合的 unitig。

**短 unitig 跨接验证跳过**（§4.5 大步长限制的落地）：u/v 短于当前 `k−1`
窗口时跳过跨接验证（保留边，merge_chains 用实际端点匹配决定）——避免
大步长（21→61）下短 unitig 的边被误删（G37 曾因此把主路径边全删）。

**性能**（G37 40× / 155k reads，release）：auto k（6 轮）~3 s / ~10 s CPU。

### 4.11 v6：被删分支回灌（2026-08-14，借鉴 megahit bubble 回灌）

**动机**（用户指出 megahit 宏基因组评价好，`references/megahit.md` §8.6）：
megahit 把被合并的气泡序列（`bubble_seq.fa`）喂回下一轮 `seq2sdbg` 建图，
metaMDBG 把上一轮 unitigs 反馈计数——**菌株/多态序列参与组装而非丢弃**。
multik 原先把渐进过滤删的分支只输出（dropped），不参与后续轮。

**实现**（`multik.rs`）：`progressive_filter` 移回迭代中（每轮 recompact 后
调用，cutoff_cap 中位 25% + 直链保护不变），删的分支返回 `Vec<Unitig>`，
下一轮循环开头 `unitigs.append(&mut carried)` 加回——分支参与后续轮的
跨接验证/压实，可能被连接（回灌成功）或再删（最终 carried 输出）。

**G37 Quast**（对照参考 580,076）：

| 指标 | v5（无回灌） | v6（回灌） |
|---|---:|---:|
| # misassemblies | 0 | **0** |
| Largest contig | 43,790 | **52,807**（+20%，分支连接） |
| N50 | 26,562 | 23,585 |
| Genome fraction | 95.99% | 95.33% |
| # contigs | 487 | 295 |

**不回归**：Lambda 46,457 / N50 46,457；20k 环状单条 100%；374 测试全绿。

**权衡**：回灌让分支参与组装（最长 +20%、contigs 487→295），但 Genome
fraction 略降（-0.66%，部分回灌分支在更大 k 下被重新评估/切分）。对宏基因
组（菌株序列重要）方向正确；阈值/回灌范围可后续调（如只回灌相似度高的
分支，或调整 cutoff）。

### 4.12 v7：megahit 清洗（tip + weak link，2026-08-14）

借鉴 megahit 清洗算法族（`references/megahit.md` §5）到 multik 最终阶段：

* **`tip_remover`**（megahit `tip_remover.cpp`）：短 unitig（≤ 2×k0）+ 是
  tip（一端无连接）+ 深度 `> 20×` 低于邻居 → 删（错误尖端）；
* **`weak_link_remover`**（megahit `weak_link_remover.cpp`）：分支点（出度
  ≥2），邻居深度 ≤ 0.05× 邻居总深度 → 断开（删边保留节点，菌株共享区
  不错连）。

**教训**：**每轮迭代做清洗会拆主路径**（G37 最长 52.8k → 32.6k——k0=21
的碎片大量是"短 tip"，8×/0.1 阈值误删/误断）；**移到最终（unitigs 压实后）
无副作用**（真实尖端少、弱连接明确）。

**G37 实测**：最终清洗后 misassemblies 0、最长 52,807（保持 v6）、Lambda
46,457、20k 100%、374 测试全绿。单菌株 G37 无弱连接/尖端问题（渐进过滤已
覆盖），tip/weak_link 的**价值在宏基因组**（多菌株弱连接、低覆盖菌株共享
区）——待真实宏基因组数据验证（todo §4）。

### 4.13 v8：严格链唯一性（2026-08-14，借鉴 SKESA "predecessor == 1"）

SKESA 扩展不变量（`references/skesa.md` §7.2）的 unitig 图对应：
`merge_chains` / `recompact_graph` 原来用**先到先得**的端点占用检查——
汇点（两个前驱）会被吞进先遍历到的链，谁被吞取决于 link 遍历顺序而非
证据。改为**严格两端唯一**：

* 两遍扫描：先按链段统计每个定向端点的出/入度，再合并时要求链段两端度数
  都恰为 1（`out_deg[ln]==1 && in_deg[rn]==1`）；占位检查保留作防御。
* **易错点（已修）**：`compute_links` 对同一 junction 从两端各发一条 link
  （u 的 out-link + v 的 in-link），方向解析后是同一条链段；逐 link 计度数
  会把对称 link 翻倍，使线性链误判为分叉。必须先用 `HashSet` 按
  `(ln, rn)` 去重再计度数。
* 语义效果：汇点/分叉处**宁断不吞**，谁当前驱交给渐进过滤（低丰度分支被
  剪后下轮自然唯一、链恢复合并）——与 SKESA "过滤后前驱恰好 1 才扩展"
  同族；`bridge_kmer ≥2` 只验"有支撑"，严格唯一性补上"唯一 + 方向一致"。

**G37 实测**：misassemblies 保持 0、mismatches 27.7/100kbp（历史最佳）、
N50 24,445（v7 26,562，-8%：宁断勿嵌合的正确性代价）、Genome fraction
95.86%（-0.13pp）、0 N；图结构递减放缓（unitigs 1345→807、edges 2316→
1864——分支不再被吞进链，故中间图更大）。Lambda/20k 合成长读与全部 374
测试不回归。

## 5. 数据结构与复用

| 组件 | 复用 | 用途 |
|---|---|---|
| pass 0 | `libs/asm/assemble.rs::assemble_unitigs_core` + `compute_links` | 初始 unitig 图 |
| 计数 | `TadpoleTable::build_supermer` / `build_streamed` | 当前 k solid 表 |
| 查表 | `TadpoleTable::get_count`（排序键二分） | 跨接/内部/triplet 验证 |
| k-mer 编解码 | `pgr::libs::kmer::key::Kmer`（FastK 字节键） | 跨接 k-mer 构造 |
| 方向 | `compute_links` 的 `Link{to, from_rc, to_rc}` | 边方向语义 |
| 合并 | `libs/olc/consensus.rs` 的精确缝合思路（或直接序列拼接） | 压实 |
| CLI | `cmd/args.rs` 标准参数 | 与 `asm unitig`/`olc` 一致 |

新增逻辑放 `libs/asm/multik.rs`（anchr 业务库），`cmd/asm/multik.rs` 薄壳；
与 `asm-olc.md` 的分层原则一致（算法在 libs，CLI 在 cmd）。

## 6. 验证计划

### 单元测试（libs/asm/multik.rs）

* 跨接验证：构造两 unitig 共享 (k1−1)-mer 端点，reads 含/不含跨接 k-mer →
  边保留/断开；
* 嵌合清理：unitig 内部混入非 reads 子串 → 整条删除；
* 小 unitig：len < k 的 unitig 经 triplet 验证后正确连入/删除；
* 压实：链式边合并后序列正确（u + v[(k−1)..]）、无重复碱基；
* 确定性：相同输入两次运行逐字节一致；
* 单调性：已确认连接在下一轮不因更大 k 被切碎。

### 集成测试（tests/cli_asm_multik.rs）

* 合成基因组（含模拟菌株 bubble：两条平行路径，一条高覆盖一条低覆盖）：
  `asm multik` 迭代后长 unitig 沿高覆盖主路径连通、低覆盖分支被剪；
* 与 `asm olc` 对照：Lambda 真实 reads 上 multik 输出 contig 数 / N50 /
  最长 contig，reads 回贴验证（`asm map`）无 gap、无嵌合；
* 零 panic：畸形输入（空文件、单 read、全 N）友好报错。

### 验收门

`cargo fmt` + `cargo clippy --all-targets -- -D warnings` + `cargo test`
全绿；命令注册与 docs 由 `cli_consistency.rs` 约束。

## 7. 风险与待决

* **k 序列选择**：v3 已实现 `--kmer auto`（读长 N50 驱动）；短读/长读
  默认序列见 §4.7。边界情况（极短 reads < 100 bp、读长混杂）待真实数据
  调参。
* **跨接验证的阈值**：`count ≥ 2` 是 metaMDBG 语义；宏基因组低丰度物种
  可能整条 unitig 的 k-mer 都只有 1–2×——"unitig 参与计数"缓解此问题，
  但阈值是否需随 k 调整（如大 k 放宽）待真实数据定。
* **渐进丰度过滤未并入**：同轮内的低丰度分支剔除（§6.3）v1 不做，先验证
  跨接验证单独的效果；若真实数据暴露菌株分支残留，再加。
* **与 OLC 的最终关系**：迭代输出是否直接是最终 contigs，还是需 OLC 收尾，
  等真实数据端到端结果定（无 N 判据：覆盖完整 + 无 gap + 无嵌合）。
* **重复区处理**：跨接验证只能区分"有/无 reads 证据"，纯重复区（两个拷贝
  序列相同）两条边都有证据——需 reads 桥接（`RepeatRemover` 语义，v1 待办）
  或覆盖度证据收尾，不指望跨接验证单独解决。

## 8. 相关文档

* `references/metaMDBG.md` §4.1.1（跨接验证机制，2026-08-14 重读补充）
* `design/asm-olc.md`（现有并行多 k OLC；v1 素材与对照）
* `design/asm-unitig.md` §8（`asm unitig` 语义与 L: 边）
* `pgr: libs/kmer/supermer.rs`（FASTA 默认计数路径）

<!-- 以下内容并入自 `asm-multik-misassembly.md`（2026-08-15 文档合并） -->

## 9. 防 misassembly 方案（bridge_filter / split_by_bridge）

> 2026-08-14。背景：G37 Quast 质检（`benchmarks/multik-g37-quast.md`）
> 显示 multik 有 8 个 misassembly（全 relocation），unitig_1 把参考相隔
> ~48 万 bp 的两段（4286-60431 与 540652-571614）连在一起。参考这两段
> **无共享 k-mer**（非简单重复区）。老流程 anchors 的目的就是防
> misassembly，但其实现（bbwrap→basecov→spanr→hnsm→contained→orient→
> merge）极冗余。本文规划精简方案：**只取 anchors 的核心思想（reads 回贴
> + 覆盖度区间过滤），复用 multik/pgr 已有机制**。

### 1. anchors 防 misassembly 的核心机制（读 `templates/anchors.tera.sh` 确认）

1. reads 映射回 unitigs（bbwrap perfectmode）→ 每碱基覆盖度；
2. 覆盖度区间 `[lower, upper]`：
   `lower = max(mincov, (median − mscale×MAD)/lscale)`、
   `upper = (median + mscale×MAD)×uscale`；
3. **只保留覆盖度在区间内的区域**（低覆盖 = 错误区、高覆盖 = 重复区，
   都不可靠 → 排除）→ "properly covered regions" 才是 anchors。

**防 misassembly 的原理**：错误连接通常发生在覆盖度异常的区域（重复区
reads 过多、错误区 reads 过少）；把这些区域排除（宁可断开、不错连），
misassembly 自然减少。

### 2. multik 的现状与差距

* multik 已有：unitig `cov=`（组成 k-mer 平均覆盖度）、渐进丰度过滤
  （只删低覆盖分支）、跨接验证（k-mer 计数 ≥2 选边）；
* **缺**：`upper`（高覆盖/重复区排除）——渐进丰度过滤只处理低覆盖，
  没有"高覆盖=重复区→断开"的机制。这是 8 个 relocation 的根源（错连
  发生在覆盖异常区）。

### 3. 方案（三个层次，从简到繁）

#### 层次 1：覆盖度一致性选边（迭代内，零额外工具）

* 位置：跨接验证后、压实/合并前；
* 机制：multik 的 unitigs 已有 `cov=`；对每条候选边 (u→v)，
  **若 v.cov 显著高于 u.cov（或高于全局 upper）→ 断开该边**（v 可能是
  重复区，连接不可靠）；
* 阈值：`upper = (median + mscale×MAD)×uscale`（anchors 同款，unitig 级
  cov 分布）；
* 优点：零额外工具、零成本（cov 已有）；
* 局限：unitig 级平均覆盖掩盖局部异常（重复区只是 unitig 的一部分时
  平均 cov 可能正常）——G37 的 misassembly contigs cov 38-109 都 < upper
  （196），**层次 1 可能不够**。

#### 层次 2：reads 回贴 + 局部覆盖剖面（复用 pgr map 工具链）

* 位置：multik 输出后（或迭代末），一个独立过滤步骤；
* 机制：**复用 anchr 已有工具**（不用 bbwrap/spanr/hnsm）：
  1. `pgr asm map`（perfect mode，已实现）把 reads 映射回 unitigs；
  2. `pgr sam to-rg` + `pgr rg coverage`（todo 已列）→ 每碱基覆盖度；
  3. 覆盖度区间 `[lower, upper]`（anchors 同款参数 mscale=3/lscale=3/
     uscale=2/mincov=5）；
  4. 覆盖度在区间外的区域 → **断开 unitig**（按 rg 边界切分，不输出
     跨异常区的连接）；
* 优点：精确复刻 anchors 核心（每碱基剖面），但只有 3 个已有命令，
  无 spanr/hnsm/contained/orient/merge 冗余链；
* 局限：perfect mode 对含错误 reads 可能漏贴（低覆盖误判）——可用
  k-mer 命中容错（multik 的 bridge_kmer 思路）替代；
* 这是**推荐路径**（层次 1 不够时的自然升级）。

#### 层次 3：reads 桥接选边（metaMDBG RepeatRemover 语义，最精确）

* 对断开点/分支点，用 reads 同时覆盖两侧的证据（`nbBridgingReads != 0`）
  选择正确连接（metaMDBG §4.2.1）；
* 之前已分析（`asm-multik.md` §4.9 附近）：纯重复区（两拷贝都有证据）
  需要 reads 桥接或覆盖度证据收尾；
* 复杂度最高，作为层次 2 不足时的最终手段。

### 4. 验证计划

1. **层次 1**：G37 上实现 cov 一致性断边，Quast 复测 misassemblies
   （预期部分减少，但 unitig 级平均覆盖可能漏）；
2. **层次 2**：复用 map/to-rg/coverage 做局部剖面，Quast 复测
   （预期 8 个 relocation → 显著减少，代价是 Genome fraction 略降——
   断开的异常区不输出）；
3. 同时验证：Lambda（46,467 零缺口）和合成长读（100% 单条）不回归
   （过滤只在覆盖异常区生效，正常区域不切）。

### 5. 与老流程的对比（为什么精简）

| 老流程 anchors | 本方案层次 2 |
|---|---|
| bbwrap（外部工具） | `pgr asm map`（已有） |
| basecov 解析 + perl 过滤 | `pgr sam to-rg` + `rg coverage`（已有） |
| spanr cover/stat/some/span/compare（外部） | 无（直接按 rg 覆盖阈值切分） |
| hnsm range 提取 | 无（在 unitig 上直接切分） |
| contained→orient→merge（4 步） | 无（multik 已有压实） |

本方案只保留 anchors 的核心判据（覆盖度区间），把"切分/去冗余/合并"
交给 multik 已有机制，去掉全部外部工具链与后处理冗余。

### 6. 待决

* 层次 1 vs 层次 2 先做哪个（若 unitig 级 cov 已能抓住错连，层次 1 足够；
  证据表明平均 cov 掩盖局部异常，倾向直接做层次 2）；
* 断开策略：异常区切分后是丢弃（像 anchors 不选）还是保留为独立 unitig
  （像 multik 的 dropped）；
* 阈值参数是否沿用 anchors 默认（mscale=3/lscale=3/uscale=2/mincov=5）。

### 7. 实现结果（2026-08-14，层次 3 落地）

**实现**（`multik.rs`）：
* `bridge_filter`：对每条 unitig 间边构造连接探针（u 尾 30bp + v 延续 30bp
  = 60-mer），查 reads 表（`TadpoleTable::build_supermer`），count ≥ 2
  保留、否则断开——metaMDBG `computeBridgingReads` 语义；每轮 recompact
  前 + split 后各跑一次；
* `split_by_bridge`：对每个 unitig 内部滑动 60-mer 窗口，无任何 reads
  支撑（count 0）的窗口为嵌合连接点 → 切分 unitig（丢弃 < k0 碎片）；
  **这是关键**——渐进丰度过滤的 recompact 会把错连边固化成单个 unitig
  （G37 的 89,411 单 unitig 含 3 个 relocation），unitig 间桥接管不到，
  必须切内部；
* **探针长度 60-mer（probe_half=30）**：100-mer 对含错误 reads（Lambda
  原始 FASTQ 0.1% 错误）误切 4.5% 窗口；60-mer 容错（0.999^60 ≈ 94% 全
  匹配），Lambda 46,457 全窗口有支撑。

**G37 Quast 复测**（对照参考 580,076）：

| 指标 | 修复前 | 层次 3 后 |
|---|---:|---:|
| # misassemblies | **8** | **0** |
| N50 | 24,527 | **26,562** |
| Largest contig | 91,246 | 43,790 |
| Genome fraction | 95.58% | 95.99% |
| mismatches / 100 kbp | 31.4 | 33.7 |
| indels / 100 kbp | 2.3 | 2.5 |
| # N's / 100 kbp | 0 | 0 |

**其他不回归**：Lambda 46,457 / N50 46,457（修复前 46,467）；20k 环状
单条 100%。374 测试全绿、fmt/clippy 干净。

**代价**：G37 最长 contig 91,246 → 43,790（错连 unitig 被切）；Genome
fraction 基本持平（95.99%）。**无错连优先于最长 contig**（符合"无 N"目标：
正确性 > 完整性）。真实数据上错连与最长 contig 的权衡可后续调
（探针长度/阈值/是否切 vs 保留）。

### 8. 2026-08-16 补充：MG1655 重复介导嵌合（mis 4→0）

`bridge_filter`/`split_by_bridge`（层次 3）把 G37 的 relocation 压到 0，
但 MG1655 5 组合并后仍有 4 个真嵌合（contig_2/24/37 + contig_7）。逐链
复核确认嵌合在 **multik 迭代轮内**形成（非 merge 阶段），两类源头：

1. **嵌合 merged reads 桥接**：`fq merge`（bbmerge 语义）在 IS 倒转重复
   （`TTGGTTTGGGAGAA` 14 bp TIR 基序）处把两条 reads 错接成 300-400 bp
   嵌合 reads；其 k-mer 在 pass 0 形成 84-122 bp 桥接 unitig（如
   `unitig_679` = [侧翼 rc][IS 起始]），通过 reads 表使 60-mer 探针 count
   ≥2，桥接链接通过验证，recompact 折返链把侧翼与 IS 连成
   `rc(flank)+IS`。
2. **重复核心片段**：多拷贝重复（如 ref 925/1097/2835 kb 三拷贝）的核心
   在 pass 0 连接 4 个侧翼 unitig（2 in + 2 out，如 `unitig_445` 122 bp）。
   严格端唯一性在 4 链接时本应断链，但 `bridge_filter` 修剪掉一个链接后
   其度降至唯一，recompact 折返链跨重复连接两拷贝 → 171 kb 缺失式嵌合
   （contig_7）。

**修复**（`multik.rs`）：

* 链连接最短 unitig 长度 `max(2×(k−1), 90)`（`oriented_segment`）：排除
  嵌合 reads 桥接片段——其 begin/end (k−1)-mer 重叠，链接方向可折返；
  短片段保持独立输出（宁断勿错）。
* pass 0 快照 ≥4 个不同链接伙伴的 unitig 为分支节点（重复核心），其链接
  永不参与链压实；标志随 `retain_graph`/`recompact_graph`/`split_by_bridge`/
  `progressive_filter` 的 unitig 重索引传播。气泡（≤3 伙伴，菌株分歧）由
  丰度过滤解析，不受影响。

**验证**：

| 数据集 | 修复前 | 修复后 |
| :--- | :--- | :--- |
| MG1655 5 组 multik51 全链 | 4 mis，N50 65.8K，GF 97.36% | **0 mis**，N50 60.3K，GF 97.22% |
| G37 MRX40P000 6 主 K 链 | 0 mis，N50 55.4K，GF 97.05% | 0 mis，N50 37.6K，GF 96.97% |

G37 的重复/分支节点链原本正确，现在被保守断开，N50 代价 ~32% 是
**宁断勿错**的正确性取舍（GF 基本持平）；MG1655 的代价仅 ~8%。397 测试
全绿、fmt/clippy 干净。待真实宏基因组/长读数据到位后评估该保守策略是否
过度（可在 `min_chain_len`/伙伴数阈值上调参）。

**unitig/bcalm 链 mis 归属修正**：旧报告（2026-08-15 生成）里 unitig/bcalm
链各 1 mis 不是它们本身的问题，而是 **旧 OLC 生成方式的遗留**——unitig
链 mis（contig_1）由两个各自干净的 unitig（contig_24 → ref 2.76Mb 连续、
contig_107 → ref 280-290kb 连续）经 **80 bp exact overlap** 在 olc 合并时
错连（80 bp ≪ min-overlap 1000）；用当前 `asm olc --unitigs
--min-overlap 1000` 重跑 unitig/bcalm 链均为 0 mis（unitig：N50 78.8K/
112 contigs/GF 97.73%；bcalm：N50 73.7K/118/GF 97.76%）。multik 链的
4 mis 才是 multik 自身（unitig 内嵌合），已由本修复归零。

**unitig / bcalm / multik 三链端到端对比**（同输入、统一口径）：multik
N50 95.5K / 107 contigs / GF 97.61%，unitig 与 bcalm 均为 N50 67.4K /
121 / GF 97.54%，三链 0 mis。multik 全面占优（N50 +42%），unitig 与
bcalm 等价（自研可替代外部依赖）。详见
`benchmarks/mg1655-unitig-bcalm-multik.md`。

<!-- 以下内容并入自 `metaMDBG-vs-multik.md`（2026-08-15 文档合并） -->

## 10. metaMDBG 与 multik 实现对比

> 基于完整源码阅读：metaMDBG 1.4（`metaMDBG-metaMDBG-1.4/`，C++）与
> anchr multik（`src/libs/asm/multik.rs`，Rust）。两者是同一条
> "multi-k 迭代 + 图验证"路线的两种实现；本文对比实现层差异。

### 1. 定位

| 维度 | metaMDBG 1.4 | anchr `asm multik` |
|---|---|---|
| 目标 | 长读（HiFi/ONT）宏基因组组装 | 短读/长读通用，无 N 染色体 |
| 语言/并行 | C++20 + OpenMP，多进程调度 | Rust + rayon，单进程内存组合 |
| 建图空间 | minimizer 空间（density 0.005 采样） | 碱基 k-mer（FastK 字节键） |
| 命令形态 | 单二进制多子命令 + checkpoint 断点续跑 | 单命令，进程内完成 |

### 2. multi-k 迭代

| 维度 | metaMDBG | multik |
|---|---|---|
| k 语义 | k′-min-mer 的 minimizer 数 | 碱基 k-mer 长度 |
| 范围 | 4 → N50×1%（≈150 minimizers），每轮 +1 | auto 21/41/61/81/101/121（碱基），大步长 |
| 单轮跨度 | k/density（每轮 +200 bp 碱基） | +20~30 bp（碱基窗口） |
| 轮数 | 150+ | 4-6 |
| 上限 | 读长驱动（N50×density×2） | `Kmer::MAX_K=128` 硬限制 |

metaMDBG 的 k 是"minimizer 窗口长度"（间隔 200 bp），每轮 +1 覆盖 +200 bp
碱基跨度；multik 的 k 是碱基窗口，大步长（20-30 bp/轮）。**metaMDBG 150
轮 vs multik 4-6 轮**——metaMDBG 的 k′-min-mer 图节点天然长（minimizer
流滑窗），multik 需要跨接验证补偿大步长。

### 3. unitig 反馈与跨接验证

| 维度 | metaMDBG | multik |
|---|---|---|
| 反馈形式 | 图结构（unitigGraph_prev 加载 + solveEdges 验证边）+ 序列参与计数 | compute_links 边 + 序列参与计数 |
| 跨接验证 | doublet（pred 尾 k-1 minimizers + succ 延续 1）查 k-min-mer ≥2；短 unitig triplet | bridge_kmer（u 尾 k-1 碱基 + v 延续 1）查碱基 k-mer ≥2；短 unitig 跳过 |
| 边更新 | solveEdges 逐边验证（createDoubletNode） | 跨接验证 + recompact 后 compute_links 重算 |
| 方向解析 | minimizer 端匹配 | 实际碱基序列匹配（from_rc/to_rc 符号不用，实测匹配） |

metaMDBG 的 solveEdges 是"逐边验证 + doublet 压实成新节点"（图结构保持）；
multik 是"跨接 k-mer 查表 + recompact 合并链 + compute_links 重算边"
（更粗粒度）。**碱基空间的跨接 k-mer 比 minimizer doublet 更细**（单碱基
窗口 vs minimizer 窗口），验证更精确但边更多。

### 4. 渐进丰度过滤

| 维度 | metaMDBG | multik |
|---|---|---|
| 机制 | removeAbundanceNoQueue：t=1.1 起步、10% 步长到 maxAbundance，删 `abundance < t` + 每轮 recompact | progressive_filter：cutoff 上限 cov 中位 25%，只删分支/孤立，直链（主路径）永不删 |
| 图简化 | simplify()：superbubble（BFS 找出口 + collapseSuperbubble2 删低丰度分支，repeatSolver 保护）+ tip | 无 superbubble（分支靠跨接验证/丰度过滤） |
| 输出 | cutoff 快照（每个 cutoff 存图，generateContigs3 从高到低倒序输出） | 被删 unitigs 进 dropped（独立输出） |
| 时机 | 每轮 contig 阶段做 | 只在最终做一次（迭代轮只 recompact，不删丰度） |

**关键差异**：metaMDBG 删到 maxAbundance（依赖"主路径是单高丰度节点"的
宏基因组假设）；multik 用中位 25% + 直链保护（单菌株覆盖波动不误删主路径
——G37 重复区 600× 曾导致误删）。metaMDBG 的 superbubble 简化（平行路径
选高丰度）multik 没有——multik 的分支由跨接验证 + 探针桥接处理。

### 5. 嵌合清理与防 misassembly

| 维度 | metaMDBG | multik |
|---|---|---|
| 嵌合清理 | removeUnsupportedUnitigs：内部 k-min-mer 缺失即删整条 | remove_unsupported：内部 k-mer 缺失 <2% 容错（覆盖波动不误删） |
| 防错连 | RepeatRemover：fragment 按 unitig 边界切 → 覆盖均值 → 2×source 判重复 → 桥接 reads（nbBridgingReads≠0）才连 | bridge_filter（unitig 间 60-mer 探针 ≥2）+ split_by_bridge（unitig 内部 60-mer 窗口切分） |
| reads 映射 | ReadVsContigMapper（minimizer 索引 + chaining，容错） | 完美 60-mer 探针（无容错 chaining） |

metaMDBG 的 RepeatRemover 在最终 contig 上按 unitig 边界切 fragment 再
验证桥接；multik 在 unitig 图合并前用探针验证（unitig 间）+ 切分
（unitig 内部）。**multik 的 split_by_bridge 是 metaMDBG 没有的**——因为
multik 的 recompact 会把错连固化成单 unitig（G37 89,411），必须切内部；
metaMDBG 的 unitig 是图路径（内部无错连），只在 contig 层处理。

### 6. 性能与规模

| 维度 | metaMDBG | multik |
|---|---|---|
| 计数 | 外部分区（nbBases/20Gb 分区，clamp [cores, 5000]），disk scale-out | supermer 两段排序（内存） |
| minimizer 提取 | 一次（convertReadsToMinimizerSpace），每轮只滑窗 | 每轮 count_at 全量 supermer（未复用） |
| removeUnsupported | 节点级查表 | O(序列总长×k) 逐窗口编码+查表（瓶颈，见 multik-complexity.md） |
| G37 实测 | —（长读数据，无直接对照） | ~4-5 s（155k reads，含探针验证） |
| 1 Mb 合成 | — | ~10 s / 816 MB |

metaMDBG 的"前期抹除计算复杂度"（minimizer 一次 + 节点级处理）在 multik
未完全实现：count_at 每轮全量、remove_unsupported 序列级扫描（基准确认，
见 `benchmarks/multik-complexity.md`）。

### 7. 借鉴与差异总结

**借鉴（metaMDBG → multik）**：
1. multi-k 迭代 + unitig 反馈（核心）；
2. 跨接验证（doublet → bridge_kmer，k-mer 计数选边）；
3. 渐进丰度过滤（removeAbundanceNoQueue → progressive_filter，含
   recompact）；
4. 嵌合清理（removeUnsupportedUnitigs → remove_unsupported）；
5. 桥接 reads 防错连（RepeatRemover → bridge_filter + split_by_bridge）。

**multik 的简化/差异**：
1. 碱基 k-mer 空间（无 minimizer 采样）——k 上限 128 硬限制 vs metaMDBG
   读长驱动；
2. 大步长 4-6 轮 vs metaMDBG 150 轮（短 unitig 大步长需跳过验证/探针
   补偿）；
3. 渐进过滤用中位 25% + 直链保护（单菌株不误删）vs metaMDBG 删到
   maxAbundance（宏基因组假设）；
4. remove_unsupported 容错 <2%（覆盖波动）；
5. split_by_bridge 切 unitig 内部（metaMDBG 无——其 unitig 内部无错连）；
6. 完美探针 vs minimap2 容错映射（multik 假设 reads 干净/unitig 精确）。

**metaMDBG 有而 multik 未做**：
1. superbubble 简化（平行路径选高丰度）——multik 靠跨接验证 + 探针；
2. cutoff 快照分级输出——multik 直接输出 + dropped；
3. 容错 reads 映射（minimap2/chaining）——multik 完美匹配；
4. checkpoint 断点续跑、外部分区计数——multik 单进程内存。

### 8. 结论

两者是同一核心思想（multi-k 迭代 + 图验证）的两种实现：metaMDBG 面向
长读宏基因组（minimizer 空间、150 轮、maxAbundance 假设、容错映射），
multik 面向通用/无 N（碱基空间、大步长、单菌株保护、完美探针）。multik
的两个独有设计（直链保护、split_by_bridge 内部切分）解决的是碱基空间/
单菌株特有的问题（覆盖波动误删、recompact 固化错连），是 metaMDBG 不
需要的。跨接验证与丰度过滤是共同骨架，实现细节因空间（minimizer vs
碱基）而异。

## 2026-08-17 all-masters 单次调用重构

### 动机

模板曾按 master 逐个调用 multik（每个 master 一条命令，validated by
larger ks）：M 个 master、K 个 k 时 reads 计数次数为 O(K²/2)（每个 master
在其验证链的每个 k 全量重计 reads），且每次调用重复解压/解析输入。G37
（K=6）串行 25.2 s；MG1655（K=9）串行 7:23（release、-p 4）。

### 设计（`assemble_all_masters`，k-major 顺序）

单次调用内每个 k 都是一个 master（`--all-masters`），按 k 升序迭代：

1. 每个 k：建该 k 的 reads 计数表一次，供 (a) master-k 自身的 pass 0，
   (b) 所有更早 master 在该 k 的验证 round 共享（master 自身 unitig 单独
   计数，`SumView` 查表求和 = 联合计数）。
2. round 之间各 master 状态独立 → `par_iter_mut` 并发跑（图遍历是串行
   代码，正好填满线程池；输出与串行逐字节一致）。
3. `--use-guide`：第一阶段 ks[0] 跑完整阶梯产出 guide（finalize 后按
   solid 阈值重复为伪 reads）；第二阶段其余 master 每 k 直接对
   reads+guide×reps 联合计数一次（与"缓存表+merge_counts"算术等价）。

### 内存教训

第一版把第一阶段每 k 的表全部缓存供第二阶段 merge 复用：K 张表同时存活，
MG1655 (K=11, auto 51..251) 峰值 26.6 GB，5 组并发直接打爆 88 GB 机器。
修正：表即建即弃（第二阶段重建联合表），峰值 ~2 张表（MG1655 单进程
实测 5.7 GB）。代价是 guided 模式 reads 计数 ~2K 次而非 K 次。

### auto 阶梯（`auto_ks`）

旧公式 `k_min=clamp(N50/3,31,51)、step=clamp(N50/100,20,30)、
k_max=min(0.8*N50,256)` 对 MG1655 (N50=339) 产出 51..251：k≥211 的
master 被残余错误打碎（一个错误杀死所有覆盖它的窗口，k≈0.6×N50 时几乎
每个窗口都覆盖某个错误），guide 也救不回（N50 9.4K、5 mis）。改为固定
验证阶梯 `31,41,51,61,71,81,101,121,128,160,192` 截断于
`clamp(N50/2, 81, 192)`（MG1655→31..160，G37→31..192，150bp→31..81），
数据集调优用显式 `--kmer`。

### 验证结论（终版，2026-08-17；release）

- **MG1655 5 组 anchor 合并：三变体全部 0 mis**（guide 31..192 / 无
  guide 31..192 / auto 默认 31..160，当时 N50 读数 ~79.6K、GF 98.85–98.88，
  注：N50 读数受 /tmp 实验脚本跨组缺 `--unitigs` 的测量错误影响，见下）。
  guide 逐指标几乎相同 → 模板定型为 `--all-masters` 无 guide。
- **G37 7 组全流程：0 mis**（auto 31..192，guide 与无 guide 输出完全
  一致；N50 55.2K 与旧 multik 链持平、GF +0.95）。单组 auto N50 121K。
- 单组对照无回归：MG1655 新旧同 k（31..128）质量完全一致（均 2 mis、
  GF 99.40、N50 112,514）——重构质量中性；单组 0 mis 不可达，0 mis 一直是
  多组 anchor 投票的产物。
- 计时/内存（单组、-p 8）：MG1655 auto 无 guide 140–172 s / 峰值
  6.6 GB（40×）、10.8 GB（80×）；guide ~360 s（233 s，run>=2 移除后）。
  G37 单组 36.5 s→28.3 s（rounds 并行，输出逐字节一致）、869 MB。
- 数字详情见 `results/model_org.md` 2026-08-17 两节。

### 2026-08-17 下午：run>=2 过度剪除修复 + N50 归因终局

1. **run>=2 过度剪除**（commit 9f5fca7 引入的连续不支撑窗口检查）：
   真实覆盖凹陷出现连续 ≥2 个不 solid 窗口，把 MG1655 5 组链的长
   unitig 整条删掉（N50 124K→79.6K）。移除该检查，回到"容忍
   `n_kmers/50` 个孤立缺失窗口"语义；嵌合交给 reads-bridge 验证。
   修复后 N50 96.4K、0 mis。
2. **N50 96.4K vs 基线 124.0K 归因**（曾误判为验证密度）：
   真因是 /tmp 实验脚本跨组 anchors 合并漏了 `--unitigs`——anchors
   被当 reads 重新走 S0 组装切碎。补上后（multik 输出不变）：
   **90 contigs / N50 118,731 / GF 99.072%（+0.55 pp）/ 0 mis**，
   达到并部分超越 08-16 基线。正式模板 `7_merge_anchors.tera.sh`
   一直带 `--unitigs`，流程无此 bug。
3. **三项 A/B 全部中性**（5 组 anchor 投票抹平 multik 层差异）：
   验证密度（every k vs every-third，仅一处 166 kb 重复区序列不同、
   QUAST 逐指标相同）、`--use-guide`、last-k cut（单轮 master 的
   remove_unsupported）。实验用 `--validate-step` 临时开关已移除，
   cut 语义保留（assemble_one 与 --all-masters 一致）。
4. 单组 anchor 层完全等价：08-16 链 N50 112,781（95 条）vs 新链
   112,781（96 条）。
