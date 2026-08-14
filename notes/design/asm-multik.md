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
* `design/asm-assemble.md` §8（`asm unitig` 语义与 L: 边）
* `pgr: libs/kmer/supermer.rs`（FASTA 默认计数路径）
