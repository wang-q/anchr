# metaMDBG（1.4）：minimizer-space de Bruijn 图宏基因组组装器（源码分析）

> 2026-08-13 整理、2026-08-17 对照源码逐项复核修订，纯源码分析（`metaMDBG-metaMDBG-1.4/`，版本 `1.4`）。metaMDBG 是
> 面向**长而准的宏基因组 reads**（PacBio HiFi、Nanopore R10.4+）的组装器，论文
> [High-quality metagenome assembly from long accurate reads with metaMDBG]
> (Nature Biotechnology 2023)，作者 Gaëtan Benoit、Rayan Chikhi、Christopher
> Quince 等。**与 pgr 的关系**：它是 rust-mdbg（minimizer-space DBG）的宏基因组
> 工程化版本，与 pgr 的 `asm unitig`（bcalm 移植）同属 k-mer/unitig 路线，但其
> 核心创新——**local progressive abundance filter（用丰度替代气泡解析处理菌株
> 多样性）**——正好回应 pgr 讨论中"气泡不如不处理"的直觉。
> **首要借鉴点（2026-08-14 重读补充）**：对 anchr 当前最有价值的是 **multi-k′
> 迭代 + 跨接验证**——每轮 k 递增后，用当前 k 的 solid k′-min-mer（丰度 ≥ 2）
> 验证上一轮 unitig 图里每条相邻 unitig 连接（doublet/triplet），不支持的边剪
> 掉、支持的压实成新 unitig。这就是**连接 unitig 时对 bubble（分叉）的选择**：
> 更大 k 的计数天然区分低丰度菌株分支与主路径。该机制是 DBG/unitig 组装通用的
> "多轮 k 精化"模式，**与 OLC 无关**；OLC RepeatRemover（§4.2.1/§9）只是
> metaMDBG 的另一个借鉴点，不是重点。
> **与 OLC 的连接（2026-08-12）**：pgr `asm olc` 已落地，metaMDBG 的
> 渐进丰度过滤与 RepeatRemover 直接映射其 v1 的"覆盖度证据 repeat
> breaking"（见 §9）。

## 1. 概况

- **定位**：从长 reads（HiFi/ONT R10）构建宏基因组 contigs，输出带
  `length= coverage= circular=` 头信息的 FASTA；可选输出 GFA 组装图。
- **两种数据模式**：`getParamsHifi()`（ONT 模式默认开 read correction、HiFi
  默认关）；nanoMDBG（ONT R10 simplex）方法已集成进 1.4。
- **与 BCALM 的关系**：论文作者含 Rayan Chikhi（bcalm 作者），但实现**不基于
  BCALM/GATB**——BCALM 的 k-mer 是碱基 k-mer，metaMDBG 的节点是 **k′-min-mer
  （连续 minimizer 序列）**，走 rust-mdbg 的 minimizer-space DBG 路线。
- **语言/构建**：C++20 + OpenMP，CMake；内嵌 `ext/`（minimap2、spoa、htslib、
  TurboPFor 整数压缩），**全部 vendored，无外部运行时依赖**。
- **命令形态**：单一二进制 + 子命令（`asm` 为主，`graph`/`contig`/`toMinspace`/
  `readSelection`/`readCorrection`/`toBasespace`/`gfa` 等为底层步，见 §3）。
- **checkpoint 断点续跑**：每步写 `tmp/checkpoints/<step>.checkpoint` 文件，
  重跑同一命令自动跳过已完成步骤（README 明示，`AssemblyPipeline.hpp` 里
  `createCheckpoint`/`isCheckpoint` 实现）。

### 1.1 核心算法与流程总览（先读这节）

**一句话**：reads → FracMinHash 式 minimizer 采样（密度 0.5%）→ 以
**k′-min-mer（连续 k 个 minimizer）** 为节点建 minimizer-space de Bruijn 图 →
**multi-k 迭代**（k 从 4 到 ≈N50×密度×2，每轮 +1）：跨轮用更大 k 的 solid
k′-min-mer 验证上一轮 unitig 连接（跨接验证，§4.1.1），同轮内渐进丰度过滤
（1.1× 起步、~10% 步长逐级删低覆盖 unitig，§6.3）→ 最终轮 minimap2+POA
重建碱基序列并抛光。

```
reads（HiFi / ONT R10）
  → readSelection：FracMinHash minimizer 采样（hash < density×maxHash）
      （HiFi 同聚体压缩；ONT 先 minimap2 read 纠错）
  ── multi-k 迭代（firstK=4 … lastK≈N50×0.005×2，每轮 k+1）──
  ① createGraph(k)：计数 solid k′-min-mer（reads + 上轮 unitig 序列）
      + computeNextUnitigGraph：加载上轮 unitig 图
          solveEdges               跨接 doublet（前 1 + 后 k−1）查表选边/压实
          removeUnsupportedUnitigs 内部 k′-min-mer 无支撑 → 整条删
          solveSmallUnitigs        单 k′-min-mer 小节点 → triplet 两两重连
  ② generateContigs(k)：ProgressiveAbundanceFilter
        simplify（superbubble/tip）⇄ removeAbundanceNoQueue（1.1×起、
        ~10% 步长、删后 recompact）→ 每个新 cutoff 存图快照
        generateContigs3：从最高 cutoff 倒序消费快照 → contigs.nodepath
  ③ toMinspace：nodepath → minimizer 序列 → unitig_data.txt 反馈下一轮
  ── 最终轮（isFinalPass）──
  derepSmall（包含去重）→ removeOverlaps（首尾重叠截断）
  → removeRepeats（reads 桥接证据断开未桥接重复）
  → toBasespace：minimap2 reads→contig + POA 抛光 → contigs.fasta.gz
```

| 核心块 | 机制 | 详见 |
|---|---|---|
| minimizer 采样 | FracMinHash 式阈值采样（非窗口最小），首尾各丢 1 个 k-mer，重复 minimizer 剔除 | §5 |
| k′-min-mer 计数 | 外部分区（nbBases/20Gb，clamp [cores,5000]）+ singleton rescue | §6.1 |
| 跨接验证 | doublet/triplet 用当前 k 的 solid 表选边——multi-k 迭代核心 | §4.1.1 |
| 渐进丰度过滤 | 1.1× 起步 ~10% 步长逐级删低覆盖 unitig + recompact；cutoff 快照倒序输出 | §6.3 |
| unitig 丰度 | 组成 k′-min-mer 丰度向量取中位数（merge 保持有序不变量） | §6.2 |
| 碱基重建 | 内嵌 minimap2 分片映射 + spoa POA 抛光；覆盖度裁两端 75 bp | §7.2 |

## 2. 仓库结构

```
metaMDBG-metaMDBG-1.4/
├── src/
│   ├── MdbgAssembler.cpp        # 主入口：子命令分发
│   ├── Commons.hpp              # 全局类型/常量/工具（8378 行）
│   ├── pipeline/AssemblyPipeline.hpp   # ★ asm 主流水线（multi-k′ 调度）
│   ├── graph/
│   │   ├── CreateMdbg.{cpp,hpp}        # ★ k′-min-mer 计数 + MDBG 建图
│   │   ├── ProgressiveAbundanceFilter.hpp  # ★ 图简化 + 渐进丰度过滤
│   │   ├── Graph.hpp              # UnitigGraph2 数据结构
│   │   ├── GenerateContigGraph.hpp # final 轮的 repeat solver（嵌合清理）
│   │   ├── GraphPOA.hpp / GraphSimplify.hpp / GenerateGfa.hpp / GfaParser.hpp
│   ├── assembly/GenerateContigs.hpp   # ★ unitig 路径 → contigs.nodepath
│   ├── toBasespace/
│   │   ├── ToMinspace.hpp        # contig nodepath → minimizer 序列（反馈下一轮）
│   │   ├── ToBasespace2.hpp      # ★ minimizer contig → 碱基序列（minimap2+POA 抛光）
│   │   ├── RepeatRemover.hpp / OverlapRemover2.hpp / DerepSmallContigs.hpp
│   │   ├── ContigPolisher.hpp / ContigTrimmer.hpp / ReadVsContigMapper.hpp
│   ├── readSelection/
│   │   ├── ReadSelection.hpp     # reads → minimizer 表示（含 MinimizerParser 入口）
│   │   ├── ReadCorrection.hpp    # ONT 的 reads 纠错（minimap2 内嵌）
│   ├── contigFeatures/           # KminmerCounter/KmerCounter（当前未参与主流程）
│   └── utils/                    # args.hxx / BooPHF / BloomFilter / edlib / MurmurHash3
└── ext/                          # vendored: minimap2, spoa, htslib, TurboPFor
```

> 源码总量约 9.7 万行（`wc -l` 覆盖 `src/` 下全部 `.cpp/.hpp/.h`，含 `utils/`
> 的 phmap/edlib/MurmurHash3），大头是 `Commons.hpp`（8378）、`CreateMdbg.cpp/hpp`、
> `ProgressiveAbundanceFilter.hpp`、`ToBasespace2.hpp`、`ReadCorrection.hpp`。
> 文件多、注释少（大量被注释的旧代码），核心算法集中在 ProgressiveAbundanceFilter
> 与 ToBasespace2。

## 3. 命令入口（`MdbgAssembler.cpp`）

| 命令 | 实现 | 作用 |
|---|---|---|
| `asm` | `AssemblyPipeline` | 完整组装（唯一面向用户的命令） |
| `graph` | `CreateMdbg` | 计数 k′-min-mer，建 unitig 图 |
| `contig` | `GenerateContigs` | 图简化（superbubble/tip/丰度过滤），导出 nodepath |
| `toMinspace` | `ToMinspace` | nodepath → minimizer 序列 |
| `readSelection` | `ReadSelection` | reads → minimizer 表示 |
| `readCorrection` | `ReadCorrection` | reads 纠错（ONT） |
| `toBasespace` | `ToBasespace2` | minimizer contig → 碱基序列 + 抛光 |
| `toBasespaceGfa` / `gfa` | `ToBasespaceGfa` / `GenerateGfa` | GFA 输出 |
| `derepSmall` / `removeOverlaps` / `removeRepeats` | 对应类 | 最终后处理 |
| `map` | `MappingContigToGraph` | contig → 图映射（调试用） |

## 4. asm 主流水线（`AssemblyPipeline::execute_pipeline`）

`execute()` 先跑 `convertReadsToMinimizerSpace()`，再跑多轮
`executePass(k, prevK, pass)`；**最终后处理（`derepSmallContigs → removeOverlaps →
removeRepeats → toBasespace`）不是循环结束后的独立步骤，而是在最后一轮
`executePass(_lastK, ...)` 的 `isFinalPass` 分支内联执行**（`AssemblyPipeline.hpp:1111`
起），输出 `contigs.fasta.gz`。首轮/中间轮只生成 unitig 序列反馈下一轮，不做碱基重建。

### 4.1 multi-k′ 迭代

```cpp
_firstK = 4;
_lastK = Commons::computeLastK(_minimizerDensityAssembly, readStats._n50ReadLength, _firstK, _maxK);
// lastK = N50ReadLength * density * 2.0（10 kb 读长、density 0.005 → lastK ≈ 100）
// 每轮 k += Commons::getMultikStep(k)  // 1.4 里恒为 1
// 最后一轮再显式补一次 executePass(_lastK, ...)
```

每轮 `executePass(k)`：

1. `createGraph(k, pass)`：`graph` 子命令。**首轮**（pass 0）从
   `read_data_corrected.txt` 计数全部 k′-min-mer（`--min-abundance`，默认 0 =
   rescue 模式）；**后续轮**加载上一轮 refined abundance
   （`loadRefinedAbundances`），并从 `read_data_corrected.txt` **和**
   `unitig_data.txt`（上一轮的 unitig 序列）一起计数——这是
   "unitig 反馈进下一轮"的**序列部分**；此外后续轮在计数完成后还会调用
   `computeNextUnitigGraph`（`CreateMdbg.cpp:3712`）加载上一轮 unitig 图做
   **跨接验证**（图结构部分，见 §4.1.1）——这才是 multi-k 迭代的核心机制。
2. `generateContigs(k, pass)`：`contig` 子命令。加载 unitig 图 →
   `ProgressiveAbundanceFilter::execute`（图简化，见 §6）→ `generateContigs3`
   从各 cutoff 快照生成 `contigs.nodepath`。
3. `toMinspaceContigs(...)`：把 nodepath 转回 minimizer 序列，写入
   `unitig_data.txt`（非 final）或 `contig_data_init.txt`（final）。
4. 每轮结束 `savePassData(k)` 把 `assembly_graph.gfa` + `parameters.gz` 备份到
   `pass_k<k>/`（`AssemblyPipeline.hpp:1435`）。"unitig 反馈"的中间状态落盘：
   `graph` 步在 k>firstK+1 时调 `computeNextUnitigGraph` 加载上一轮
   `unitigGraph_prev.*`（`CreateMdbg.cpp:3744`）；`contig` 步写
   `unitigGraph.nodes.refined_abundances.bin`（`GenerateContigs.hpp:764`），
   供下一轮 `graph` 的 `loadRefinedAbundances`（`CreateMdbg.cpp:391`）复用。

**k′ 的语义与长度换算**（`AssemblyPipeline::writeParameters`,
`AssemblyPipeline.hpp:1479`）：multi-k 里的 `k` 是**一个 k′-min-mer 包含的
minimizer 个数**，不是碱基长度。换算关系写在 `parameters.gz` 里，供各子命令
`Parameters::load` 复读：

```cpp
minimizerSpacingMean = 1 / assemblyDensity;   // 相邻 minimizer 的平均间距(碱基)
kminmerLengthMean   = minimizerSpacingMean * (k-1);   // 一个 k-min-mer 的期望碱基跨度
kminmerOverlapMean  = kminmerLengthMean - minimizerSpacingMean; // 相邻 k-min-mer 重叠
```

故 `k` 每 +1，k′-min-mer 期望碱基长度约 +`1/density`（assembly density 0.005 →
每轮约 +200 bp）。这就是"unitig 反馈 + k 递增"实现**渐进长单元化**
（longer k-min-mer → 更多直链、更少分支）的机制。

**assembly graph 导出节奏**：`--gen-graph` 默认在第 11 轮（`_nextGenGraphIteration=11`，
之后每 +10）导出一次 GFA（`doesGenerateAssemblyGraph`，
`AssemblyPipeline.hpp:831`；首轮过大不导出）——用"隔轮导出"控制磁盘/内存，
pgr 的 `pl` 管道若多轮组装可参考。

### 4.1.1 跨接验证：连接 unitig 时对 bubble 的选择 ★

`CreateMdbg::computeNextUnitigGraph`（`CreateMdbg.cpp:3712`，仅
`k > firstK+1` 时执行）是 multi-k 迭代的核心：**上一轮的 unitig 图结构本身
被加载进当前轮**（`unitigGraph_prev.*`），然后按当前（更大的）k 重新验证每条
相邻 unitig 连接。三个子步骤：

```cpp
solveEdges(unitigGraph);            // 1. 逐条验证 unitig 间连接
removeUnsupportedUnitigs(unitigGraph); // 2. 删内部 k′-min-mer 无支撑的 unitig
solveSmallUnitigs(unitigGraph);     // 3. 处理长度恰为单个 k′-min-mer 的小 unitig
```

1. **`solveEdges`**（`CreateMdbg.cpp:3903`）：对每条相邻 unitig 边
   `pred → succ`，用 `getDoublet2` 构造**跨接 k′-min-mer（doublet）**：
   前驱末尾 1 个 minimizer + 后继开头 k−1 个 minimizer（长度恰为**当前** k）。
   查当前轮计数表 `_mdbgNodesLight`（丰度 ≤ 1 的不进表，见下）：
   - `isEdgeSupported`（`CreateMdbg.hpp:3229`）：doublet 存在 → `createDoubletNode`
     （`CreateMdbg.cpp:4063`）把 doublet 压实成新 unitig 节点，替换原直接边
     （`pred → edgeNode → succ`）；
   - doublet 不存在 → `removeSuccessor` 删除该边。
   无论支持与否原边都会被移除（`CreateMdbg.cpp:3997-3998`），区别只在于
   支持时经由新 edge node 重建连接、不支持时直接断开。
2. **`removeUnsupportedUnitigs`**（`CreateMdbg.cpp:4138`）：unitig 内部所有
   当前 k 的 k′-min-mer 必须在 `_mdbgNodesLight` 中存在，否则整条 unitig 删除
   ——上一轮拼错/嵌合的 unitig 在更大 k 下内部 k′-min-mer 不 solid，被剔除。
3. **`solveSmallUnitigsSub2`**（`CreateMdbg.cpp:4489`）：对上一轮长度恰为一个
   k′-min-mer 的小 unitig（`_nbMinimizers == _kminmerSizePrev`），构造前向
   triplet（前驱末尾 1 + 小 unitig 全部）与后向 triplet（小 unitig 全部 +
   后继开头 1），各自也是当前 k 长度；分别查计数表得 `supportedPredecessors` /
   `supportedSuccessors`，支持的前驱 × 后继**两两建边**，最后删掉小节点本身。

**`_mdbgNodesLight` 的丰度传播**（`IndexKminmerFunctor`，
`CreateMdbg.hpp:951`）：后续轮从 reads 和 `unitig_data.txt` 解析当前 k 的
k′-min-mer 时，丰度取"其组成的前一轮 k−1 k′-min-mer 丰度"相邻两者的**最小值**
（`getAbundance`，`CreateMdbg.hpp:1006`，`nbSubKminmers=2`），≤ 1 的直接
`continue` 不进表（`CreateMdbg.hpp:1445`）——即表里只存 **solid（丰度 ≥ 2）**
的当前 k 的 k′-min-mer。

**这就是"连接 unitig 时对 bubble 的选择"**：上一轮 k 的 unitig 图在分歧点有
多个候选后继（bubble）；进入更大 k 后，跨接每个候选的 k′-min-mer 互不相同，
只有被 reads（+ unitig 序列）以丰度 ≥ 2 支撑的才保留连接，低丰度分支的边被
剪掉、高丰度主路径保留并被压实成更长 unitig。**与 §6.3 渐进丰度过滤互补**：
后者是同一轮内按 unitig 丰度删低覆盖节点（丰度阈值过滤），前者是**跨轮按
更大 k 的 solid k′-min-mer 选边**（结构验证）——两者共同实现"渐进长单元化
+ 丰度驱动消歧"。这也是"unitig 反馈进下一轮"的图结构部分：上一轮 unitig
不仅作为序列参与计数，其图结构还决定下一轮哪些连接候选被验证。

### 4.2 最终后处理（isFinalPass）

调用链（`AssemblyPipeline.hpp:1111-1145`）：`generateContigs → toMinspaceContigs →
derepSmallContigs → removeOverlaps → removeRepeats → toBasespace`。后三步是**纯碱基
空间**的 contig 后处理，不依赖 minimizer 图：

- **`derepSmallContigs`**（`DerepSmallContigs.hpp`）：候选小 contig 与其余 contig
  做 minimap2 比对，**若被更高覆盖的 contig 完全包含则丢弃**。宏基因组里低丰度
  物种常被高丰度物种的 contig 片段包含，这一步是"包含去重"。
- **`removeOverlaps`**（`OverlapRemover2.hpp`）：contig 两两 minimap2 比对
  （`countContigs → indexContigs → mapContigs`），找首尾互相覆盖的 overlap，**截断
  较长 contig 被覆盖的重叠部分**（不删整条）。粒度粗、不需要 reads 证据。
- **`removeRepeats`**（`RepeatRemover.hpp`）★：`ReadVsContigMapper` 把 reads 映射
  回 contig，找未桥接的重复位点并断开。完整流水线见下：

#### 4.2.1 `removeRepeats`（`RepeatRemover`）—— 桥接 reads 证据 ★

核心流水线（`RepeatRemover.hpp:220-233`）：

```
alignReads → fragmentContigs → computeFragmentsCoverage → breakUnbridgedRepeats
```

1. **`alignReads`**：`ReadVsContigMapper` 把 reads 映射回 contig（见 4.2.2），产出
   `(contigIndex, contigStart, contigEnd)` 对齐列表。
2. **`fragmentContigs`**（`RepeatRemover.hpp:557-580`）：把每条 contig 按"组成
   k-min-mer 所属 unitig 是否变化"切成**片段**（`FragmentFunctor`，
   `RepeatRemover.hpp:682-703`）——片段边界 = 图里 unitig 的边界，等价于把 contig
   拆回"unitig 串"。
3. **`computeFragmentsCoverage`**（`RepeatRemover.hpp:761-785`）：每片段 coverage =
   组成 k-min-mer 丰度的**均值**（注意：这里用均值，unitig 丰度用中位数，两处
   语义不同，见 §10.14）。
4. **`breakUnbridgedRepeats`**（`RepeatRemover.hpp:950-`）：
   - **`computeBridgingReads`**（`RepeatRemover.hpp:1328-1370`）：对每条 reads 对齐，
     若该对齐同时覆盖 **≥2 个相邻 fragment**（三种情况：contained / 左叠到右 /
     右叠到左），则这些 fragment 两两之间 `_nbBridgingReads++`——"这条 read 把两个
     片段连起来"的证据。
   - **`getCovPath` 双向延伸**（`RepeatRemover.hpp:1374-1462`）：从每个 source
     fragment 出发，`minRepeatCoverage = sourceCoverage * 2.0`；延伸规则：下一片段
     coverage ≥ 2×source 视作重复片段跳过；相邻片段直接相连；否则**必须有桥接
     reads（`nbBridgingReads != 0`）才连**，否则路径在此断开。
   - 路径分配 `_finalContigIndex`，`nbContigsFinal > 1` → `isCircular=0` 并沿
     finalContigIndex 边界切分 contig（`RepeatRemover.hpp:1256-1323`）。

#### 4.2.2 `ReadVsContigMapper`（桥接证据的来源，`ReadVsContigMapper.hpp`）

- contig 建 **minimizer 索引**：每个 contig 的 minimizer 哈希 → `MinimizerPosition`
  （含 contig 内坐标），排序后建 `_minimizerLookupTable`（哈希 → 位置区间）。
- `mapRead2`：对 read 的每个 k-min-mer 查表得候选锚点 → 按 `(contigIndex, 方向)`
  分组 → `chainAnchors` 做 chaining（按坐标排序、重叠过滤）→ 输出
  `ReadMapping2(contigIndex, contigStart, contigEnd)`。
- **内存/精度取舍**：低密度 minimizer（0.5%）锚点少、链短，比全 k-mer 索引省内存，
  适合"证据回放"这类不要求逐碱基精度的用途。

> 对 pgr 的意义：`removeRepeats` 的"read 同时覆盖多片段 → 连接证据；无证据且覆盖
> 异常 → 断开"正是 `design/asm-olc.md` 待决的 repeat breaking 成熟实现（见 §9.2）；
> `ReadVsContigMapper` 的 minimizer chaining 则可作为 pgr `asm map`（全 k-mer 完美
> 匹配）在容错/长读场景下的对照（pgr 不引新依赖，仅作语义参考）。

### 4.3 关键参数与默认值（`AssemblyPipeline.hpp:100-284`）

`asm` 的可调参数集中在 `AssemblyPipeline` 构造函数里，均带 `args` 默认值：

| 参数 | 默认 | 说明（源码行） |
|---|---|---|
| `--kmer-size` | 15 | minimizer 长度，**cap ≤ 16**（`:202`） |
| `--density-assembly` | 0.005 | 组装 minimizer 采样密度（`:117`） |
| `--density-correction` | 0.025 | 纠错 minimizer 密度（`:125`） |
| `--max-k` | 0(=自动) | 覆盖 `computeLastK` 的结果（`:118,211-214`） |
| `--min-abundance` | 0(=rescue) | k′-min-mer 最低丰度（`:119`） |
| `--max-bubble-length` | 50000 | superbubble 弹出上限（`:120`） |
| `--max-tip-length` | 50000 | tip 剪除上限（`:121`） |
| `--min-read-quality` | 0.0 | read 平均质量过滤（`:106`） |
| `--min-contig-length` | 50 | 输出 contig 最短长度（`:107`，floor 50 `:274`） |
| `--min-contig-coverage` | 1 | 输出 contig 最低覆盖（`:108`，floor 1 `:275`） |
| `--skip-correction` | off | 跳过 read 纠错（`:128`） |
| `--max-memory` | 8GB | toBasespace 内存预算（见 §7.2，`ToBasespace2.hpp:226`） |
| `-t/--threads` | 默认 | OpenMP 线程数 |

**数据模式差异**（`getParamsHifi:292` / `getParamsNanopore:309`）：

| 属性 | HiFi | Nanopore(ONT R10) |
|---|---|---|
| read 纠错 | 关 | **开** |
| homopolymer 压缩 | **开** | 关 |
| readCorrectionMinIdentity | 0.99 | 0.96 |
| 抛光覆盖上限 `_usedCoverageForContigPolishing` | 50 | 100 |
| minimap2 预设 | `-x map-hifi` | `-x map-ont` |

**`--density-correction ≥ --density-assembly` 校验**：仅当 `_useReadCorrection`
（即 ONT 模式且未 `--skip-correction`）时执行，违反则打印 parser 并 `exit(0)`
（`:278-284`）。HiFi 默认不开纠错，故不受此限。

## 5. minimizer 提取（`ReadSelection` + `Kmer.hpp::MinimizerParser`）

### 5.1 编码与同聚体压缩

- 2-bit 编码（`DnaBitset`），`KmerModel` 滚动 k-mer；`EncoderRLE` 对 HiFi 做
  homopolymer 压缩（ONT 关）。
- `MinimizerParser(_minimizerSize, _minimizerDensity, ...)`：默认
  `--kmer-size 15`（cap ≤ 16），assembly density 0.005、correction 0.025
  （`AssemblyPipeline.hpp:117,125`）。**约束**：`--density-correction ≥
  --density-assembly`，违反即打印 parser 并 `exit(0)`（`AssemblyPipeline.hpp:278-284`，
  仅在 `--read-correction` 开启时校验）。

### 5.2 采样规则（FracMinHash 式"通用 minimizer"）

```cpp
_minimizerBound = minimizerDensity * maxHashValue;   // u_int64
// 对每个 k-mer：
u_int64_t kmerHash = MurmurHash3_x64_128(&kmerValue, sizeof(kmerValue), 42);
if(kmerHash < _minimizerBound){ minimizers.push_back(kmerValue); ... }
```

即**不是窗口内取最小 k-mer**，而是对每个 k-mer 哈希、保留低于阈值者——等价于
rust-mdbg 的 universal minimizer / FracMinHash 采样。这样 minimizer 是无序空间
均匀采样，密度 ≈ 0.5%（assembly）/2.5%（correction）。

> 与 pgr 的对照：pgr `asm map` 用全量 k-mer 索引（讨论过 minimizer/syncmer 但
> 结论是不优化）；metaMDBG 的 minimizer 密度极低（0.5%），是**组装图节点**，
> 不是比对种子，两者目的不同。

### 5.3 读取边界与"可重复 minimizer"过滤

- **首尾各丢弃一个 k-mer**（`Kmer.hpp:1395`）：`MinimizerParser::parse` 遍历
  `pos=_trimBps(=1)` 到 `kmers.size()-_trimBps`，即**每条 read 的第一个和最后
  一个 k-mer 不参与采样**（`_trimBps=1`，`:1362`），避免边界截断误差。
- **可重复 minimizer 剔除**（`Kmer.hpp:1434-1437`）：哈希命中阈值后，若该
  minimizer 出现在 `_isRepetitiveMinimizers`（全数据集上出现次数过高的低信息量
  minimizer）集合中则跳过。该集合由上游构建，用于抑制高重复区域导致的图污染。

## 6. 建图 + 图简化（`CreateMdbg` + `ProgressiveAbundanceFilter`）

### 6.1 k′-min-mer 计数（`CreateMdbg::createMDBG`）

- 分区数 `_nbPartitions = nbBases / 20Gb`，clamp 到 `[nbCores, 5000]`。
- `KminmerCounter`：把每条 read 的连续 minimizer 切成 k′-min-mer（`KmerVec`，
  canonical normalize，`hash128()` = MurmurHash3 128-bit），按 `hash % nbPartitions`
  写分区文件 → 分区内去重计数 → 合并写 `kminmerData_min.txt` +
  `kminmerData_abundance.txt`。
- **rescue 机制**（`--min-abundance 0` 默认）：首轮计数后把 `abundance==1` 的
  singleton k′-min-mer 单独"营救"一遍——对每条 read，若其 k-min-mer 中**多数是
  solid（丰度 >1）、少数是 singleton**，则这些 singleton 很可能是低覆盖基因组
  的真实 k-mer 而非测序错误，予以保留（`RescueKminmerFunctor`：read 内丰度中位
  数的 10% 为阈值，`median*0.1 <= 1` 才营救）。
- 节点表 `MdbgNodeMapLight`：`phmap::parallel_flat_hash_map<u_int128_t,
  DbgNodeLight>`（10 分区 + mutex）。
- 之后 `indexEdges → computeUnitigNodes → computeDeterministicUnitigs →
  indexUnitigEdges`，输出 `unitigGraph.nodes.bin` / `edges.successors.bin` /
  `nodes.abundances.bin` / `stats.bin`。

### 6.2 UnitigGraph2（`Graph.hpp`）

- `UnitigNode`：`_unitigName`、`_successors_forward/reverse`、`_nbMinimizers`、
  `_abundance`（float）、`_abundances`（每个组成 k-min-mer 的丰度向量）、
  `_unitigMerge`。
- **丰度语义**：unitig 丰度 = 组成 k-min-mer 丰度向量的**中位数**
  （`computeMedianAbundance`，`Graph.hpp:252-287`）；`recompact` 合并两个 unitig
  时把两个丰度向量 merge 后再取中位（`mergeWith`，`Graph.hpp:293-321`）。
  - **"向量始终有序"是硬不变量**：`computeMedianAbundance` 内部的 `sort` 被注释
    掉，中位数正确性依赖 `_abundances` 始终升序。三个保证点：`createEdgeNode`
    先 `std::sort(_abundances)` 再取中位（`CreateMdbg.cpp:5026-5027`）；图加载后
    `AbundanceSortFunctor` 统一排序（`Graph.hpp:651-653`）；`mergeWith` 用
    `std::merge`（要求两输入已排序）合并后重新中位数。**Rust 移植时须显式维护
    排序，或改用无需排序的顺序统计量**（详见 §10.13）。
  - 初始值：每个 k-min-mer 若不在 solid 表（`abundance==1`，可能是测序错误）
    记 1，否则记实际丰度——中位数对少数低覆盖读的污染不敏感。
- 长度估计 `getLength = (nbMinimizers-1) * _minimizerSpacingMean`（minimizer
  空间里没有碱基坐标，长度是期望值）。

### 6.3 渐进丰度过滤（`ProgressiveAbundanceFilter`）★

`execute` → `simplifyProgressive(functor)` 主循环：

```cpp
maxAbundance = min(图内最大丰度, 10000);
currentCutoff = 0;
while(true){
    isModification = simplify();          // superbubble + tip + repeat solver
    checkSaveState(currentCutoff);        // 每个新 cutoff 存一张图快照
    nbErrorRemoved = removeAbundanceNoQueue(maxAbundance, currentCutoff);
    if(!isModification && !nbErrorRemoved) break;
}
```

**`simplify()`**（图结构简化）：

- `SuperbubbleRemoverOld`：找 `nbSuccessors>1` 的节点做 BFS 找出口，`isSuperbubble`
  判定，`collapseSuperbubble2` 收集并删除低丰度分支（丰度高于
  `currentCutoff/0.25` 的受 repeat solver 保护不移除），邻接 unitig `recompact`。
- `TipRemover`：按 `_nbMinimizers` 升序队列删 tip。
- final 轮额外挂 `_repeatSolver`（`GenerateContigGraph`）做嵌合 unitig 清理。

**`removeAbundanceNoQueue`**（丰度渐进过滤，核心，`ProgressiveAbundanceFilter.hpp:2183-2343`）：

```cpp
float t = 1.1;                       // 阈值从 1.1x 起步
while(t < abundanceCutoff_min){
    currentCutoff = t;
    _maxAbundance = currentCutoff*2; // 每轮重置，供下游"高于 2×cutoff 受保护"判定
    // 遍历 _validNodes2（当前全部有效 unitig）：
    for(node : _validNodes2){
        if(node->_abundance >= t) continue;
        // 收集受影响邻居：前驱(正向) 与 后继(取反向) → recompactNodes
        recompactNodes.insert(前驱节点);
        recompactNodes.insert(后继节点的反向索引);
        removeNode(node);            // 真删除（标记而非物理释放）
    }
    // recompactNodes 按 BubbleSideComparatorRev 排序后逐个 recompact
    //   → 单前驱+单后继的邻居合并成更长 unitig，丰度向量 std::merge 后重取中位
    float newT = t * (1 + 0.1);            // 每次约 +10%
    float increaseStep = min(newT - t, 10);// 步长封顶 10
    t += increaseStep;
    if(_usingFunctor){               // functor 模式(contig 生成路径)
        if(isCutoffProcessed2 未含 currentCutoff){
            记录 currentCutoff; return 1;  // 每个新 cutoff 只删一轮，回外层再 simplify
        }
    }
    if(nbErrorsRemoved > 0) break;   // 有删除即退出，等外层 simplifyProgressive 循环
}
```

即**从丰度 1.1x 起步、按 ~10% 的步长逐步抬升阈值，每次移除低于阈值的 unitig
并对受影响邻居重新压实**。关键设计：

- **删除是"标记 + 收集邻居"而非就地重连**：一次遍历收集所有受波及的前驱/后继，
  排序后统一 recompact，避免一边删一边改邻接表造成的顺序依赖。
- **合并方向**：删掉的 unitig 的前驱直接参与压实，后继取**反向索引**参与——因为
  压实总是在节点"另一端"进行，反向表示保证合并方向一致。
- **交替收敛**：外层 `simplifyProgressive` 每轮先 `simplify()`（结构）再
  `removeAbundanceNoQueue`（丰度），两者都无改动才退出；functor 模式每新 cutoff
  只删一轮就返回，让外层重新 simplify——保证"删一批 → 简化一批"的节奏。

这是论文里 "local progressive abundance filter" 的实现：它不试图区分菌株气泡里的
正确路径，而是用丰度把低覆盖分支逐级删掉，让高丰度物种的主路径自然收敛——**与
pgr 讨论中"气泡经常引入不确定性、不如不处理"的直觉一致**。

**cutoff 快照（`dumpUnitigs`）**：每个新 cutoff 把当前 unitig 图导出到
`filter/unitigs_<idx>.bin`（node path + circular/repeat 标记 + 丰度），
`_cutoffIndexes` 记录 `{idx, cutoff}`。供 `generateContigs3` 从高 cutoff 到低
cutoff 倒序消费（见下）。

## 7. contig 生成 + 碱基空间重建

### 7.1 `GenerateContigs::generateContigs3`（`GenerateContigs.hpp:549-757`）

- 从 `_cutoffIndexes` **最后一个（最高 cutoff）倒序**读快照
  `filter/unitigs_<idx>.bin`：每个快照条目 = `{size, isCircular, isRepeatSide,
  contigAbundance(float), nbMinimizers, nodePath[]}`。
- 逐条判定是否接受（`GenerateContigs.hpp:575-641`）：
  - `_minUnitigAbundance = cutoff / 0.5`，`contigAbundance < _minUnitigAbundance`
    的路径跳过（`isValidContigAbundance`）；
  - `isContigAssembled(nodePath)`（`GenerateContigs.hpp:484-495`）：路径上**任一**
    unitig 已在 `_processedNodeNames`（被更高 cutoff 快照组装过）→ 整条路径跳过，
    避免重复输出；
  - final 轮 `_repeatedUnitigNames`（repeat solver 判定的重复 unitig）跳过。
- 圆形 contig 特殊处理：`isCircular` 时 `nbMinimizers += 1` 以在长度估算上闭合环
  （`GenerateContigs.hpp:642-644`），但 nodepath 本身不重复追加首节点。
- 组装成功后将路径上所有 unitig 写入 `_processedNodeNames`，并记录
  `_nodeNameAbundances[unitigName] = {contigAbundance, nbMinimizers}`
  （`GenerateContigs.hpp:725-733`）。
- 输出 `contigs.nodepath` + 写 `unitigGraph.nodes.refined_abundances.bin`
  （`dumpUnitigAbundances`，`GenerateContigs.hpp:759-777`）：refined abundance 取
  `ceil`，长 unitig（`nbNodes-kminmerSize+1 > kminmerSize`）保底 ≥ 2 —— 供下一轮
  `graph` 的 `loadRefinedAbundances`（`CreateMdbg.cpp:391`）复用。

### 7.2 `ToBasespace2`（minimizer → 碱基）

- 内嵌 minimap2（`mm_dbg_flag |= MM_DBG_NO_KALLOC`），HiFi `map-hifi`/`ava-pb`、
  ONT `map-ont`/`ava-ont`。
- 流程：reads 映射到 minimizer contigs → `partitionReads`（按内存分片，
  `minimapBatchSize = peakMemory/8`）→ `createBaseContigs` 每片读回 reads
  用 POA（内嵌 spoa）抛光 → 输出 `contigs.fasta.gz`（header
  `ctg<id> length= coverage= circular=yes|no`）。
- `--skip-correction` / `--min-contig-length 50` / `--min-contig-coverage 1`
  可过滤输出。

**输出 FASTA 头格式**（`Commons.hpp:2212-2222` `createContigHeader`）：
`>ctg<id> length=<len> coverage=<cov:两位小数> circular=yes|no`——与 pgr
`asm unitig` 的 `cov=` 头同形，可直接被解析复用。coverage 由 `ContigPolisher`
把 reads 映射回 contig、对每列算覆盖并**裁掉两端各 75 bp 后再平均**
（`ContigPolisher.hpp:658-684`），是稳健覆盖度而非简单平均。

**输出过滤**（`ContigPolisher.hpp:2785-2793`，三类之一即丢弃该 contig）：
coverage ≤ `--min-contig-coverage`；或 length < `--min-contig-length`；或
**length < 7500 且 coverage < 4**（一个额外的经验性低质量小 contig 剔除规则，
`_minContigLength=50` 之外的隐式硬阈值）。

**圆形闭合**：`GenerateContigs::generateContigs3` 对 `isCircular` 的 unitig
路径做 `nbMinimizers += 1` 以在长度估算上闭合环
（`GenerateContigs.hpp:642-644`），但 nodepath 本身不重复追加首节点。

## 8. 与 pgr 的对应/借鉴点

1. **multi-k 跨接验证 = 连接 unitig 时 bubble 的选择**（**首要借鉴**，
   2026-08-14 重读补充，见 §4.1.1）：每轮 k 递增后用当前 k 的 solid
   k′-min-mer（丰度 ≥ 2）验证上一轮 unitig 图每条相邻连接（doublet/triplet），
   支持的压实成新 unitig、不支持的剪边；同时 `removeUnsupportedUnitigs` 剔除
   内部 k′-min-mer 无支撑的嵌合 unitig。这给 pgr 的多 k 迭代（`asm unitig`
   多趟、SKESA clean_reads 式反馈）提供了一个**结构验证骨架**：每轮之间
   回放上一轮 unitig 边界处的 k-mer 计数即可选边，不需要额外启发式。与 OLC
   无关，是"多重多次 Kmer"的核心价值。
2. **丰度过滤替代气泡解析**（同轮机制，次优先）：`ProgressiveAbundanceFilter` 的
   "1.1x 起步、~10% 步长、边删边压实"策略可以直接映射到 pgr `asm contig`/
   `unitig` 的 `--min-coverage` 语义：不是单阈值一刀切，而是**多轮渐进 + 每轮
   重压实**，低覆盖菌株分支逐步被吞并。pgr 目前只有全局 `--min-coverage`，
   可考虑加"渐进模式"。
3. **cutoff 快照倒序输出**：metaMDBG 存多个 cutoff 的图快照、生成 contig 时
   从高丰度往低丰度补——天然适合"先出高置信 contig、再补低丰度"的宏基因组
   输出策略，pgr 的 `asm contig --min-coverage` 是单值，可借鉴快照思路。
4. **unitig 丰度 = 中位数向量**：与 pgr `asm unitig`（bcalm 移植）的
   `km:f:` 平均丰度不同，metaMDBG 保留每个 unitig 的丰度向量并取中位数，
   merge 时合并向量。若 pgr unitig 要输出稳健丰度（宏基因组场景），可参考
   中位数语义。
5. **k′-min-mer = minimizer 序列**：pgr 的 kmer 表是碱基 k-mer（u128 ≤ 64），
   metaMDBG 的节点是"minimizer 序列"（k′ 个 minimizer 的向量，hash128）。
   两者维度不同：minimizer 空间天然支持长读（HiFi/ONT），pgr 目前是短读工具，
   这一块暂不对齐，但知道差距在哪。
6. **内嵌 minimap2+spoa 抛光**：pgr 的 `asm map` 是完美匹配、无 gap；metaMDBG
   的 toBasespace 用 minimap2 容忍错误 + POA 抛光。若 pgr 未来要支持长读纠错
   或容错比对，metaMDBG 是"内嵌依赖"的参考，但 pgr 目前不引新依赖（用户约束）。
7. **断点续跑**：checkpoint 文件机制简单实用，pgr 的 `pl` 管道若做多步任务
   可参考（不过 pgr 目前坚持原语路线，优先级低）。
8. **multi-k 反馈 = 迭代式参数精化**（新增，2026-08-12）：把"组装"重构成
   **"参数化子命令 + 磁盘中间文件 + 循环调度"**——同一 `graph`/`contig`/
   `toMinspace` 子命令被 `AssemblyPipeline` 以不同 `k` 反复调用，跨进程只通过
   `parameters.gz`（gzip 二进制参数 blob，`Parameters::load/save`）和
   `unitig_data.txt`/`refined_abundances` 传递状态。pgr 的单进程 `libs/` 路线
   不必照搬子进程，但**"迭代长度参数 + 反馈 unitig"** 的骨架可直接映射到
   `asm` 的多趟 OLC/unitig 循环（跨接验证语义见本列表第 1 条）；`parameters.gz`
   可类比 pgr 用 struct 传参。
9. **外部分区计数（scale-out）**：k′-min-mer 计数不把全量 k-mer 塞内存，而是
   `hash128 % nbPartitions`（`nbBases/20Gb`，clamp `[nbCores, 5000]`）写分区文件
   → 分区内去重计数 → 合并（`KminmerCounter::partitionKminmers`，
   `CreateMdbg.hpp:3652`）。这是典型的"外排序式"大数据手法，pgr 若做超大
   数据集（如 `kmer count` 溢出内存）可参考分区+归并，而非一味加大内存。
10. **内存驱动的批量分片**：toBasespace 用 `--max-memory`（默认 8 GB，
   `_maxMemoryGB/8`，clamp `[1,100]`）决定 minimap2 一次读入多少 reads
   （`ToBasespace2.hpp:337`）——峰值内存预算显式控制批大小。pgr 若加长读抛光，
   可把内存预算作为一等参数。
11. **fragment = unitig 边界切分**（新增，2026-08-14）：把 contig 按"组成
   k-min-mer 所属 unitig"切回片段再逐片算覆盖——"先按图结构分片、再回放 reads"
   的模式与 pgr `asm map` 输出天然对齐（`to-rg` 后按 rg 边界即片段），是
   RepeatRemover 桥接证据的地基（见 §4.2.1）。
12. **read 桥接证据 = 图类型无关的 repeat breaking**（新增，2026-08-14）：
   "read 同时覆盖 ≥2 片段 → 连接证据；无证据且覆盖异常 → 断开" 可复用到任何
   "用 reads 验证连接"的场景，不限于 OLC（见 §4.2.1、§9.2）。
13. **minimizer chaining 比对的取舍**（新增，2026-08-14）：低密度 minimizer 锚点 +
   按 contig 分组 + chaining 过滤（`ReadVsContigMapper`），内存友好、容忍错误；
   与 pgr `asm map` 的全 k-mer 完美匹配是不同取舍——pgr 精度优先，metaMDBG
   内存/速度优先（见 §4.2.2）。

> 结论：metaMDBG 对 pgr 的价值分三层（2026-08-14 更新）：
> - **直接可移植（首要）**：§4.1.1 multi-k 跨接验证（doublet/triplet 用当前 k
>   的 solid k′-min-mer 选边 + `removeUnsupportedUnitigs` 嵌合清理）——连接
>   unitig 时 bubble 的选择骨架；§6.3 渐进丰度过滤（同轮内丰度过滤，含交替
>   收敛节奏、`2×cutoff` 保护、删后压实）、unitig 丰度中位数、多 cutoff 快照
>   倒序输出、RepeatRemover 桥接 reads 阈值（`2×source` 与 `nbBridgingReads != 0`）。
> - **架构参考**：multi-k 迭代调度、外部分区计数、内存驱动批分片、checkpoint 断点续跑。
> - **语义对照**：minimizer-space 节点、minimap2+POA 抛光，与 pgr 短读+完美匹配
>   路线距离较远，暂不借鉴但已知差距。

## 9. OLC v1 借鉴映射（2026-08-12）

承接 `design/asm-olc.md` 的 v1 待决项：

1. **渐进丰度过滤 → unitig 覆盖度驱动的布局前过滤**：
   `ProgressiveAbundanceFilter::removeAbundanceNoQueue`
   （`ProgressiveAbundanceFilter.hpp:2183`）：`t=1.1` 起步、`~10%` 步长、
   每轮删 `abundance < t` 的 unitig 并 recompact 邻接——不是单阈值一刀切。
   pgr `asm unitig` 头部已带 `cov=`，v1 可在 `asm olc` 布局前按 unitig
   丰度多轮剔除（或给 `asm unitig` 加渐进 `--min-coverage` 模式）。
2. **RepeatRemover 的桥接 reads 证据 → OLC repeat breaking 的实现路径**：
   `RepeatRemover.hpp:283` 把 reads 映射回 contig（`ReadVsContigMapper`，
   minimizer 索引）→ 按比对边界分片算覆盖度与 `_nbBridgingReads`
   （`:1195`）→ 无桥接的片段边界断开（判定在 `RepeatRemover.hpp:1254`，
   `nbContigsFinal > 1` 即把该 read 判为跨重复、`isCircular=0` 并拆分，`:1254-1257`）。
   这正是 `canu.md` §8.3 预言的"pgr
   `asm map` + `sam to-rg` + `rg coverage` 回放"的成熟实现——pgr 设施
   齐全，v1 可直接照搬语义
   （桥接 reads = 覆盖度证据，对应 Celera 的 6/15 阈值）。
   **具体阈值**（2026-08-14 补充，见 §4.2.1）：`minRepeatCoverage =
   sourceCoverage * 2.0` 判重复片段；`nbBridgingReads != 0` 才有连接证据；
   `isCircular=0` + 沿 `_finalContigIndex` 切分。
3. **cutoff 快照倒序输出 → 宏基因组 contig 分级输出**：多个 cutoff 的图
   快照、从高丰度往低丰度补——适合"先出高置信 contig、再补低丰度"策略。
4. **unitig 丰度中位数语义**：pgr `asm unitig` 的 `cov=` 是平均丰度；
   宏基因组场景若要稳健丰度，可参考"组成 k-min-mer 丰度向量取中位数、
   merge 时合并向量再取中位"（`Graph.hpp` `computeMedianAbundance`）。

## 10. 源码 quirks / 边界行为汇总（2026-08-13 校对补充）

这些是阅读源码时发现的、对理解行为或移植时有影响的细节，多数是隐式/硬编码：

1. **`getMultikStep` 恒为 1，后面是死代码**（`Commons.hpp:1984-1996`）：
   函数体第一行直接 `return 1;`，后续"k<20 → 1 / k<40 → 2 / else → 5"的分支
   永不执行。故 1.4 每轮 `k += 1` 是硬编码事实，不要被注释里残留的多步逻辑误导。
2. **`computeLastK` 用的是 N50 读长，不是平均读长**（`Commons.hpp:1726-1741`）：
   `lastK = n50ReadLength * assemblyDensity * 2.0`，再 `max(lastK, firstK+2)`。
   `--max-k` 直接覆盖该结果。注意 density 是组装密度 0.005（非纠错密度）。
3. **`_snpmerSize = 21` 是死配置**（`AssemblyPipeline.hpp:207`）：初始化后仅写入
   `parameters.gz`，从未被实际消费（snpmer 相关子命令调用全被注释，`:727` 等）。
4. **minimizer 长度硬 cap ≤ 16**（`AssemblyPipeline.hpp:202`）：`--kmer-size` 超
   过 16 会被静默截断，不报错。
5. **`--density-correction` 校验只在纠错开启时生效**（`:278-284`）：HiFi 默认
   关纠错，即使违反 `density-assembly > density-correction` 也不会报错。
6. **每条 read 首尾各丢一个 k-mer**（`Kmer.hpp:1395`，`_trimBps=1`）：影响最
   末端 minimizer 的采样，属边界截断设计而非 bug。
7. **`.checkpoint` 断点文件不校验内容**（`AssemblyPipeline.hpp` `isCheckpoint`）：
   只要文件存在即跳过该步，若上一步被用户中断但 checkpoint 已写入（或参数变更
   后重跑），会得到不完整的中间结果。重跑同目录需手动清 `tmp/checkpoints/`。
8. **`--min-abundance > 1` 时关闭 rescue**（`CreateMdbg.cpp:317`）：rescue 只在
   `_minAbundance <= 1` 时执行；且 rescue 内部 `median*0.1 > 1` 或 read 全为
   singleton（`allAbundanceOne`）时直接跳过（`CreateMdbg.hpp:4611-4613`）。
9. **unitig 长度是期望值**（`Graph.hpp:222-224`）：`(nbMinimizers-1) *
   _minimizerSpacingMean`，minimizer 空间没有真实碱基坐标，所有"长度"（含
   `length=` 头）都是 minimizer 数乘期望间距的估计。
10. **toBasespace 的 `--max-memory` floor 在 4GB**（`ToBasespace2.hpp:268`），
    批大小 `= _maxMemoryGB/8` clamp 到 `[1,100]`（`:337-349`）——内存预算不是
    上限而是"每批 reads 数"的缩放因子。
11. **`removeAbundanceNoQueue` 里 `_maxAbundance = currentCutoff*2`**（
    `ProgressiveAbundanceFilter.hpp:2203`）在每次抬升 cutoff 时重置，供下游
    "高于 2×cutoff 的受保护"判定使用（对应 repeat solver 保护）。
12. **大量 `exit(0)` 而非抛异常**：CLI 解析失败、参数校验失败均直接
    `exit(0)`（`:152,159,282`），退出码为 0 而非非零——脚本里判断"成功"时
    需以输出文件存在性为准，不能只看退出码。
13. **`computeMedianAbundance` 内部 `sort` 被注释**（`Graph.hpp:252-287`）：
    中位数正确性完全依赖 `_abundances` 向量"始终有序"不变量——建图时
    `createEdgeNode` 排序、图加载后 `AbundanceSortFunctor` 排序、`mergeWith`
    用 `std::merge` 保持有序。**Rust 移植时若直接用 `Vec<u32>` 需显式维护
    排序**，或改用无需排序的 quickselect（见 §6.2）。
14. **丰度统计两处语义不同**：unitig 丰度 = **中位数**（`Graph.hpp`
    `computeMedianAbundance`，稳健、抗低覆盖污染）；**片段覆盖度 = 均值**
    （`RepeatRemover.hpp:761-785`，对覆盖突变敏感）——metaMDBG 用均值做 repeat
    判定、用中位数做过滤阈值，移植时不要混用（见 §4.2.1）。
