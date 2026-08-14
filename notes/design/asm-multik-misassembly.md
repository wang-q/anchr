# multik 防 misassembly 方案：覆盖度区间过滤（anchors 核心，精简实现）

> 2026-08-14。背景：G37 Quast 质检（`benchmarks/multik-g37-quast.md`）
> 显示 multik 有 8 个 misassembly（全 relocation），unitig_1 把参考相隔
> ~48 万 bp 的两段（4286-60431 与 540652-571614）连在一起。参考这两段
> **无共享 k-mer**（非简单重复区）。老流程 anchors 的目的就是防
> misassembly，但其实现（bbwrap→basecov→spanr→hnsm→contained→orient→
> merge）极冗余。本文规划精简方案：**只取 anchors 的核心思想（reads 回贴
> + 覆盖度区间过滤），复用 multik/pgr 已有机制**。

## 1. anchors 防 misassembly 的核心机制（读 `templates/anchors.tera.sh` 确认）

1. reads 映射回 unitigs（bbwrap perfectmode）→ 每碱基覆盖度；
2. 覆盖度区间 `[lower, upper]`：
   `lower = max(mincov, (median − mscale×MAD)/lscale)`、
   `upper = (median + mscale×MAD)×uscale`；
3. **只保留覆盖度在区间内的区域**（低覆盖 = 错误区、高覆盖 = 重复区，
   都不可靠 → 排除）→ "properly covered regions" 才是 anchors。

**防 misassembly 的原理**：错误连接通常发生在覆盖度异常的区域（重复区
reads 过多、错误区 reads 过少）；把这些区域排除（宁可断开、不错连），
misassembly 自然减少。

## 2. multik 的现状与差距

* multik 已有：unitig `cov=`（组成 k-mer 平均覆盖度）、渐进丰度过滤
  （只删低覆盖分支）、跨接验证（k-mer 计数 ≥2 选边）；
* **缺**：`upper`（高覆盖/重复区排除）——渐进丰度过滤只处理低覆盖，
  没有"高覆盖=重复区→断开"的机制。这是 8 个 relocation 的根源（错连
  发生在覆盖异常区）。

## 3. 方案（三个层次，从简到繁）

### 层次 1：覆盖度一致性选边（迭代内，零额外工具）

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

### 层次 2：reads 回贴 + 局部覆盖剖面（复用 pgr map 工具链）

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

### 层次 3：reads 桥接选边（metaMDBG RepeatRemover 语义，最精确）

* 对断开点/分支点，用 reads 同时覆盖两侧的证据（`nbBridgingReads != 0`）
  选择正确连接（metaMDBG §4.2.1）；
* 之前已分析（`asm-multik.md` §4.9 附近）：纯重复区（两拷贝都有证据）
  需要 reads 桥接或覆盖度证据收尾；
* 复杂度最高，作为层次 2 不足时的最终手段。

## 4. 验证计划

1. **层次 1**：G37 上实现 cov 一致性断边，Quast 复测 misassemblies
   （预期部分减少，但 unitig 级平均覆盖可能漏）；
2. **层次 2**：复用 map/to-rg/coverage 做局部剖面，Quast 复测
   （预期 8 个 relocation → 显著减少，代价是 Genome fraction 略降——
   断开的异常区不输出）；
3. 同时验证：Lambda（46,467 零缺口）和合成长读（100% 单条）不回归
   （过滤只在覆盖异常区生效，正常区域不切）。

## 5. 与老流程的对比（为什么精简）

| 老流程 anchors | 本方案层次 2 |
|---|---|
| bbwrap（外部工具） | `pgr asm map`（已有） |
| basecov 解析 + perl 过滤 | `pgr sam to-rg` + `rg coverage`（已有） |
| spanr cover/stat/some/span/compare（外部） | 无（直接按 rg 覆盖阈值切分） |
| hnsm range 提取 | 无（在 unitig 上直接切分） |
| contained→orient→merge（4 步） | 无（multik 已有压实） |

本方案只保留 anchors 的核心判据（覆盖度区间），把"切分/去冗余/合并"
交给 multik 已有机制，去掉全部外部工具链与后处理冗余。

## 6. 待决

* 层次 1 vs 层次 2 先做哪个（若 unitig 级 cov 已能抓住错连，层次 1 足够；
  证据表明平均 cov 掩盖局部异常，倾向直接做层次 2）；
* 断开策略：异常区切分后是丢弃（像 anchors 不选）还是保留为独立 unitig
  （像 multik 的 dropped）；
* 阈值参数是否沿用 anchors 默认（mscale=3/lscale=3/uscale=2/mincov=5）。

## 7. 实现结果（2026-08-14，层次 3 落地）

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
