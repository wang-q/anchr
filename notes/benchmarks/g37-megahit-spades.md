# G37 组装策略对比：anchr vs MEGAHIT/SPAdes（2026-08-16）

> 数据来源：`results/model_org.md`（2026-08-16 全流程）、`/tmp/g37_run`
> （2026-08-16 实际运行产物）、`megahit-1.2.9/src/` 源码精读。本笔记回答
> 三个问题：为什么 megahit/spades 的 N50 更高？为什么它们引入错误组装？
> 我们该改进哪些步骤（附已验证的实验证据）。

## 1. 结论速览

* **我们的 GF 缺口 ≈ 低复杂度/重复区**：G37 参考 580,076 bp 中 `pgr rept`
  标注 53 个重复/低复杂度区（16,875 bp，占 2.9%）。QUAST 对齐显示我们的
  merge 链未覆盖 17,419 bp，其中 **13,299 bp 落在这些重复区**（76%），仅
  4,120 bp 是普通覆盖缺口。megahit/spades 的 N50 优势主要来自"敢穿过这些
  区"，代价就是错误组装全部发生在这些区的边界。
* **前提修正**：在原始 trim reads 上，megahit 的 N50（39.7K）其实**低于**
  我们（55.1K）；真正拉开差距的是 spades（163.8K）与 merged-read 链
  （mr_megahit 319.2K、mr_spades 580.5K，megahit/spades 跑 merge+纠错后的
  长 reads）。
* **已验证的最大单项改进 = 提高 MR 链的 k 上限**：模板把主 K 硬编码为
  31..81，而 MR reads N50≈363 bp（merge 后），megahit mr 用 k 45..225。
  把 MR 链主 K 扩到 31..121（7 组全链实验）：

| 指标 | 现链 merge_mr_multik | 高 k 链（31..121） |
| --- | ---: | ---: |
| # contigs | 17 | **13** |
| N50 | 55,049 | **83,618**（+52%） |
| GF% | 96.963 | **98.633**（+1.67 pp） |
| # misassemblies | 0 | **0** |
| mm/100 kbp | 77.46 | 219.47* |

  *mm/100 kbp 上升是**reads 与参考在低复杂度区的真实差异**（见 §7.1），
  不是组装错误：新覆盖区的 consensus 与 reads 一致率 99.963%（37/100 kbp）。

## 2. 三方 G37 数据（QUAST，results/model_org.md）

| Assembly | # contigs | Largest | N50 | GF% | Dup | # mis | mm/100 kbp | indel/100 kbp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| merge_multik（我们） | 15 | 179,712 | 55,098 | 96.997 | 1.000 | 0 | 76.44 | 19.20 |
| merge_mr_multik（我们 MR） | 17 | 114,442 | 55,049 | 96.963 | 1.003 | 0 | 77.46 | 19.14 |
| merge_anchors（我们） | 15 | 179,712 | 55,098 | 97.004 | 1.000 | 0 | 76.78 | 19.55 |
| spades | 35 | 236,302 | 163,847 | 99.056 | 1.002 | **1** | 222.74 | 66.49 |
| mr_spades | 2 | 580,506 | 580,506 | 100.000 | 1.001 | 0 | 300.77 | 91.64 |
| megahit | 56 | 89,673 | 39,730 | 97.507 | 1.004 | **4** | 85.23 | 18.31 |
| mr_megahit | 44 | 319,186 | 319,186 | 99.893 | 1.025 | **3** | 311.14 | 88.39 |

关键读数：
1. **错误组装的坐标全部验证在重复区边界**（QUAST all_alignments vs
   `1_genome/repetitive/repetitive.json`）：megahit k81_70 在 174,996（重复区
   174,996-175,164 起点）、k81_38 在 225,052/226,031 附近、k81_26 在
   224,299/230,404（两端都贴着重复区）、k225_12 在 168,544/428,366（两个
   重复区边缘）；spades NODE_2 断点在 87,631（85,550-87,693 簇尾）与
   314,937（314,726-314,928 区内）。K31/K51 单主 unitigs 的 1 个 mis 也在
   224,299/230,404——与 megahit k81_26 **同一位点**。
2. **我们 0 mis 的代价**：anchor 上界过滤 + multik branch/bridge 防御 +
   `--min-overlap 1000` 精确 overlap，使链在重复区边界"宁断勿错"。
3. **序列质量**：我们的 mm/indel 全链最低；mr_* 链 mm 300/100 kbp 主要是
   覆盖了低复杂度区（那里 reads 与 NC_000908 真实差异高达数千/100 kbp，
   bwa mpileup 实测重复区 5,006/100 kbp vs 非重复区 96/100 kbp）。

## 3. 我们的瓶颈（为什么 N50 只有 55K）

按对最终结果的贡献排序：

1. **主 K 上限 81 太低**（模板 `4_unitigs`/`6_unitigs` 硬编码 `31 41 51 61
   71 81`）。MR reads merge 后 N50≈363 bp，megahit mr 用 k 45..225；我们的
   `Kmer::MAX_K=128` 且 `auto_ks` 对 363 bp reads 会给 51/71/91/111。
   实测单组 MRX40P000 各主 K unitig N50：K31 12.3K → K81 20.4K → K101
   35.1K → K121/128 **44.9K**。
2. **multik 是"单主骨架 + 从 k 验证"，从 k 不重组装**（`asm-multik.md`
   §4）：later ks 只验证/压实 k0 图，所以主 K 直接决定能解析的重复规模。
   megahit 每轮 `seq2sdbg` **用上一轮 contigs/bubbles 重建新 k 图** +
   `iterate` 从 reads 收集跨 contig 端点的 (k+step+1)-mer 迭代边，因此
   高 k 轮是真正的新图（megahit.md §4.3/§4.6）。
3. **anchor 上界过滤把重复区整段排除**（`anchor.rs`：`upper=(median+3·MAD)
   ·uscale`，G37 Q25L60X40 upper≈110）：低复杂度区多定位 reads 堆高覆盖 →
   过滤 → 链在这些区断开（17,419 bp 未覆盖中的 13,299 bp）。
4. **OLC 合并只认 ≥1000 bp 精确 overlap**（`--min-overlap 1000`，
   `olc/overlap.rs`）：跨组/跨主 K unitigs 若都在同一重复边界附近断开
   （断点差 <1000 bp）就无 overlap 可并。
5. **没有 megahit 式气泡合并**：变异/错误造成的小气泡我们靠覆盖度剪枝
   （`progressive_filter` 丢低丰度分支），主路径不"吞并"备选路径；megahit
   的 naive/complex bubble 把 ≥95% 相似、≤20k 的备选路径合并进主路径。
6. **没有本地组装**：megahit `local_assemble` 用双端 reads 延伸 contig 端点
   （局部区 ≤650 bp，IDBA 内核），补小缺口；我们有 `fq extend` 但模板不
   用于 unitig 端点。

## 4. megahit 为什么更长（源码机制）

精读 `megahit-1.2.9/src/`，其 N50 优势来自四件事：

1. **多 k 重建 + 引导 + 迭代边**（`main_sdbg_build.cpp` seq2sdbg、
   `main_iterate.cpp`）：每轮 k 用上一轮 contigs/bubbles/local contigs 做
   种子建新图；`ContigFlankIndex` 存 contig 端点 k-mer 及后续碱基，`iterate`
   把 reads 里跨端点、连续 `step+1` 个已匹配位置的 (k+step+1)-mer 收进
   新图（`contig_flank_index.h` `FindNextKmersFromReads`）。高 k 的新图真正
   越过了低 k 的断点。
2. **清洗轮**（`main_assemble.cpp`，默认 5 轮）：tip 移除（倍增阈值、下一
   节点深度 >8× 才删，`tip_remover.cpp`）→ naive bubble → complex bubble
   （分叉汇合两端唯一、长度差 ≤ (1-sim)、编辑距离相似度 ≥0.95 才合并，
   `bubble_remover.cpp` `GetSimilarity` 是 banded DP）→ weak link 断开
   （分支深度 ≤0.1× 总深度才断，`weak_link_remover.cpp`）→ 低深度修剪
   （局部窗口深度均值，`min_depth×1.1` 迭代，`low_depth_remover.cpp`）。
3. **本地组装**（`localasm/local_assemble.cpp`）：reads 回帖 contig 端点
   （`hash_mapper`：seed 31、sparsity 8、唯一回帖、相似度 ≥0.8、
   `min_mapping_len=75`），端点局部区间（`min(2·mean, mean+3·sd)` 封顶
   650）内收集 reads，IDBA 内核逐 k 组装延伸。
4. **长 reads + 大 k**：mr 链 merge reads ~450 bp，k 到 225，短重复在
   k-mer 层面直接变唯一。

## 5. megahit/spades 为什么错（机制 + 证据）

* **复杂气泡合并过宽**：95% 相似度 + 20k 长度上限，对"两个拷贝相似但不同
  基因组位置"的重复副本也会合并 → relocation。G37 的 mis 全是 relocation。
* **weak link 断不开近等深重复**：断链条件是分支深度 ≤10% 总深度；单倍体
  基因组的两个重复副本深度约 50/50，条件不触发 → 图保持错连。低复杂度区
  多定位 reads 深度反而更高，更不会被断。
* **多 k 重建传播低 k 错误**：低 k 图里的错连若在更高 k 仍被 reads 支持
  （重复区 reads 本身就多定位），高 k 重建也消不掉。
* 证据：全部 4+3+1 个 mis 的断点坐标都落在 `pgr rept` 重复区边界（§2），
  且我们 K31/K51 单主 unitigs 在同一位点也有 1 个 mis——说明这不是
  megahit 独有，而是该区本身的图结构难题，我们的下游防御把它拦住了。

## 6. 我们的 0 mis 从哪来（防御清单，勿动）

* `bridge_filter`：每个连接必须被 reads 完整覆盖 60-mer 探针
  （`multik.rs` `probe_kmer`）——嵌合连接无 bridging reads → 丢。
* `split_by_bridge`：unitig 内部无 reads 支持的 100-mer 窗口 → 切断。
* branch 标记：≥4 伙伴（重复片段扇形展开）的节点不参与 recompact。
* anchor 覆盖窗口：上界排重复区、下界排低覆盖/错配区。
* `asm olc`：精确 seed-and-verify + mutual-best 贪心布局 +
  `--min-overlap 1000`（legacy 对齐 daligner ≥1000 bp / ≥99.9% 口径）。
  实测：去掉 anchor 直接合并 7 组高 k unitigs，N50 83,639 / GF 98.679 但
  **1 mis**；加 anchor 后 83,618 / 98.633 / **0 mis**——anchor 以 ~0.05 pp
  GF 换掉 1 个 mis，方向正确。

## 7. 已验证的实验：高 k MR 链（改进方向 P0）

在 `/tmp/g37_mr128` 用与模板相同的 7 个 MR 组（MRX40P000-004 +
MRX80P000-001），唯一改动：主 K 从 31..81 扩到 31..121（8 个主 K，
31/41/51/61/71/81/101/121，-p 2/主，8 并发），下游 olc/anchor 参数不变。

| 阶段 | 现链（31..81） | 高 k 链（31..121） |
| --- | ---: | ---: |
| 单组 unitigs（olc 跨主 K 后）N50 | ~31.3K | 37.9-83.6K（X40 组） |
| 单组锚后 N50 | ~31.2K | 37.8-83.6K |
| 7 组合并 N50 | 55,049 | **83,618** |
| 合并 # contigs | 17 | **13** |
| 合并 GF% | 96.963 | **98.633** |
| 合并 # misassemblies | 0 | **0** |
| Dup | 1.003 | 1.005 |

剩余未覆盖 7,928 bp（7 个 run：85612-85702、167294-169516、226286-227199、
229535-230731、350368-351594、389356-390534、428327-429423）reads 覆盖
93-100%（mean 7.5-41.3×）——reads 在，但需要重复区解析/本地组装机制才能
接上（P2）。

### 7.1 mm/100 kbp 上升的归因

高 k 链 mm 219/100 kbp（现链 77）看似劣化，但：
* 新覆盖 ~9.5K bp 集中在低复杂度区，这些区 reads 与 NC_000908 差异
  **5,006/100 kbp**（bwa mpileup 实测），consensus 忠实反映 reads 差异；
* 新 contigs 与 reads 回帖一致率 99.963%（bwa mem，37.2/100 kbp）；
* 即 mm 上升 = QUAST 对比参考的口径问题，不是 consensus 变差。报告时可
  加一列"reads-vs-contigs 一致率"，或按覆盖 bp 归一化 mm。

## 8. 改进建议（按优先级）

### P0 模板 k 上限（✅ 2026-08-16 已实现并验证）
* `6_unitigs`（MR 链）：主 K `31 41 51 61 71 81 101 121`（或按组 reads
  N50 用 `--kmer auto`，MR reads 会得 51/71/91/111）。预期：G37 N50
  55K→84K、GF +1.7 pp、0 mis 保持。
* `4_unitigs`（150 bp reads）：主 K 上限提到 91（K101 单组 N50 回落
  19.9K vs K91 30.7K，k 逼近读长后反而碎），即
  `31 41 51 61 71 81 91`。
* 实现：`templates/{4,6}_unitigs.tera.sh` 三个 unitigger 分支统一用
  `KS` 变量（6_ 含 101/121、4_ 含 91）；MR 链 7 组全链实测 N50 55.0K→
  **83.6K**、GF 96.96→**98.63%**、contigs 17→13、**0 mis**（§7）。
* 注意：高主 K 的 K31/K51 单主自身有 1 个 mis（224,299/230,404），保持
  现有 bridge/anchor/olc 防御即可，最终链 0 mis（实验已验证）。

### P1 multik 增加 megahit 式气泡合并（✅ 2026-08-16 已实现并验证）
在 multik 图内实现 naive/complex bubble：分叉-汇合两端唯一 + 长度差 ≤
(1-sim) + banded 编辑距离 ≥0.95 + 长度 ≤20k → 主路径吸收高深度支路，
低深度支路 carry 输出（不删）。这样变异/错误气泡不打断主路径，unitig 在
变体位点连续。成功判据：小变异数据集上 unitig N50 提升且 0 mis；复用
`bridge_filter` 保证合并后连接仍有 reads 支持。
* 实现：`multik.rs` `bubble_merge`（megahit `SearchAndPopBubble` 语义：
  分叉汇合两端唯一、`in_deg` 对齐、banded DP 相似度、`merge_len*k`
  长度上限；branch 标记节点不参与；被吸收的备选路径作为独立 unitigs
  输出）+ CLI `--merge-similar`（默认 0.95）/ `--merge-len`（默认 20）。
* G37 MR 链实测（主 K 31..121 + 气泡合并）：contigs 13→**12**、Dup
  1.005→**1.001**、mm 219→**211**、N50 83.6K 持平、GF 98.63→98.59
  （-0.05 pp，低覆盖备选路径被 anchor 滤掉，多为 ≥95% 冗余序列）；
  **0 mis 保持**。单元测试 5 个（相似度/主导路径/异源拒绝/超长拒绝/
  异汇合拒绝）+ CLI 冒烟测试。

### P1 本地组装 / 缺口闭合（✅ 2026-08-16 已实现并验证）
对最终 unitigs 端点做 megahit `local_assemble` 式延伸：`asm map` 回帖 →
端点局部区（insert 估计）收集 reads → `asm unitig` 组装延伸。目标：剩余
4.1K 普通缺口 + 部分重复区缺口。legacy 的 fill 步骤是同类思路的旧实现
（`anchr-legacy-pipeline.md`），可对照。
* 实现：新命令 `anchr asm extend CONTIGS READS...`（`libs/asm/extend.rs` +
  `cmd/asm/extend.rs`）——对每个 contig 端点做 base-by-base 的 reads
  k-mer 图游走（`RefineTable::fill_right_counts/fill_left_counts`），只有
  当前碱基有严格多数支持（`>= --min-support` 且 `>= 2x` 次高）才延伸，
  交汇/重复上下文自动停（不会跨位点错连）；参数 `--kmer 31`、
  `--max-extend 500`、`--min-support 2`、`--min-extend 0`。模板
  `{4,6}_unitigs` 三 unitigger 分支在 olc 跨主 K 合并后自动加一步
  `asm extend`。
* G37 MR 链实测（主 K 31..121 + 气泡合并 + extend）：GF 98.585→
  **98.809%**、N50 83.6K→**83.8K**、未覆盖 8,210→**6,907 bp**、
  **0 mis 保持**（延伸进低复杂度区的碱基均被 reads 支撑）。

### P2 重复区解析（高风险，大工作量，GF 100% 的关键）
剩余 7,928 bp 全部是低复杂度区且 reads 覆盖充足，需要：
1. 高 k 单主（121/128）unitigs 若已穿过该区，直接取用（本实验已自动完成
   一部分：17,419→7,928 bp）；
2. 对仍未接上的重复边界用双端/merge reads 的连接证据裁决：两侧唯一锚 +
   insert 一致 + 无竞争连接才接，否则保持断开。这是 spades 重复解析的
   核心逻辑，也是"0 mis"红线的边界，必须做歧义检测。
（待做，见 `todo.md`。）

### P2 指标口径
QUAST mm/indel 按覆盖 bp 归一化 + 补 reads-vs-consensus 一致率列，避免
把真实菌株差异当组装错误（§7.1）。
（待做，见 `todo.md`。）

## 9. 关联

* `notes/references/megahit.md`（源码全景分析，本笔记 §4/§5 是其机制聚焦）；
* `notes/design/asm-multik.md`（multik 设计，主 K/从 K 机制与 §3.2 对应）；
* `notes/benchmarks/multik-g37-quast.md`（2026-08-14 multik 初版对比）；
* `notes/todo.md`（P0/P1/P2 已同步）。

## 10. 2026-08-16 实现总账（G37 MR 链，7 组）

| 阶段 | N50 | # contigs | GF% | Dup | 未覆盖 bp | # mis |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 原模板（主 K 31..81） | 55,049 | 17 | 96.963 | 1.003 | 17,419 | 0 |
| + 主 K 31..121（P0） | 83,618 | 13 | 98.633 | 1.005 | 7,928 | 0 |
| + 气泡合并（P1） | 83,618 | 12 | 98.585 | 1.001 | 8,210 | 0 |
| + 端点延伸（P1） | 83,756 | 12 | 98.809 | 1.003 | 6,907 | 0 |

改动文件：`templates/{4,6}_unitigs.tera.sh`（KS 上限 + extend 步）、
`src/libs/asm/multik.rs`（bubble_merge）、`src/libs/asm/extend.rs`（新）、
`src/cmd/asm/{multik,extend}.rs`、`src/cmd/asm/mod.rs`、
`tests/cli_asm_multik.rs`、`tests/cli_asm_extend.rs`（新）。
`cargo test -- --test-threads=1` 408 全绿；fmt/clippy `-D warnings` 干净。
