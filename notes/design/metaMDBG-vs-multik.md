# metaMDBG 与 anchr `asm multik` 实现对比（2026-08-14）

> 基于完整源码阅读：metaMDBG 1.4（`metaMDBG-metaMDBG-1.4/`，C++）与
> anchr multik（`src/libs/asm/multik.rs`，Rust）。两者是同一条
> "multi-k 迭代 + 图验证"路线的两种实现；本文对比实现层差异。

## 1. 定位

| 维度 | metaMDBG 1.4 | anchr `asm multik` |
|---|---|---|
| 目标 | 长读（HiFi/ONT）宏基因组组装 | 短读/长读通用，无 N 染色体 |
| 语言/并行 | C++20 + OpenMP，多进程调度 | Rust + rayon，单进程内存组合 |
| 建图空间 | minimizer 空间（density 0.005 采样） | 碱基 k-mer（FastK 字节键） |
| 命令形态 | 单二进制多子命令 + checkpoint 断点续跑 | 单命令，进程内完成 |

## 2. multi-k 迭代

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

## 3. unitig 反馈与跨接验证

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

## 4. 渐进丰度过滤

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

## 5. 嵌合清理与防 misassembly

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

## 6. 性能与规模

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

## 7. 借鉴与差异总结

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

## 8. 结论

两者是同一核心思想（multi-k 迭代 + 图验证）的两种实现：metaMDBG 面向
长读宏基因组（minimizer 空间、150 轮、maxAbundance 假设、容错映射），
multik 面向通用/无 N（碱基空间、大步长、单菌株保护、完美探针）。multik
的两个独有设计（直链保护、split_by_bridge 内部切分）解决的是碱基空间/
单菌株特有的问题（覆盖波动误删、recompact 固化错连），是 metaMDBG 不
需要的。跨接验证与丰度过滤是共同骨架，实现细节因空间（minimizer vs
碱基）而异。
