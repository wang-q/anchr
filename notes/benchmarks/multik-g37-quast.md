# G37 Quast 对比：multik vs 老流程（2026-08-14）

> 用 Quast（`~/.cbp/bin/quast.py`）对 multik 的 G37 输出做最终质量确认，
> 对照参考基因组 `~/data/anchr/ref/g37/genome.fa`（580,076 bp）。multik
> 输入：`4_down_sampling/Q25L60X40P000/pe.cor.fa`（155k reads，40×）；
> 老流程结果来自 `~/data/anchr/g37/9_quast/report.tsv`。

## 指标对比（contigs ≥ 500 bp）

| 指标 | **multik** | merge_anchors | mr_spades | mr_megahit |
|---|---:|---:|---:|---:|
| # contigs | 31 | 17 | 4 | 7 |
| Largest contig | 91,246 | 179,685 | 395,878 | 357,750 |
| N50 | 24,527 | 55,041 | 395,878 | 357,750 |
| Genome fraction (%) | 95.58 | 96.86 | 99.26 | 99.69 |
| **# mismatches / 100 kbp** | **31.4** | 70.8 | 249.1 | 293.4 |
| **# indels / 100 kbp** | **2.3** | 17.3 | 75.4 | 90.2 |
| # misassemblies | **8** | **0** | 1 | 0 |
| # N's / 100 kbp | **0** | — | — | — |
| Duplication ratio | 1.001 | 1.000 | 1.000 | 1.004 |

## 解读

**multik 的优势（序列质量）**：
* mismatches/indels 全流程最低（31/2.3 vs 老流程 70-293/15-90）——multik
  的 contigs 是 reads 的 solid k-mer 路径（跨接验证 + 渐进过滤），不经过
  外部组装器的错误累积；
* 0 N（无 N 目标达成）；Duplication 1.001（无冗余）。

**multik 的劣势（连接正确性 + 完整性）**：
* **# misassemblies = 8**（全是 relocation，`contigs_reports/...mis_contigs.info`：
  unitig_1 内部两个、unitig_2/3/7/8/12/15 各一个）——multik 在重复区/
  复杂区的连接把参考上不同位置的序列连在一起（跨接验证对"两个拷贝都有
  reads 证据"的重复区无法区分，见设计 §7 风险）；
* Genome fraction 95.58%（略低于老流程 96.9-99.7%）——395 条 contigs 里
  366 条 <500 bp（渐进过滤 dropped 的低丰度碎片），未参与 ≥500 bp 统计；
* N50 24.5K 低于老流程（55K-395K）——主路径碎片化（dropped 多）。

## 结论

真实数据暴露 multik 的两个待改进点（与设计 §7 风险一致）：
1. **重复区连接嵌合**（8 个 relocation）——需要 anchors 式的"高覆盖
   （重复区）排除"或 reads 桥接（v5 方向）；
2. **dropped 碎片多**（366 条 <500 bp）——渐进过滤/输出策略可调
   （如按 Quast 口径只输出 ≥500 bp，或收紧 dropped 条件）。

序列质量（mismatch/indel/0 N）已是全流程最优——"无 N 染色体"的序列
正确性目标达成，剩余是连接正确性与完整性的调优。

## 更新（2026-08-14，层次 3 reads 桥接落地）

实现 `bridge_filter`（unitig 间探针）+ `split_by_bridge`（unitig 内部
60-mer 窗口切分，设计 `asm-multik-misassembly.md` §7）后：

| 指标 | 修复前 | 层次 3 后 |
|---|---:|---:|
| # misassemblies | 8 | **0** |
| N50 | 24,527 | **26,562** |
| Genome fraction | 95.58% | 95.99% |

**8 个 relocation 全部消除**（错连 unitig 在合并前被切分），N50 反升、
Lambda/20k 环状不回归。剩余 366 条 <500 bp 碎片（dropped）是输出策略
问题，与错连无关。
