# `fq s-filter` vs quorum 丢弃判定对照（G37 Q25L60 全量，2026-08-15）

> P1-1a（`design/fq-validation.md` #16）：验证 `fq s-filter` 是否复现老流程
> quorum 步骤的丢弃判定。老流程语义：quorum_error_correct_reads 输出头部
> 带 `:sub:`/`trunc` 标记的 reads 全部丢弃（`R.discard.lst`），pe.cor =
> 未修正原始序列（见 `references/anchr-legacy-pipeline.md` §2.1）。

## 方法

* 输入：G37 `2_illumina/trim/Q25L60/` 全量 reads（R1 320,845 + R2 322,417 +
  Rs 17,220），用与老流程相同的 `faops interleave -q -p pe/se` 重建
  `pe.renamed.fastq`（643,262 PE）+ `se.renamed.fastq`（34,335 SE）——
  共 675,254 条记录；
* quorum 侧：直接使用 G37 老产物 `R.discard.lst`（36,390 行，unique
  35,711；本机 quorum 二进制与新系统不兼容、段错误，无法重跑）；
* s-filter 侧：参数对齐 quorum（`-k 24 --good 3 --anchor-count 4
  --min-count 3`，bits 7、质量阈值 detect+5，对应 quorum `-m 24 -m 3 -g 3
  -a 4 -b 7 -q 38`），`target/release/anchr fq s-filter` 全量 5.8 s /
  峰值 5.0 GB；`--skip 1` 变体（对齐 quorum `-s 1`）结果相同。

## 结果

| 集合 | 条数 | 占总 reads 比例 |
|---|---:|---:|
| 总 reads | 675,254 | — |
| quorum 丢弃（老产物） | 35,711 | 5.29% |
| `fq s-filter` 丢弃 | 52,823 | 7.82% |
| **交集（两边都丢）** | **35,536** | — |
| 仅 quorum 丢弃 | 175 | 0.026% |
| 仅 s-filter 丢弃 | 17,287 | 2.56% |

**quorum 丢弃集合的 99.5%（35,536/35,711）被 s-filter 复现**；s-filter
额外多丢 2.56% 的 reads。

## 差异机制（对照 `quorum-1.1.2/src/error_correct_reads.cc`）

`fq s-filter` 是 `pgr::libs::kmer::qcheck` 的封装，逐碱基判定、**遇到第一个
错误事件立即丢弃**。quorum 的判定更宽容，四处实现差异均使 s-filter 更严格：

1. **窗口容忍未实现**：quorum `err_log::check_nb_error` 用 `-w 10 -e 1`
   （窗口内错误数 ≥2 才触发丢弃）；s-filter 无窗口概念，首错即丢；
2. **多候选替代不标记**：quorum 在"多个候选且无唯一解"时原样输出、不记录
   sub/trunc（`error_correct_reads.cc:500-552`）；s-filter 直接判
   `Substitution` 丢弃；
3. **高质量碱基豁免未实现**：quorum 多候选时 `*qual >= qual_cutoff`（-q 73
   ≈ Q40）直接保留；s-filter 的 `CheckParams` 无该参数；
4. **无 `--no-discard` 路径差异**：老流程本来就把带标记 reads 全丢，等价于
   s-filter 的丢弃语义（此处不构成差异）。

仅 quorum 丢弃的 175 条（0.026%）：**根因已定位（2026-08-15 补充）**——
`get_best_alternatives` 的计数表清理差异：quorum（`mer_database.hpp:303`）
在更高 quality level 出现时**只清空当前位置之前的 counts**（`for j<i`），
当前位置之后的旧 level 计数残留；pgr `qcheck::best_alternatives` 清空
全部。残留计数让 quorum 在部分 reads 上读到低计数 → 误判 substitution/
truncation，而 s-filter 判定保留。补充验证：175 条中 167 条提取（8 条在
SE），`asm map` perfect 回贴参考仅 93 条命中、每条 1 个位置（**非参考
重复区多匹配**）；样本 reads 彼此高度相似（pe650/1 与 pe651/1 移位重复），
属 reads 级冗余 + quorum 判定噪声的组合，非系统性差异。`--skip` 变体
不能复现。影响 0.026%，不修（对齐 quorum 需复刻其部分清理行为，代价 >
收益）。

## 结论

* **语义对齐成立**：s-filter 复现 quorum 丢弃判定的 99.5%，差异集中在
  quorum 的纠错器容错细节（窗口/质量豁免/多候选），属"筛选器 vs 纠错器"
  的定义差异，不是实现缺陷；
* **参数定稿**：`fq s-filter -k 24 --good 3 --anchor-count 4 --min-count 3`
  对齐 G37 quorum 参数；若需与 quorum 行为完全一致（消除多丢的 2.56%），
  需在 qcheck 补窗口容忍与质量豁免（暂缓，等真实宏基因组数据评估影响）；
* **后续**：`fq-validation.md` #16 状态 ⚠️ → ✅（对照完成，差异已记录）。
