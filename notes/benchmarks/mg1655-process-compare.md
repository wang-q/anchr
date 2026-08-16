# MG1655 新链处理与对比（2026-08-16）

> 目标：用当前改进后的 multik 链（主 K 31..128 + 气泡合并 + extend
> `--min-len 1000` + `--min-contig-len 200`）处理 MG1655，与旧链及
> megahit/spades 对比。5 组 reads：`6_down_sampling/MRX40P000/P001/P002` +
> `MRX80P000/P001`（同 `mg1655-unitig-bcalm-multik.md` 口径），参考
> `1_genome/genome.fa`（NC_000913，4,641,652 bp）。

## 结果（QUAST，contigs ≥ 0 bp）

| 组装 | N50 | # contigs | GF% | Dup | # mis | 备注 |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| **新 multik 链（本次）** | **123,988** | 88 | **98.364** | 1.002 | **0** | 同 5 组 reads |
| 旧 multik 链（2026-08-16 文档） | 95,478 | 107 | 97.61 | ~1.00 | 0 | K31..81，同 5 组 |
| 旧 merge_anchors（全组） | 95,484 | 103 | 97.669 | 1.000 | 0 | 旧运行 |
| 旧 fill_anchors（全组） | 105,614 | 88 | 97.673 | 1.000 | 0 | 旧运行 |
| mr_spades（旧，全量 merge reads） | 148,607 | 75 | 98.405 | 1.000 | 0 | 输入比 5 组多 |
| mr_megahit（旧，全量 merge reads） | 132,896 | 64 | 98.685 | 1.000 | **1** | 输入比 5 组多 |
| spades（旧） | 125,603 | 76 | 98.169 | 1.002 | 0 | 4_ 链输入 |
| megahit（旧） | 82,825 | 105 | 98.071 | 1.004 | **1** | 4_ 链输入 |

关键读数：
1. **新链相对旧 multik 链大幅提升**：N50 95.5K→**124.0K**（+30%）、GF
   97.61→**98.364%**（+0.75 pp）、0 mis 保持——改进来自 K128 主 K、气泡
   合并与 min-contig-len 200 保留短碎片。
2. **与外部组装器同口径（5 组）对比**：我们的 N50 123,988 已超过 megahit
   （82,825）与 spades（125,603 持平），且 0 mis vs megahit 1 mis；
   mr_spades 148.6K 仍更高（其输入为全量 reads，非同一口径）。
3. **mm/indel 极低**：mismatches 0.22/100 kbp、indels 0.09/100 kbp——
   MG1655 与参考同株，consensus 高度准确。

## extend 安全门槛（本次发现并修复）

extend 对**短碎片**的延伸会把重复元件拷贝接成嵌合体：238 bp 的 K128
碎片被 extend 长成 1,238 bp，两段 627 bp 分别对齐 1.2 Mb 处两个重复拷贝
（QUAST 报 3 mis：contig_83 inversion、contig_84/85 relocation）。

| 配置（K31..128 + bubble + min200） | N50 | GF% | # mis |
| :--- | ---: | ---: | ---: |
| extend 全部 contig（旧模板） | 125,790 | 98.575 | **3** |
| 不 extend | 124,007 | 98.166 | 0 |
| **extend 仅 ≥1000 bp（`--min-len 1000`，本次修复）** | **123,988** | **98.364** | **0** |

修复：`asm extend` 新增 `--min-len`（默认 1000），短于阈值的 contig 原样
输出；`{4,6}_unitigs.tera.sh` 三 unitigger 分支显式传 `--min-len 1000`。
G37 代价：GF 99.083→98.983%（短碎片延伸收益约 0.1 pp），0 mis 保持。

## olc 性能修复（支撑 MG1655 规模）

`consensus_with_ratio` 的近似去重/边界拼接原是 O(n²×L)（每对重建 k-mer
索引/窗口扫描）：G37 池约 1.5 min 可接受，MG1655 单组 9 主 K 池（约 6K
contigs）需数小时。改为：
* `dedup_contained_ratio`：全局 31-mer→contig 索引，候选按共享种子预筛；
* `merge_overlapping_contigs`：per-kept 种子集 + 头种子，O(1) 预筛拒绝
  大多数对（必要条件是 query 头/尾种子在 keeper 中，反之亦然）。
* 输出与旧实现在 G37 池逐位一致（41 contigs 完全相同），单测 31+14 全绿；
  G37 池 1m34s→21s，MG1655 9 主池 3m15s（旧版无法完成）。

## 数据与复现

* 运行目录：`/tmp/mg1655_new/`（主 K 产物、per-group unitigs/extend/anchor、
  `merge_b.fasta`、`quast_b/`）。
* 命令：与模板 6_unitigs multik 分支一致（K_LIST 验证、`--min-overlap 1000`
  `--min-contig-len 200`、extend `--min-len 1000`、anchor、最终合并
  `--min-contig-len 200`）。
