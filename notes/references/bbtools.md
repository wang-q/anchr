# BBTools: reads 处理与组装的主参考工具包

> 整理于 2026-08，源自对仓库根 `BBTools-40.01/`（2026-02-11）源码的分析。
> BBTools 是 fq/asm 迁移的主参考：`fq` 家族以 BBTools 39.38 逐字节核对
> （golden 见 [audit-fq.md](../audit/audit-fq.md)），`asm` 家族对照 BBTools
> tadpole（contig）与 bbmap `perfectmode`（map，见
> [audit-asm.md](../audit/audit-asm.md)）。本文汇总与 anchr 相关的模块与
> 语义要点，作为已随迁设计文档的索引与补全。

## 1. 简介

`BBTools`（Brian Bushnell, JGI/DOE）是 Java 实现的快速序列处理工具包，
以单机大内存、多线程为设计目标。本仓库内的 `BBTools-40.01/` 是其完整源码
发行（脚本壳 + Java 源码 + 资源），目录自带 `.gitignore`（`*`），仅作参考，
不被 git 跟踪。

- **版本基线**：本地实际安装 39.38（`~/.cbp/libexec/bbtools/`，`bbduk.sh`
  2024-11-18，入口 `jgi.BBDuk`）；仓库根 40.01 是更新版（入口
  `bbduk.BBDukS` 家族）。运行对照以 39.38 为准，40.01 只作行为参考。
- **运行方式**：`*.sh` 脚本壳调用 `java -cp BBTools.jar <入口类>`；
  `current/` 下为按包组织的 Java 源码，`jni/` 下为原生扩展（两版本均禁用）。

## 2. 目录结构

| 路径 | 内容 |
| :--- | :--- |
| `*.sh`（顶层 ~320 个） | 命令入口脚本，每个工具一个 `.sh` |
| `current/` | Java 源码，按包组织（`clump/`、`jgi/`、`assemble/`、`map/`、`bbduk/` 等） |
| `resources/` | 参考数据（如 `lambda.fa.gz`，OLC 冒烟用） |
| `jni/` | JNI 原生扩展（`BBMergeOverlapper.c`、`BandedAlignerJNI.c` 等，已禁用） |

## 3. 与 anchr 相关的模块映射

BBTools 入口 → anchr/pgr 目标（映射出处 [anchr-trim-replace.md](../design/anchr-trim-replace.md) §4.2、
[anchr-merge-replace.md](../design/anchr-merge-replace.md)、
[fq-assemble.md](../design/fq-assemble.md)）：

| BBTools 工具 | 入口类 | 源码位置 | anchr 目标 |
| :--- | :--- | :--- | :--- |
| clumpify | `clump.Clumpify` | `current/clump/Clumpify.java` | `fq clump` |
| bbnorm | `jgi.KmerNormalize` | `current/jgi/KmerNormalize.java` | `fq norm` |
| reformat 降采样 | `jgi.ReformatReads` | `current/jgi/ReformatReads.java` | `fq sample` |
| bbduk | `jgi.BBDuk`（39.38）/ `bbduk.BBDukS`（40.01） | `current/jgi/BBDuk.java`、`current/bbduk/` | `fq trim-adapter` / `fq trim` |
| bbmerge / tbo/tpe | `jgi.BBMergeOverlapper.mateByOverlapRatio` | `current/jgi/BBMergeOverlapper.java` | `fq merge` / `fq overlap` |
| repair | `jgi.SplitPairsAndSingles rp` | `current/jgi/SplitPairsAndSingles.java` | `fq split` |
| kmercountexact | `jgi.KmerCountExact` | `current/jgi/KmerCountExact.java` | k-mer 基础（留 pgr） |
| tadpole | `assemble.Tadpole` | `current/assemble/Tadpole.java` | `asm contig`/`asm unitig`、`fq extend`/`fq ec-kmer` |
| bbmap | `map.AbstractMapThread` | `current/map/AbstractMapThread.java` | `asm map`（`perfectmode`） |

## 4. 已确认的关键语义要点

以下结论来自已随迁的 audit/design 文档，实现与 golden 核对时直接引用：

1. **质量值存储**：`Read.quality` 存 phred 值（输入时减 ASCII_OFFSET），输出再加回
   （[anchr-trim-replace.md](../design/anchr-trim-replace.md) §4.1）。
2. **clumpify 桶序**：外部 bucket 路径输出为"按桶拼接"序，与 BBTools 大数据行为
   一致，但与内存路径不同（[audit-fq.md](../audit/audit-fq.md)）；`--dedupe
   --dupesubs 0` 整对去重，N 通配精确匹配，保留期望错误更少的一对。
3. **tadpole 内存模型**：多字 `Kmer`（`Vec<u64>` 共 2k 位，镜像 BBTools long
   array），k 无上限；逐字节一致需复刻其 `-Xmx` 相关内存语义，已文档化偏差
   （[fq-assemble.md](../design/fq-assemble.md)）。
4. **bbmap perfectmode**：种子-验证、完美匹配（无错配无缺口）、`ambiguous=all`
   语义在 `AbstractMapThread.java:1371` 确认（[asm-map.md](../design/asm-map.md)）。
5. **golden 数据**：`tests/bbtools/Lambda/`（`R1.2k.fq.gz`、`golden/` 等）是
   fq 逐字节对照的基准，随 fq 测试迁移批次从 pgr 迁入。

## 5. 版本差异与注意事项

- 39.38 → 40.01 的主要变化：`bbduk` 入口从 `jgi.BBDuk` 改为 `bbduk.BBDukS`；
  其余 7 个 trim 流水线工具的入口类不变（[anchr-trim-replace.md](../design/anchr-trim-replace.md) §4.2）。
- `BBMergeOverlapper` 的 Java fallback 有 JNI 版本，anchr 侧用纯 Rust 移植，
  两版 JNI 均不采用。
- 性能对照见 [bbtools-vs-anchr.md](../benchmarks/bbtools-vs-anchr.md)
  （hyperfine 实测 BBTools 39.38 8 线程 vs anchr release 单线程）。
