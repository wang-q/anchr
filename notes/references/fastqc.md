# FastQC: reads 质控报告参考工具

> 整理于 2026-08-13，源自对仓库根 `FastQC-master/` 源码的分析。
> anchr 的 `templates/2_fastqc.tera.sh` 当前调用外部 `fastqc` 命令做 reads
> 质控；本文记录 FastQC 的结构与核心模块，作为后续用 Rust 迁移到 anchr
> （替换外部 fastqc 调用）的参考。实现层面的第二参考是
> [falco.md](falco.md)（FastQC 的 C++ 仿制实现，输出格式兼容、组件边界清晰）。

## 1. 简介

`FastQC`（Babraham Bioinformatics，Simon Andrews，Java）是 reads 质控的事实
标准工具：读取 FASTQ（或 BAM/SAM），按 14 个 QC 模块逐项统计，输出 HTML 报告
（+ 可选 zip/txt）。本仓库内的 `FastQC-master/` 是其完整源码发行（Java 源码 +
  jar + Perl 启动脚本），目录自带 `.gitignore`（`*`），仅作参考，不被 git 跟踪。

## 2. 目录结构

| 路径 | 内容 |
| :--- | :--- |
| `fastqc.pl` / `fastqc` / `run_fastqc.bat` | 启动脚本（Perl/shell/bat） |
| `uk/ac/babraham/FastQC/` | Java 源码（Modules/、Sequences/、Graphs/、Report/ 等） |
| `Templates/` | HTML 报告模板（ftl 风格，含 Icons/） |
| `Configuration/` | 模块配置（adapter 序列等） |
| `Help/` | 用户文档（3 Analysis Modules 等） |
| `test/` | 单元/集成/UI 测试（含测试数据） |
| `*.jar`（htsjdk/commons-*/cisd-jhdf5） | 依赖库（HTSeq、HDF5 等） |

## 3. 与 anchr 的关系

- **流程位置**：`templates/2_fastqc.tera.sh` 的质控步骤（`fastqc -t
  {{opt.parallel}}`，输出 `1_fastqc.html`），属于 0/2 模板链的 reads 质控；
- **迁移动机**：与 BBTools 替换同理——外部 Java 工具慢、依赖 JVM，目标用
  Rust 实现 QC 统计，输出与 FastQC golden 对照；
- **可复用基础**：质量编码检测可用 `pgr::libs::fq::qual::detect_quality_base`；
  FASTQ 读入用 `pgr::libs::fmt`；GC/k-mer 统计可参考 anchr 已有的
  `fq` 命令（kmer 表等）。

## 4. 核心模块分析（Modules/，14 个 QC 模块）

| 模块 | 统计内容 | 算法要点 |
| :--- | :--- | :--- |
| `BasicStats` | reads 数、序列长度、%GC | 基本统计 |
| `PerBaseQualityScores` | 按位置质量分布 | `QualityBoxPlot` 箱线图（mean/median/下上四分位），`BaseGroup` 位置分组（长序列自动分箱） |
| `PerSequenceQualityScores` | 每条 reads 平均质量分布 | 直方图 |
| `PerBaseSequenceContent` | 按位置 A/C/G/T 组成 | 位置碱基比例 |
| `PerSequenceGCContent` / `GCModel` | GC 含量分布 vs 理论正态 | GC 模型拟合 |
| `NContent` | 按位置 N 比例 | 位置 N 统计 |
| `SequenceLengthDistribution` | 长度分布 | 直方图 |
| `DuplicationLevel` | 重复 reads 水平 | 序列去重/降采样估计 |
| `OverRepresentedSeqs` | 过度代表序列 | k-mer/文库污染检测（报告 top 序列） |
| `AdapterContent` | 接头含量 | k-mer 前缀匹配（`Configuration/` 的 adapter 库） |
| `KmerContent` | k-mer 富集 | 观察 vs 期望计数（位置特异） |
| `PerTileQualityScores` | 按 tile 位置质量 | 测序 tile 偏差 |

输出模型：每个模块产出"通过/警告/失败"三态 + 图/表数据，报告由
`Report/HTMLReportArchive` 用 `Templates/` 渲染。

## 5. 输入与编码

- 输入 FASTQ/FASTQ.gz（htsjdk 读 BGZF/普通 gzip）、BAM/SAM、Casava 拆分名；
- 质量编码检测（Sanger/Illumina 1.3+/1.5/1.8）在读取时推断——对应
  `pgr::libs::fq::qual` 的 `detect_quality_base` 语义；
- 长 reads（如 PacBio/ONT）下 per-base 模块按 `BaseGroup` 分箱，避免
  每位置一行。

## 6. 迁移考量（设计提案，未定稿）

- **输出形态**：FastQC 的 HTML 报告是用户习惯的交付物；Rust 迁移可先对齐
  txt 统计（FastQC 的 `--format txt` 或 zip 内数据），HTML 用 anchr 的
  `templates/`（tera）渲染，golden 对照参照 BBTools 替换的既有流程
  （`anchr-trim-replace.md` 的方法论）；
- **模块裁剪**：reads 质控主用 BasicStats/PerBaseQuality/
  PerSequenceQuality/GC/NContent/AdapterContent/OverRepresented；
  PerTile（需 tile 坐标）与 DuplicationLevel（大内存）可后置；
- **与现有命令的关系**：`fq` 命令组的统计原语（trim-qual 的质控、
  `fq s-filter` 的 k-mer 计数）可复用；`detect_quality_base` 已在 pgr；
- **golden 数据**：FastQC 自带的 `test/data/`（fastq 样例）可作对照输入，
  输出与 FastQC-master 实跑逐字节/数值核对。

## 7. 注意事项

- FastQC 的 `GCModel` 有正态拟合与"测序文库污染"判定，语义细节需在迁移时
  对照源码确认（`uk/ac/babraham/FastQC/Modules/GCModel/`）；
- `OverRepresentedSeqs`/`AdapterContent`/`KmerContent` 的阈值与 k-mer 算法
  各自独立，迁移优先级建议：统计型模块（BasicStats/质量/GC/N/长度）先行，
  富集型（k-mer/adapter）后置；
- 版本基线：仓库 `FastQC-master`（2026-08 下载）；本地 PATH 的 `fastqc`
  版本需记录，golden 以固定版本生成。
