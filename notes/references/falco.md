# Falco: FastQC 的 C++ 仿制实现（QC 方案第二参考）

> 整理于 2026-08-13、2026-08-17 对照源码逐项复核修订，源自对仓库根 `falco-2.0.1-Source/`（v2.0.1，2026-08-03）源码的分析。
> 与 [fastqc.md](fastqc.md) 互补：FastQC 是"语义/输出标准"，Falco 是
> "高效实现"参考——anchr 设计自有 QC 方案时同时参考两者。

## 1. 简介

`Falco`（Smith Lab，Andrew D. Smith）是 FastQC 的开源仿制实现（C++，
MIT），v2.0 重写后结构与 FastQC 不再一一对应，但输出保持 FastQC 兼容。
本仓库内的 `falco-2.0.1-Source/` 是源码发行（6.3M），自带 `.gitignore`
（`*`），仅作参考。

**输出**（`falco -o output input.fq` 生成 3 个文件）：
- `fastqc_data.txt`：QC 指标文本摘要（FastQC 兼容）；
- `fastqc_report.html`：可视化 HTML 报告；
- `summary.txt`：每个分析的 pass/warn/fail 判定（制表符分隔）。

### 1.1 核心算法与流程总览（先读这节）

**一句话**：按后缀/内容识别输入格式（FASTQ plain/gz/BGZF、SAM、BAM）→
`task_queue` 把 reads 切成任务块，`analyze(n_threads, mode)` 起 n 个 worker
**线程级并行**消费，每个 worker 持有 per-file 的 `results_collector` 累积
各项统计 → 合并 → `falco_grade` 对各模块做 pass/warn/fail 判定 →
三路输出（`fastqc_data.txt` 文本 / `summary.txt` 判定 / `html.cpp` 报告）。

```
输入文件（FASTQ/FASTQ.gz/FASTQ.BGZF/SAM/BAM）
  → 读入层（fastq_file / fastq_gz_file / fastq_bgzf_file / sam_file / bam_file）
  → task_queue 任务块 → analyze() n_threads 个 worker 并行
      每 worker × 每 file 一份 results_collector：
      basic stats（长度/GC/N/编码检测 unknown/sanger/solexa）
      per-base/per-tile 质量（tile_processor）
      adapter 含量（adapter_matcher 前缀 k-mer）、contaminants、kmer 富集
      duplication（original_duplicates + duplication_results）
      长 reads 位置分箱（base_groups）
  → 合并各 worker 的 collector → falco_grade（pass/warn/fail）
  → fastqc_data.txt + summary.txt + fastqc_report.html
```

| 核心块 | 机制 | 文件 |
|---|---|---|
| 读入层 | 按格式分层（plain/gz/BGZF/SAM/BAM），BGZF 走块解压 | `fastq_*_file.cpp`、`sam/bam_file.cpp` |
| 并行模型 | `task_queue` 生产者-消费者 + 每 worker 独立 `results_collector`，免锁累积后合并 | `falco_analyzer.cpp`、`task_queue.hpp` |
| 编码检测 | 只分 sanger(+33)/solexa(+64)/unknown 三类（FastQC 有更细分档） | `quality_score.cpp` |
| 位置分箱 | 长 reads 的 per-base 统计按 `base_groups` 分组 | `base_groups.cpp` |
| 富集类 | adapter 前缀匹配 / contaminants / k-mer 计数 / duplication | `adapter_*.cpp`、`kmer_counter.cpp` 等 |
| 判定与输出 | `falco_grade` 阈值判定；文本 + HTML 分节报告 | `falco_grade.cpp`、`results_summary.cpp`、`html.cpp` |

## 2. 源码结构（src/，29 个 .cpp）

| 文件 | 职责 |
| :--- | :--- |
| `falco_analyzer.cpp/hpp` | 分析入口（`analyze(n_threads, mode, ...)` 并行） |
| `fastq_file.cpp` / `fastq_gz_file.cpp` / `fastq_bgzf_file.cpp` | FASTQ 读入（plain/gzip/BGZF） |
| `sam_file.cpp` / `bam_file.cpp` | SAM/BAM 输入 |
| `quality_score.cpp/hpp` | 质量编码检测（`encoding`：unknown/sanger/solexa） |
| `base_groups.cpp/hpp` | 位置分箱（长 reads 的 per-base 分组） |
| `adapter_matcher.cpp` / `adapter_set.cpp` | 接头含量（k-mer 前缀匹配） |
| `contaminants.cpp` | 过度代表序列 |
| `kmer_counter.cpp` | k-mer 富集 |
| `original_duplicates.cpp` / `duplication_results.cpp` | 重复水平 |
| `tile_processor.cpp` | per-tile 质量 |
| `results_summary.cpp` / `report.cpp` / `html.cpp` | 摘要/报告生成（HTML 分节） |
| `falco_grade.cpp` | pass/warn/fail 判定 |
| `task_queue.hpp` / `run_mode.cpp` | 并行任务调度 |

## 3. 与 FastQC 的差异（对移植的启示）

- **编码检测简化**：只分 Sanger（Phred+33）/Solexa（+64）两类，`unknown`
  兜底；FastQC 有更细的 Illumina 1.3/1.5/1.8 细分——anchr 用 pgr
  `detect_quality_base` 可覆盖更全；
- **模块平铺**：FastQC 的 `Modules/` 14 个类，Falco v2 是平铺的 analyzer
  + 各统计组件，无统一 module 接口——anchr 可自行设计模块 trait；
- **HTML 分节**：`results_summary.cpp` 按 basic_stats/adapter/duplication
  等节组装——报告结构与 FastQC 兼容但实现更直接；
- **并行**：`task_queue` + `n_threads`，reads 级并行——对应 anchr 的
  rayon 模式（`fq` 命令已用）。

## 4. 对 anchr QC 方案的价值

1. **golden 兼容**：Falco 的 `fastqc_data.txt` 与 FastQC 兼容，且 Falco 是
   C++（无 JVM 依赖、可重复编译）——可作为 anchr Rust 实现的**中间对照**
   （FastQC 语义太宽、Falco 数值确定，双 golden 更稳）；
2. **模块裁剪参考**：Falco 的组件边界（adapter/contaminants/kmer/
   duplication/tile）比 FastQC 更清晰，便于挑出 anchr 先做的统计型模块；
3. **输出管线**：HTML 分节 + summary.txt 的组装方式可直接借鉴到
   `templates/`（tera）渲染；
4. **编码/读入**：确认 Phred 检测与多格式输入（plain/gz/BGZF/SAM/BAM）
   的分层，anchr 复用 pgr `fmt` + `fq::qual` 即可覆盖。

## 5. 落地状态

设计提案，未定稿——与 [fastqc.md](fastqc.md) §6 合并评估：
- 优先级：BasicStats/质量/GC/N/长度（统计型）先行，adapter/kmer/
  duplication/tile（富集/大内存）后置；
- golden 对照：FastQC-master 实跑 + Falco 实跑（同一输入，双基准），
  anchr 输出与两者数值/格式对齐；
- 报告：`fastqc_data.txt`/`summary.txt` 文本先行，HTML 用 tera 渲染
  （参考 Falco 的 `results_summary` 分节）。
