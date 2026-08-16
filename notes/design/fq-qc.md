# anchr QC 方案设计（FastQC / Falco 双参考）

> 2026-08-13 设计稿（实施细节版）。目标：用 Rust 实现自有 reads 质控
> （QC），替换 `templates/2_fastqc.tera.sh` 的外部 `fastqc` 调用。
> 参考：[fastqc.md](../references/fastqc.md)（语义/输出标准）、
> [falco.md](../references/falco.md)（高效实现参考）。

## 1. 背景与目标

`2_fastqc.tera.sh` 对 R/S/T 三组 reads 各跑一次外部 `fastqc`（JVM 启动慢、
多线程受 Java 限制、输出 HTML 依赖 Java 模板）。目标：

- **输出兼容**：`fastqc_data.txt`（FastQC 格式，Falco 兼容）+ `summary.txt`
  （pass/warn/fail）+ `fastqc_report.html`（FastQC 风格）；
- **性能**：单遍流式为主，rayon 并行（对齐 Falco 的 `task_queue`）；
- **零 Java 依赖**：Rust 实现，复用 pgr 基础层（fmt/fq::qual/kmer）；
- **可替换**：`2_fastqc.tera.sh` 的 `fastqc -t N ...` → `anchr fq qc -p N ...`，
  输出命名兼容模板期望。

## 2. 命令与接口设计

### 2.1 命令：`anchr fq qc`

`fq` 命令组的子命令（与 `trim-qual`/`s-filter` 并列）——QC 是 reads
质控统计，归属 FASTQ 处理组；对应模板步骤 `2_fastqc`。

```
anchr fq qc <infiles>... [OPTIONS]
  -o, --outdir <dir>       输出目录（默认 "."，不存在则创建；
                           helper 参照 pgr `outdir_arg_with_default`）
  -p, --parallel <n|auto>  reads 级并行（复用 `parallel_arg_with_default`
                           + `parse_parallel_auto`，支持 auto/1..=1024）
  -f, --format <fmt>       输出格式：html（默认）/ txt / both
      --modules <list>     限定模块（逗号分隔，默认全部启用）
      --adapter-file <f>   自定义 adapter 库（默认内置 Illumina 常用）
      --max-kmer <n>       富集类模块的 k-mer 表上限（内存控制）
      --no-gcmodel         跳过 GCModel 正态拟合（只输出 GC 直方图）
```

> **参数风格对齐**：并行一律 `--parallel`（-p，`auto` 支持），与
> `fq`/`asm`/顶层流程命令一致（不复用 FastQC 的 `-t/--threads`）；
> `-o/--outdir` 沿用 pgr `outdir_arg_with_default` 的风格（QC 输出目录，
> 需把该 helper 补进 anchr `src/cmd/args.rs`）；`--format`/`--modules` 等
> 长参数命名与现有命令的 `--kmer`/`--len`/`--idt` 风格一致（kebab-case、
> 短参数仅常用项）。

### 2.2 输出文件（对齐 Falco/FastQC）

输出到 `outdir/`（与输入 basename 无关，Falco 风格；模板适配）：

| 文件 | 内容 |
| :--- | :--- |
| `fastqc_data.txt` | `>>模块名\tPASS/WARN/FAIL` + `#列头` + 数据行 + `>>END_MODULE` |
| `summary.txt` | `模块名\t状态`（制表符分隔，可 grep） |
| `fastqc_report.html` | tera 渲染的分节报告（basic_stats/adapter/...） |

> 模板替换时 `2_fastqc.tera.sh` 的 `fastqc ... -o .` 改为 `anchr fq qc
> ../${PREFIX}1.fq.gz ... -o .`，HTML 命名从 `${PREFIX}1_fastqc.html`
> 改为 `fastqc_report.html`（或让 qc 支持 `--report-name`）。

### 2.3 输入

- FASTQ / FASTQ.gz（普通 gzip 与 BGZF，pgr `fmt::seq` 已支持）；
- 单端/双端（多文件各自输出报告，模板逐组调用）；
- stdin 支持（`anchr fq qc stdin`）；
- 质量编码：`pgr::libs::fq::qual::detect_quality_base`（比 Falco 的两类
  检测更细，兼容 FastQC 的 Illumina 1.3/1.5/1.8 判定）。

## 3. 模块架构（Rust）

### 3.1 模块 trait

```rust
pub trait QcModule: Send + Sync {
    fn name(&self) -> &'static str;          // FastQC 模块名（报告节名）
    /// 每读一次（流式）；内部累积状态
    fn consume(&mut self, rec: &SeqRecord, stats: &SharedStats) -> Result<()>;
    /// 输出模块数据（fastqc_data.txt 段 + grade）
    fn report(&self, w: &mut dyn Write) -> Result<QcGrade>;
}
```

- 大部分模块**单遍流式**（`consume` 逐 read 累积）；需要两遍的模块
  （OverRep/Kmer 的 k-mer 表）在 `consume` 里同步建表；
- `SharedStats`：BasicStats 的共享计数（reads 数、长度、GC、质量编码），
  供其他模块（如 PercentOfPairs、长度归一）引用；
- 并行：reads 分块（rayon `par_chunks`），每块内模块状态独立，结束时
  `merge`（各模块实现 `merge(&mut self, other)`；归并顺序用
  `pgr::libs::par::ordered_map` 保证确定性）——对齐 Falco
  `task_queue` + `results_collector` 的思路。

### 3.2 模块清单与优先级

| # | 模块 | 实现方式 | 内存 | 优先级 |
| :--- | :--- | :--- | :--- | :--- |
| 1 | BasicStats | 流式计数（reads/长度 min-max-median/%GC） | O(1) | M1 |
| 2 | PerBaseQualityScores | 位置分箱（BaseGroup）+ 每箱质量直方图 | O(n_groups×42) | M1 |
| 3 | PerSequenceQualityScores | 每条平均质量直方图 | O(1) | M1 |
| 4 | PerBaseSequenceContent | 位置 A/C/G/T/N 计数 | O(n_groups×5) | M1 |
| 5 | PerSequenceGCContent | GC 比例直方图（0-100） | O(1) | M1 |
| 6 | NContent | 位置 N 计数 | O(n_groups) | M1 |
| 7 | SequenceLengthDistribution | 长度直方图 | O(1) | M1 |
| 8 | AdapterContent | adapter k-mer 前缀匹配 | O(adapter 表) | M2 |
| 9 | OverRepresentedSeqs | k-mer 计数 → top 序列 | O(kmer 表) | M2 |
| 10 | KmerContent | 观察 vs 期望 k-mer 富集 | O(kmer 表) | M2 |
| 11 | DuplicationLevel | reads 去重/降采样估计 | 内存受限路径 | M2 |
| 12 | PerTileQualityScores | tile 坐标质量（仅 Illumina 命名） | O(tile×group) | M3 |
| 13 | GCModel（附于 #5） | GC 分布正态拟合 + 污染判定 | O(1) | M3 |

### 3.3 位置分箱（BaseGroup，对齐 Falco `base_groups.cpp`）

```
make_base_groups(n_bases, n_initial, n_groups_target):
  若 n_bases <= n_initial + n_groups_target: 每组 1 个位置
  否则: 前 n_initial 个位置各自一组；其余按 get_linear_interval
        (n_bases - n_initial, n_groups_target) 均分为约 n_groups_target 组
```

- 默认 `n_initial`/`n_groups_target` 对照 Falco 的 `get_default_base_groups`
  （长 reads 自动分箱，避免每位置一行）；
- 短 reads（≤ 阈值）每位置一组，输出 `#Base` 为 `1..n` 或 `1-50` 组区间
  （FastQC 的 `#Base` 列格式：`1`、`2`、... 或 `1-50`）。

## 4. 算法细节（逐模块）

### 4.1 BasicStats

- 计数：n_reads、总碱基；每 read 长度 → min/max/median（保留长度直方图
  或按需排序流）；GC 计数（G/C 碱基 / 总碱基，含 N 不计）；
- 输出 `##Basic Statistics` 节：Filename/File type/Encoding/Total Sequences/
  Sequences flagged as poor quality/Sequence length/%GC（对齐 FastQC）；
- 编码：`detect_quality_base` 结果 → `Sanger / Illumina 1.9` 等文本。

### 4.2 PerBaseQualityScores

- 每箱（位置或组）统计：`min/max/sum/sumsq` + 质量直方图
  （bins 0..=41，u64）；报告期计算 mean/median/lower/upper quartile：
  - mean = sum / count；quartile 从直方图累积（lower = 25%、median = 50%、
    upper = 75%，对齐 FastQC 的 `QualityBoxPlot`）；
- 输出 `#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile`；
- **grade**（对齐 Falco/FastQC）：median < 20 → FAIL，< 25 → WARN；
  lower quartile < 5 → FAIL，< 10 → WARN。

### 4.3 PerSequenceQualityScores

- 每条 read 平均质量（质量值求和 / 长度）→ 直方图（0..=41）；
- 输出 `#Quality\tCount`；grade：整体均值 < 20 → FAIL，< 27 → WARN。

### 4.4 PerBaseSequenceContent

- 每箱 A/C/G/T/N 计数（u64×5）；输出 `#Base\tG\tA\tT\tC`；
- grade：任意位置 G/C 或 A/T 偏差 > 20% → FAIL，> 10% → WARN
  （偏差 = |比例 - 50%|）。

### 4.5 PerSequenceGCContent + GCModel

- 每条 read GC 比例（0-100%）→ 直方图；输出 `#GC Content\tCount`；
- GCModel：对 GC 直方图做正态拟合（mean/σ，`--no-gcmodel` 跳过），
  输出理论正态曲线对照 + 文库污染判定（偏离 > 15% → WARN，> 30% → FAIL）；
- 拟合算法对照 FastQC `GCModel/` 源码（正态参数估计 + 峰检测），
  初期可只输出直方图 + 简化均值/σ，污染判定标注"待对照"。

### 4.6 NContent / SequenceLengthDistribution

- NContent：每箱 N 计数；输出 `#Base\tN-Count`；grade：任意位置 N > 20%
  → FAIL，> 5% → WARN；
- 长度分布：长度直方图；输出 `#Length\tCount`；grade：所有 read 同长
  → WARN/FAIL（对齐 FastQC 的 sequence length 判定，warn=1/fail=1 语义
  对照 Falco：全部同长 → WARN，混合 → PASS）。

### 4.7 AdapterContent（M2）

- adapter 库（内置 Illumina TruSeq/通用 adapter，`Configuration/` 对照）；
- 匹配：reads 前缀/内部 k-mer 命中 + 延伸（Falco `adapter_matcher` 思路：
  k-mer seed → 允许错配的延伸比对），统计含 adapter 的 reads 比例；
- 输出 `#Sequence\tCount\tPercentage\tPossible Source`；
- grade：> 10% → FAIL，> 5% → WARN。

### 4.8 OverRepresentedSeqs（M2）

- k-mer 计数表（复用 `pgr::libs::kmer`，`--max-kmer` 控制内存；超限转
  外部桶——参考 `fq norm` 的 `--mem` 路径）；
- 找高频 k-mer → 组装候选序列 → 与已知污染库/adapter 比对标注来源；
- 输出 top N（FastQC 默认 10）：`#Sequence\tCount\tPercentage\tPossible Source`；
- grade：单条序列占 reads > 1% → FAIL，> 0.1% → WARN。

### 4.9 KmerContent（M2）

- 全 k-mer 表（k=7 默认，FastQC 语义）观察计数 vs 期望
  （位置特异：每位置的 k-mer 期望 = 位置修正），富集 z-score；
- 输出 top 富集 k-mer 及位置；grade：z-score 超阈（FastQC 5%/2% 判定，
  Falco warn 2/fail 5）。

### 4.10 DuplicationLevel（M2）

- reads 精确去重（质量字符参与）→ 去重后 reads 比例；
- 内存受限：降采样/外部桶（对齐 `fq clump` 的 `--mem` 路径）；
- 输出 dedup 百分比曲线；grade：> 70% 重复 → WARN，> 50% → FAIL
  （对齐 Falco 0.70/0.50 的 FastQC 语义——FastQC 的 duplication 模块
  判据为 deduplicated percentage）。

### 4.11 PerTileQualityScores（M3）

- 仅 Illumina tile 命名（`instrument:run:flowcell:lane:tile:x:y`）启用；
- 按 tile × BaseGroup 的质量均值 vs 整体基线偏差；输出热图数据；
- grade：tile 偏差 > 10 → FAIL，> 5 → WARN（Falco 阈值）。

## 5. 判定阈值总表（对齐 Falco `falco_grade.cpp`）

| grader | WARN | FAIL |
| :--- | :--- | :--- |
| quality_base_median | < 25 | < 20 |
| quality_base_lower | < 10 | < 5 |
| quality_sequence | < 27 | < 20 |
| n_content | > 0.05 | > 0.20 |
| sequence（A/T/G/C 偏差） | > 10% | > 20% |
| sequence_length | 全部同长 | 全部同长 |
| gc_sequence | 偏离 > 15% | 偏离 > 30% |
| duplication | > 70% | > 50% |
| overrepresented | > 0.1% | > 1% |
| kmer | z > 2 | z > 5 |
| tile | 偏差 > 5 | 偏差 > 10 |
| adapter | > 5% | > 10% |

> 阈值以 Falco 2.0.1 为基线，FastQC-master 语义复核；`graders` 做成可配置
> （`--grade-file` 或模块常量），版本差异记录在 golden 对照。

## 6. 输出格式细节

### 6.1 fastqc_data.txt（对齐 FastQC/Falco）

```
##FastQC\t0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Filename\tin.fq.gz
...
>>END_MODULE
>>Per base sequence quality\twarn
#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile
1\t32.1\t33\t...
...
>>END_MODULE
```

### 6.2 summary.txt

```
Module\tStatus
Basic Statistics\tpass
Per base sequence quality\twarn
...
```

### 6.3 HTML（M3，tera）

- 分节模板（`templates/qc/`）：每模块一节（统计表 + 图），参考 Falco
  `results_summary.cpp` 的分节组装；
- 图：质量箱线图/直方图用内联 SVG（Rust 生成或纯 CSS 条形图，避免
  引入绘图依赖）；`project-understanding` 的 `plot` 路线在 pgr，anchr
  先用简单 SVG；
- 不追求与 FastQC 像素级一致，结构兼容（标题、模块名、表列一致）；
  **HTML 与图不作为验收标准**——验收只看数据输出。

## 7. 复用与依赖

| 需求 | 来源 |
| :--- | :--- |
| FASTQ 读入（plain/gz/BGZF/stdin） | `pgr::libs::fmt::seq` |
| 质量编码检测 | `pgr::libs::fq::qual::detect_quality_base` |
| k-mer 表（OverRep/Kmer/Adapter） | `pgr::libs::kmer`（`--max-kmer` 控制） |
| 并行（线程数/`auto`） | `pgr::libs::sys::logical_cpus`（`parse_parallel_auto` 已用） |
| 并行归并 | `pgr::libs::par::ordered_map`（模块分块结果 merge） |
| 直方图数据层 | `pgr::libs::plot::histogram::{calc_hist, calc_density, create_table}` |
| 热图数据层（PerTile） | `pgr::libs::plot::heat::{heatmap, GcHeatmap}` |
| 渲染（HTML） | tera（已有，`templates/qc/*.html`） |
| 临时文件/外部桶 | tempfile（已有，`fq clump/norm` 模式） |

> **pgr plot 复用边界**：pgr 的 `plot` 渲染输出是 **LaTeX（pgfplots，
> tectonic 编译）**，与 QC 的 HTML/SVG 目标不同——**数据计算层复用**
> （histogram 的 bin/密度/表格、heat 的 GcHeatmap 数据结构），**渲染层
> 不用 pgr 的 LaTeX**——QC 报告统一用 tera + 内联 SVG。

**不新增**外部依赖（对比 FastQC 的 htsjdk/commons、Falco 的 C++ 依赖）；
复用优先顺序：pgr 基础层（fmt/qual/kmer/io/sys/par/plot 数据层）→
anchr 已有（rayon/tera/tempfile）→ 自实现。

## 8. 里程碑（分阶段，每阶段验证）

### M1：统计型模块（文本输出）

- 实现模块 1-7（BasicStats/PerBase/PerSeqQuality/Content/GC/N/长度）；
- `anchr fq qc` 命令 + `fastqc_data.txt` + `summary.txt` 输出；
- 并行：reads 分块 + 模块 merge；
- 验证：`cargo test` 单元/集成；与 FastQC/Falco 实跑数值对照
  （Lambda 数据 + FastQC test/data），文本数据逐字段数值一致
  （小数点位数不强制）。

> **2026-08-13 完成**：`anchr fq qc`（M1）已落地——模块 1-7、文本输出、
> `tests/cli_fq_qc.rs`（6 集成）+ libs 单元测试全绿；与 fastqc 0.12.1
> golden（`tests/qc/golden/`）逐模块数值对照：per-base/per-seq quality、
> GC/N/长度**零数值差异**；per-base content 仅 Java double 精度噪声。
> 分箱复刻已装 0.12.1 的 `BaseGroup`（与 FastQC-master 源码不同，已反汇编
> 确认）。剩余：`--parallel` 暂单线程、GC grade 简化（GCModel 在 M3）、
> M2 富集模块未开始。

> **2026-08-13 M1-M4 全部完成**：`--parallel` 并行生效（分块 + merge，
> p4 与单线程输出零差异）；M2 富集模块（Adapter/Overrep/Duplication/
> Kmer）实现，Kmer 复刻 FastQC 的 2% 采样 + binomial 检验（不完全 beta）；
> M3 GCModel grade（summary 与 golden 逐字节一致）、PerTile（无 tile
> 不输出）、HTML 报告（tera + SVG，结构兼容）；M4 `2_fastqc.tera.sh`
> 改为 `anchr fq qc`。剩余：KmerContent 显著性待真实富集数据验证、
> HTML 视觉打磨、Falco 双 golden 对照。

> **2026-08-13 补充**：KmerContent 富集路径已验证（合成全 adapter 数据，
> 2% 采样命中时输出 adapter 7-mer，obs/exp 正确）；并行分块的 kmer 采样
> 改为全局序号（修复分块破坏 skip 语义）；边界测试补齐（空输入/单 read/
> 变长/超长 reads）。

> **2026-08-13 Falco 双 golden 完成**（htslib 1.21 装好后编译 falco
> 2.0.1）：交叉验证确认 anchr 与两个参考的**底层数值一致**——Falco 的
> 单位置值反推即 FastQC 组均值（如位置 10/11 = 36.5315/36.574 → 组均值
> 36.55275）。Falco 自身与 FastQC 的实现变体（不分箱、GC 整数计数、
> duplication 舍入、grade 阈值）不是 anchr 的差异；anchr 严格对齐
> FastQC 0.12.1（主参考，已完整验证）。Falco golden 存于
> `tests/qc/golden/falco/`。

### M2：富集与重复模块

- AdapterContent / OverRepresentedSeqs / KmerContent / DuplicationLevel；
- k-mer 表复用 + 内存受限路径（`--max-kmer`、外部桶）；
- 验证：adapter/overrep 在含污染数据上的 golden 对照（构造 +
  真实数据）。

### M3：报告与剩余模块

- HTML（tera 分节 + 内联 SVG）、PerTileQuality、GCModel 拟合；
- `--modules` 限定、`--format html/txt` 等参数完善；
- 验证：HTML 结构对照（模块名/表列），tile 用 Illumina 命名数据。

### M4：模板替换与双 golden 收尾

- `2_fastqc.tera.sh` 改为 `anchr fq qc`（输出命名适配）；
- 双 golden：FastQC-master 实跑 + Falco 实跑（同一输入），anchr 输出
  与两者数据逐字段数值对齐（允许记录的已知差异：编码细分、HTML/图
  仅结构对照、小数点位数）；
- 阈值表定稿（可配置 grader）；文档更新（fastqc.md/falco.md 落地状态、
  todo 销账、project-understanding）。

## 9. 测试与 golden 策略

- **单元**：模块级合成数据（构造已知统计值）；
- **集成**：`tests/cli_qc.rs`（AnchrCmd 风格）：Lambda reads、FastQC
  `test/data/`、Falco 测试数据；断言 fastqc_data.txt/summary.txt 关键字段；
- **golden**：`tests/qc/golden/` 存放 FastQC/Falco 实跑输出（生成脚本 +
  版本记录）。**验收口径**：`fastqc_data.txt`/`summary.txt` 的**数据
  逐字段数值一致**（如 Mean/Median/Count/百分比；小数点位数不强制，
  数值本身对齐即可）；HTML 与图只做结构对照（标题/模块名/表列），
  不逐像素/逐字节要求一致；
- **边界**：空输入、单 read、长度 0、质量编码未知、超长 reads
  （BaseGroup 分箱）、stdin、gz/BGZF。

## 10. 风险与决策点

1. **per-base quartile 的内存**：长 reads（Mb 级）逐位置直方图爆内存——
   用 BaseGroup 分箱（Falco 思路）限制为 O(n_groups×42)，`n_groups_target`
   默认对齐 Falco；
2. **GCModel 拟合语义**：FastQC 的污染判定有版本细节，M3 前先对照
   `GCModel/` 源码确认，初期输出直方图 + 均值/σ；
3. **Duplication 内存**：精确去重只对短 reads；长 reads/大数据走降采样
   或外部桶（复用 `fq clump` 模式）；
4. **HTML 兼容范围**：决定"结构兼容不追求像素一致"，在 M3 定稿；
5. **阈值版本差异**：FastQC-master 与已装 fastqc 版本阈值可能有出入，
  以 Falco 2.0.1 为基线 + FastQC-master 复核，记录差异表；
6. **k-mer 模块内存**：OverRep/Kmer 的 k-mer 表是主要内存源，`--max-kmer`
  上限 + 外部桶兜底。

## 11. 落地状态

设计稿（本文）。实施从 M1 开始，模块化推进；每个里程碑对照
`notes/references/fastqc.md` + `falco.md` 的源码细节（BaseGroup、
adapter_matcher、grade 阈值、报告分节）。
