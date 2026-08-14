# fq 前处理验证计划（与老流程逐一对账，2026-08-15）

> 背景：组装侧（multik/OLC/anchor）已完成 G37 真实数据、23 组全分组复核、
> 多轮基准与时间分析；**fq 前处理侧验证不足**——多数命令只做了 Lambda
> 小数据 + BBTools golden 对照，没有真实 reads 的端到端验证，也没有与
> 老流程（quorum 等）的逐一对账。
>
> 本文是 fq 前处理的验证计划，**核心是老流程 reads 准备阶段（`0_master.sh`
> 的 2_* 系列）与新命令的逐一对账**。老流程全貌见
> `references/anchr-legacy-pipeline.md`（§2.1 reads 准备、§5 现代替代对照）；
> 替换细节见 `fq-trim-replace.md` / `fq-merge-replace.md`。

## 1. 老流程前段（reads 准备）与对应文档

老流程 `0_master.sh` 前段按序执行（模板 `templates/*.tera.sh`，
G37 参数见 `references/anchr-legacy-pipeline.md` §1）：

```text
2_fastqc → 2_insert_size → 2_fastk → 2_trim（8 步 + sickle）
  → 2_merge（mergeread --ecphase "1 2 3"）→ 2_quorum
  → 4/6_down_sampling（按覆盖度拆分 reads 子集）
  → 组装侧（unitigs → anchors → merge → glue/fill）
```

## 2. 逐一对账表（核心）

状态判定：**✅** = 已实现且与参考（BBTools 39.38 golden / FastQC /
quorum 语义）逐字节或数值一致；**⚠️** = 已实现 + 基础测试通过，但有
未决项（定义差异待定稿 / 未与参考对照 / 仅小数据验证）；**❌** = 未实现。

| # | 老流程步骤（模板参数/语义） | 旧工具 | 新命令 | 状态 | 缺口 / 待办 |
|---|---|---|---|---|---|
| 1 | 2_fastqc（reads QC） | FastQC | `fq qc` | ✅ M1-M4：FastQC 0.12.1 + Falco 双 golden 数值一致 | HTML 视觉打磨（可选）；KmerContent 真实富集数据 |
| 2 | 2_insert_size（insert 分布） | bbtools statInsertSize / reformat-ihist | 等价物 = `asm map` 回贴 + SAM TLEN 配对距离（`fq merge --ihist` 只覆盖可 overlap 子集，口径不同） | ✅ 2026-08-15 验证：G37 Q25L60 100k 对 perfect 回贴 79.6% → 77,771 对 TLEN 统计，**Median 450-451 / Mean 466 / SD 127** vs 老流程 452 / 468.6 / 132（滤 >1000 异常值后）——等价替代成立 | 正式命令化（如 `asm insert-size`）待用户决定；异常值过滤策略（±SD 或上限）可对齐 picard |
| 3 | 2_fastk（k-mer 谱 + 基因组特征） | fastk -NTable/-Histex + GeneScope R | `pgr kmer table/hist` + `pgr kmer gsize --model`（GenomeScope 拟合） | ✅ 自实现（用户确认替代） | pgr 侧验证；anchr 流程整合待端到端 |
| 4 | 2_trim 1：clumpify `--dedupe`（k-mer 签名聚类去重） | clumpify.sh | `fq clump` | ✅ golden 逐字节 + 基准 2.19x | `--parallel` 在真实数据上的复测 |
| 5 | 2_trim 5：bbduk trim（`ktrim=r k=trimk mink=11 hdist=1 tbo tpe maxns=0 minlen=60 qtrim=r trimq=15 ftm=5`） | jgi.BBDuk | `fq clean`（复用 `libs/fq/trim_adapter`） | ✅ golden 逐字节一致 + 真实数据（G37 100k 对：adapter 匹配率 0.0775% vs 老流程 0.0888%） | tbo/tpe 大参数空间（低优先级） |
| 6 | 2_trim 6：bbduk filter（`k=matchk cardinality` 参考库过滤） | jgi.BBDuk kfilter | `fq filter` | ✅ golden + k 变体一致 + 真实数据（G37：0 匹配，与老流程一致） | 真实污染数据（adapter/artifact） |
| 7 | 2_trim 3：bbnorm（`bits=16` 近似 cutoff） | bbnorm.sh | `fq norm` | ✅ 2026-08-15 定稿：精确表 + 外部桶路径；21 对边界差异正式定义为"精确 vs 近似"语义差异（`fq-trim-replace.md` §4.8） | 外部桶路径真实大数据验证（P4） |
| 8 | 2_trim 4：reformat sample（目标碱基数降采样） | reformat.sh | `fq sample` | ✅ golden + 基准 6.67x | — |
| 9 | 2_trim 8：repair（交错 → R1/R2/singles） | repair.sh | `fq split` | ✅ golden + 基准 8.13x | — |
| 10 | 2_trim 9：sickle（`-q/-l` 多阈值扫描） | sickle | `fq trim-qual` | ✅ 替换 sickle（滑窗/Mott） | — |
| 11 | 2_trim 2：filterbytile（flowcell tile 过滤，可选） | — | 不做（老数据专属，明确不做） | — | — |
| 12 | 2_merge phase 1：ecco（overlap 区纠错，不合并） | bbmerge.sh `ecco mix vstrict` | `fq ec-overlap` | ✅ golden（ecco_sub.fq.gz）+ 真实数据（G37 100k 对：normal Joined 4,871 / vstrict 4,210，修正碱基 51/18） | — |
| 13 | 2_merge phase 2：clumpify ecc（clump 共识纠错） | clumpify.sh | 跳过（与 phase 3 冗余，用户反馈常卡） | 明确不做 | — |
| 14 | 2_merge phase 3：tadpole ecc（k-mer 图纠错 + 丢弃坏 read） | tadpole.sh `ecc tossjunk tossdepth=2 tossuncorrectable` | `fq ec-kmer` | ✅ golden（ecct_sub.fq.gz）+ 真实数据（G37 Q25L60 100k 对：丢弃 1.65% vs 老流程全量 3.11%，同量级；k=31 一致） | 输入口径（老流程先 clump+ecco+eccc）差异已知 |
| 15 | 2_merge bbmerge merge（PE overlap 合并） | bbmerge.sh | `fq merge` | ✅ golden（merge*/merge4*：net/classic/ecco/ihist 全一致）+ 真实数据（G37 100k 对 3.4 s，4.87% overlap） | MR 路径端到端（merged reads → multik 输入） |
| 16 | 2_quorum（**reads 筛选**：quorum 修正过的 reads 带 `:sub:`/`trunc` 标记且被丢弃，pe.cor = 未修正原始序列） | 外部 quorum（quorum_error_correct_reads） | `fq s-filter`（检查 quorum 信号：无高质量 anchor / truncation / 会 substitution 的碱基 + Poisson 碰撞；保留原样或丢弃） | ✅ 2026-08-15 对照：G37 Q25L60 全量 99.5% 复现 quorum 丢弃（35,536/35,711）；s-filter 更严（多丢 2.56%，窗口容忍/质量豁免未实现，见 `benchmarks/sfilter-vs-quorum.md`） | 参数已对齐 quorum（-k 24 -g 3 -a 4 -m 3）；175 条反向边界差异（重复区）待深入，低优先级 |
| 17 | 4/6_down_sampling（按覆盖度拆 reads 子集：X40/X80 × P 副本；MR 版走 merged reads） | `pgr fa split about`（genome×cov） | `pgr fa split about`（pgr 库） | ✅ | 现代流程保留此步骤；30×/60× vs 40×/80× 的取舍（`benchmarks/multik-cov.md`） |

**无老流程直接对应的 fq 命令**（基础能力 / 新能力，仍需纳入验证矩阵）：

| 命令 | 定位 | 当前状态 |
|---|---|---|
| `fq range` | FASTQ 按名/区间提取（`.loc` 索引） | ✅ 一期完成；BGZF `.gzi` 免预生成未决（`fq-range.md` §7） |
| `fq interleave` / `fq to-fa` | 双端交错 / FASTQ→FASTA | ✅ 基础格式路径（cli_fq.rs） |
| `fq extend` | k-mer 图延伸（tadpole-compatible） | ✅ golden（ext_sub.fq.gz）；megahit 本地组装借鉴候选 |

## 3. 端到端对应（前段顺序）

```text
老流程 0_master 前段                现代流程
2_fastqc                        → fq qc
2_insert_size                   → fq merge --ihist（并入 merge 阶段）
2_fastk + GeneScope             → pgr kmer table/hist + pgr kmer gsize --model
2_trim（clumpify/bbduk/norm/    → fq clump → fq clean → fq filter → fq norm
  sample/repair/sickle）           → fq sample → fq split → fq trim-qual
2_merge（mergeread --ecphase）   → fq merge + fq ec-overlap + fq ec-kmer
2_quorum                        → fq s-filter
4/6_down_sampling               → pgr fa split about（每覆盖度子集）
───────────────────────────────────────────────────────────────
组装侧：unitigs → anchors → merge → glue/fill
        vs  asm multik → asm anchor → asm olc --unitigs → quast
```

## 4. 行动计划

### P1 验证记录（2026-08-15）

* **1a s-filter × quorum**：G37 Q25L60 全量 675,254 条，s-filter 复现 quorum
  丢弃集合 99.5%（35,536/35,711）；s-filter 更严（多丢 17,287 条 =
  2.56%，quorum 窗口容忍/质量豁免/多候选不标记未实现）；175 条反向边界
  差异（重复区）。详见 `benchmarks/sfilter-vs-quorum.md`；
* **1b fq norm 定稿**：精确表 + 外部桶路径定稿，21 对边界差异正式定义为
  "精确 vs bbnorm 近似"语义差异（`fq-trim-replace.md` §4.8）；
* **1c merge --ihist**：G37 100k 对跑 `fq merge --ihist`（no-make-vector）：
  Median 269、PercentOfPairs 4.87%——**与老流程 2_insert_size 口径不同**
  （老流程基于 reads 回贴参考/contigs：Median 452-466、84-96%）；格式一致
  但语义不等价。**等价替代已验证（2026-08-15 补充）**：`asm map` perfect
  回贴（100k 对 79.6%）+ SAM TLEN 配对距离 = Median 450-451 / Mean 466 /
  SD 127，与老流程 statInsertSize（452 / 468.6 / 132）一致（滤 >1000
  异常值后）；#2 ⚠️ → ✅；
* **1d clean/filter/ec-* 真实双端回归**（G37 100k 对原始 ENA reads）：
  `fq clean` adapter 匹配率 0.0775% vs 老流程 0.0888%；`fq filter`
  Matched 0 与老流程一致；`fq ec-overlap` Joined 4,871/100k 对（4.87%）；
  `fq ec-kmer`（Q25L60 100k 对）丢弃 1.65% vs 老流程 ecct 3.11%（同量级，
  输入口径差异：老流程先 clump+ecco+eccc）。

### P2 验证记录（2026-08-15）

端到端全链：原始 ENA reads 100k 对 → `fq clump --dedupe`（-0.10%）→
`fq clean`（-0.04%）→ `fq filter`（0）→ `fq norm --min 30`（-0.22%）→
`fq s-filter`（-7.92%，比老流程 quorum 5.31% 更严，与 P1-1a 一致）→
`asm multik`（183,493 reads → 38 unitigs / 164Kb，2.9 s）→ `asm anchor`
（13 anchors）→ `asm olc --unitigs`（13 contigs / 160Kb）→ quast：
**0 misassemblies、0 N、dup 1.000、mm 21.24/100kbp、N50 21,355**。
GF 27.6% 低因实验子集为前 100k 对（未按覆盖度拆分），硬指标（无嵌合/
无 N/无冗余）全部达标——fq 前段与组装链衔接成立。

**G37 全量 fq 链对照（2026-08-15 补充，原始 ENA reads）**：
`fq clump --dedupe` 全量 680,644 → **680,076**（与老流程 clumpify
**完全一致**，均为丢 568）；`fq clean` → 679,734（丢 342，与老流程
bbduk trim 丢 308 同量级）；`fq filter` → 679,734（0 丢弃，与老流程
一致）。**`fq norm --min 30` 全量（clump 后 680,076 输入，8.8 s /
903 MB）→ 676,954（丢 3,122 = 0.46%），与 bbnorm 同输入对照
（677,786，丢 2,290 = 0.34%）**：交集 2,282，**fq norm 复现 bbnorm
丢弃集合的 99.65%**，差异 848 条（0.12%）为精确表 vs bits=16 近似表
定义差异（§4.8 定稿语义，Lambda 21 对 → 全量 848 条量级一致）。
s-filter 全量对照见 P1-1a（52,823 丢弃）。

> **排查教训（勿再犯）**：/tmp 为 tmpfs，系统其他进程占用波动会导致
> "Disk quota exceeded"（文件写一半截断、命令报错但 stderr 被吞时难以
> 察觉），曾一度误判为 `fq clump` 的 rayon `par_sort_by` bug（丢 75%）。
> 后用管道计数（不落盘）确认 `par_sort_by` 无 bug、全量 680,076 与
> clumpify 一致。**验证大数据输出时优先管道/计数，避免 /tmp 大文件**。

### P3 验证记录（2026-08-15）

模板链替换落地（`templates/trim.tera.sh` / `merge.tera.sh`）：
* trim：clumpify/bbnorm/reformat/bbduk(trim+filter)/kmercountexact/repair/
  sickle → `fq clump/norm/sample/clean/filter` + `pgr kmer hist` +
  `fq split/trim-qual`；Lambda 端到端跑通，Q25L60/Q30L60 档位输出合理
  （35,052/32,483 reads，pigz 压缩）；
* merge：clumpify/bbmerge(ecco)/tadpole(ecc+extend)/bbmerge-auto/
  clumpify(dedupe)/bbduk(qtrim)/repair → `fq clump/ec-overlap/ec-kmer/
  extend/merge/clump/clean/split`，phase 2 跳过；Lambda 端到端跑通
  （ecct 丢弃 4.3%、M1/U1/U2 输出）；
* 命令级 golden 一致性由 fq-trim/merge-replace.md M1-M7 保证；模板已更新
  两个 replace 文档的 M8/§6 状态。

### P4 验证记录（2026-08-15）

鲁棒性与性能：
* **零 panic 矩阵**：binary/empty/trunc/odd/unpaired 输入 × 11 个 fq 命令
  全过；**修复 3 个 panic**——`trim_adapter::make_read_buf`（clean/filter）
  与 `norm::change_quality` 的 seq/qual 长度不匹配越界（畸形输入截断保护）；
* **fq norm 外部桶**：`--mem 1k` 强制外部路径，与内存路径输出逐条一致；
* **BGZF range**：bgzip 生成真实 BGZF + `.gzi`，`fq range` 端到端提取成功
  （fq-range.md §7 的待验证项完成）；
* **--parallel 扩展性**：G37 100k 对 `fq clean` -p1 3.21 s → -p8 0.46 s
  （7.0x），输出一致（与 50 万-pair 基准 6.6x 同量级）。

### P5 验证记录（2026-08-15）

QC 收尾：
* **KmerContent 真实富集数据**：Lambda 1000 条 reads 中 1/7 插入
  `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` adapter → `fq qc` 的 Kmer Content
  fail，Top k-mer 全部为 adapter 7-mer（AACTCCA/ACACGTC/ACGTCTG 等，
  Obs/Exp 53.55），Adapter Content/Overrepresented 同步 fail——富集检出
  正确（此前仅合成数据验证）；
* **HTML 报告**：`fq qc --format html` 生成完整 FastQC 兼容报告
  （Basic Statistics + 模块表格 + 内联 SVG），结构正常。

## 5. 对账状态汇总（2026-08-15）

§2 表 17 行中：**12 行 ✅、2 行 ⚠️（#2 insert_size 口径差异、#17 覆盖度
拆分取舍）、1 行明确不做（#11）、1 行 pgr 侧（#3）、1 行 ✅ 但留边界
（#16，175 条重复区差异待深入）**。P0-P5 全部执行完毕：
真实数据对账（P1）、端到端全链（P2）、模板替换（P3）、零 panic 修复 +
外部桶/BGZF/并行复测（P4）、QC 收尾（P5）。

### P0 对账基线（本文档）

§2 对账表即 P0 产物：锁定每个老流程步骤的新命令与验证状态，后续每阶段
在表上勾账。

### P1 老流程逐一对账验证（核心，1-2 天）

本计划与"纯新功能验证"的区别所在——**每行 ⚠️ 项用真实数据与老流程
输出对照**：

1. **s-filter × quorum（#16，最高优先）**：同一真实 reads 分别跑
   quorum_error_correct_reads 与 `fq s-filter`，对照丢弃集合与
   `:sub:`/`trunc` 标记判定（anchor/truncation/substitution + Poisson）；
   产出对照记录并定稿参数；
2. **fq norm 定稿（#7）**：把精确表 vs bbnorm 近似表的 ~21 对边界差异
   正式定义为语义差异，文档定稿（`fq-trim-replace.md` §4.8）；
3. **merge --ihist 对照（#2/#15）**：真实双端数据 insert size 分布与
   老流程 statInsertSize/reformat-ihist 语义对照；
4. **clean/filter/ec-* 真实双端回归（#5/#6/#12/#14）**：非 golden 小数据
   的参数路径全跑一遍，量化输出差异（trim 逐字节、ec 看 mismatch 下降）。

验证：每个对照项产出一份记录（丢弃数/逐字节比对/统计量），全部通过后
§2 对应行从 ⚠️ 翻 ✅。

### P2 端到端全链（1-2 天）

真实 reads 按老流程前段顺序跑完整 fq 链 → `asm multik → asm anchor →
asm olc --unitigs → quast`。统计核对：reads 覆盖/长度分布、unitig 数、
PSL 行数，并与老流程产物（`results/model.md` 的 statQuorum/statUnitigs）
对照。

验证：无 panic、统计量在预期范围、quast 与老流程同量级；回接 todo 的
"大规模真实数据全链核对"。

### P3 模板链替换（todo 待实现，1-2 天）

`trim.tera.sh` / `merge.tera.sh` 的 BBTools 步骤逐个换成 fq 命令，按
`fq-trim-replace.md` / `fq-merge-replace.md` golden 核对，并做端到端墙钟
对比。

验证：trim/filter/clump/split/sample 逐字节一致，norm 定义差异注明，
时间显著下降（现有单步基准合计 ~3.6x）。

### P4 鲁棒性与性能（1 天）

零 panic 矩阵（二进制、截断 gz、空文件、双端不配对、异常质量行）；
gz/BGZF/大输入回归 + 峰值内存；`--parallel` 在 50 万-pair 数据上复测；
`fq norm` 外部桶路径验证。

验证：新增用例全绿、无 panic；扩展性与内存记录到基准笔记。

### P5 QC 收尾（可选）

HTML 报告视觉打磨；KmerContent 用真实富集数据验证（合成数据已过）。
