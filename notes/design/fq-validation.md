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
| 2 | 2_insert_size（insert 分布） | bbtools statInsertSize / reformat-ihist | `fq merge --ihist`（bbmerge ihist 格式） | ⚠️ ihist 格式已对齐 golden；真实数据未对照 | 真实双端 insert size 分布对照；reformat-ihist 仍属独立缺口（todo 挂账） |
| 3 | 2_fastk（k-mer 谱 + 基因组特征） | fastk -NTable/-Histex + GeneScope R | `pgr kmer table/hist` + `pgr kmer gsize --model`（GenomeScope 拟合） | ✅ 自实现（用户确认替代） | pgr 侧验证；anchr 流程整合待端到端 |
| 4 | 2_trim 1：clumpify `--dedupe`（k-mer 签名聚类去重） | clumpify.sh | `fq clump` | ✅ golden 逐字节 + 基准 2.19x | `--parallel` 在真实数据上的复测 |
| 5 | 2_trim 5：bbduk trim（`ktrim=r k=trimk mink=11 hdist=1 tbo tpe maxns=0 minlen=60 qtrim=r trimq=15 ftm=5`） | jgi.BBDuk | `fq clean`（复用 `libs/fq/trim_adapter`） | ✅ 19 组 trim 变体逐字节一致（39.38 `ordered=t`） | 真实双端 reads；tbo/tpe 大参数空间 |
| 6 | 2_trim 6：bbduk filter（`k=matchk cardinality` 参考库过滤） | jgi.BBDuk kfilter | `fq filter` | ✅ golden + k 变体一致 | 真实污染数据（adapter/artifact） |
| 7 | 2_trim 3：bbnorm（`bits=16` 近似 cutoff） | bbnorm.sh | `fq norm` | ⚠️ 精确 canonical 表；与 bbnorm 有 ~21 对边界差异，定义未定稿 | 用户裁定精确 vs 近似并文档定稿（`fq-trim-replace.md` §4.8）；外部桶路径真实大数据 |
| 8 | 2_trim 4：reformat sample（目标碱基数降采样） | reformat.sh | `fq sample` | ✅ golden + 基准 6.67x | — |
| 9 | 2_trim 8：repair（交错 → R1/R2/singles） | repair.sh | `fq split` | ✅ golden + 基准 8.13x | — |
| 10 | 2_trim 9：sickle（`-q/-l` 多阈值扫描） | sickle | `fq trim-qual` | ✅ 替换 sickle（滑窗/Mott） | — |
| 11 | 2_trim 2：filterbytile（flowcell tile 过滤，可选） | — | 不做（老数据专属，明确不做） | — | — |
| 12 | 2_merge phase 1：ecco（overlap 区纠错，不合并） | bbmerge.sh `ecco mix vstrict` | `fq ec-overlap` | ✅ golden（ecco_sub.fq.gz） | 真实双端 |
| 13 | 2_merge phase 2：clumpify ecc（clump 共识纠错） | clumpify.sh | 跳过（与 phase 3 冗余，用户反馈常卡） | 明确不做 | — |
| 14 | 2_merge phase 3：tadpole ecc（k-mer 图纠错 + 丢弃坏 read） | tadpole.sh `ecc tossjunk tossdepth=2 tossuncorrectable` | `fq ec-kmer` | ✅ golden（ecct_sub.fq.gz） | 真实数据纠错效果（mismatch 下降量化） |
| 15 | 2_merge bbmerge merge（PE overlap 合并） | bbmerge.sh | `fq merge` | ✅ golden（merge*/merge4*：net/classic/ecco/ihist 全一致） | MR 路径端到端（merged reads → multik 输入） |
| 16 | 2_quorum（**reads 筛选**：quorum 修正过的 reads 带 `:sub:`/`trunc` 标记且被丢弃，pe.cor = 未修正原始序列） | 外部 quorum（quorum_error_correct_reads） | `fq s-filter`（检查 quorum 信号：无高质量 anchor / truncation / 会 substitution 的碱基 + Poisson 碰撞；保留原样或丢弃） | ⚠️ 自实现 + 语义对齐（用户确认）；**仅 4 个测试，从未与 quorum 输出对照** | **真实 reads 上跑 quorum 与 `fq s-filter`，对照丢弃集合与标记判定**；参数（k/anchor-count/prior/poisson）定稿 |
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
