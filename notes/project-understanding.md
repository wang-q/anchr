# anchr 项目理解

> 整理于 2026-08-13，fq/asm 迁移（`notes/design/fq-asm-migrate.md`）完成后。
> 本文是项目整体理解（架构、命令模块、核心库、现状评估、设计模式），
> 含各文档的索引与定位，需要时查阅。

## 1. 项目定位

`anchr`（**A**ssembler of **N**-free **CHR**omosomes）是染色体级组装流程编排器。
它不重复造格式 I/O 与基础算法——FASTA/FASTQ 读取、Phred 编码、k-mer 表、
PAF 解析等基础层来自 **pgr** crate（`pgr::libs::*`）；anchr 负责 reads 处理
（`fq` 命令组）、组装（`asm` 命令组）以及染色体流程编排（模板 + 顶层命令）。

### 1.1 与 pgr 的分工

| 层 | 内容 | 归属 |
| :--- | :--- | :--- |
| 基础层 | 格式 I/O（`fmt/`）、Phred 编码（`fq::qual`）、`PairReader`、k-mer、PAF、io/ds/loc/sys | pgr |
| 业务层 | `fq` 清洗/合并/纠错/归一/采样/分块、`asm` 组装（unitig/contig/OLC/map） | anchr |
| 命令壳 | `fq`/`asm` 的 clap 解析与执行 | anchr |
| 编排 | 模板渲染（`templates/`）、流水线顶层命令 | anchr |

### 1.2 与外部生态的对照

- `fq` 家族以 **BBTools 39.38** 为主参考，逐字节核对（golden 见
  `tests/bbtools/Lambda/golden/`，分析见 `notes/references/bbtools.md`）；
- `asm` 家族对照 BBTools tadpole（contig）、BCALM2/GATB（unitig）、
  bbmap `perfectmode`（map），OLC 三段参考 canu/celera/metaMDBG/skesa；
- 流程模板原使用 BBTools/dazzler 工具链，逐步替换为 anchr 自己的命令。

## 2. 架构全景

### 2.1 双层结构

- **`src/cmd/`**：命令壳。每个命令 `make_subcommand()`（clap 定义）+
  `execute()`（参数解析 → 调用 libs → 输出）。顶层命令平铺
  （`anchors`/`contained`/…），fq/asm 按子命令组组织（`cmd/fq/`、`cmd/asm/`）；
- **`src/libs/`**：业务逻辑（`asm/`、`olc/`、`fq/`、`map.rs`、`overlap.rs`），
  与 pgr 的 libs 风格一致（复杂逻辑进 libs，命令壳保持薄）。
- **`src/utils.rs`**：项目内共享工具（`read_fasta`、`write_lines`、`ucfirst` 等），
  `pub mod utils` 导出。

### 2.2 命令分发模式

`src/anchr.rs` 的 `main` 注册全部子命令：

```rust
env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
let app = Command::new("anchr")...;
match app.get_matches().subcommand() {
    Some(("fq", sub)) => cmd::fq::execute(sub)?,
    Some(("asm", sub)) => cmd::asm::execute(sub)?,
    Some(("overlap", sub)) => cmd::overlap::execute(sub)?,
    ...
}
```

`cmd/fq/mod.rs` / `cmd/asm/mod.rs` 各自聚合子命令并二次分发。

### 2.3 依赖策略

- **pgr**：git 依赖 + rev 锁定（发布形态可复现）；本地开发用
  `.cargo/config.toml` 的 `[patch."https://github.com/wang-q/pgr"]` 覆盖为
  `../pgr`（该文件被 gitignore，不随发布）；pgr 基础 API 变更时显式 bump rev。
- 其余依赖与 pgr 对齐（clap/itertools/petgraph/indexmap/cmd_lib/rayon/bstr/
  log/lru/env_logger/flate2 等），避免两边因外部 crate 版本出现行为差异。
- 曾用 `bio`/`intspan` 外部 crate，已移除：序列 I/O 用 `pgr::libs::fmt`，
  区间用 `pgr::libs::ds::IntSpan`，分层覆盖改用 `pgr::libs::runlist`
  （`covered` 命令及 `covered_benchmark` 已废弃，覆盖基准在 pgr 侧）。

### 2.4 构建配置

- `Cargo.toml`：`publish = false`；release `lto = true`；
  `[[bench]]` 2 个（`fq_assemble_benchmark`、`asm_map_benchmark`，
  criterion + rand dev 依赖）；
- `rust-toolchain.toml`：stable + rustfmt/clippy/rust-src；
- CI：`.github/workflows/`（build/codecov/publish；publish 用 zigbuild 交叉编译）。

## 3. 命令模块全景

### 3.1 顶层流程命令（Overlaps / Assembling）

| 分组 | 命令 | 说明 |
| :--- | :--- | :--- |
| Dependence | `dep check/install` | 依赖检查 |
| Download | `ena meta/manifest` | ENA 元数据抓取（JSON）/ 下载清单（tsv/ftp/md5/ascp） |
| Overlaps | `contained`/`dazzname` | 包含去冗余 / 重命名（8_*/2_insert_size 用） |
| Assembling | `trim`/`quorum`/`mergeread`/`template` | 流程原语与模板编排 |

### 3.2 `fq` 命令组（15 个，随迁自 pgr）

| 命令 | BBTools 对照 | 说明 |
| :--- | :--- | :--- |
| `to-fa` / `interleave` | — | 格式转换 / 双端交错 |
| `clump` | clumpify | k-mer 签名聚类、去重 |
| `split` / `sample` | repair / reformat | 拆分交错文件 / 降采样 |
| `clean` / `filter` | bbduk | 接头与质量清洗 / 参考 k-mer 过滤 |
| `s-filter` | quorum | 按自身 k-mer 计数过滤 |
| `merge` / `ec-overlap` | bbmerge | 双端 overlap 合并 / 纠错（ecco） |
| `ec-kmer` / `extend` | tadpole | k-mer 重装纠错 / 沿图延伸 |
| `norm` | bbnorm | 按 k-mer 深度过滤 |
| `range` | — | 按名字/区间提取（.loc 索引） |
| `trim-qual` | sickle/cutadapt | 质量修剪（滑窗/Mott） |

### 3.3 `asm` 命令组（7 个，随迁自 pgr）

| 命令 | 对照 | 说明 |
| :--- | :--- | :--- |
| `contig` / `unitig` | tadpole / BCALM2 | k-mer 图组装 |
| `multik` | metaMDBG | multi-k 迭代组装（unitig 图跨轮验证，无 N 核心方向，见 `design/asm-multik.md`） |
| `ovlp` / `layout` / `cns` | OLC 三段 | 精确 overlap → 布局 → 一致序列 |
| `olc` | — | 多 k unitig + OLC 组合流水线 |
| `map` | bbmap perfectmode | reads 回帖（完美匹配） |

## 4. 核心库层详解

### 4.1 `libs/asm/` — k-mer 图组装与 read 精修

- `refine.rs`：read 级 k-mer 图精修——质量门控计数表、局部重装纠错、
  保守延伸、junk/low-depth 过滤（最初移植自 BBTools `tadpole.sh`，后续
  已发展出 long-k 路径/打包表/流式计数，仅保留 CLI 兼容对照）；
- `assemble.rs`：`assemble`（contig）/`assemble_unitigs`（unitig）入口，
  种子深度阈值、分支处理。
- `multik.rs`：multi-k 迭代组装（unitig 图跨轮验证：桥接 k-mer 选边 +
  嵌合清理 + 渐进丰度过滤/主路径保护 + recompact；设计
  `design/asm-multik.md`）。

### 4.2 `libs/olc/` — OLC 三阶段

- `overlap.rs`：unitig 伪 reads 的精确后缀/前缀 overlap（radix sort + canonical k-mer）；
- `layout.rs`：链式布局（坐标、strand）；`consensus.rs`：布局缝合为 contig。

### 4.3 `libs/fq/` — reads 处理业务

- `trim.rs` / `trim_adapter.rs`：质量修剪（滑窗/Mott）与接头去除（bbduk 语义）；
- `merge.rs` / `bbnet.rs`：bbmerge overlap 合并 + CellNet 神经网络推理；
- `norm.rs` / `clump.rs`：bbnorm 归一与 clumpify 聚类（`temp_dir_for` 外部桶路径）；
- `sample.rs` / `split.rs` / `overlap.rs`：降采样、拆分、overlap 记录。

> 边界：`qual`（Phred 转换/检测）、`pairs`（PairReader）留在 pgr，anchr
> 一律 `pgr::libs::fq::qual/pairs`。

### 4.4 `libs/map.rs` — perfect-mode 回帖

种子-验证、无错配无缺口（`AbstractMapThread.java:1371` 语义），输出 SAM
（`--outm`）与 unitig 覆盖（`--outu`）。

## 5. 设计模式与约定

### 5.1 命令模式

每个命令模块两个公开函数：`make_subcommand() -> Command`（`.about` 第三人称、
`after_help` 用 `r###"..."###`）+ `execute(&ArgMatches) -> anyhow::Result<()>`。
共享 clap helper 集中在 `src/cmd/args.rs`（`outfile_arg`/`get_outfile`/
`ensure_outfile_distinct`/`infiles_arg_with_numargs`/`parse_parallel_auto` 等 16 个）。

### 5.2 零 Panic 策略

任何用户输入（畸形数据、二进制文件）不应 panic；`main` 的 match 用 `?`
传播错误（`Error: ...` + 非零退出码）。测试覆盖错误路径（`run_fail`）。

### 5.3 测试约定

- 集成测试 `tests/cli_*.rs`，用 `common::AnchrCmd`（模仿 pgr 的 `PgrCmd`：
  `cargo_bin("anchr")` + stdin/args/current_dir/run/assert/run_fail）；
- 依赖 pgr 命令的流程测试（如 `sam to-rg`、`rg coverage`）用
  `common::pgr_cmd()` 定位 pgr 二进制，找不到则跳过；
- 测试数据 `tests/bbtools/Lambda/`（golden 来自 BBTools 39.38，双轨逐字节对照）；
- 性能敏感改动先写 criterion 基准（`benches/`）。

### 5.4 流程命令约定

`contained`/`dazzname` 为自实现 overlap 工具（无外部依赖）；dazzler
流水线（overlap/orient/merge 及外部 LAshow/daligner）已废弃
（2026-08-15，现代 OLC 用 `asm ovlp/layout/cns`）。

## 6. 项目现状评估

### 6.1 已完成的（成熟）

- fq/asm 从 pgr 迁移完成（阶段 1-4）：25 命令 + 业务 libs + 21 测试文件
  + golden 数据；双轨核对 22/22 逐字节一致（`scripts/verify-migrate.sh`）；
- pgr 依赖规范化：git rev 锁定 + 本地 patch + `.cargo/config.toml`；
- `bio`/`intspan` 外部依赖移除，改用 pgr 基础层（含 Coverage→runlist 迁移，
  基准验证 sweep 路径快 2-30×）；
- 基础设施与 pgr 对齐：AGENTS.md、rust-toolchain、CI（zigbuild）、
  notes 结构（design/references/audit/benchmarks）、`docs/` 用户文档。

### 6.2 进行中的（活跃开发）

- 染色体级流程编排（templates 模板链的替换与验证）；
- `fq range` 对 BGZF 输入依赖 pgr 的 `.gzi` 索引生成（基础层机制，未封装 CLI）。

### 6.3 待补全的（TODO / 设计阶段）

- ~~迁移收尾~~（2026-08-15 销账：`fq-asm-migrate.md` 已标注"阶段 1-4 均已
  落地"（阶段 4 走"只写文档"路线，pgr 代码保留）、todo/索引已同步）；
- 既有 warning 清理（`Overlap` dead_code（bin/lib 双份 libs 结构问题）、
  `refine` 未读字段（2026-08-15 已删 2 处））；
- ~~部分流程命令（`dep`/`ena`/`template` 链）的外部工具版本核对~~
  （2026-08-15 完成：现代依赖清单见 §7.1，`check_dep.sh` 已更新，必需项
  全在位、已替代工具标注 legacy）。

### 6.4 不做 / 不适合做的

- 不重复基础层：格式 I/O、Phred、k-mer 表、PAF 解析一律走 pgr；
- 不引入新的外部算法库替代已迁移实现（BBTools/BCALM 语义已本地化）；
- 不内置 pgr 的比对/索引/遮蔽命令（那些留在 pgr）。
- 不搬 pgr 专属基准：`benches/` 33 个 .rs 与 `notes/benchmarks/` 12 篇
  均服务 pgr 剩余命令/基础层，留在 pgr；
- 不做并行 walk（`--parallel-walk` 已撤：多阶段线程池叠加会卡死系统，
  教训见 `design/asm-assemble.md` §12.2）。

## 7. 与周边项目的关系

- **pgr**：基础库提供者（git rev 锁定 + 本地 patch 双轨）；
- **BBTools 39.38/40.01**：fq/asm 语义与 golden 的主参考（源码在
  `BBTools-40.01/`，仅作参考，被 gitignore）；
- **bcalm/canu/celera/metaMDBG/megahit/skesa/quorum**：asm 与纠错的源码参考
  （笔记见 `notes/references/`）。

### 7.1 外部依赖清单（2026-08-15 版本核对，现代流程）

| 类别 | 工具 | 用途 | 备注 |
|---|---|---|---|
| 必需 | `anchr` / `pgr` | 流程主命令 / 基础库命令（`kmer hist`、`fa split/range` 等） | 自建（cargo） |
| 必需 | `ureq`（Rust crate） | `ena meta` 的 ENA portal API 客户端（替换 Perl LWP::Simple） | 2.12.1，`features=["json"]` |
| 必需 | `parallel` `jq` `pigz` | 模板编排、JSON env、压缩 | cbp |
| 必需 | `quast` | 最终组装质控（用户明确保留） | cbp |
| 可选 | `fastqc` | QC 对照（`fq qc` 主参考） | 已替代但保留对照 |
| 可选 | `spades` / `megahit` | 组装对照（用户保留参考） | 现代主路线 = multik+OLC |
| 可选 | `bwa` `samtools` `gatk` `mosdepth` | 变异检测模板（`3_*`） | 非组装主流程 |
| 已替代 | BBTools（bbduk/clumpify/bbnorm/reformat/repair） | trim 流水线 | → `fq clean/filter/clump/norm/sample/split`（P3 模板替换后不再调用） |
| 已替代 | `sickle` | 质量修剪 | → `fq trim-qual` |
| 已替代 | `kmercountexact`/`jellyfish` | k-mer 直方图 | → `pgr kmer hist` |
| 已替代 | `quorum` | reads 筛选 | → `fq s-filter` |
| 已替代 | `bcalm`/`bifrost` | unitigs | → `asm unitig`/`asm multik` |
| 已替代 | `picard` | insert size | → `asm map` + SAM TLEN |
| 已替代 | `masurca` | quorum 脚本 PATH 依赖 | quorum 不再用时不需要 |

本机 cbp 环境 2026-08-15 核对：必需项全部在位；`bifrost` 未装（现代流程
不需要）；dazzler 工具链（daligner/fasta2DB/LAshow/dazz）已随流水线废弃
不再需要。

## 8. 关键风险与技术债

1. **pgr rev 漂移**：pgr 基础 API 演进需显式 bump rev 并重跑
   `scripts/verify-migrate.sh`；本地 patch 掩盖了"发布版与本地版"差异，
   发布前需确认 patch 未影响构建；
2. **双轨遗留**：pgr 删除后部分笔记/索引仍引用旧路径（anchr 副本待同步）；
3. **外部工具依赖**：流程命令依赖外部工具（bwa/picard/quast 等），CI/
   容器环境需预装；测试中缺失时跳过（可能掩盖回归）；
4. **`fq range` BGZF 路径**：`.gzi` 索引生成未封装 CLI，属基础层缺口；
5. **既有 warning**：8 处 dead_code 等（迁移代码带入 + 原有），clippy
   `-D warnings` 前需清理。

## 9. 主题链路索引（按技术线，跨目录）

- **组装流程链**：`fq`（trim/ec/merge/norm）→ `asm unitig/contig` →
  `asm ovlp/layout/cns/olc` → `asm map` → `template` 编排；
  文档：`docs/fq.md`、`docs/asm.md`、`notes/design/asm-assemble.md`、
  `notes/design/asm-olc.md`、`notes/design/asm-map.md`；
- **BBTools 替换链**：`fq clump/split/sample/clean/filter/merge/ec-*/norm`
  对照 39.38 golden；文档：`notes/design/fq-trim-replace.md`、
  `notes/design/fq-merge-replace.md`、`notes/references/bbtools.md`；
- **Overlap 工具**：`contained`/`dazzname`（8_*/2_insert_size）；`covered`
  命令已废弃（pgr runlist 替代）；基准 `benches/covered_benchmark.rs`。

## 10. 设计笔记索引（notes/design/）

| 文档 | 内容 |
| :--- | :--- |
| `fq-asm-migrate.md` | 迁移方案档案（阶段 1-4 + 批次 + 核对清单） |
| `fq-trim-replace.md` | trim 流水线 BBTools 替换（8 步 + 入口映射） |
| `fq-merge-replace.md` | merge/ec 系列 BBTools 替换 |
| `asm-assemble.md` | contig/unitig/olc 组装设计 |
| `asm-map.md` | perfect-mode map 移植 |
| `fq-range.md` | FASTQ `.loc` 索引（range 命令） |
| `asm-olc.md` | OLC 三段设计 |
| `asm-multik.md` | multi-k 迭代组装（unitig 图跨轮验证，借鉴 metaMDBG；无 N 染色体核心方向；§9 防 misassembly（bridge_filter/split_by_bridge）、§10 metaMDBG 对比，2026-08-15 并入） |
| `asm-olc-modern-flow.md` | 现代组装流程总结（2026-08-15 会话收官：fq→multik→anchor→olc--unitigs→quast；完整总结与用户裁定） |
| `fq-validation.md` | fq 前处理验证计划（2026-08-15：与老流程 reads 准备阶段逐一对账 + P0-P5 行动计划） |
| `qc.md` | 自有 QC 方案设计（FastQC/Falco 双参考，M1-M4 里程碑） |

## 11. 外部工具参考索引（notes/references/）

| 文档 | 服务对象 |
| :--- | :--- |
| `bbtools.md` | fq/asm 主参考（tadpole/bbduk/bbmap/clumpify 等） |
| `bcalm.md` | `asm unitig` 移植来源 |
| `canu.md` / `celera.md` | `asm olc` 参考 |
| `metaMDBG.md` / `skesa.md` / `megahit.md` | multi-k 迭代（metaMDBG 首要借鉴，见 `asm-multik.md`）/ 保守扩展 + 多 k 迭代本体与 Rust 移植（含与 multik 的对比，见 skesa.md §7.2）/ 图清洗参考 |
| `cutadapt.md` / `sickle.md` | `fq trim-qual` 算法来源 |
| `fairy.md` | `fq norm` 大数据方案调研 |
| `quorum.md` | read 纠错参考（`fq ec-kmer`） |
| `fastqc.md` | reads 质控参考（`2_fastqc.tera.sh`，迁移候选） |
| `falco.md` | FastQC 的 C++ 仿制实现（QC 方案第二参考，输出兼容） |
| `anchr-legacy-pipeline.md` | 老手工流程理解（G37 例子：多覆盖度拆分 → 各部分 unitigs/anchors → merge → OLC 拼装（glue/fill）→ spades/megahit 对照；2026-08-15 修正） |

## 12. 笔记根 / audit / benchmarks 索引

- `notes/audit/audit-fq.md`：fq 命令族审计记录（BBTools 对照结论）；
- `notes/audit/audit-asm.md`：asm 命令族审计记录；
- `notes/benchmarks/bbtools-vs-anchr.md`：fq 命令 vs BBTools 39.38 CLI 基准
  （随迁自 pgr 并适配，`kmer hist` 行仍属 pgr）；
- `notes/benchmarks/qc-bench.md`：`fq qc` 端到端基准（anchr vs FastQC vs
  Falco，hyperfine；fastqc 1.1s 启动主导、falco 最快、anchr 慢 3× 且
  绝对值小，多线程小数据负优化）；
- `notes/benchmarks/multik.md`：`asm multik` 吞吐 sanity check（1 Mb 合成
  基因组 → 单条 100%，9.4 s / 816 MB）；
- `notes/benchmarks/README.md`：benchmarks 目录索引；
- `notes/todo.md`：待办清单（仅 actionable 项，历史结论与细节见各
  design/benchmark 文档，会话交接见 `design/asm-olc-modern-flow.md`）；
- 本文档为项目理解与索引入口；`AGENTS.md` 为行为准则。
