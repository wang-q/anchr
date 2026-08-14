# asm 命令族代码审核记录（2026-08-13）

> 维护约定：本文是迁移时的审计快照；后续 anchr 侧对 asm 的修改应更新本
> 审计记录（而非 pgr 侧）。

对 `anchr asm` 全部 7 个子命令（contig/unitig/ovlp/layout/cns/olc/map）及相关库
文件（`libs/asm/assemble.rs`、`libs/asm/tadpole.rs`、`libs/olc/overlap.rs`、
`libs/olc/layout.rs`、`libs/olc/consensus.rs`、`libs/map.rs`，以及
`cmd/asm/common.rs`）和全部测试/文档进行审核。以下仅保留有借鉴意义的结论；
逐轮验证过程已精简（第一轮发现并修复 5 类问题，第二、三轮复核未再发现新问题，
第四轮为迁移后代码路径复核 + 文档准确性修正）。

> 注：本报告初版（2026-08-12）审核的是迁移前的 `src/cmd_pgr/asm/`；`asm` 已迁移至
> `src/cmd/asm/`，本版路径已全部更新为迁移后的位置。

| 子命令 | 命令文件 | 主要库 |
| ------ | -------- | ------ |
| contig | `src/cmd/asm/contig.rs` | `libs/asm/assemble.rs`, `libs/asm/tadpole.rs` |
| unitig | `src/cmd/asm/unitig.rs` | `libs/asm/assemble.rs`, `libs/asm/tadpole.rs` |
| ovlp   | `src/cmd/asm/ovlp.rs`   | `libs/olc/overlap.rs` |
| layout | `src/cmd/asm/layout.rs` | `libs/olc/layout.rs` |
| cns    | `src/cmd/asm/cns.rs`    | `libs/olc/consensus.rs` |
| olc    | `src/cmd/asm/olc.rs`    | 上述全部 |
| map    | `src/cmd/asm/map.rs`    | `libs/map.rs` |

共享：`src/cmd/asm/common.rs`。

## 与外部参考实现的语义一致性核对

asm 家族对照 BBTools tadpole（contig）、BCALM2/GATB `ograph.cpp graph3`
（unitig）、BBTools bbmap `perfectmode`（map）语义，OLC 三段各自对照参考流程，
并逐字节/逐行对拍验证。已确认的一致/有意行为：
- `contig`：k-mer 计数质量门控、种子深度阈值、双向贪心扩展、泡泡消除，复现
  tadpole 默认 contig 模式（`popbubbles=t`）。
- `unitig`：仅当唯一后继且其前驱唯一时扩展，停于分支/接点/覆盖缺口/环，与 k
  扫描顺序无关（BCALM2 压缩语义）。
- `map`：种子-验证、完美匹配（无错配无缺口）、`ambiguous=all`（单端）、
  `--paired` 携带配对标志/伴侣坐标/TLEN，作为 bbmap `perfectmode maxindel=0`
  的替代。
- OLC：unitig → ovlp（精确重叠）→ layout（贪心互惠接点）→ cns（精确重叠拼接）
  三段闭环。
- 确定性与复现：种子扫描按规范化 k-mer 排序快照线性迭代（不依赖哈希表内存序），
  多线程输出保序。

## 排除的疑点（安全不变量，经核验无需修复）

- `find_overlaps`：seed_k 先 `min(最短 unitig)` 再 `ensure!(seed_k >= 1)`，
  空/超短 unitig 先归 0 被友好错误拦截，不会在 `query_unitig` 的
  `seq[0..seed_k]` 越界 panic；seed_k 已钳到最短长度，所有 unitig 窗口切片有效。
- `verify_seed`/`extend`：四个方向扩展均有边界守卫，`q[p..p+k]`/`t[tpos..tpos+k]`
  不越界；`query_unitig` 的 tid/tpos 来自排序后 payload，索引安全。
- `build_layouts`/`extend`：`placed` 保证不死循环；`2*u+end` 索引在
  `qid/tid < n` 保证下安全；`is_repeat` 的 `length*10` 在 64 位 usize 下不溢出。
- `consensus`/`dedup_contained`：`overlap <= seq.len() && overlap <= piece.len()`
  防切片越界；`contains` 对空 needle 有守卫。
- `make_contig`：`bb.len() == k` 守卫在种子处即拒绝 contig；扩展循环
  （`extend_to_right`/`right_counts_of`/`left_counts_of`/`calc_coverage`）与
  tadpole 移植一致；`resolved_min_contig_len` 的 `2*k` 因 k ≤ 128 不溢出。
- `ovlp` 的 `--overlap-k 0`/`--min-overlap 0`：`--overlap-k 0` 使 `seed_k` 钳为 0，
  `ensure!(seed_k >= 1)` 报友好错误（帮助范围是 `1..=128`，不 panic）；`min_overlap=0`
  仅放行所有 found overlap，亦不 panic。
- `map`：单端/配对统计、SAM FLAG/TLEN/RNEXT 计算正确；`--max-reads` 单端跨文件
  与配对成对语义均正确。
- `common.rs`：`to_paf`/`read_unitigs`/`format_cov` 边界安全。

## 记录项（未改，低风险 / 待决策）

- `contig`/`unitig` 对**单个文件**按"交错配对"解释（1 交错或 2 配对），单端 reads
  需用 2 个文件或显式交错文件；奇数条 reads 时报 `unpaired trailing read`。与文档
  一致，属设计语义，未改。
- `map` 配对模式仅报告两端**首个**命中（`ambiguous=all` 语义仅单端成立），属配对
  模式固有约束，未改。

## 已知限制（有意保留）

- `map --max-reads` 奇数时按"成对计两条"取整到偶数（配对原子性）。
- `contig`/`unitig` 的 `-p/--parallel` 仅为兼容而接受，实际单遍确定性处理
  （文档已注明）。

## 修复的缺陷（根因模式）

### Zero-Panic / 溢出

- **`-k` 超上限 panic**（contig/unitig/olc）：只校验 `k >= 1` 未校验上界，
  `-k 129`（> `key::Kmer::MAX_K = 128`）在 `count_read_kmers` → `Kmer::new(k)
  .expect()` 处 panic。修复：`assemble()` 与 `assemble_unitigs_core()` 加
  `ensure!(k <= MAX_K)` 友好报错。回归 `command_asm_contig_rejects_kmer_above_
  limit` / `command_asm_unitig_rejects_kmer_above_limit`。
- **`cns` 布局 `contig_0` 下溢 panic**：`parse_index(...)? - 1` 对 1-based 布局
  id 直接减一，`contig_0` 使 usize 下溢（debug panic / release 环绕成
  `usize::MAX` 触发 OOM）。修复：先校验 `>= 1` 再减一。回归
  `command_asm_cns_rejects_contig_zero`。

### 数据安全（`-o` 同输入保护 / 输出截断）

- **`contig`/`unitig` 的 `-o` 可覆盖输入**：writer 在读取 reads 前打开，`-o`
  指向输入会先截断 reads。修复：打开 writer 前 `ensure_outfile_distinct`。回归
  `command_asm_contig_outfile_not_input` / `command_asm_unitig_outfile_not_input`。
- **`map` 的 `--outm`/`--outu` 可覆盖输入**：SAM writer 在流式读取 reads
  **之前**打开并截断，`--outm`/`--outu` 指向 reads 文件会丢数据。修复：对
  `outm`/`outu` 与 ref+reads 全部 `ensure_outfile_distinct`。回归
  `command_asm_map_outm_not_input`。
- **一致性补充**：`ovlp`/`layout`/`cns`/`olc` 四个写 `-o` 的子命令虽在打开
  writer 前已把输入读入内存（无中途丢失），为对齐项目统一约定也补上
  `ensure_outfile_distinct`（分别对 unitig FASTA / PAF+FASTA / TSV+FASTA /
  reads）。

### 文档一致性 / 帮助文本

- **`-k` 上界描述错误**：`docs/asm.md` contig/unitig 段写 "no upper bound"，
  map 段与 `map` help 写 "1..=64"，与代码实际上限 128 不符。修复：contig/unitig
  段改 "up to 128, the k-mer key table limit — k > 64 uses multi-word k-mers"；
  map 段与 help 改 "1..=128"。

## 实机边界验证（均无 panic，返回友好错误或正常退出）

- `contig` reads 长度 < k → 正常退出，0 contigs
- `unitig` 空输入文件 → 正常退出，0 unitigs
- `ovlp` 空 unitigs → 友好错误 `cannot overlap: unitigs are empty...`
- `map -k 129`（超上限） → 友好错误 `k-mer length must be in 1..=128`
- `map --paired` 两端 reads 数目不等 → 友好错误 `paired files have different read counts`
- `map` 空 reads → 正常退出，0 mapped

## 第四轮复核（迁移后代码路径 + 文档准确性）

对迁移后位于 `src/cmd/asm/`、`src/libs/asm/`、`src/libs/olc/`、`src/libs/map.rs`
的现行代码做了完整的纵深复核：

- 逐文件重读了全部 7 个命令文件与 4 个库文件，未发现新的逻辑缺陷；此前的
  Zero-Panic、数据安全（`-o` 同输入保护）、`-k` 上界校验等修复在迁移后代码中
  均完整保留并继续生效。
- `map` 的 `--parallel` 文档与代码一致（`1..=1024`，默认 8；`parallel_arg_with_default`
  与 `docs/asm.md` 同为 1024 上界），无文档-代码偏差。
- 复核了 `olc.rs` 的 `--keep-dir` 语义，发现**文档准确性缺陷**：原帮助/文档称
  中间文件可用于 "re-running the stage commands separately"，但 `--keep-dir`
  落盘的 `unitigs.fa`/`ovlp.paf`/`layout.tsv` 名称（`k<k>:unitig_<id>`，无
  `stem:` 前缀）与独立 `ovlp`/`layout`/`cns` 通过 `read_unitigs` 派生的
  `stem:name` 前缀不一致，直接回灌会报 "PAF query ... not found in unitigs"。
  修复：`olc.rs` 帮助与 `docs/asm.md` 改为明确仅用于调试/检查，并注明不可直接
  经独立阶段命令回灌。

第四轮结论：迁移后代码逻辑复核通过（无新缺陷）；仅修正 `--keep-dir` 的文档表述。

## 第五轮复核（`map` 输出数据安全）

- 复核 `map` 的 `--outm`/`--outu` 输出安全时发现一个**数据安全缺陷**：`--outm`
  与 `--outu` 指向同一路径时未加防护。`map_files_inner` 会为两者各开一个 writer，
  后开者截断先写者已写入的头部/记录，静默产生损坏的 SAM（与既有 `-o`/`--outm`
  同输入防护属同一数据安全类别，但此前只校验了输出 vs 输入，未校验两输出之间）。
  修复：在 `cmd/asm/map.rs` 的 `execute` 用 `pgr::libs::io::same_path` 校验
  `--outm`/`--outu` 不同，回归 `command_asm_map_outm_neq_outu`。

## 第六轮复核（`ovlp`/`olc` 的 `--overlap-k` 上限）

- 复核 `find_overlaps` 的 seed 钳制时发现一个**静默错误结果缺陷**：`seed_k` 只
  钳到最短 unitig 长度并校验 `>= 1`，却未校验 `Kmer::MAX_K` 上界。当最短 unitig
  > 128 bp 且用户传 `--overlap-k > 128` 时，`seed_k` 落在 `(128, 最短长度]`，
  `canonical_keys` 因 `k > MAX_K` 直接返回空（不 panic），导致**静默输出 0 条
  重叠**而非友好报错——与 `contig`/`unitig`/`map` 对 `-k > 128` 一律友好报错的
  约定不一致。修复：在 `libs/olc/overlap.rs::find_overlaps` 钳制后追加
  `ensure!(seed_k <= Kmer::MAX_K)` 友好报错（同时覆盖 `ovlp` 与 `olc` 两个入口），
  并把 `--overlap-k` 的帮助文本/`docs/asm.md` 注明 `1..=128`。回归
  `command_asm_ovlp_rejects_overlap_k_above_limit`。

## 第七轮复核（文档一致性 + `cns` 布局解析健壮性）

- **`map` 帮助文本误用复数**：`cmd/asm/map.rs` 的 `ref` 参数 help 写 "Reference
  FASTA file(s)"，暗示可传多个引用文件，但该参数为单值（`num_args` 默认 1），且
  `after_help` 与 `docs/asm.md`（`<ref.fa>` 单数）均明确引用是单个。改为
  "Reference FASTA file to map against"。
- **`contig`/`unitig` 的 `-k` 帮助未标上限**：家族内 `map`/`ovlp`/`olc` 的 `-k`/`--overlap-k`
  帮助都注明 `1..=128`，而 `contig`/`unitig` 只写 "K-mer length"，与 `docs/asm.md`
  明确的上限 128 不一致。补为 "K-mer length (1..=128)"。
- **`cns::parse_layouts` 畸形超大 `contig_N` 触发巨型分配 abort**：原逻辑
  `layouts.resize(ci + 1, ...)` 对畸形布局行（如 `contig_999999999`）会一次性扩容到
  巨大容量，导致容量溢出/分配失败（进程 abort），违反 Zero-Panic 约定。修复：改为
  强制 contig id 连续（首步 `si==0` 时要求 `ci == layouts.len()`，后续步要求
  `ci < layouts.len()`），既消除了巨型分配，也移除了 `resize`。回归
  `command_asm_cns_rejects_noncontiguous_contig_id`（同时保留 `contig_0` 与正常
  拼接测试）。

第七轮结论：`cargo fmt`/`clippy` 干净（无新增告警），`cargo test asm` 全部通过。

## 第八轮复核（报告准确性自检）

- 复核本报告"排除的疑点"一节对 `ovlp --overlap-k 0` 的描述时发现**报告描述与代码
  不符**：报告称 seed_k 经 `.max(1)` 归 1，而现行 `find_overlaps` 是
  `seed_k = opts.seed_k.min(最短长度)` 后 `ensure!(seed_k >= 1)` 报友好错误（帮助
  范围 `1..=128`，0 属越界，报错不 panic）。已修正该条目为与代码一致的行为描述。
- 第八轮另对 map 的单端/配对命中、`cns` 连续 id 约束、`common.rs` 标签去重与
  `format_cov` 边界等做了复核，未发现新缺陷。

第八轮结论：仅修正报告自身的一条行为描述不准确，无代码/文档缺陷，审核收敛。

## 第九轮复核（`--links` 悬空引用 + 基准编译）

本轮对 `unitig` 的零悬空策略做了 GFA / FASTA 两条输出路径的对称复核，并全量编译
`--all-targets`（含 benches），发现两处新缺陷：

- **`unitig --links`（bcalm FASTA 模式）悬空引用缺陷**：上一轮为 GFA 输出加
  零悬空策略（`L` 边不引用被 `--min-contig-len` 丢弃的 segment），但同一循环里
  FASTA `--links` 分支仍把 `links[i]` 原样写入 `L:+:<to>:<ori>` 头部条目，未过滤
  被丢弃的 unitig id，导致头部引用不存在的 `unitig_<id>`——与 GFA 属同一类
  悬空引用缺陷，只因 `command_asm_unitig_links_header` 测试用 `--min-contig-len 1`
  （不触发丢弃）而未被覆盖。修复：在 `libs/asm/assemble.rs` 的 FASTA `--links`
  分支用与 GFA 相同的 `kept` 过滤 `links[i]` 后再写入头部。回归
  `command_asm_unitig_links_no_dangling`（bubble 数据 + 默认 `--min-contig-len`，
  校验每个 `L:` 引用的 id 都存在于输出头部；已证实在修复前失败）。
- **`benches/asm_map_benchmark.rs` 编译错误**：`MapOptions` 初始化缺 `parallel`
  字段（`libs/map.rs` 的 `MapOptions` 新增了 `parallel` 而基准未补），`cargo
  check --all-targets` 报 E0063。修复：补 `parallel: 8`（与 `map` 命令默认一致）。
  编译通过。

第九轮结论：`cargo fmt` 干净，`cargo clippy --all-targets` 无新增告警（仅存的
`bases_in` 字段未读等为既有死代码告警，非本轮引入），`cargo test` asm 全部通过
（cli_asm 11 + cli_asm_contig 9 + cli_asm_map 11 + cli_asm_olc 12 + cli_asm_unitig 11）。

## 第十轮复核（零悬空策略闭环 + 全量回归）

本轮对第九轮的两处修复做收敛复核，并对整个 asm 命令族做再扫描，未发现新缺陷：

- **零悬空策略闭环确认**：`emit_links`/`emit_gfa` 仅在 `unitig` 命令暴露
  （`cmd/asm/unitig.rs`），且两条输出分支（GFA `L` 边、FASTA `L:` 头部）均已用
  同一 `kept` 过滤，无遗漏的悬空引用路径。`olc` 管线用 `min_contig_len: 0`
  （自动 max(124, 2k)）做中间 unitig 过滤，overlap/layout/consensus 都在该已过滤
  集合上运行，内部一致、无悬空风险。
- **文档一致性**：`docs/asm.md` 中 `-k`/`--overlap-k` 的 `1..=128` 边界与
  `contig`/`unitig`/`map`/`ovlp`/`olc` 帮助文本一致，无残留旧边界措辞。
- **并行边界**：`map` 用 `parallel_arg_with_default("8")`（`value_parser`
  1..=1024 强制）；`contig`/`unitig` 用 `parse_parallel_auto` 校验（auto/1..=1024）。
  均无越界风险。
- **边界/溢出复核**：`overlap.rs` 空 unitig 集合时 `seed_k` 归 0 → `ensure` 友好报错
  不 panic；`consensus` 除零（seq.len()>=k）不可能；`layout` 坐标 `overlap_len<=prev_end`
  防下溢；`tadpole.rs` 的 `min_prob` 滑动乘积对 >127 质量值做了 `min(127)` 钳制，
  `extension_rollback` 用 `saturating_add`/`min` 防溢出。

第十轮结论：经纵深复核（含第九轮修复的回归验证、全量 `cargo check --all-targets`
编译与 54 个 asm 用例通过），未发现新的代码/文档缺陷，零悬空策略闭环、边界与
溢出不变量健全，审核收敛。

## 结论

`asm` 命令族审核完成（累计修复 5 类问题 + 4 处 `-o` 防护统一 + 1 处 `--keep-dir`
文档修正 + 1 处 `--outm`/`--outu` 冲突防护 + 1 处 `--overlap-k` 上限校验 + 第七轮
3 处：`map` 帮助单复数、`contig`/`unitig` 的 `-k` 上限帮助、`cns` 布局 id 连续性与
巨型分配防护 + 第九轮 2 处：`unitig --links` 悬空引用、`asm_map_benchmark.rs` 编译
错误），经纵深复核收敛；与 BBTools/BCALM
语义对拍、边界输入验证零 panic，`cargo fmt`/`clippy` 干净（asm 相关无新增告警），
相关集成测试与 `cargo test --lib` 全部通过。
