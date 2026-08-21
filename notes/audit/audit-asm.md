# asm 命令族代码审核记录（2026-08-13）

> 维护约定：本文是迁移时的审计快照；后续 anchr 侧对 asm 的修改应更新本
> 审计记录（而非 pgr 侧）。

对 `anchr asm` 全部 10 个子命令（anchor/contig/unitig/extend/ovlp/layout/cns/multik/olc/map）及相关库
文件（`libs/asm/assemble/`、`libs/asm/table.rs`、`libs/asm/refine.rs`、`libs/asm/extend.rs`、
`libs/asm/multik/`、`libs/olc/overlap.rs`、`libs/olc/layout.rs`、`libs/olc/consensus.rs`、
`libs/olc/anchor.rs`、`libs/map.rs`，以及 `cmd/asm/common.rs`）和全部测试/文档进行审核。
以下仅保留有借鉴意义的结论；逐轮验证过程已精简（第一轮发现并修复 5 类问题，第二、三轮复核未再发现新问题，
第四轮为迁移后代码路径复核 + 文档准确性修正；第十一、十二轮为新增子命令 anchor/extend/multik
的纵深审核与文档同步）。

> 注：本报告初版（2026-08-12）审核的是迁移前的 `src/cmd_pgr/asm/`；`asm` 已迁移至
> `src/cmd/asm/`，本版路径已全部更新为迁移后的位置。

| 子命令 | 命令文件 | 主要库 |
| ------ | -------- | ------ |
| contig | `src/cmd/asm/contig.rs` | `libs/asm/assemble/contig.rs`, `libs/asm/assemble/bubble.rs`, `libs/asm/table.rs`, `libs/asm/refine.rs` |
| unitig | `src/cmd/asm/unitig.rs` | `libs/asm/assemble/unitig.rs`, `libs/asm/table.rs`, `libs/asm/refine.rs` |
| ovlp   | `src/cmd/asm/ovlp.rs`   | `libs/olc/overlap.rs` |
| layout | `src/cmd/asm/layout.rs` | `libs/olc/layout.rs` |
| cns    | `src/cmd/asm/cns.rs`    | `libs/olc/consensus.rs` |
| olc    | `src/cmd/asm/olc.rs`    | 上述全部 |
| map    | `src/cmd/asm/map.rs`    | `libs/map.rs` |
| anchor | `src/cmd/asm/anchor.rs` | `libs/olc/anchor.rs`, `libs/map.rs` |
| extend | `src/cmd/asm/extend.rs` | `libs/asm/extend.rs` |
| multik | `src/cmd/asm/multik.rs` | `libs/asm/multik/` |

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
  （不触发丢弃）而未被覆盖。修复：在 `libs/asm/assemble/unitig.rs` 的 FASTA `--links`
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

## 第十一轮复核（新增子命令 anchor/extend/multik 纵深审核）

对 `asm` 迁移后新增的三个子命令（anchor/extend/multik）及其库做了逐行纵深审核，发现并修复以下缺陷：

- **`extend` 的 `-o` 数据安全缺陷**：`cmd/asm/extend.rs` 只对 `-o` 校验了 contigs 文件
  （`ensure_outfile_distinct(outfile, [contigs])`），未校验 reads 输入文件。由于 writer
  在 walk 前打开并截断目标，`-o` 指向任一 reads 文件会直接销毁用户输入。修复：校验
  扩展到 `contigs + 全部 infiles`，并前移到读取输入之前（fail-fast）。回归
  `command_asm_extend_outfile_not_reads`（校验失败后 reads 文件内容保持不变）。
- **`anchor --stats` 冲突校验延迟**：`cmd/asm/anchor.rs` 的 `--stats` 与 `-o`/输入的
  冲突校验原先在 `map_files` + SAM 解析之后才执行，无效路径会白跑一轮映射。修复：移到
  映射前，与 `-o` 的 `ensure_outfile_distinct` 并列执行。
- **`docs/asm.md` 缺失 anchor/extend 子命令**：`cmd/asm/mod.rs` 注册 10 个子命令，而
  `docs/asm.md` 只列了 8 个（缺 anchor/extend）。修复：子命令列表补两行，并新增完整的
  `## anchor`、`## extend` 文档段（参数、语义、示例，与帮助文本一致）。
- **`docs/asm.md` k-mer 上限过时（128→256）**：pgr `key::Kmer::MAX_K` 已升至 256
  （对标 megahit，2026-08-16），全部 asm CLI help 也已同步为 `1..=256`（contig/unitig/
  map/ovlp/olc/multik），但 `docs/asm.md` 仍写 `up to 128`/`1..=128`。修复：统一为
  256。注：multik 的 `auto_ks` 阶梯独立封顶 192（读长驱动、残留错误经验值），不受此影响，
  文档本已与 `auto_ks` 代码一致。
- **`extend` 的 `-k` 帮助未标范围**：家族内其他命令的 `-k`/`--overlap-k` 帮助都注明
  `1..=256`，`extend` 只写 "(default 31)"，而库层 `extend_contigs` 要求 `k in 2..=MAX_K`。
  补为 "(default 31; 2..=256, the k-mer key limit)"。
- **`docs/asm.md` multik 段缺 `--merge-similar`/`--merge-len`**：命令定义了这两个泡泡
  合并选项且做范围校验，但 multik 文档的 Options 列表漏列。修复：补两行（参数与默认值，
  与帮助文本一致）。
- **`multik --print-ks` 忽略 `--list-files`**：`--print-ks` 分支原先在解析 infiles
  之前判断 `is_list`（默认 false），传 `--list-files` 时把列表文件当普通 FASTA 读，
  导致报错。修复：`execute` 开头即读 `list_files` flag，`--print-ks` 与主路径共用
  `resolve_paths(f, is_list)`。回归 `command_asm_multik_print_ks_list_files`。
- **`multik`/`anchor` 的 `--parallel` 无范围校验**：`multik` 用
  `RangedU64ValueParser 0..=1024`（0 = 全部核心），`anchor` 用 `1..=1024`，防越界值
  破坏 rayon 线程池构建。

复核确认无需修改的路径：

- `multik`：`--kmer` 解析在库层 `assemble_multik` 校验 `1..=256` 并对 `--list-files`
  空解析报友好错误；`--print-ks` 已尊重 `--list-files`；`--merge-similar`/`--merge-len`/
  `--parallel` 均已在前置校验。k 序列内部 `sort + dedup`，`auto_ks` 空 reads 返回空序列
  报友好错误，均无 panic。
- `anchor`：空映射（read_len=0 / 全 0 覆盖）时 `anchor_stats` 返回全 0 边界，
  `anchor_regions` 的 `half=0` 走非边缘分支、无除零（`read_len` 作除数仅在 `half > 0`
  时进入），`extract_anchors` 切片边界由 regions 生成规则保证，无越界。
- `extend`：`extend_contigs` 校验 `k in 2..=MAX_K`；`min_len`/`min_extend` 语义与文档
  一致（短 contig 原样通过）；`cross_contig_kmers`/`is_foreign` 的跨 contig 护栏
  （`MIN_FOREIGN_RUN=5` 回滚）与低覆盖 seam 护栏（`MIN_LOW_RUN=5` 回滚）均有单元测试
  覆盖，无 panic 路径。

第十一轮结论：`cargo fmt`/`clippy --all-targets` 干净，extend/anchor/multik 集成测试
（4+4+10）与全部 asm 集成测试、107 个 lib 测试通过。

## 第十二轮复核（第十一轮修复的复核 + 未覆盖路径）

对第十一轮的修复逐项复核，并补扫此前未覆盖的路径，发现并修复 3 处新缺陷：

- **`bridge_filter` 文档与行为不一致**：`libs/asm/multik/bridge.rs` 的 `bridge_filter`
  注释称"短 unitig 的 link 保守保留"，而实现实际是**丢弃**无法构建 junction 探针的
  link（短 unitig 探针构建失败即删除）。注释与代码矛盾会误导后续维护者。修复：注释改
  为与实现一致的行为描述（未经验证的 junction 不压缩，碎片保持独立输出）。
- **`extend --min-support 0` 静默错误扩展**：`extend_contigs` 未校验 `--min-support`
  下界，`0` 时扩展计数全 0 仍能通过 `>= 2x runner-up` 多数判定，导致无 read 支持地
  追加碱基（静默扩展出无支撑序列）。修复：库层 `ensure!(opts.min_support >= 1)` 友好
  报错。回归 `command_asm_extend_rejects_zero_min_support`。
- **调试 trace 函数无效目录 panic**：`master.rs` 的 `trace_chains`/`trace_graph` 在
  `ANCHR_MULTIK_TRACE_DIR` 指向无效/不可写路径时 `File::create(...).unwrap()` panic。
  trace 是调试辅助，不应因环境变量错误而中断主流程。修复：改用 `if let Ok(f) = ...`
  静默跳过。复核时确认：`trace_graph` 的 `links[i]` 索引与 `unitigs` 一一对应，循环内
  `i` 均 < `unitigs.len()`，无越界。

第十二轮结论：`cargo fmt`/`clippy --all-targets` 干净，extend/anchor/multik 集成测试
（5+4+12）与全部 asm 集成测试通过。

## 第十三轮复核（multik 内部 + olc 驱动 + 文档一致性纵深）

对 multik 调度/图/桥/精化内部与 olc 驱动做纵深复核，发现并修复 3 处文档-实现不一致与
1 处默认值来源不统一：

- **`--supermer-m` 默认 minimizer 长度与 pgr 不一致（两处内联）**：
  `cmd/asm/unitig.rs` 的默认值计算 `(12).min((5).max(k / 4))` 与
  `libs/asm/multik/schedule.rs::pass0_opts` 的内联式重复了 pgr
  `supermer::minimizer_len`，但缺 `k-1` 上界、且用截断除法 `k / 4` 而非
  `ceil(k/4)`（k%4≠0 时差 1），与 pgr 单点事实来源漂移。修复：两处统一调用
  `pgr::libs::kmer::supermer::minimizer_len(k)`，同时修正 `--supermer-m` 帮助文本
  为 `min(12, max(5, ceil(k/4)), k-1)`。
- **`docs/asm.md` 的 `--supermer` 描述缺 `k-1` 上界**：帮助/文档与实现不一致，
  补为 `min(12, max(5, ceil(k/4)), k-1)`。
- **`progressive_filter` 注释把 cutoff 上限写成 "graph maximum"**：
  `libs/asm/multik/graph.rs` 的实现实际是 **25% 的 unitig 覆盖度中位数**（metaMDBG
  语义；重复区会把图最大丰度顶到几百 x，若按最大会误删主链）。修正注释与实现一致。
- **`docs/asm.md` cns 输出头格式笔误**：`>contig_<id>,len=...` 逗号应为空格
  `>contig_<id> len=...,cov=...`，与实现输出对齐。

第十三轮结论：`cargo fmt`/`clippy` 干净；`--supermer` 各 k 输出逐字节不变（minimizer
默认值对齐 pgr 后 G37 full 验证一致），全部 asm 测试通过。

## 第十四轮复核（olc 内部纵深 + 输入校验闭环）

对 OLC 的 consensus/coverage 数值内核与 unitig/contig 输入校验做纵深复核，发现并修复
3 处问题：

- **`identity`（banded Levenshtein）首碱基漏计 bug**：`i==0 && j==0` 的基例写死
  `best = 0`，不比较首碱基，导致**首碱基不匹配**的近重复 contig 编辑距离少算 1、
  identity 虚高（如 `AC` vs `GC` 被算成 1.0 而非 0.5），可能把边界差异的近重复
  误判为 dedup。修复：`best = usize::from(ai != bj)`。新增回归测试
  `identity_counts_first_base_mismatch`（首碱基失配计 1；单碱基对失配 identity=0）、
  `identity_scores_substitutions`（纯替换打分 `1 - mm/len`）。
- **`coverage` 的 `dominant_offset` 与 `overlap_geometry` 重复实现**：两处都是
  "31-mer 索引 + offset 直方图"，但 `dominant_offset` 缺 `overlap_geometry` 的
  `EXACT_WORK_CAP` 有界回退——homopolymer 富集 contig 上会二次方发散（20k×20k 对）。
  修复：`coverage` 复用 `overlap_geometry`（含 heavy/light 分治 + 精确验证回退），
  删除 `dominant_offset`；`bounded_geometry_matches_exact_histogram` 与
  `identity_*` 等测试全绿。
- **`--min-count-seed 0` 未校验（unitig/contig）**：`0` 会把所有 k-mer 视为 solid、
  抹掉错误过滤，静默改变组装语义。修复：`assemble_unitigs_core` 与 `contig::assemble`
  加 `anyhow::ensure!(min_count_seed >= 1)` 友好报错（multik/olc 已通过共享参数
  `MultikOptions` 校验覆盖）。

第十四轮结论：`cargo fmt`/`clippy` 干净；`cargo test --lib`（110）+ cli_asm_unitig
（14）+ cli_asm_olc（18）+ cli_asm_multik（12）+ cli_asm_extend（5）+ cli_asm_anchor
（4）全部通过。

## 第十五轮复核（olc 驱动 + overlap/layout 库纵深 + 输入解析健壮性）

对 olc 驱动（`cmd/asm/olc.rs`）与 overlap/layout/consensus 三个 OLC 库（
`libs/olc/overlap.rs`、`libs/olc/layout.rs`、`libs/olc/consensus.rs`）及独立阶段
命令（ovlp/layout/cns）的输入解析做了纵深复核，发现并修复 1 处校验缺口：

- **`cns --dedup-ratio` 缺范围校验**：`consensus_with_ratio` 直接消费 `ratio`
  未做边界校验。`ratio <= 0` 时 `dedup_contained_ratio` 的
  `coverage(...) >= ratio` 恒真，任何与已保留 contig 共享一个 31-mer 的 contig
  都会被丢弃（静默删掉几乎全部 contig）；`ratio > 1` 则静默退化为
  `ratio >= 1.0` 的精确子串语义，与帮助文本"<1.0 合并边界差异近重复"不符。
  修复：库层 `anyhow::ensure!(ratio > 0.0 && ratio <= 1.0)` 友好报错，与 multik
  驱动对 `--merge-similar` 的既有范围校验保持一致。回归
  `rejects_out_of_range_dedup_ratio`（0/-0.1/1.1 报错，1.0/0.5 正常）；`--dedup-ratio`
  帮助文本补 `(0.0, 1.0]` 范围。

复核确认无需修改（经逐例推导验证）的路径：

- **`overlap.rs::verify_seed`/`extend` 反向链延伸边界**：对反向 dovetail
  （a 后缀 = rc(b 前缀)）做了完整的手工推导，`minus` 分支
  `q[qs-1] == complement(t[te])` / `q[qe] == complement(t[ts-1])` 的方向与四端
  边界守卫（`qs>0 && te<m`、`qe<n && ts>0`）均正确，能精确还原 10 bp 反向重叠，
  无越界。
- **`layout.rs` 边界 contain 边**：contain 重叠位于长 contig **端部**时（q 是 t 的
  前缀/后缀）生成的扩展边可正确重建长 contig（短前缀 + t 剩余 = t）；仅中间 contain
  不生成边（`contain_overlaps_do_not_chain`）。`is_repeat`/`mutual-best`/`placed`
  防环与多义停止逻辑经推演正确。
- **`filter_contained` 等长与传递 contain**：等长重复取 id 小者、传递 contain
  （a⊃b⊃c）全部收敛到最长者，无双向误删。
- **`drop_cross_chimeras`**：头/尾 `flank` 覆盖、junction `span` 跨盖、`min_groups=2`
  的判定逻辑正确；单文件（全部同 tag）时 `min_groups` 不满足，静默 no-op 属设计语义。
- **`consensus`/`merge_overlapping_contigs`/`dedup_contained_ratio`**：首碱基 identity
  修复（第十四轮）、`overlap_geometry` 有界回退、合并后 `extend_seed_set` 的
  `kept_heads` 陈旧仅为少一次预筛（`kept_seeds` 已更新，无漏检）。
- **输入解析健壮性**：`cns` 布局 TSV 的 `contig_0` 下溢、超大 `contig_N` 巨型分配、
  非连续 id、缺名、非法 strand 均友好报错；`layout` PAF 缺名友好报错；`ovlp`/`olc`
  的 `--overlap-k`/空 unitigs 均友好报错；`read_unitigs` 的 tag 净化与冲突 `.i` 后缀
  确定。

记录项（未改，低风险 / 设计语义）：

- `--unitigs` 模式下 `--kmer` 仍被无条件解析，传无效值（如 `abc`）会报错——与帮助
  "ignored in this mode"（指无输出效果）略有出入，但 fail-fast 对拼写错误更友好，保留。
- `--kmer 21,21` 重复 k 值产生重名 unitig（`k21:unitig_<id>` ×2），下游
  `filter_contained` 收敛、输出正确，仅浪费计算；`--keep-dir` 的 unitigs.fa 会有重名
  头（仅调试用途）。

第十五轮结论：`cargo fmt`/`clippy --all-targets` 干净（仅既有依赖告警），`cargo test`
全量（lib 111 + 全部集成测试）通过。

## 第十六轮复核（组装内核纵深：extend/refine/unitig/contig/bubble/dfa）

对组装内核做了逐行纵深复核：`libs/asm/extend.rs`（交叉 contig 归属护栏）、
`libs/asm/refine.rs`（tadpole correct/extend/discard 移植）、`libs/asm/assemble/`
（unitig/contig/bubble/mod）与 `libs/asm/dfa.rs`（cuttlefish 风格 DFA 状态分类）。
发现并修复 1 处注释-实现不一致：

- **`unitig.rs` 覆盖度重建注释把"种子"误写成"最右 k-mer"**：
  `build_unitigs` 把左右两次走步的每条 k-mer 覆盖度重排成输出序列顺序。经
  逐例推导（含 R>0 具体数值例证），被丢弃的 `left_counts[0]` 是**种子 k-mer**
  的计数（RC 走步首窗口 = rc(seed)，canonical 计数与种子相同），与
  `right_counts[0]`（右走步首条目 = 种子）重复——是**首个**右条目，而非注释所写
  "last right entry"。代码正确，注释误导（若后续维护者按注释去改会把正确代码改坏）。
  修复：注释改为 "the seed k-mer / the first right entry"。

复核确认无需修改（经逐例推导验证）的路径：

- **extend.rs 的低覆盖 seam 回滚**：`low_run` 达 5 后 `truncate(MIN_LOW_RUN-1)`
  恰好裁掉该段推入的全部低覆盖碱基（第 5 步在 push 前 break），回滚方向对；
  `foreign_run` 达 5 后 `truncate(foreign_run)` 裁掉进入他 contig 领地的全部碱基；
  `saturating_sub` 防下溢；`top_two` 平局按索引序确定。
- **refine.rs 的 tadpole 移植**：`is_junction` 语义、`extend_to_right2` 的
  junction 基碱序（先捕 `b=right_max_pos` 再推进窗口）、`reassemble_inner` 的
  `ca=a+1-k` 由 `len>=k` 保证不下溢、`regenerate_counts` 重算窗口 `ca..=ca+k`
  （含一个无害的左侧多余窗口）、`clear_window2` 滑窗与 `count_errors` 的
  `i+k-1` 边界均安全。
- **contig.rs 的多遍播种**：`(1..contig_passes).rev()` + 末遍 `min_count_seed`
  共 16 遍，与注释一致；`pass_threshold`/`max_len=1e9` 上限不可由 CLI 触发
  （`min_extension`/`contig_passes`/`contig_pass_mult` 未暴露，仅默认值 2/16/1.7）。
- **bubble.rs 的 BubblePopper**：`merge`/`pop` 的 k-1 重叠裁剪、`is_loop`、
  `find_mutual_dest` 汇合判定、`mid_nodes_concur` 一致性检查均忠实移植。
- **dfa.rs**：`step`/`in_count_at` 的 RC 走向翻转（`3-in_base`、in/out 互换）、
  `canon_idx_pair` 二分、`succ_in/succ_out` 的 `NO_SUCC` 哨兵降级（非 panic）均正确。

第十六轮结论：`cargo fmt`/`clippy` 干净；asm 相关 `cargo test` 通过。

## 第十七轮复核（multik 内部纵深：bridge/master/graph/schedule/mod）

对 `libs/asm/multik/` 剩余的 bridge.rs / master.rs / graph.rs / schedule.rs / mod.rs
及 CLI（`cmd/asm/multik.rs`）做了逐行纵深复核。发现并修复 3 处问题：

- **`sequence_similarity` 零带宽退化**（真 bug，可复现）：`max_indel` 由
  `max(n,m) * (1 - min_similarity)` 截断为 0 时（即 `--merge-similar 1.0`，或
  短序列 + 近 1.0 相似度），原代码直接 `return 0.0`——**完全相同的两条序列也返回
  0.0**，导致 `bubble_merge` 判定 `0.0 < merge_similar` 而永不合并泡泡，
  `--merge-similar 1.0`（合法 CLI 范围 (0.0, 1.0]）静默失效。修复：`max_indel == 0`
  时按离散结果返回 `a == b ? 1.0 : 0.0`（任何替换/插入缺失都会使相似度 <
  `1 - 1/max(n,m)` < min_similarity，故 0.0 仍正确落在阈值之下）。新增回归测试
  `sequence_similarity_zero_width_band`。
- **`Master::round` 计时标签错误**（调试仅，`ANCHR_MULTIK_TIMING`）：`prog` 槽位误用
  `t5.elapsed()`（recompact 之后到计时点的总时长，混入 progressive_filter + split +
  trace），并非 progressive_filter 耗时。修复：`progressive_filter` 返回后记录
  `t6`，`prog` 改为 `t6.duration_since(t5)`。
- **`bubble_merge` 文档-实现不符**：doc 写替代路径长度上限 "`merge_len * k`"，
  实现为 `(merge_len * k / merge_similar).round()`（`max_len`，`bubble_merge_rejects_long_middles`
  测试亦依赖此式）。修复：doc 改为 "`merge_len * k / merge_similar`"。

复核确认安全（经推导/边界核对）的关键路径：

- **bridge.rs**：`bridge_kmer` 四个方向（u 左右端 × 伙伴正向/RC）的接点窗口
  `upstream 尾 (k-1) + downstream 续基` 均与 `u_ext == v_begin / rc(v_end)` 的
  实际端点匹配一致，续基下标（`v[k_build-1]`、`rc(v)[k_build-1]=comp(v[len-1-k_build])`）
  逐例验证正确；`probe_kmer` 的 `probe_half*2` 窗口两侧各取 `probe_half` 碱基，
  边界 `v_dir.len() >= km1+probe_half` 守卫完备；`split_by_bridge` 的 `RollCanon::new`
  由 `n < probe_len` 早退保护，`cut=0` 只作 run 跟踪标记不产生零长切片。
- **master.rs**：`SumView` 的 `is_solid` 短接（reads 单侧达标即 solid）与
  `count` 的 `saturating_add` 正确；round 第一步 link 校验用 `view.count`（`get_count`
  内部 canonical），`remove_unsupported` 用 `view.is_solid(km.canon())`
  （`get_count_canonical` 假定已 canonical），二者与 table.rs 的
  `get_count`/`get_count_canonical` 语义逐一对应；carried 重喂时
  `links.resize_with` 保持与 unitigs 索引对齐。
- **graph.rs**：`tip_remover` 的 `in_src`（去重计数）与 `in_cov`（最大源覆盖）O(n+E)
  预计算与旧 `any(l.to == i)` 语义一致；`weak_link_remover` 的
  `to = usize::MAX` 哨兵在两次 `from_rc` 扫描后才统一 retain，无相互干扰；
  `probe_stats` 的 `hi` 分位在 `left`/`right` 各至少含 1 个窗口时安全
  （`n >= window_len` 早退保证）；`recompact_graph`/`merge_chains` 的链头回退
  （`seen` 防环）+ 右向走步（`visited` 防环）+ 严格度 1 不变量一致。
- **schedule.rs**：`assemble_one`/`assemble_all_masters` 的表生命周期
  （bound-lookahead 窗口、`OnceLock` 共享、cap-1 channel 反压）保证峰值内存
  约两表 + probe；`auto_ks` 阶梯与 `docs/asm.md` 一致。

**记录在案、未改行为的设计观察**（改动需质量门禁 A/B 验证）：

- **`weak_link_remover` 只删物理边的单侧镜像，实际近乎 no-op**：
  `compute_links` 把同一物理接点从两端各发一条记录（u 的出边 `from_rc=false` +
  v 的入边 `from_rc=true`），段图按 `seen` 去重后只认物理边。`weak_link_remover`
  仅标记 `links[i]` 侧 `to = usize::MAX`，镜像记录仍存活 → 物理弱边在
  `merge_chains`/`recompact_graph` 的段图中依然存在。又因弱边总出现在分支端
  （`depths.len() > 1` 才删除），而链压缩要求两端严格度 1（`out_deg == in_deg == 1`），
  分支端本就永不压缩 → 该函数在最终链路中**不影响输出**（近似无效而非有害）。
  完整修复（双侧镜像删除）会允许主路径穿透原分支节点、改变 N50/GF，属
  启发式行为变更，须按 `results/model_org.md` 质量门禁先 A/B 验证再合入；
  本轮仅记录，未改动。

第十七轮结论：`cargo fmt`/`clippy` 干净；multik 20 个单元测试 + 全量
`cargo test -- --test-threads=1` 通过（含新增回归）。

## 第十八轮复核（libs/asm/table.rs 整文件纵深）

对 `libs/asm/table.rs` 全文（build 四路径 / prefix_index / locate / get_count /
sorted/solid 快照 / count_keys_seq / merge_tables / Kmer / base_code 族 /
count_read_kmers_packed）逐行纵深复核，重点核对与调用方（contig/unitig/refine/
fq merge/multik）的契约。发现并修复 3 处问题：

- **空 reads 构建的 k=0 空表可被查询 panic（真 bug，可复现）**：
  `build_threaded` 的 `reads.is_empty()` 分支返回 `KmerTable::default()`（`k=0`、
  `keys/counts` 空）。此后任何查询都会 panic：`prefix_index` 的 `keys.len()/kb`
  除零、`locate` 对空 `q` 取 `q[0]/q[1]` 越界、`sorted_entries`/`solid_entries` 的
  `chunks_exact(0)`。实际触发路径：`anchr asm contig`/`unitig` 对**空输入文件**
  （0 reads）在 parallel=0 路径建出 k=0 表，随后 `scan_table`→`sorted_entries`、
  `assemble_unitigs_from_table`→`solid_entries` 直接 panic（violates Zero-Panic）；
  `fq merge`/`refine` 目前因处理循环由同一份空 reads 驱动、不会查询空表而侥幸
  安全，但属于脆弱隐式不变量。修复：空分支改为保留 `k` 的
  `KmerTable { k, keys: vec![], counts: vec![] }`，使空表与其它 k 集空表一致
  （n=0，所有查询安全返回空/0）。新增回归 `empty_reads_table_stays_queryable`
  （get_count / find_row / sorted_entries / solid_entries / solid_row_ranks /
  fill_left_counts / fill_right_counts 全部安全且为空）。
- **`build_threaded` 的 `threads` 参数文档与实现不符**：doc 声称
  "`threads == 0` uses the rayon global pool; otherwise a private pool of exactly
  `threads` workers"，实现完全忽略 `_threads`，始终用环境 rayon 池（body 注释
  已说明这是有意为之：调用方已用 `--parallel` 包一层池，再建私有池会超订）。
  修复：doc 改为说明 `threads` 仅为调用侧对称保留的 advisory 参数，不建私有池。
- **`base_code`/`base_comp_code` 文档与实现不符**：doc 写非 ACGTU 返回 -1，
  实现返回 0（返回类型 `u8` 亦不可能为 -1）。调用方均先 `base_defined` 门控
  （`count_read_kmers_packed`、refine 各窗口、contig/unitig/extend），`_ => 0`
  分支实际不可达。修复：doc 改"其它碱基映射为 0，调用方必须先经 `base_defined`
  门控"。

**排除的疑点（经核验无需修复）**：

- `count_read_kmers_packed` 的 `fw[..half] <= rc[..half]` 半字节规范判定对
  `k % 4 != 0` 同样成立（边界字节零填充对称性 + revcomp 对合性），已由
  `k81_table_counts_match_bruteforce`（k=81 非 4 的倍数）实测背书。
- `to_bytes()` 只返回 `key_bytes` 字节，`extend_from_slice(fw)` 不会多 emit。
- supermer / direct / streamed 三条路径空输入均保留 `k`（pgr `count_keys` 与
  `build_impl` 的 n_records==0 分支），仅 `build_threaded` 空分支是 k=0 唯一来源。
- `Kmer::default()`（derive 出 k=0）无生产调用点；生产只用校验过的 `Kmer::new(k)`。
- `count_keys_seq`/`merge_tables` 的计数 `u32` 上限与 `(j-i).min(u32::MAX)` 截断：
  需 >4G 次相同 k-mer 才可能溢出，现实不可达，未改。

第十八轮结论：`cargo fmt`/`clippy` 干净；`cargo test --lib libs::asm::table`
6 个测试通过（含新增回归 `empty_reads_table_stays_queryable`）。

## 第十九轮复核（libs/asm/dfa.rs + assemble/bubble.rs 整文件纵深）

对 `libs/asm/dfa.rs`（unitig 行走 DFA 状态机）和 `libs/asm/assemble/bubble.rs`
（contig 泡泡消除）整文件逐行纵深复核，重点验证与 BBTools Tadpole
`BubblePopper.java`（Debian 源 39.01+dfsg-2）参考实现的语义一致性。发现并修复
1 处真 bug，其余疑点经核验排除：

- **`bubble.rs` 间接泡泡路径恒失效（真 bug，可复现）**：`expand_right` 中
  `self.dest = dest_id` 赋值被错误放在 `mid_nodes_concur` 调用**之后**（上一版
  结构）。而 `mid_nodes_concur` 用 `self.dest`（此时仍为 L244 重置的
  `usize::MAX`）与每个 mid 的 `right_dest` 比较（L448 `right_dest != dest_id`），
  条件恒成立 → 恒返回 false → 间接泡泡（indirect bubble）合并路径完全失效，
  仅简单（direct）泡泡合并有效，contig 图只简化掉最外层分支。参考实现
  `BubblePopper.expandRight` 在 `midNodesConcur` 调用前已明确赋值 `dest`。修复：
  将 `self.dest = dest_id` 移至 `mid_nodes_concur` 调用前（L294-299，含注释说明
  与参考实现的对齐依据）。

**排除的疑点（经核验无需修复）**：

- `dfa.rs` 状态分类：`VertexStates` 的 in/out 度、唯一后继（`succ_out`/`succ_in`）
  与 `step()` 的朝向翻转（`fw_is_canon` → `out_base`，否则 `3 - in_base`）与
  unitig 行走的 canonical 约定一致；并行/串行两类 build 路径状态向量长度对齐
  `sorted_entries`，无越界。
- `bubble.rs` `expand` 的翻转（`center.flip(destMap.get(center.id))`）顺序与
  参考一致；`find_representative_mid_edge`/`fetch_mid_nodes` 的 `min_len`
  （`2*k-1`）门控与 `BubblePopper` 构造参数一致。
- `bubble.rs` `pop` 的 right_mid_edge 兜底 `get_right_edge(dest_id, Some(1))`
  在 `mid_nodes_concur` 已保证各 mid 右 dest 一致的前提下，退化为取首个边，
  语义等价于参考的 `findRightmostEdge`，无风险。

第十九轮结论：`cargo fmt`/`clippy` 干净；`cargo test --lib libs::asm::assemble`
通过（含 `contig`/`bubble` 相关测试）；L1 smoke 字节级 diff 无回归。

## 第二十轮复核（cmd/asm 命令层 + assemble 库入口 + 文档一致性）

对 `cmd/asm` 命令层（mod.rs / common.rs / contig.rs / unitig.rs / layout.rs /
olc.rs）与 `assemble` 库入口（contig.rs / unitig.rs）及 `docs/asm.md` 纵深
复核。发现并修复 3 处问题：

- **`is_junction` 跨模块重复实现（漂移风险）**：`refine.rs` 与
  `assemble/mod.rs` 各有一份逐字相同的 `is_junction(max, second, opts)`，
  仅 `opts` 类型不同（`RefineOptions`/`AssembleOptions`）。公式若在任一处
  修改，`asm contig`/`asm unitig` 与 `fq` 校正/延伸行为将静默分叉。参照第
  十四轮删除 `dominant_offset` 重复的先例，上提为 `libs/asm/mod.rs` 的标量
  参数共享函数（单一事实来源），两处调用点（refine 5 处、contig 9 处）
  改为透传 4 个标量，`assemble/mod.rs` 改为 `pub(crate) use super::is_junction`
  保持 contig.rs 导入不变。
- **`read_unitigs` 标签去重边界缺陷（真 bug，可复现）**：`.i` 后缀合成后
  不复查唯一性。特定文件茎组合（如 `["foo.fa", "foo.2.fa", "foo.fa"]`，
  `foo.2` 文件夹在两个 `foo` 之间）下，第三个文件被分配与第二个相同的
  tag `foo.2`，`used` 集合不去重合成结果 → 输出重复 unitig 名。影响
  ovlp/olc/cns/layout 四个共用 `read_unitigs` 的命令（重复名破坏 PAF/布局
  的 name→id 引用）。修复：冲突时从文件下标起循环递增直至 tag 唯一。
- **`assemble_unitigs_core` 空输入防护缺失（潜在 panic）**：`--supermer`
  的 FASTA 探测 `SeqReader::new(&infiles[0])` 位于任何空输入检查之前，
  `read_records` 的空 `ensure!` 在其后才触发。所有 CLI 调用方（unitig/olc/
  multik）目前先校验非空而侥幸安全，但 `pub(crate)` 函数缺前置守卫属脆弱
  隐式不变量（与第十八轮空表 panic 同类）。修复：入口首行加
  `ensure!(!infiles.is_empty())`。

**排除的疑点（经核验无需修复）**：

- k 上界 256 在六处库入口（contig/unitig/refine/multik/extend/table）均以
  `ensure!` 兜底，`-k` 传超界返回友好错误而非 panic（零 panic 满足）。
- `--supermer-m` 由 pgr `build_table_slices_with_m` 校验 `m ∈ [2, k-1]`，
  非法值友好报错。
- `unitig` 的 FASTQ→直连计数回退（docs 声称"FASTQ 自动回退"）在
  `assemble_unitigs_core` 真实实现（`fasta_input` 探测 + `use_supermer`
  `&& fasta_input` 门控）。
- `format_cov` 的负输入分支不可达（coverage 为 k-mer 计数均值，恒非负）。
- `cmd/asm/mod.rs` 十个子命令注册/分发完整；`docs/asm.md` 与各命令 CLI
  帮助文本一致（含 multik 固定阶梯、olc `--keep-dir`/`--cross-validate`）。

第二十轮结论：`cargo fmt`/`clippy` 干净；`cargo test --lib libs::asm` 34 个
通过；`cli_asm_contig`/`cli_asm_unitig`/`cli_asm_olc` 集成测试 10/14/18 通过；
L1 smoke 字节级 diff 无回归。

## 第二十一轮复核（剩余 cmd 命令层 + olc 库入口 + 参数解析）

对前几轮未深覆盖的命令层与库入口做收尾纵深复核：`cmd/asm/multik.rs`、
`cmd/asm/cns.rs`、`cmd/asm/ovlp.rs`、`cmd/asm/extend.rs`、`cmd/asm/anchor.rs`、
`cmd/args.rs::parse_parallel`、`libs/olc/consensus.rs`、`libs/olc/overlap.rs`、
`libs/olc/layout.rs`。**未发现新问题**，核验到的既有防护如下：

- `multik`：`--merge-similar ∈ (0.0, 1.0]`、`--merge-len ∈ 1..=1024`、
  `--parallel ∈ 0..=1024` 均 cmd 层校验；`--use-guide` 依赖 `--all-masters`
  （clap requires）；`--print-ks` 经 `auto_ks_for_reads` 派生；guide 伪 reads
  临时文件由 `_guide_keep` 保活至组装结束（生命周期正确）。
- `cns`：`--dedup-ratio` 范围校验在库层 `consensus_with_ratio`
  （`(0.0, 1.0]`）；`parse_layouts` 对 7 字段、`contig_0` 拒收、contig/step
  连续性、布局 id 巨大值（防巨型分配）全部友好报错（零 panic）。
- `ovlp`：`find_overlaps` 对 seed_k 同时校验空/过短 unitigs（`seed_k >= 1`）
  与上界（`<= MAX_K`），并 clamp 到最短 unitig。
- `extend`：库层 `extend_contigs` 校验 `k ∈ 2..=MAX_K`、
  `min_support >= 1`（第 12 轮修复确认）。
- `anchor`：`--lscale > 0`/`--uscale > 0` 除零防护、`--stats` 与 `-o` 及输入
  的三方冲突前置校验（第 11 轮修复确认）；SAM 解析对未知参考/坏 CIGAR/
  坏 POS 均 bail，`pos + mlen - 1` 因 SAM POS 1-based（pos≥1）无下溢。
- `parse_parallel`：`auto`=全核、`half`=`div_ceil(2).clamp(1,8)`、整数
  `1..=1024`，非法值友好报错。
- `build_layouts`：overlap 超过前一步末端有 `ensure!` 防护（零 panic）。

第二十一轮结论：纵深复核未发现新问题，收尾轮干净；`cargo fmt`/`clippy`
干净，`cargo test --lib` 与相关集成测试通过，L1 smoke 字节级一致。

## 结论

`asm` 命令族审核完成（累计修复 5 类问题 + 4 处 `-o` 防护统一 + 1 处 `--keep-dir`
文档修正 + 1 处 `--outm`/`--outu` 冲突防护 + 1 处 `--overlap-k` 上限校验 + 第七轮
3 处：`map` 帮助单复数、`contig`/`unitig` 的 `-k` 上限帮助、`cns` 布局 id 连续性与
巨型分配防护 + 第九轮 2 处：`unitig --links` 悬空引用、`asm_map_benchmark.rs` 编译
错误 + 第十一轮 6 处：`extend` 的 `-o` 数据安全、`anchor --stats` 校验 fail-fast、
`docs/asm.md` 补 anchor/extend 文档、`docs/asm.md` k 上限 128→256、`extend` 的 `-k`
帮助范围、`docs/asm.md` multik 段补 `--merge-similar`/`--merge-len` + 第十二轮 3 处：
`bridge_filter` 文档-行为一致性、`extend --min-support 0` 校验、multik trace 无效目录
容错 + 第十三轮 4 处：`--supermer-m` minimizer 长度与 pgr 单点事实来源对齐
（`minimizer_len(k)`，含 `k-1` 上界与 `ceil(k/4)`）、`docs/asm.md` `--supermer` 补
`k-1` 上界、`progressive_filter` 注释 cutoff 上限改 25% 中位数、`docs/asm.md` cns
输出头逗号改空格 + 第十四轮 3 处：banded Levenshtein 首碱基漏计修复、
`dominant_offset` 重复实现删除（复用 `overlap_geometry`）、`--min-count-seed 0`
校验（unitig/contig）+ 第十五轮 1 处：`cns --dedup-ratio` 范围校验 (0.0, 1.0] +
第十六轮 1 处：`unitig.rs` 覆盖度重建注释"last right entry"→"the seed / first
right entry"（代码正确，注释误导）+ 第十七轮 3 处：`sequence_similarity` 零带宽
退化修复（`--merge-similar 1.0` 时完全相同的序列不再返回 0.0）、`Master::round`
计时 `prog` 槽位改用 `t6.duration_since(t5)`、`bubble_merge` 长度上限文档改
`merge_len * k / merge_similar`，+ 第十八轮 3 处：空 reads 构建的 k=0 空表
查询 panic 修复（空输入文件不再 `sorted_entries`/`solid_entries`/`prefix_index`
panic）、`build_threaded` 的 `threads` 参数文档修正（advisory、不建私有池）、
`base_code`/`base_comp_code` 文档改"其它碱基映射 0、须 `base_defined` 门控"，
第十九轮 1 处：`bubble.rs` `expand_right` 的 `self.dest = dest_id` 赋值移到
`mid_nodes_concur` 之前（修复间接泡泡路径恒失效，对齐 BBTools
`BubblePopper.expandRight`），
第二十轮 3 处：`is_junction` 重复实现上提为共享标量函数（消除 refine 与
assemble 的公式漂移风险）、`read_unitigs` 标签去重循环递增至唯一（消除
重复 unitig 名数据完整性缺陷）、`assemble_unitigs_core` 入口补空输入
`ensure!`（消除 `--supermer` FASTA 探测的 `infiles[0]` 潜在 panic），
第二十一轮为干净收尾轮（未发现新问题，复核 multik/cns/ovlp/extend/anchor
命令层、olc 库入口与 parse_parallel，确认既有防护完备），
经纵深复核收敛；与 BBTools/BCALM
语义对拍、边界输入验证零 panic，`cargo fmt`/`clippy` 干净（asm 相关无新增告警），
相关集成测试与 `cargo test --lib` 全部通过。
