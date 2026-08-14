# pgr fq 命令族代码审核记录（2026-08-11）

> 维护约定：本文是迁移时的审计快照；fq 命令族已迁至 anchr，后续 anchr 侧
> 对 fq 的修改应更新本审计记录（而非 pgr 侧）。

对 `pgr fq`（10 个子命令）命令族及相关库文件（`libs/fq/{clump,norm,pairs,
sample,split,trim,trim_adapter}`、`libs/fmt/fq`、`libs/fmt/seq`、`libs/loc`、
`libs/par`）和全部测试/文档进行审核。本报告由 `audit-fa-fq-2bit.md` 中的 fq
部分拆分并扩写，覆盖 fq 全命令族。以下仅保留有借鉴意义的结论；验证过程已精简。

## 命令族

`fq` 现有 10 个子命令：`interleave`(`il`)、`to-fa`、`clump`、`split`、
`sample`、`clean`、`filter`、`norm`、`range`、`trim-qual`。

## 与外部参考实现的语义一致性核对

fq 家族以 BBTools 39.38 为主参考，逐字节核对（golden 数据见
`tests/bbtools/Lambda/README.md`），另有 `fa range` 复用的 `.loc` 索引与
`sickle`/`cutadapt` 风格的质量修剪。均已核验一致：

- `clump` → `clumpify.sh`（k=31, seed=1 默认输出顺序，逐字节）。
- `split` → `repair.sh` `rp` 模式（R1/R2/singles 逐字节）。
- `sample` → `reformat.sh samplebasestarget sampleseed`（精确模式，去重下无
  上行采样）。
- `clean`/`filter` → `bbduk.sh`（anchr 管线两趟调用，`ordered=t` 确定性）。
- `norm` → `bbnorm.sh passes=1 bits=16 min=<n> target=9999999` 的读决策逻辑
  （计数用精确 canonical 表，非近似 `bits=16` 哈希表）。
- `trim-qual` 的 sliding/mott → `sickle`/`cutadapt` 的质量修剪算法。
- `range` → 复用 `fa range` 的 `.loc` 索引语义（name、plain offset、record
  size），BGZF/plain 支持一致。

有意差异（均已记录）：
- `norm` 计数用精确表而非 `bits=16` 近似计数，深度边界附近的读可能因计数精度
  不同而被取舍（文档已注明）。
- `clump` 外部 bucket 路径输出为"按桶拼接"序，与 BBTools 大数据行为一致但与
  全局内存序不同（文档已注明）。

## 排除的疑点（安全不变量，经核验无需修复）

- `fq is_fq` 对目录输入：`File::open` 成功但 `read_exact` 失败（EISDIR），返回
  友好错误而非 panic（Zero-Panic）。
- `sample` 的 `score_cov`/`read_stats` 索引：`cov` 非空后，`above_limit` 回退到
  -1 时分支短路、`depth_al` 保持 -1，不会越界索引。
- `norm` 外部路径 `load_table` 置 `k: 0`：`read_stats`/`canonical_keys` 均用
  `opts.k`，不依赖表的 `k` 字段。
- `trim.rs` `sliding_cut`/`mott_cut`：窗口尺寸至少 1、空序列提前返回，无除零/
  越界。
- `trim_adapter.rs` 的 `make_codes`/`antialias`/`JavaRandom`：与 JVM 序列化
  逐项有单元测试锁定。
- `fq range` 的 LRU 缓存：`LruCache<String, Vec<u8>>` 用 `&str` 借钥命中，
  `Borrow` 关系成立。

## 记录项（未改，低风险 / 待决策）

- `fq range` 的 `name:0-N`（start=0）被当作"整条记录"而非子序列（与 `fa range`
  的 `start==0` 约定一致）；未文档化该边界，但不返回错误数据（仅取整条）。
- `norm` 输出会对保留的读施加 `changequality` 归一化（N 质 0、ACGT 质最低 2），
  与 bbnorm 默认一致；文档已说明 changequality 被应用。

## 已知限制（有意保留）

- `interleave` 双文件格式只按 `infiles[0]` 探测（`is_fq` 只看第一个文件），若
  第二个文件格式不同会在实际读取时报错（非静默）。
- `trim-qual` 配对质量编码只从 R1 自动检测（文档已注明）。

## 修复的缺陷（根因模式）

按根因模式分组；标注"承袭"者自 `audit-fa-fq-2bit.md` 沿袭，其余为本次审核新发现。

### Zero-Panic / clap 参数缺失

- **`split` 缺少 `--outfile-2` 时 panic**：`args.get_one::<String>("outfile_2")
  .unwrap()` 对 `None` 解包崩溃。修复：`outfile_2` 参数加 `.required(true)`，
  clap 先于执行校验，缺失时输出友好用法错误。新增回归测试
  `command_fq_split_missing_outfile2_is_clap_error`。
- **`sample` 缺少 `--bases` 时 panic**：`args.get_one::<i64>("bases").unwrap()`
  对 `None` 解包崩溃。修复：`bases` 参数加 `.required(true)`。新增回归测试
  `command_fq_sample_missing_bases_is_clap_error`。
- **`clean`/`filter`/`norm`/`clump` 的 `--parallel` 未做 1..=1024 范围校验**：
  四个命令的 `--parallel` 均接受 `auto` 或整数，但只 `parse::<usize>()` 未限上
  界，数值直接进入 `rayon::ThreadPoolBuilder::num_threads` 创建线程池
  （clean/filter 经 `par::ordered_map`、norm/clump 直接建池），越界值（如
  `--parallel 1000000`）会尝试创建海量线程，导致系统资源耗尽。违反全局硬约束
  "`--parallel` 必须经 clap 校验 1..=1024"。修复：新增共享助手
  `cmd_pgr::args::parse_parallel_auto`（`auto` 取逻辑 CPU 数，整数须在
  1..=1024，越界友好报错），四处统一改用。新增回归测试
  `command_fq_{clean,filter,norm,clump}_parallel_out_of_range_*`。
- **`clump` 的 `--buckets` 未校验范围，`--buckets 0` 触发除零 panic**：
  `clump_buckets` 中 `key.kmer as u64 % buckets as u64`（`libs/fq/clump.rs`）
  对 `--buckets 0` 除零崩溃。修复：`execute` 中校验 `--buckets` 须在
  1..=4096，越界友好报错（与 clump 内部 `MAX_BUCKETS=4096` 一致）。新增回归
  测试 `command_fq_clump_buckets_out_of_range_is_friendly_error`。
- **`sample` 输入以空（0 碱基）记录结尾时除零 panic**：`sample` 循环里
  `remaining` 每轮按记录碱基数递减，当输入末尾存在空记录时 `remaining` 先减到
  0，下一轮 `target / remaining`（`libs/fq/sample.rs`）除零崩溃。修复：在除法
  前若 `remaining == 0`（说明剩余全为空记录）直接 `break`。新增回归测试
  `command_fq_sample_trailing_empty_records_do_not_panic`。

### 数据安全（`-o` 同输入保护）

- **流式命令允许 `-o` 覆盖输入文件**：`fq to-fa`/`fq interleave` 已统一加入
  `ensure_outfile_distinct`。

### 输入校验 / 静默错误

- **`interleave` 双文件交错对读取计数不匹配静默截断**（`zip` 取较短者）。修复：
  `interleave_read` 中任一文件先读完而另一未读完即 `bail!`。

### 行为一致性 / 算法

- **`interleave` 单文件虚拟 R2 两路径不一致**：单 FQ→FA 为 `"\n"`（空序列）、单
  FA→FA 为 `"N"`；帮助与 `docs/fq.md` 均声明 "N's"。修复：统一为 `b"N"`。
- **`interleave` 双文件路径返回的最终索引错误**：更新后的 `idx` 被丢弃，最终返回
  未递增的 `start`，违背 pub fn 契约。修复：两文件分支改为
  `idx = interleave_read(..)?`。

### 文档一致性

- **`trim-qual` 命令名错写为 `trim-q`**：`trim_qual.rs` 帮助文本/示例、
  `libs/fq/trim.rs` 的 `TrimOptions` doc、`docs/fq.md` 均改为 `trim-qual`。
- **`docs/fq.md` 子命令清单不完整**：补全 10 个子命令列表。
- **`docs/fq.md` 缺失 `trim-qual` 小节**：新增完整 Options/Examples。
- **`clean` 文档 gzip 输出示例误导**（`io::writer` 写端不压缩）：
  `-o unmerged.trim.fq.gz` 改为 `-o unmerged.trim.fq` 并仅指明输入可为 gzipped。
- **`to-fa` 文档误置于 `norm` 小节**：移回其所属小节。

### 死代码 / 功能不可达

- **`clean` 的 `--mask-kmers` 静默失效（死代码）**：`trim_adapter.rs` 中按
  `ktrim_right` 分派"ktrim / kmask / filter"三个分支，但 `clean.rs` 把
  `ktrim_right` 硬编码为 `true`，使 kmask 分支永远不可达。于是文档化的
  `--mask-kmers`、`--mask-fully-covered`、`--trim-pad` 三个选项被静默忽略，
  掩码功能完全失效（`filter` 用 `ktrim_right: false` 走 filter 分支，不受影响）。
  修复：
  - `ktrim_right` 改为 `kmask.is_none()`——默认保持 ktrim=right（与 bbduk
    逐字节 golden 行为不变），指定 `--mask-kmers` 时切到 kmask 掩码分支。
  - 新增守卫：`--mask-fully-covered` / `--trim-pad`（仅掩码语义）在未给
    `--mask-kmers` 时报友好错误，避免静默无效；`--mask-kmers` 在未给 `--ref`
    时报友好错误（与文档"requires --ref"一致），避免静默无操作。
  - 新增回归测试 `command_fq_clean_kmask_masks_instead_of_trims`（掩码为 N、
    全长保留 vs 默认 ktrim 截短）、
    `command_fq_clean_kmask_mask_only_options_require_mask_kmers`（友好错误）、
    `command_fq_clean_kmask_requires_ref`（缺 --ref 友好错误）。
  `docs/fq.md` 中 `--mask-kmers`/`--mask-fully-covered`/`--trim-pad` 的描述
  现与行为一致。

## 第二轮：merge 家族（merge / ec-kmer / ec-overlap / extend / s-filter）审核

首轮报告完成后，代码库新增了 anchr merge 管线的 5 个子命令及其库（`libs/fq/{merge,
overlap,bbnet}`、`libs/asm/{tadpole,assemble}`、`libs/kmer/{quality,qcheck}`）。
本报告拆分后这部分未被首轮覆盖，故新增第二轮审核。逐命令/库核对 Zero-Panic、clap
参数、数据安全（`-o` 同输入保护）、算法边界与文档一致性。

### 排除的疑点（安全不变量，经核验无需修复）

- `merge`/`ec-overlap` 的 `process_pair` 对 `seq<2` 短读提前返回 `Ambiguous`，
  `join_reads` 的 `overlap==0` 分支切片起点恒 ≤ 终点（`insert` 已由
  `min_insert>0` 保证），`corrected_pair` 的 `rc2.len()-c2len` 由
  `c2len ≤ b.len()` 保证不欠位；`expected_tip_errors` 中 `bases`/`quals` 同长。
- `bbnet::parse` 对 23 维输入层、非空层、逐层权重维度做了校验，
  `feed_forward` 的最后 `last()/first()` 均有保障，畸形 `.bbnet` 返回友好错误
  而非越界。
- `tadpole::run` 对 `k∈[1,128]` 校验；`count_read_kmers` 对 phred 索引
  `min(127)` 限位；`extend_read` 的 rollback 用 `saturating_add(1)` 防 `%0` 与
  `usize::MAX+1` 溢出。
- `libs/kmer/quality.rs` 对 `k > MAX_K` 提前返回空表；`qcheck::check_read` 对
  短读/无锚点返回 `NoAnchor` 而非 panic。

### 修复的缺陷

- **`s-filter` `-k ≥ 65` 触发 u128 移位越界 panic**：`qcheck` 的锚点/延伸扫描用
  u128 滚动 k-mer（每碱基 2 bit，最多 64 碱基），`masks` 返回 `rc_top=2k-2`，对
  `k≥65` 时 `1u128 << rc_top` 等移位在 debug 下 panic
  （实测 `pgr fq s-filter -k 65` 崩）。修复：`execute` 校验 `k∈1..=64` 并报友好
  错误。新增回归测试 `command_fq_s_filter_kmer_out_of_range_is_friendly_error`。
- **`s-filter` 的 `--discard-file` 可覆盖输入文件**（数据安全）：仅 `-o` 受
  `ensure_outfile_distinct` 保护，`--discard-file` 指向输入会破坏原文件。修复：
  对 `discard_file` 也加 `ensure_outfile_distinct`。新增回归测试
  `command_fq_s_filter_discard_file_same_as_input_rejected`。
- **`ec-overlap` 的 `--efilter 0` 未真正禁用**（行为不一致）：`merge` 用
  `(*x>0.0).then_some(*x)`（0 置 `None` 禁用），`ec-overlap` 却 `Some(*x)`，
  使 `--efilter 0` 在 `process_pair` 中 `(expected+offset)*0 < bad` 恒真，意外
  抑制 pfilter。修复：`ec-overlap` 与 `merge` 对齐。
- **`bbnet` 单元测试与新增 23 维校验冲突**：`parse_and_forward_minimal` 用旧的
  2 维输入网，被 `dims[0]==23` 校验拒绝后 `unwrap()` panic。修复：测试改用 23 维
  输入（隐藏层补足 23 权重），并新增 `parse_rejects_wrong_input_dims` 用例。

### 文档一致性

- **`docs/fq.md` 的 `ec-kmer`/`extend` `-k` 标注 "no upper bound" 不实**：实际
  `tadpole::run` 限制 `k≤128`。改为 "up to 128, the k-mer key table limit"。
- **`docs/fq.md` 的 `s-filter` `-k` 未注明上限**：标注 `1..=64, the u128
  rolling-key limit`，与新校验一致。
- **`docs/fq.md` `ec-overlap` 与 `extend` 之间缺 `---` 分隔**：补全，保持小节间
  分隔一致。

## 第三轮：迁移重构后复核 + 未覆盖路径补查

期间外部并行完成了 anchr 迁移阶段 1（`notes/design/fq-asm-migrate.md`）：
`detect_quality_base`/`PHRED33`/`PHRED64` 从 `libs/fq/trim.rs` 抽到基础层
`libs/fq/qual.rs`，`kmer qhist`/`fq s-filter`/`trim` 改引用新位置；kmer 的
`base_codes`、`count::count_keys` 由 `pub(crate)` 改 `pub`；新增
`tests/migrate_api.rs`（外部 crate 视角编译证明 anchr 依赖的基础符号可达）。

复核结论（均绿）：
- 重构未引入新问题：`cargo build`/`clippy --all-targets -D warnings`/`fmt --check`
  clean，全量 `cargo test` 通过（含 `migrate_api` 与既有 fq/asm 测试）。
- 无残留旧路径引用：`grep trim::detect_quality_base|trim::PHRED` 无命中；
  `detect_quality_base` 的 `seq[i]` 因 `SeqReader` 强制 seq/qual 等长而安全。
- 补查此前未逐行覆盖的路径，确认校验齐全：
  - `clump`：`-k` 2..=31（`fill_max` 的 `x2 << (2k-2)` i64 移位上限）、
    `--buckets` 1..=4096、`--parallel` 走 `parse_parallel_auto`；
  - `norm`：`-k` 2..=31、`--parallel` 走 `parse_parallel_auto`；
  - `sample`：`remaining==0` 提前 break 防除零，`--bases` 必填；
  - `s-filter`：`-k` 1..=64（第二轮修复）、`--discard-file` 防覆盖输入。
  未发现新缺陷。

## 结论

首轮：`fq`（10 子命令）命令族合计修复本次缺陷（Zero-Panic / clap 参数缺失 4：
`split --outfile-2`、`sample --bases`、`clump --buckets 0` 除零、
`sample` 空记录结尾除零；`--parallel` 范围校验 1 [clean/filter/norm/clump
四处]；死代码/功能不可达 1 [clean `--mask-kmers`]），均含回归测试；另有承袭自
`audit-fa-fq-2bit.md` 的数据安全 1、输入校验 1、行为一致性/算法 2、文档一致性 5。
全部 fq CLI 集成测试与 fq 库单元测试通过，`cargo clippy -- -D warnings` 与
`cargo fmt --check` clean。

第二轮：merge 家族（`merge`/`ec-kmer`/`ec-overlap`/`extend`/`s-filter` 5 个子命令
及 `libs/fq/{merge,overlap,bbnet}`、`libs/asm/refine`、`libs/kmer/{quality,qcheck}`
）修复缺陷 4（Zero-Panic 1 [s-filter `-k≥65` 移位越界]、数据安全 1 [s-filter
`--discard-file` 覆盖输入]、行为一致性 1 [ec-overlap `--efilter 0`]、测试修复 1
[bbnet 23 维校验与单元测试对齐]），文档一致性 3，均含回归测试或已修正文档。

第三轮：anchr 迁移阶段 1 重构（`detect_quality_base` 抽到 `qual.rs`、kmer pub 化）
复核全绿，无残留旧引用；补查 `clump`/`norm`/`sample`/`s-filter` 校验齐全，未发现
新缺陷。

经多轮纵深复审收敛，未再发现新问题。

---

# 第四轮（2026-08-13）：anchr 迁移后 `fq` 全命令族复核

`fq` 已迁移至 `src/cmd/fq/`（16 个子命令：to-fa / clump / interleave / merge /
norm / range / sample / split / clean / ec-kmer / ec-overlap / extend / filter /
s-filter / trim-qual），相关库位于 `src/libs/fq/`（trim_adapter / trim / clump /
norm / sample / split / merge / overlap / bbnet）与 `src/libs/asm/`（tadpole /
assemble，与 asm 共享）。本节对迁移后的现行代码做全命令族纵深复核。

## 复核范围与方法

- 逐行重读全部 16 个命令文件、`cmd/args.rs` 共享参数助手，以及
  `libs/fq/`、`libs/asm/refine.rs`、`libs/asm/assemble.rs` 的关键路径。
- 逐命令核对 Zero-Panic（索引/溢出/除零/unwrap）、clap 参数校验、数据安全
  （`-o` 同输入保护 + 多输出互斥）、算法边界（k 上界、滑动窗口、桶除零、
  质量修剪、merge/join 切片）、文档一致性（`docs/fq.md` ↔ 帮助 ↔ 行为）。
- 对首轮此前未覆盖/未复核的边界逐项核验（见下）。

## 复核确认（安全不变量，经核验无需修复）

- **`trim.rs` `sliding_cut` / `mott_cut`**：`n==0` 提前返回 `(0,0)`，
  `window_size = max(1, n/10) ≥ 1`，循环条件 `window_start + window_size <= n`
  保证 `q(i)` 索引不越界；`mott_cut` 的 `start/stop` 满足 `start >= stop` 时返回
  `(0,0)`，不会产生欠位切片。空序列、极短读均安全。
- **`norm.rs` `score_cov`**：`above_limit` 递减到 -1 后循环条件 `>= 0` 短路，
  `depth_al` 保持 -1，`cov[(-1)]` 不会被求值；`load_table` 的
  `chunk[kb..].try_into().unwrap()` 因 `bytes.len() % (kb+4) == 0` +
  `chunks_exact` 保证每块恰 `kb+4` 字节而不会 panic。
- **`clump.rs` `fill_max` 的 `1 << (2k-2)`**：`k` 在 clap 层已限 `2..=31`，
  `2k-2 ≤ 60`，i64/u128 移位安全；`clump_buckets` 的 `% buckets` 由 `--buckets`
  `1..=4096` 校验覆盖。
- **`merge.rs` `join_reads` / `corrected_pair`**：`vec![0; insert]` + 各处
  `.min()`/`saturating_sub()` 保证切片终点 ≤ 长度；重叠分支的
  `i = insert as isize - 1` 从 0 起步即 `-1` 不进入循环；`rc2.len()-c2len` 由
  `c2len ≤ b.len()` 保证不欠位。
- **`tadpole.rs::run`**：`k` 校验 `1..=128`（`Kmer::MAX_K`），`ec-kmer`/`extend`
  的 `-k > 128` 会友好报错而非 panic；`extend_read` 的 rollback 用
  `saturating_add(1)` 防 `% 0` 与 `usize::MAX+1` 溢出。
- **`bbnet` `feed_forward` 23 维**：`debug_assert_eq!(v.len(), 23)` 为 debug 断言，
  构造上 `v` 恒 23 元素，release 无影响；`parse` 已校验输入层 23 维。
- **merge / ec-overlap 的 `--net` 必填**：make-vector 模式下缺 `--net` 由
  `merge()` 库入口（`libs/fq/merge.rs`）统一友好报错；`merge` 命令额外在 CLI 层
  复检，`ec-overlap` 依赖库层守卫，行为一致（友好错误，非 panic）。

## 修复的缺陷（根因模式）

### 数据安全（多输出互斥缺失）

- **`split` / `merge` / `ec-overlap` / `s-filter` 的多个输出可指向同一路径**：
  这四个命令会打开多个 writer（split 的 R1/R2/singles；merge、ec-overlap 的
  merged/outu/ihist；s-filter 的 kept/discard），但此前只对"输出 vs 输入"做了
  `ensure_outfile_distinct`，未校验"输出 vs 输出"。当两个输出指向同一文件时，
  后开的 writer 会截断先写者已写内容，静默产生损坏的输出。这与 `trim-qual`、
  `range` 已有的输出互斥守卫及 asm 审核第五轮（`--outm`/`--outu` 冲突）属同一
  数据安全类别，但 fq 这四个命令此前遗漏。
  修复：在 `cmd/args.rs` 新增共享助手 `ensure_outfiles_distinct`（跳过 `stdout`
  哨兵，两两 `same_path` 校验），四个命令在打开任何 writer 之前调用；并保留
  `trim-qual`/`range` 原有内联守卫（行为不变）。新增回归测试：
  - `command_fq_split_rejects_same_r1_r2_outfile`
  - `command_fq_merge_rejects_outu_same_as_outfile`
  - `command_fq_ec_overlap_rejects_outu_same_as_outfile`
  - `command_fq_s_filter_rejects_outfile_same_as_discard_file`

## 第四轮结论

迁移后的 `fq` 全命令族复核发现 1 类数据安全缺陷（多输出互斥缺失，涉及
`split`/`merge`/`ec-overlap`/`s-filter` 四命令），已修复并含回归测试；其余
Zero-Panic、clap 校验、算法边界、文档一致性经复核均无新问题。新增 4 个集成测试
全部通过，`cargo clippy` 对本轮改动无新增告警。

---

# 第五轮（2026-08-13）：`clean`/`filter` 的 `--stats` 输出与库边界复核

第四轮之后，对 `fq` 全命令族再做一轮纵深复核。本轮逐行复查了
`libs/fq/trim_adapter.rs` 的 `kmask`/`ktrim`/`qtrim`/`trim_by_amount`/
`detect_poly_*` 边界、`libs/fq/sample.rs`/`split.rs` 的库内部边界，以及
`clean`/`filter`/`norm`/`clump`/`range`/`interleave`/`ec-kmer`/`extend`/
`s-filter`/`sample`/`trim-qual` 等命令行文件的参数校验与输出互斥，并核对
`docs/fq.md` ↔ 帮助文本 ↔ 行为的一致性。

## 复核确认（安全不变量，经核验无需修复）

- **`trim_adapter.rs` `kmask` 的 `marked[lo..hi]` 切片**：`lo = i.saturating_sub(k-1-
  trim_pad)`、`hi = (i + trim_pad + 1).min(n)`，因 `i ≥ k-1` 且 `hi ≥ lo+1` 恒成立，
  `fill` 不会越界；`maskfullycovered` 的 `fill(true)`/`fill(false)` 同界安全。
- **`Masks::new(k)` 的 i64 移位**：`1i64 << (2k)` 与 `x2 << (2k-2)` 需 `k ≤ 31`；
  `clean`/`filter`/`norm`/`clump` 的 clap/execute 均已限 `k ∈ 2..=31`，`ec-kmer`/
  `extend` 由 `tadpole::run` 校验 `k ≤ 128`，均无移位越界。
- **`trim_adapter.rs` `test_right_window`/`test_optimal`/`avg_quality` 的质量索引**：
  phred 索引均 `min(127)` 或直接对 128 表取值；`change_quality` 的 `qual[i]` 依赖
  `SeqReader` 强制 seq/qual 等长（与第三轮结论一致）。
- **`sample.rs` 库**：`--bases` 负数由 execute 拒绝；`remaining==0` 提前 break 防
  除零；`target` 为负时 `prob<0` 恒不选中，无越界。
- **`range.rs` 区域切片**：`end < start || end > len` 时 `bail!`，`seq[start-1..end]`
  不会越界；`--cache` 用 `NonZeroUsize` 排除 0。
- **`interleave`/`to-fa`/`norm`/`clump`/`ec-kmer`/`extend`/`sample`/`s-filter`** 的
  `-o` 同输入保护均已在第四轮及此前就位。

## 修复的缺陷（根因模式）

### 数据安全（辅助输出路径未校验）

- **`clean` / `filter` 的 `--stats` 可覆盖输出或输入文件**：`trim_adapter` 在全部
  记录写完后才调用 `write_stats` 落盘 stats 文件，但 `execute` 只对 `-o` 做了
  `ensure_outfile_distinct`，未校验 `--stats` 路径。当 `--stats` 与 `-o` 或输入同
  路径时，会先用 `-o` 写入修剪/过滤结果再被 stats 文本覆盖（或直接破坏输入），
  静默产生数据损坏。这与第四轮"多输出互斥"属同一数据安全类别，但 `--stats` 这一
  辅助输出此前遗漏。
  修复：`clean.rs`/`filter.rs` 在打开任何 writer 前，对 `--stats` 依次做
  `ensure_outfile_distinct`（防覆盖输入）与 `ensure_outfiles_distinct`（防与 `-o`
  互斥）。新增回归测试：
  - `command_fq_clean_rejects_stats_same_as_outfile`
  - `command_fq_filter_rejects_stats_same_as_outfile`

## 第五轮结论

本轮纵深复核确认了 `trim_adapter` 的 `kmask`/`Masks` 移位、质量修剪索引、`sample`/
`range` 边界等安全不变量无需修复；新发现 `clean`/`filter` 的 `--stats` 辅助输出
路径未做互斥校验（数据安全类别），已修复并含 2 个回归测试。其余命令的 Zero-Panic、
clap 校验、算法边界、文档一致性经复核均无新问题。新增测试全部通过，`cargo clippy`
对本轮改动无新增告警。

---

# 第六轮（2026-08-13）：`clean`/`filter` 的 `--hamming-distance` 上界

第五轮之后对 `fq` 全命令族再复查。重点核对了 `trim_adapter` 参考表构建
`add_kmer` 的 hdist 递归复杂度，以及各命令参数之间的关系校验
（`k`/`mink`/`hdist`/`trim_pad` 等）。

## 复核确认（安全不变量，经核验无需修复）

- **`add_kmer` 的 hdist 递归**：`k` 由 clap 限 `2..=31`；`mink > k` 时短路 short
  k-mer 分支，`min_len = 1.max(k.min(mink))` 不会越界。
- **`--qtrim-window 0` + `--qtrim w`**：`execute` 将 `qtrim_window=0` 视为禁用窗口
  模式，`process_pair` 走 `test_optimal` 分支，不进入 `test_right_window`，无除零。
- **`--trim-pad`/`--mask-fully-covered` 依赖 `--mask-kmers`**、`--mask-kmers` 依赖
  `--ref`：第四轮已加守卫并含回归测试。

## 修复的缺陷（根因模式）

### 参数校验缺失（指数级资源耗尽）

- **`clean` / `filter` 的 `--hamming-distance`（hdist）无上界**：`add_kmer`
  （`libs/fq/trim_adapter.rs`）对参考序列按 `(4*k)^hdist` 枚举单碱基替换变体，
  递归深度即 `hdist`。`hdist` 无上限时（如 `--hamming-distance 5`，k=31 约 290
  亿次调用）参考表构建呈指数增长，导致资源耗尽/近乎死循环（非 panic 但不可接受）。
  修复：`clean.rs`/`filter.rs` 在构建表前校验 `hdist ∈ 0..=3`，越界友好报错；
  帮助文本与 `docs/fq.md` 标注 `0..=3`。新增回归测试：
  - `command_fq_clean_rejects_hamming_distance_above_limit`
  - `command_fq_filter_rejects_hamming_distance_above_limit`

## 第六轮结论

本轮新发现 `clean`/`filter` 的 `--hamming-distance` 无上界导致的指数级资源耗尽
风险（参数校验缺失类别），已修复并含 2 个回归测试，帮助文本与文档同步。其余命令
的 Zero-Panic、clap 校验、算法边界、文档一致性经复核均无新问题。新增测试全部通过，
`cargo clippy` 对本轮改动无新增告警。

---

# 第七轮（2026-08-13）：`kmask` 的 `--trim-pad` 极大值溢出

第六轮之后对 `fq` 全命令族再复查。重点核对了 `trim_adapter` 中 `kmask` 对
`--trim-pad`（usize，无上界）的算术边界。

## 复核确认（安全不变量，经核验无需修复）

- `kmask` 的 `minus = k.saturating_sub(1).saturating_sub(trim_pad)` 与右端
  `lo = i.saturating_sub(trim_pad)` 本就用 saturating 运算，安全。
- `ktrim`/`qtrim`/`trim_by_amount` 的 `min`/`saturating_sub` 对极大长度阈值
  （`minlen`/`maxlength`/`force_trim_*`/poly 阈值）均钳制，无越界。

## 修复的缺陷（根因模式）

### Zero-Panic（usize 溢出）

- **`kmask` 对 `--trim-pad` 极大值（接近 `usize::MAX`）溢出 panic**：
  `libs/fq/trim_adapter.rs` 的 `kmask` 用普通加法 `plus = opts.trim_pad + 1` 与
  `i + opts.trim_pad + 1`（左端短 k-mer），`--trim-pad` 为 usize 且 clap 无上界，
  传 `18446744073709551615` 时 `trim_pad + 1` 在 debug 下溢出 panic。修复：改为
  `saturating_add`（`plus`、主循环 `hi`、左端 `hi` 三处），极大 `--trim-pad` 变为
  掩码整条读的无崩溃行为。新增回归测试
  `command_fq_clean_kmask_huge_trim_pad_no_overflow`。

## 第七轮结论

本轮新发现 `kmask` 对极大 `--trim-pad` 的 usize 加法溢出（Zero-Panic 类别），已
修复并含回归测试。其余命令的 Zero-Panic、clap 校验、算法边界、文档一致性经复核
均无新问题。新增测试全部通过，`cargo clippy` 对本轮改动无新增告警。

---

# 第八轮（2026-08-13）：收敛复核

第七轮之后对 `fq` 全命令族做最终收敛复核：复查 `trim_adapter` 其余用户可控
usize 参数（`minlen`/`maxlength`/`force_trim_*`/`qtrim-window`/poly 阈值）的
算术与索引、`ktrim`/`count_set_kmers`/`qtrim` 路径，以及各命令的文档 ↔ 帮助 ↔
行为一致性；并全量运行测试与 clippy。

## 复核确认（安全不变量，经核验无需修复）

- `trim_to_position`/`trim_by_amount` 的 `saturating_sub`/`min` 对极大长度阈值
  均钳制；`process_pair` 的 forceTrim 用 `i64` 运算（`as i64` 回绕但不 panic），
  `right = (len-b-1).max(0)` 恒非负。
- `test_right_window` 的 `window as u32` 截断仅在 `window ≤ qual.len()`（小窗口）
  时到达，极大 `--qtrim-window` 走 `qual.len() < window` 分支返回 0，无截断误判。
- `ktrim`/`count_set_kmers` 不依赖 `trim_pad`，`k ∈ 2..=31` 限位下 i64 滚动安全。
- `ec-kmer`/`extend` 的 `-k` 由 `tadpole::run` 校验 `1..=128`；`clean`/`filter`/
  `norm`/`clump` 的 `-k` 由 execute 校验 `2..=31`。
- 全部 16 个子命令的 `-o` 同输入保护、多输出互斥、`--stats`/`--discard-file`/
  `--outu`/`--ihist`/`--outfile-2` 等辅助输出互斥均在第四至七轮校验并含回归测试。

## 第八轮结论

本轮收敛复核未发现新问题。`cargo fmt` clean；`cargo clippy --all-targets` 对本次
改动涉及文件（`clean.rs`/`filter.rs`/`trim_adapter.rs` 及测试）无新增告警；全量
`cargo test` 26 个测试二进制全部通过（含本轮新增回归测试）。`fq` 全命令族审核至此
收敛，无再发现问题。

---

# 第九轮（2026-08-13）：`test_optimal`/`avg_quality` 质量索引越界

第八轮之后对 `fq` 全命令族做下一轮纵深复核，重点复查 `clean`/`filter` 质量修剪
路径中所有以质量值为索引访问固定长度概率表的位置，以及 `merge`/`norm`/`overlap`
的质量索引。

## 复核确认（安全不变量，经核验无需修复）

- **`merge.rs` `probability` / `expected_mismatches`**：`aqual[i]`/`bqual[i]` 均
  `.min(PROB_CORRECT4.len()-1)`（钳到 59）再索引 60 项表，越界安全；
  `expected_tip_errors` 的 `limit` 由 `max_bases.min(quals.len())` 界定，索引
  `bases[i]/quals[i]` 不越界。
- **`norm.rs` `codes[b as usize]`**：`b` 为碱基字节 0..=255，`codes` 表按字节寻址，
  安全；`cov[...]` 索引由 `covlast`/`above_limit` 边界保证。
- **`trim_adapter.rs` `table[(kmer & 0xFF) as usize]`**：`0xFF` 掩码保证索引 0..=255，
  不越界。
- **`merge`/`ec-overlap` 的 `to_phred` 转换质量**：为 0..=93 的 phred 值，即便如
  此 `probability`/`expected_mismatches` 仍钳到 59，双保险。

## 修复的缺陷（根因模式）

### Zero-Panic（质量索引越界）

- **`test_optimal` / `avg_quality` 对 ≥128 的 phred 质量值索引越界 panic**：
  `libs/fq/trim_adapter.rs` 中 `test_optimal`（`--qtrim rl` 的 bbduk testOptimal
  路径）与 `avg_quality`（`--min-avg-quality` 路径）直接
  `prob[q as usize]` / `prob[qual[i] as usize]` 索引 128 项错误概率表。`make_read_buf`
  用 `q.saturating_sub(quality_base)` 得 phred 值，畸形/二进制输入的质量字节
  （如 0xe1=225，减 33 后为 192）会越界 panic（debug 下崩溃），违反 Zero-Panic。
  修复：两处改为 `q.min(127)` / `qual[i].min(127)`，与同文件 `expected_errors`
  （line 906）及 `clump.rs:620` 的既有 `min(127)` 钳制一致；≥127 的 phred 质量
  取表末项（约 0 错误概率），行为无崩溃。新增回归测试：
  - `command_fq_clean_qtrim_high_quality_byte_no_panic`（`--qtrim rl` 走
    `test_optimal`）
  - `command_fq_clean_min_avg_quality_high_quality_byte_no_panic`（
    `--min-avg-quality` 走 `avg_quality`）
  两测试均以 0xe1 质量字节的原始 FASTQ 输入验证不再 panic。

## 第九轮结论

本轮新发现 `clean`/`filter` 质量修剪路径 `test_optimal`/`avg_quality` 对 ≥128
phred 质量值的概率表索引越界（Zero-Panic 类别），已修复并含 2 个回归测试；其余
`merge`/`norm`/`overlap`/`trim` 的质量索引经核验均已钳制或表寻址安全。全量 fq
测试通过，`cargo fmt` clean，本次改动文件（`trim_adapter.rs`、`cli_fq_clean.rs`）
在 `cargo clippy --all-targets` 下无新增告警。

> 旁注（非 fq 范畴）：`tests/cli_asm_unitig.rs::command_asm_unitig_gfa_no_dangling_links`
> 在全量串行运行下偶发失败（单独运行稳定通过），属 `asm` 命令 GFA 输出的非确定性
> 悬空链接问题，应在 `asm` 审核中跟踪，与 `fq` 无关。

---

# 第十轮（2026-08-13）：独立交叉复核，收敛

第九轮修复后，本轮用全新视角独立复核 `fq` 全命令族，并交叉核验外部审核提出的
疑点。

## 复核确认（疑点均排除，经核验无需修复）

对独立审核提出的若干"疑似缺陷"逐条对照现行代码核验，全部为误报或已覆盖：

- **`merge` 命令 `--buckets` 除零**：`merge` 命令无 `--buckets` 参数（分桶仅属
  `clump`，其 `--buckets` 已在 execute 校验 `1..=4096`，见 `clump.rs:132-135`）。
- **`clump` `--max-isoforms` / `sample` `--target` / `norm` `--min-count`**：这些参数
  名在 fq 命令族中不存在（`sample` 用 `--bases` 且必填；`--min-count` 属 `s-filter`，
  类型 u64 无除法）；grep 确认无此类参数、无 `.expect()` 调用。
- **`ec-kmer` `-k` 无上界**：`tadpole::run` 统一校验 `1..=128`（`Kmer::MAX_K`），
  超限友好报错，已覆盖。
- **`bbnet::feed_forward` 索引越界 / i32→usize 溢出**：`parse` 校验 `dims.len()>=2`、
  输入层 23、各层非空、权重维度一致；`feed_forward` 的 `values` 恒有 ≥2 层且末层非
  空，`last()/first()` 均有保障。无 i32→usize 越界路径。
- **`merge.rs` 库 `join_reads`/`corrected_pair` 空数组访问**：`process_pair` 对
  `seq<2` 短读提前返回，切片均有 `.min()`/`saturating_sub()` 钳制，与第四轮结论一致。
- **`extend` `.expect()` 非必填参数**：`el`/`er` 用 `unwrap_or(100)`，无 `expect()`。
- **`to-fa`/`interleave` `-o` 覆盖输入**：两命令 `-o` 均有 `ensure_outfile_distinct`
  保护（首轮/第四轮已就位）。

## 结论

本轮独立交叉复核未发现任何新的确定缺陷；外部审核提出的全部疑点均核实为误报或已
被先前各轮覆盖修复。全部 fq 集成测试二进制（`cli_fq`、`cli_fq_clean`、
`cli_fq_clump`、`cli_fq_ec_kmer`、`cli_fq_ec_overlap`、`cli_fq_extend`、
`cli_fq_filter`、`cli_fq_merge`、`cli_fq_norm`、`cli_fq_range`、`cli_fq_sample`、
`cli_fq_s_filter`、`cli_fq_split`、`cli_fq_trim_qual`）全部通过（含第九轮新增 2 个
回归测试）；`cargo fmt` clean；本次改动文件（`trim_adapter.rs`、`cli_fq_clean.rs`）
在 `cargo clippy --all-targets -- -D warnings` 下无新增告警。

`fq` 全命令族审核至此收敛：第九轮修复 Zero-Panic（质量索引越界）后，第十轮未再
发现新问题。
