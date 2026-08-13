# fairy（fairy-prime）：FracMinHash 稀疏采样 + 宏基因组 coverage（源码分析）

> 2026-08-13 整理，基于本地 `fairy-prime/`（v0.5.8，2024 Microbiome，
> bluenote-1577/fairy，从 sylph fork）。功能：多样本宏基因组 MAG binning
> 的 contig coverage 计算，替代 all-to-all read alignment（BWA/minimap2），
> 声称快 100-1000×。对应 pgr 语境：`pgr fq norm` 大数据量方案调研中
> "reads 采样 vs k-mer 采样"的参考——fairy 是典型的 **k-mer 采样**路线。

## 1. 概况

- **定位**：`fairy sketch`（reads → 每样本 `.bcsp` 索引）+
  `fairy coverage`（contigs × 样本 → coverage 矩阵，MetaBAT2/MaxBin2/
  SemiBin2 兼容）。作者明确 caveat：不适于单样本 binning、PacBio HiFi、
  菌株分辨组装。
- **依赖**：needletail（FASTA/Q 解析）、rayon（并行）、fxhash、
  serde+bincode（`.bcsp` 落盘）、scalable_cuckoo_filter（近似去重）、
  statrs（泊松/伽马）、fastrand（bootstrap）、memory-stats（ram-barrier）、
  simple_logger（stderr 日志）、human-sort（输出样本列自然排序，
  contain.rs:37）；musl 静态编译时切 jemalloc（main.rs:13-15）。
  Cargo.toml 还声明了 flate2/zlib-ng、nalgebra、rand、regex，但 10 个源码
  文件中未见直接使用，疑为历史遗留。
- **入口**：`main.rs` 把 `coverage` 恒以 pseudotax=true 调用
  （`contain(args, true)`）→ 实际走 pseudotax 分支；sketch 的 dedup
  默认关闭（见 §8 quirk）。

## 2. 架构

| 模块 | 行数 | 内容 |
|---|---|---|
| `sketch.rs` | 853 | sketch 主流程、read dedup、序列化 |
| `contain.rs` | 1202 | coverage 主流程、ANI/λ 估计、矩阵输出 |
| `seeding.rs` | 209 | FracMinHash 标量实现（滚动 2-bit） |
| `avx2_seeding.rs` | 266 | 同上 AVX2 4 通道版（仅 x86_64） |
| `types.rs` | 189 | SequencesSketch / GenomeSketch / AniResult |
| `inference.rs` | 121 | λ 的矩估计/二分（nb 路径） |
| `cmdline.rs` | 126 | clap 参数 |
| `constants.rs` | 14 | 共享常量：`.bcsp`/`.bcdb` 后缀、CUTOFF、MIN_ANI、MAX_DEDUP_*、VAR_CUTOFF 等 |
| `main.rs` | 28 | 入口；coverage 恒以 `contain(args, true)` 走 pseudotax |

## 3. FracMinHash 采样（seeding.rs）

- **2-bit 编码**：A=0、C=1、G=2、T=3；反向互补 `nuc_r = 3 - nuc_f`
  （A↔T、C↔G）；canonical = `min(f, r)`。注意编码与 khmer
  （A=0,T=1,C=2,G=3）不同。`BYTE_TO_SEQ` 表里 **'U' 也映射为 3**（U 当 T 处理），
  小写字母一并覆盖；**非 ACGT（含 N）一律得 0(A)**——含 N 的 read 会生成
  人为的 A-kmer。
- **哈希**：`mm_hash64` = murmur64 最终化变体（源码带
  `//TODO this is bugged. Fix after release` 注释）。标量版
  `!key.wrapping_add(key << 21)` 因 `!` 优先级低于方法调用，实为
  `!(key + (key<<21))`，与 AVX2 版 `key = key+(key<<21); key = key ^ cmpeq(key,key)`（即取反）
  结果**一致**——"bug"注释只是作者对哈希质量的顾虑，不造成标量/AVX2 分歧；
  若有问题应指"对 canonical min(f,r) 直接哈希"这一层。
  **重要：fairy 直接对 canonical 2-bit 整数做 murmur64**（`mm_hash64(min(f,r))`），
  而不是先编码成字节再哈希——这正是 pgr `libs/hash.rs::seq_fracminhash`（L650-656）
  刻意避免的（raw 2-bit 值结构化、需先 `.to_le_bytes()` 再 rapidhash）。
- **哈希族混用**：fairy 有两套哈希：`mm_hash64`（u64→u64，采样用）+ 可逆版
  `rev_hash_64`（seeding.rs:18-52，逐步骤逆推 murmur 最终化以反解原始 kmer，
  未被调用，属死代码）；另有 `mm_hash`（usize 版）+ `MMHasher`/`MMHashSet`
  （types.rs:74-99，类型为 `HashMap<K,V,MMBuildHasher>`）。**genome sketch 的
  去重 set 用 murmur 系 `MMHashSet`（sketch.rs:401）**，而 **read 计数用
  `FxHashMap`**（字段声明的类型）——两类结构哈希族不一致，是 fairy 的既有实现
  细节。注意：字段虽声明为 `FxHashMap`，但 `SequencesSketch::new`（types.rs:137）
  与 `sketch_sequences_needle`（sketch.rs:798）新建表时用的是 `std HashMap::default()`
  （默认 RandomState 哈希器），`from_enc`（types.rs:140）才真正构造 `FxHashMap`。
  另有一处易混淆的死代码 `_get_kmer_identity`（contain.rs:1180），它按
  `count==1` 的比例估计 read 的 kmer identity，与「反解 kmer」无关。
- **`mm_hash`（usize 版，types.rs:62-72）**：与 `mm_hash64` 同构，先把字节
  按 `usize::from_ne_bytes` 重解释再套同一串 murmur 变换；经 `MMHasher`/
  `MMBuildHasher`（types.rs:74-99）供 `MMHashSet<u64>` 作 genome 去重
  （sketch.rs:401）。注意 `MMHasher::write` 是**覆盖式**（`self.hash =
  mm_hash(bytes)`），若被多次 write 会互相覆盖——u64 key 恰好只 write 一次
  故无碍。表与哈希算法来自 Heng Li 的 miniprot（types.rs:1 头注释含
  MIT 声明）。另有 `decode`/`print_string`（seeding.rs:54-76）调试函数，
  `decode` 遇非 0-3 值直接 `panic!`，属死代码。
- **采样**：`threshold = u64::MAX / c`；`hash < threshold` 才保留
  → 采样率 ≈ 1/c。默认 `c=50`（约 1/50，sylph 为 1/200）。
- **滚动**：f 左移 2 位累积、r 右移 + 顶部补补链，与 pgr 现有滚动同构；
  AVX2 版把 read 切成 4 段重叠窗口并行滚动（`extract_markers_avx2`），
  只支持 k=21/31（`2(k-1)=40/60` 硬编码）。

## 4. sketch：read → .bcsp

- **存储**：`SequencesSketch.kmer_counts: FxHashMap<u64, u32>` 全内存；
  落盘前转 `Vec<(u64, u32)>`（`SequencesSketchEncode`，序列化快一个量级）
  + 元数据（c、k、file_name、sample_name、paired、mean_read_length），
  bincode 写 `.bcsp`。每样本一张表，`threads`（默认 3）个样本并行。
  输出文件名基于输入文件 basename（sketch.rs:296-303/352-355）：双端样本
  **带 `.paired` 中缀**（`<name>.paired.bcsp`），单端为 `<name>.bcsp`；
  有 sample_names 时以 sample name 取代文件名。
- **去重（关键）**：
  - **pair marker**：固定 k=16（`Marker=u32`，`k = size_of::<Marker>()*4`），
    `pair_kmer` 对双端 read 各自取**偶数位**碱基拼成一个 16-mer、**奇数位**
    碱基拼成另一个 16-mer，形成两个 `[Marker;2]`
    （`doublepairs.0=[read1偶数,read2偶数]`、
    `doublepairs.1=[read1奇数,read2奇数]`）；`pair_kmer_single` 把单端 read
    按**前半/后半**各取偶数位、奇数位拼成同样的两对。长度不足
    （双端**任一端** <33bp，即 `< 2k+1 = 33`；单端 <66bp，即 `< 4k+2 = 66`）
    或单端 >400bp → 无 marker、不去重。
  - **规则**：对每个采样的 kmer，若 `(km, marker)` 已见过 → 该 kmer
    不计数（`num_dup_removed++`）；否则插入并 `c += 1`。
  - **效果**：`c` ≈ 该 kmer 出现过的**不同 read-pair 数**——完全相同的
    read pair 只贡献一次；部分重叠的 pair 只对未见过 marker 的 kmer 计数。
  - 精确模式用 `FxHashSet`；`--fpr>0`（默认 0.0001，隐藏参数）切
    `ScalableCuckooFilter`。注意 cmdline 默认 `fpr=0.0001`（非 0），故即便
    dedup 开启，默认走的也是**近似 cuckoo filter** 而非精确 `FxHashSet`；
    而 `--no-dedup` 默认 `true` 使 dedup 整体关闭，二者叠加=默认完全不去重。
  - **fpr 只作用于双端路径**：单端路径 `sketch_sequences_needle` 恒走精确
    `FxHashSet`（`dup_removal_lsh_full_exact`），忽略 `--fpr`。
  - 双端路径中 read2 的 kmer 若已在 read1 中被采样（`temp_vec1.contains(km)`）
    则跳过，避免同对 read 内部重复计数。
  - **coverage 对裸 fastq 的内部 sketch 与 sketch 命令相反**：`get_seq_sketch`
    对 raw reads 调 `sketch_sequences_needle(..., no_dedup=false)`——去重**开启**
    （单端精确路径），而 `fairy sketch` 默认关闭。
- **单端上限**：`MAX_DEDUP_COUNT=4`，`c < 4` 才查去重（高拷贝序列的
  计数不再被门控）。双端无上限。
- 其余：`mean_read_length` 逐条 moving average；`MAX_DEDUP_LEN` 常量
  未使用。

## 5. coverage：contig × 样本

- **contig sketch**（`sketch_genome_individual`，每 contig 独立，contain.rs
  恒走此路径）：FracMinHash 采样 → **重复 kmer 去重被 `|| true` 短路禁用**
  （见 §8 quirk，contig 上重复出现的 kmer 并不会被丢弃）→ 仅
  `min_spacing=30` 过滤（代码 `pos - last_pos > min_spacing` 才保留，即
  **间距 ≤ min_spacing 的相邻 kmer 被丢弃**，`--min-spacing` 可调）→
  `genome_kmers: Vec<u64>`；`gn_size` = contig 长度。可预存 `.bcdb`。
  （但本版本 `sketch` 命令**只生成 `.bcsp`，从不生成 `.bcdb`**；coverage
  仅**读取**预存的 `.bcdb`/`.sylqueries` contig sketch——`.syl*` 后缀表明可与
  sylph 的 sketch 互通。同文件的 `sketch_genome` 合并版才真正去重，
  但 coverage 不用它。）
- **查询**：对每个 contig 遍历 `genome_kmers` 在样本表中查 multiplicity
  → `covs` 向量；`contain_count` = 命中数。过滤：`gn_kmers.len() ≥
  min_number_kmers`（默认 8）且 `covs` 非空。
- **ANI**：`naive_ani = (contain/total)^(1/k)`；随后用 λ 校正
  （`ani_from_lambda`：`contain/(1-e^-λ)/total` 再开 k 次方）。
  **输出阈值 0.95**（pseudotax 分支，main 恒走；普通分支为 0.9）
  ——与 wiki 的说法一致，但 0.9 分支实际不可达；阈值可用
  `-m/--minimum-ani`（0-100，短选项由 clap 从字段名 `minimum_ani` 自动派生
  为 `-m`，实测 `./fairy coverage -h` 确认）覆盖。
- **覆盖度估计**（`get_stats`）：
  1. `median_cov` = covs 中位数；median<30 时按 `Poisson(median)` CDF
     < 0.9999999999 剪掉高倍噪声（`max_cov`）；
  2. `full_covs` = 未命中补 0 + 命中的 covs（≤max_cov）；
  3. λ：默认 `ratio_lambda` = `count(mode+1)/count(mode) × (mode+1)`，
     要求 ≥25 个命中、mode+1 存在、两侧计数 ≥ `min_count_correct`（默认 3）；
     备选 `mme`（矩估计）/ `nb`（named `binary_search_lambda`）/ `mle`（零位 +
     Newton-Raphson）；四种估计分别由**隐藏 CLI 开关** `--ratio` / `--mme` /
     `--nb` / `--mle` 选择（cmdline.rs:103-110），并另有隐藏 `--no-ci`（关
     bootstrap CI）与 `--no-adjust`（关 λ 校正，直接用 naive_ani，
     cmdline.rs:111-114、contain.rs:868）；
  4. `final_est_cov` = λ（可估）| median<15 时 `geq1_mean_cov` | median。
- **方差**：对 `full_covs` 前 95% 窗口算（`VAR_CUTOFF=10` 以下不剪）。
- **CI**：100 次 bootstrap（`fastrand::seed(7)` 固定），5-95 分位，
  <50 次成功则输出 NA。
- **pseudotax**（恒生效）：`winner_table`（contain.rs:489-512）把共享 kmer
  按 `final_est_ani` 最高的 genome 重新分配（共享/重复 kmer 只归属一个
  genome，消除跨物种重复计数），再对每个 genome 二次 `get_stats`
  （contain.rs:410-418）。两点关键：
  - **round-1 校正不作用于 round-2 输出**：第一次 `estimate_true_cov`
    （contain.rs:395，恒 `estimate_unknown=true`）把 `final_est_cov` 乘
    `read_length/(read_length-k+1)` 与 `1/((seq_id/100)^k)`；但 winner
    reassign 后的二次结果（`stats_vec_seq_2`）里，`estimate_true_cov` 被
    注释（contain.rs:417）——**最终输出的 cov 矩阵是未经这两项校正的原始
    λ/median**，而 `winner_table` 用的 ANI 又来自校正前的第一轮。短读下
    这约低估 `read_length/(read_length-k+1) × 1/(seq_id/100)^k` 倍
    （如 150bp/k=31/seq_id=99.5：≈1.25×1.167≈1.46×）。
  - `winner_table` 的 bool 标记（`(ani, gn, changed)` 第三项）在
    `get_stats` 里被注释（`// || map[kmer].2`，contain.rs:772），未使用。
  - 被 `min-spacing` 过滤掉的 kmer 若开了 pseudotax，也会被收集进
    `pseudotax_tracked_nonused_kmers`（sketch.rs:420-422）并参与二次分配。
- **`-I/--read-seq-id`（默认 99.5，隐藏）**：kmer identity 按
  `(seq_id/100)^k` 估计（contain.rs:375），`seq_id=99.5,k=31` 时 ≈0.857；
  这是**固定假设**而非从数据估计（`_get_kmer_identity` 那个估计函数是
  死代码），长读/高错配需手动调低。
- **输出**：默认 MetaBAT2 格式（contigName/contigLen/totalAvgDepth + 每样本
  cov/var）。`--maxbin-format` 与 `--aemb-format` 均去掉 contigLen、
  totalAvgDepth、每样本 var 三列：`--maxbin-format`（注意：cmdline 里
  clap 的 long 是 `maxbin-format`，但对应**结构体字段名**是
  `concoct_format`，故源码内写的是 `args.concoct_format`——对外 CLI 是
  `--maxbin-format`，兼容 MaxBin2）保留表头（contigName + 各样本名），
  `--aemb-format` **无表头**只输出每样本 cov。`--full-contig-name` 保留
  空格后的全名。
- **输入分类与一致性检查**：coverage 的每个输入按后缀分类
  （`.bcsp`/`.sylsample` → 样本 sketch；`.bcdb`/`.sylqueries` → contig
  sketch；fasta/fastq → 内部即时 sketch）。要求至少 1 个 contig 源 + 1 个
  read 源，否则报错退出；多个 contig 文件时**每个文件输出一段独立 TSV 表头**
  （warn "not a valid TSV file"）。一致性检查：所有 genome sketch 的 `k`
  必须一致；read sketch 的 `c` 不得大于 genome 的最小 `c`（`get_seq_sketch` 与
  `get_genome_sketches` 双处校验）；sketch 的 `k` 必须与 `-k` 匹配。
  - **`.sylsample` 名不副实**：主循环判定 `is_sketch` 只认 `.bcsp`
    （contain.rs:339），`.sylsample` 虽在初筛被分进样本 sketch 桶，运行时
    却被当 raw fastq 走 `sketch_sequences_needle` 解析 bincode 二进制
    → 大概率 warn "not a valid fasta/fastq" 或被 needletail 误解析。
    相对地，`.sylqueries` 的 genome sketch 能被正确加载。这是
    「接受 `.sylsample`」与「运行时只认 `.bcsp`」的不对称。
  - **`lowest_genome_c` 命名误导**：`get_genome_sketches`（contain.rs:558-562）
    实际取的是各 genome sketch `c` 的**最大值**（`existing < c → 替换`），
    并非"最小"。随后用 `max < args.c` 判错——这比"min < args.c"更宽松，
    多 genome 文件 `c` 不一致时可能漏报。主循环传给 `get_seq_sketch` 的
    `genome_c` 是 `genome_sketches[0].c`（contain.rs:353，该文件第一个
    contig 的 c）。
  - **输出列排序**：contig 行按输入 fasta 顺序（`contig_list_sorted`
    未排序）；样本列用 human-sort **自然排序**（`sort(&mut read_list_sorted)`,
    contain.rs:37，非字典序）；`totalAvgDepth` 是该 contig 各样本 cov 的
    算术平均（未命中记 0，contain.rs:69-77）。
- **reads 检出比例日志**（`estimate_covered_bases`，仅 log）：估算
  `tentative_bases = c × Σkmer_counts × read_len/(read_len-k+1)` 与
  `covered = Σ gn_size × final_est_cov`，输出 `min(covered/tentative, 1)`
  的百分比，只在低错误 reads 下近似准确。
- **rel/seq abundance 未实现**：`AniResult.rel_abund / seq_abund` 恒为
  `None`（`get_stats` 从不赋值）；论文/README 描述的 taxonomic / sequence
  abundance 两列只存在于**死代码** `_print_ani_result`（下划线前缀、未被调用），
  实际输出恒走 `print_cov_matrix`（binning 矩阵）。

## 6. 内存与并行

- sketch：每样本一张 FxHashMap（论文：土壤样本约 4GB/样本；1/50 采样下
  内存 ≈ 采样 kmer 数 × (8+4) 字节 + 哈希表开销）；`--ram-barrier`
  （隐藏）是**软限制**：虚拟内存超限时 `sleep` 阻塞等回收，不保证上限。
- coverage：`step` 决定同一时刻驻留内存的样本 sketch 数（genome sketches
  全量常驻）。非 pseudotax 分支 `step=1`；pseudotax 分支
  `step = threads/2 + 1`。因 main 恒走 pseudotax，**CLI 下默认实际是
  `step=threads/2+1`**，`step=1` 分支不可达；`-s/--sample-threads` 可覆盖。
  论文数字：10 样本索引 9min + 查询 7min vs BWA 40h。

## 7. 对 pgr 的启示

1. **fairy = "k-mer 采样"，不是 "reads 采样"**：它保留约 1/c 的 kmer 并
   只统计这些 kmer 的 multiplicity。pgr norm 的 bbnorm 语义（highpass
   filter，按 read 内 kmer 深度分位判定）**不需要采样**——采样会直接改变
   深度分位的语义；此前讨论的"对 reads 采样"也不是 fairy 路线。
2. **dedup 思路**（pair marker 指纹 + 按 kmer 门控计数）与 pgr 现有
   `fq clump` 精确整对去重是不同抽象层级；pgr 不需要引入。
3. **pgr 已有等价实现，不必新造**：fairy 的 FracMinHash 稀疏采样路线
   在 pgr 里已有现成对应物——`src/libs/hash.rs`：
   - `seq_fracminhash`（L630，DNA canonical k-mer 以 `threshold = u64::MAX/scale`
     采样，与 fairy §3 完全同构；蛋白质走 `is_protein` 分支）、
   - `load_fracminhash`（L664，FASTA→逐记录 sketch）、
   - `mash_sketch_distances` / `set_distances`（L327/L401，含 containment
     与 ANI CI `ani_ci_from_jaccard` L441）、
   - `bottom_k_min_hashes`（L260，Mash bottom-k 累积器）。
   命令入口：`pgr dist frac`（FracMinHash，`--scale`）、`pgr dist mash`
   （Mash bottom-k）、`pgr dist mini`（minimizer）。fairy 对 pgr 的增量价值
   不在"要不要做 FracMinHash"（已做），而在**它把 sketch 用于多样本 coverage/
   丰度矩阵而非距离**这一用例，以及 read-level 去重门控计数。
4. 若 pgr 未来做 coverage/丰度类工具（fairy 的 `coverage` 输出即 MetaBAT2 兼容
   contig×样本矩阵），FracMinHash + FxHashMap + `Vec<(u64,u32)>` bincode
   序列化（`SequencesSketchEncode`，注释称快一个量级）是现成的最小实现模板；
   `c` 与内存线性反比。注意 fairy 的覆盖度估计**不是**简单的 kmer 计数——
   它叠加了泊松剪枝（median<30 时 CDF 剪高倍噪声）、λ 的 ratio/mme/nb/mle
   四种估计（默认 ratio，`ratio_lambda`：`count(mode+1)/count(mode)×(mode+1)`，
   要求 ≥25 命中、mode+1 存在、两侧计数≥`min_count_correct`=3）、100 次
   bootstrap CI、pseudotax 二次分配（`winner_table` + 二次 `get_stats`）、
   `read_length/(read_length-k+1)` 与 `1/(seq_id/100)^k` 校正——这套是 fairy
   独有、pgr 完全没有的统计层，比采样本身更值得参考。
5. **与 pgr 现有设施对照**：

| 项 | pgr 现有 | khmer | fairy |
|---|---|---|---|
| 精确计数 | `KmerTable`（u128+u32）、`count.rs` .pkt sort-merge | 无 | 无 |
| 近似计数 | 无（bbnorm 精确表语义） | CMS（u8/u16 饱和） | FracMinHash 稀疏采样 + u32 multiplicity |
| 判定 | truedepth/depthAL 分位 + toss | median ≥ cutoff（在线） | 中位数 + 泊松剪枝 + λ 校正 |
| 内存 | 外部桶（mem_cap 约束） | 固定但随装载率失真 | 与采样率反比（1/c） |

6. **fairy 不含 graph 构建、对 pgr `asm` 无直接借鉴**：它完全是「稀疏采样
   sketch → 计数 → 查询」这条丰度/覆盖度估计路线，**没有 de Bruijn 或任何图
   结构**（全库无图构建代码），也没有 error-correction。因此对 pgr 的
   `asm`（OLC overlap-layout-consensus，`cmd/asm/{olc,ovlp,layout,contig}.rs`）
   与 `fq` error-correction（`cmd/fq/ec_kmer.rs`、`ec_overlap.rs`）**没有
   可直接迁移的算法**——它们需要的是 kmer 深度 / 图 / 比对，而非稀疏采样。
   对 pgr 的借鉴价值在**稀疏采样 + 查询表**这一范式（对应 `dist frac/mash`），
   而非图算法。若要为 `fq ec` 找"重复去除/计数门控"的现成参考，fairy 的
   `dup_removal_lsh_full_exact`（sketch.rs:583）那套（`(km, [Marker;2])` 组合键
   + 按 kmer 计数门控，避免同一 read-pair 重复计数）是唯一可借鉴的 read-level
   工程技巧。

## 8. 源码 quirks

- `--no-dedup` 默认 `true` 且无反开关（clap SetTrue 语义）→ **去重实际
  默认关闭**，与 README/论文描述的 illumina 去重不符，疑似 0.5.8 有意/无意
  改动。
- `sketch_genome_individual` 里去重条件写成 `if !duplicate_set.contains(&km)
  || true`——`|| true` 使该条件恒真，**contig 级重复 kmer 去重实际被禁用**；
  同文件未使用的 `sketch_genome` 版本才是正确的 `if !duplicate_set.contains(...)`。
  coverage 恒走前者，故 contig sketch 实际只做 min_spacing 过滤。
- `mm_hash64` 带 "TODO this is bugged" 注释；AVX2 与 scalar 等价，
  问题在 canonical 编码层。
- 非 ACGT 一律编码为 A，含 N 的序列会产生人为 kmer。
- `pair_kmer` 的 k=16 与主 k=21/31 无关，仅作 read 指纹。
- coverage 经 main.rs 恒走 pseudotax 分支 → 默认 ANI 阈值实际是 0.95、
  `step=threads/2+1`；`contain()` 里 0.9 阈值、`step=1` 的分支在 CLI 下均不可达。
- **冗余去除失效**：`derep_if_reassign_threshold` 调用被注释，`-R/
  --redundancy-threshold` 目前不生效。
- **AVX2 只支持 k=21/31**：`extract_markers` 一旦检测到 AVX2 即无条件调用
  avx2 版，k 非 21/31 时 `panic!()`（`use_40` 分支）；非 x86/无 AVX2 走标量
  `fmh_seeds` 则无此限制。
- **AVX2 对短序列的行为分歧**：read 路径的 `extract_markers_avx2` 要求
  `len ≥ k+1`（avx2_seeding.rs:42），标量版只要 `len ≥ k`（seeding.rs:93）——
  长度为 k 的 read 在 AVX2 机器上不产出 kmer；contig 路径
  `extract_markers_avx2_positions` 更需 `len ≥ 2k`（avx2_seeding.rs:160，
  标量版 `fmh_seeds_positions` 只要 `≥k`，seeding.rs:156），故长度在
  [k, 2k) 的 contig 在 AVX2 机器上会产出**空 sketch**（coverage 直接跳过），
  标量机器则正常。这会造成同一输入在不同机器上结果不同。
- **coverage 最终 cov 未做 round-2 校正**：见 §5 pseudotax——winner
  reassign 后的 `estimate_true_cov` 被注释（contain.rs:417），输出 cov
  未经 `read_length/(read_length-k+1)` 与 `1/(seq_id/100)^k` 校正。
- **`.bcdb` 只读不写**：本版本 `sketch` 仅生成 `.bcsp`，contig sketch 需
  外部提供（兼容 sylph 的 `.sylqueries`）。
- **`--maxbin-format` 命名错位**：clap 的对外 long 是 `maxbin-format`（README
  亦如此），但结构体字段名却叫 `concoct_format`——源码里到处是
  `args.concoct_format`，容易误以为 CLI 是 `--concoct-format`（实际不是）。
- **大量下划线前缀死代码**：`_print_ani_result`（含 rel/seq abundance 列的
  ANI 明细表）、`_print_header`、`_get_sketches_rewrite`、
  `_derep_if_reassign_threshold`、`_get_kmer_identity`、`_ani_from_lambda_moment`
  等均未被调用；本版本实际只走 `print_cov_matrix`（binning 矩阵）路径，
  abundance 输出功能整体未落地。
- **`binary_search_lambda` 名为二分实为线性扫描**（inference.rs:26-99）：注释
  掉的二分代码都在，实际用的却是 0..10000 步的均匀线性扫描找最优点（还带
  `dbg!` 调试输出）；`nb` 路径即走此函数。
- **更多死代码**：`MultGenomeSketch` 类型（types.rs:163）整文件无引用；
  `PAIR_REGEX`（constants.rs:1）从未使用（与 `MAX_DEDUP_LEN` 同为死常量）；
  `rev_hash_64`（seeding.rs:18）在 crate 内无调用；`sketch_genome` 合并版
  （sketch.rs:443）亦无调用者（coverage 只用 `sketch_genome_individual`）。
- **`contain` 恒走 pseudotax 的另一后果**：`_print_ani_result` 的
  pseudotax 分支会 `.unwrap()` 取 `rel_abund`/`seq_abund`，但它们恒为 `None`；
  若该死代码被启用，pseudotax 路径会直接 panic——目前它不被调用故无碍。
