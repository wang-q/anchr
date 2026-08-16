# SKESA / skesa-rs：de Bruijn 图短读组装器（源码分析）

> 2026-08-13 整理、2026-08-17 对照源码逐项复核修订，纯源码分析。SKESA（`SKESA-master/`，C++）是 NCBI 的微生物
> 基因组 de-novo 短读组装器，论文 [SKESA: strategic k-mer extension for
> scrupulous assemblies](https://doi.org/10.1186/s13059-018-1540-z)（Genome
> Biology 2018）；`skesa-rs/` 是基于 SKESA v2.4.0 / SAUTE v1.3.0 快照（commit
> `27caba2`，2024-10-11）的**逐位忠实 Rust 移植**（henriksson-lab/rustification
> 项目），并追求与 C++ 输出字节级一致。两者共享同一算法族，是 pgr 做
> k-mer 计数 / de Bruijn 图遍历的**首选参考**。
> **与 OLC 的连接（2026-08-12）**：pgr 的 `asm olc`（多 k unitig 层 OLC，
> `design/asm-olc.md`）已落地，SKESA 的 fork 过滤 / 可逆性 / 迭代多 k 语义
> 直接映射其 v1 待决项（见 §7.1）。
> **与 asm multik 的连接（2026-08-14）**：`asm multik`（metaMDBG 式 unitig
> 反馈，`design/asm-multik.md`）已实现并过 G37。SKESA 与它是"多 k 迭代"的
> 两个不同实现族——SKESA 每轮**重建图 + contig 引导 + 移除已用 reads**，
> multik 是**全量 reads + unitig 序列反馈**；两者的 fork/bubble 语义对比
> 见 §7.2（本次重读补充）。

## 1. 概况

- **定位**：SKESA 面向 Illumina 短读（单/双端），用**保守启发式**在重复区
  断裂，换取序列质量；k 从 mate 长度一直增到 insert size，兼顾 N50。
- **确定性**：同输入（含 read 顺序）下输出 contig 的顺序/方向**确定**——
  依赖排序 + 稳定启发式，不依赖多线程调度（这点对 pgr 做字节级一致性很关键）。
- **版本澄清**：本地 `SKESA-master/`（C++）实际是 **SKESA 2.5.1**
  （`skesa.cpp:397` 的 `--version` 打印 `SKESA 2.5.1`）；而 Rust 移植
  `skesa-rs-main/` 的 README 声称基于 **SKESA v2.4.0 / SAUTE v1.3.0 快照**
  （commit `27caba2ed...`，2024-10-11，见 `README.md:5,141`），crate 自身版本
  `0.2.2`（`Cargo.toml`）。即：仓库内 C++ 源码比 Rust 移植所依据的快照更新一档；
  但 Rust CLI 的 about 文本（`cli.rs:58`，写 "Original: SKESA 2.5.1"）与
  `--version` 分支（`cli.rs:755` 打印 `SKESA 2.5.1`）都**自报 2.5.1**，与更新的
  C++ 对齐——即 Rust 版"CLI 版本号跟新 C++、内部算法语义跟旧快照 2.4.0"。
  对比语义时以 **2.4.0 快照**为准（README 明示），C++ 2.5.1 仅作行为参考。
- **语言/构建**：
  - C++：Boost + gcc，`make`（NGS 版）/ `make -f Makefile.nongs`（文件版）；
    variant（`boost::variant<LargeInt<1>..LargeInt<16>>`）做运行时多精度 k-mer。
  - Rust：Cargo，依赖仅 `clap`（optional）、`noodles`、`flate2`、`rayon`；
    默认不编译 CLI（库形式），`--features cli` 才引 clap。
- **两者关系**：skesa-rs 是**逐函数移植**（含复刻 bug 以求可复现），C++ 与
  Rust 文件基本一一对应；Rust 在 `--cores 4` 基准上与 C++ 时间/RSS **基本持平**
  （README 快照 wall 0.972x、RSS 1.000x，contigs 输出 SHA-256 一致）。
- **未移植**：SAUTE / saute-prot / gfa-connector 完整对等、SRA 输入（Rust 显式拒绝）。

### 1.1 核心算法与流程总览（先读这节）

**一句话**：reads 计数（排序外部归并 或 bloom+分块哈希，§3.2/§3.3）滤出可信 k-mer
建 canonical de Bruijn 图（§3.4）→ **保守扩展**：四层后继过滤 + 可逆性检查 +
"前驱唯一且可回退"不变量，在 fork 处宁可断裂（§4）→ **迭代多 k**：k 从 21 增到
max_kmer，每轮用旧 contig 标 visited 引导新种子、并清理已用 reads，渐进解重复
（§5）→ paired-end 连接补长插 → 输出。

```
reads → 计数（排序归并 | --hash_count: bloom+分块表）→ 首图 @ k=21
      （高覆盖时自动抬 min_count = coverage/50）
  → 保守组装 seeds（fork 处断，宁断勿错）
  loop k: 21 → max_kmer（steps=11 轮等距、一律取奇）:
    ConverToSContigAndMarkVisited   # 旧 contig 用到的 k-mer 标 visited
    GenerateNewSeeds                # 只从未 visited 区组装新种子
    ExtendContigsJob                # 从旧 contig 边缘 k-mer 继续扩展
    CleanReads                      # 剔除完全落在 contig 内的 reads（逐轮递减）
  若 insert N50 > 1.5×max_kmer: 长插三轮 [1.25×max, 中点, insert_N50]
  connect_pairs（3×N50 限内连 mate）→ 输出（--min-contig 200）
```

| 核心块 | 机制 | 详见 |
|---|---|---|
| k-mer 计数 | 排序计数（canonical、外部归并、确定性）或 blocked bloom filter（128B 块）+ 分块哈希表（`--hash_count`，省内存） | §3.2、§3.3 |
| de Bruijn 图 | canonical k-mer 排序数组；Node 偶正奇负 + `Index()=m_node/2-1`；count 打包 total(32)+branch(8)+plus-fraction(16) 一个 u64 | §3.4 |
| 保守扩展 | 四层后继过滤（低丰度 fork `≤0.1×Σ`、GGT/ACC 链特异噪音、不可扩展 fork <100bp、strand 平衡）+ 可逆性检查 + 前驱唯一不变量；fork 处断 | §4 |
| 迭代多 k | k: 21→max_kmer（steps=11）；旧 contig visited 引导 + `CleanReads` 移除已用 reads（read 集逐轮递减） | §5 |
| paired 连接 | 抽样 10000 对估 insert N50，`3×N50` 限内连 mate；N50 > 1.5×max_kmer 时另加长插三轮 | §5 |

## 2. C++ 仓库结构（SKESA-master/）

| 文件 | 作用 |
|---|---|
| `skesa.cpp` | 主入口，参数解析（boost::program_options）、流程编排 |
| `Integer.hpp` / `LargeInt.hpp` / `LargeInt1/2.hpp` | 大整数 / 变长 k-mer 编码 |
| `KmerInit.hpp` / `Model.hpp` | k-mer 初始化、2-bit 编码表 |
| `concurrenthash.hpp` | **并发 k-mer 计数**：blocked bloom filter + 分块哈希表 |
| `counter.hpp` | **排序计数** `CKmerCount`（variant of 有序 vector） |
| `kmercounter.cpp` | 独立 k-mer 计数工具（`kmercounter` 子命令） |
| `DBGraph.hpp` | de Bruijn 图：`CDBGraph`（排序计数版）、`CDBHashGraph` |
| `graphdigger.hpp` | **图遍历 + contig 组装**（保守扩展、fork 解析） |
| `assembler.hpp` | **迭代组装编排**（多轮 k、paired-end 连接、insert 估计） |
| `guidedassembler.hpp` / `guidedgraph.hpp` / `guidedpath_naa.hpp` | SAUTE 目标富集组装 |
| `saute.cpp` / `saute_prot.cpp` | SAUTE 入口（核酸/蛋白） |
| `gfa.hpp` / `gfa_connector.cpp` | GFA 图输出 / contig 连接成图 |
| `glb_align.cpp/hpp` / `nuc_prot_align.hpp` | 全局/核酸蛋白比对 |
| `readsgetter.hpp` | FASTA/FASTQ/gzip/SRA 读入 + 适配器裁剪 |

Rust `src/` 同名对应：`concurrent_hash.rs`、`sorted_counter.rs`/
`counter.rs`/`flat_counter.rs`、`db_graph.rs`、`graph_digger.rs`、
`assembler.rs`、`guided_*.rs`、`snp_discovery.rs`、`linked_contig.rs`、
`paired_reads.rs`、`clean_reads.rs`、`glb_align.rs`、`gfa.rs`。

## 3. 核心数据结构

### 3.1 大整数 k-mer（`LargeInt.hpp` / Rust `large_int.rs` + `kmer.rs`）

- 2-bit 编码，k-mer 长度 → 精度 `precision = (kmer_len+31)/32`（用几个 u64
  存储），最大 16×64=1024 bit = 512 nt。
- C++ 用 `boost::variant<LargeInt<1>..LargeInt<16>>` 做运行时多精度；
  Rust 用 `enum Kmer { K1(LargeInt<1>), ..., K16(LargeInt<16>) }` +
  `macro_rules! define_kmer_enum` 展开 16 个变体的全操作分派（`kmer.rs`）。
- 关键操作：`revcomp`、`shl/shr`、`oahash`（SKESA 自定义哈希，`KmerOaHasher`
  复刻为 Rust `Hasher`）、`resize`（换精度，左截断/补零 + 顶字掩码）。
- 变长 k-mer 的 `LargeInt<N>` 一律内联存储为 `[u64; N]`（Rust `large_int.rs`）；
  `Flat((u64,u64)) / 内联数组 / boxed` 的存储分级属于 `KmerCount::Storage`
  （见 §3.3），并非 k-mer 编码本身。

### 3.2 并发 k-mer 计数（`concurrenthash.hpp` / Rust `concurrent_hash.rs`）

两条计数路线（`--hash_count` 切换）：

**A. blocked counting bloom filter（`CConcurrentBlockedBloomFilter<128>`）**
- 每 `SBloomBlock` 128 字节、`alignas(64)` 缓存行对齐（`concurrenthash.hpp:145`）；
  计数元素按块内位偏移打包（每计数器 2/4/8 bit）。
- 由**两个**哈希值 `(hashp, hashm)` 生成 `hash_num` 个哈希位：`hashp += hashm`
  迭代（`concurrenthash.hpp:77-85`），块内取 `hashp & (elements_in_block-1)`。
- 每块一个 `SAtomic<uint8_t>` 自旋锁；`Insert` 返回
  `eNewKmer / eAboveThresholdKmer / eExistingKmer`，用于**只把达到 min_count
  的 k-mer 灌入真实哈希表**——bloom 过滤是内存控制的关键。
- 计数封顶 `m_max_element = (1<<counter_bit_size)-1`（饱和计数）。

**B. 分块并发哈希表（`SHashBlock<Key,V,BucketBlock=32>`）**
- 每桶一个**小数组（≤32 槽）+ 溢出 forward list**（`concurrenthash.hpp:395-445`）：
  先试哈希指定位置，再线性扫小数组，最后溢出链表——`BucketBlock=32` 折中缓存
  命中与溢出量。
- 每槽状态原子 `eAssigned / eKeyExists`，`Lock/Wait` 实现无锁读 + 自旋写；
  `CDeque` 分块并行初始化大表。
- k-mer 的 (key, count) 打包存储；`count` 在计数期低 32 位 total、高 32 位
  plus-strand（见 §3.3），`(plusf<<48)+(branches<<32)+total`（`concurrenthash.hpp:1380`）。

### 3.3 排序计数（`counter.hpp` / Rust `counter.rs`/`sorted_counter.rs`/`flat_counter.rs`）

- `CKmerCount` 是 `vector<pair<LargeInt<N>, size_t>>` 的 variant，**只存 canonical**
  （kmer 与其 revcomp 中较小的）；排序后二分查找（`lower_bound`）。
- **count 打包**（`counter.hpp:42-44` + `DBGraph.hpp:180`）：
  ```
  低 32 bit: total count（self+revcomp）
  高 32 bit: 计数期=plus-strand count；进 CDBGraph 后重排为:
             [0:31]=total | [32:39]=8bit 分支信息 | [40:47]=未用 | [48:63]=16bit plus-fraction
  ```
- `--memory`（GB）决定**多轮外部归并**：内存预算 → 每轮可装多少元素 →
  分块排序落盘再归并（`counter.hpp` 的 `SortAndExtractUniq`/`MergeTwoSorted`）。
- Rust `KmerCount` 用 `enum Storage` 区分 `Flat/Words2..8/General`，1..8 字内联、
  其余 boxed 兜底（`counter.rs`）；排序用 `rayon::par_sort_unstable_by_key`
  （>10000 并行、否则串行，保证小规模确定性）。`flat_counter.rs`/`sorted_counter.rs`
  是不同精度/计数策略的变体。

### 3.4 de Bruijn 图（`DBGraph.hpp` / Rust `db_graph.rs`）

- **节点编码**：`Node` 包一个 `size_t m_node`；偶数=正链、奇数=负链、0=无效；
  `Index() = m_node/2 - 1` 映射回数组下标（`DBGraph.hpp:102-124`）。
  图里只存 canonical k-mer，`GetNode` 对 `kmer<revcomp` 返回正链节点、否则
  revcomp 的负链节点（`DBGraph.hpp:155-164`）。
- **分支信息**：`GetNodeSuccessors` 读 count 高 32 位的 8bit 分支掩码，负链取
  高 4 位、正链取低 4 位（`DBGraph.hpp:249-267`）；`shifted=(kmer<<2)&max_kmer`
  + nt 查后继。**用打包位快速跳过无后继的碱基**，避免逐个探测。
- **strand 信息**：`PlusFraction = (count>>48)/65535`（`DBGraph.hpp:185-190`），
  供图遍历区分正负链计数。
- visited 用原子 uint8（1=永久占用、2=临时、3=多 contig），多线程安全标色
  （`DBGraph.hpp:204-222`）。
- **哈希版 visited 复用了 count 的高位**（`CDBHashGraph`，`DBGraph.hpp:432,454`）：
  排序版 `CDBGraph` 用独立 `m_visited` 数组（`DBGraph.hpp:301`），而
  `CDBHashGraph` 把 count 的 `[40:47]`（对应 `DBGraph.hpp:180` 的"未用" 8 bit）
  改作 visited 控制位（`eVisited/eTemp/eMulti = 1<<40 / 1<<41 / 1<<42`，
  `SetColor` 用 `mask<<40` 打色、`GetColor` 取 `(count>>40)&0xFF`）——省下
  一个 per-node 数组，是"状态打包进计数"的另一处实例。
- Rust `SortedDbGraph` 完全对应（`db_graph.rs`），另有 `HashNode` 对应哈希计数版 `CDBHashGraph`；两个具体图共用 `DBGraph` trait。

## 4. 图遍历与 contig 组装（`graphdigger.hpp` / Rust `graph_digger.rs`）

核心是**保守扩展 + 在重复区断裂**，用"只沿唯一、可信路径延伸"换质量：

- **fork 类型**（`graphdigger.hpp:93`）：`eNoFork/eLeftFork/eRightFork/
  eLeftBranch/eRightBranch/eSecondaryKmer`——记录左右分支与次生 k-mer。
- **后继过滤**（`FilterNeighbors` / `FilterLowAbundanceNeighbors`，
  `graphdigger.hpp:1769-1887`），按序：
  1. **低丰度 fork 剔除**：`abundance(后继) <= fraction × Σabundance` 的删除
     （`fraction` 即 `--fraction` 默认 0.1，噪音/信号比上限）；`LowCount()==1`
     且首后继丰度>5 时，把丰度==1 的尾巴删掉。
  2. **strand 特异的 Illumina 噪音**（`GGT→GG[ACG]` 现象）：对以 `GGT` 结尾的
     后继，用 `abundance×(1-PlusFraction)` 与 `fraction×am` 比较剔噪——正负链
     两处分别处理（`graphdigger.hpp:1793-1815, 1837-1859`）。
  3. **不可扩展 fork**：首后继丰度>5 时，剔除 `ExtendableSuccessor` 为假的。
  4. **strand 平衡问题**：存在 `min(plusf,minusf)>0.25` 的双链好节点时，剔掉
     `min(plusf,minusf) < 0.1×fraction×max(...)` 的偏链后继
     （`graphdigger.hpp:1861-1884`）。
- **可逆性检查**（`GetReversibleNodeSuccessors`，`graphdigger.hpp:1740-1763`）：
  扩展前验证每个后继再回退（对后继的 revcomp 求后继）能回到原节点，否则该
  fork 不可逆、断裂。另有变体 `GetReversibleNodeSuccessorsF`（`:1710`）在过滤时
  额外记录 `eRightFork/eLeftFork` 位（左右两向分别看 step_back.size()>1）。
- **core 扩展 `ExtendToRight`**（`graphdigger.hpp:2273-…`）：从初始节点向右延申，
  对当前节点的后继做 `FilterNeighbors` 后分三种情形：
  1. `successors.empty()` → 无延伸，断；
  2. `size()==1`（简单扩展，`:2285-2305`）→ 新节点须 `GoodNode`（丰度≥low_count），
     且其前驱过滤后**恰好 1 个**、且该前驱的 revcomp 能**回到当前节点**（`:2289-2294`）
     才步进——这是"只沿唯一、可回退路径"的强不变量；
  3. `size()>1` 且未开 SNP → 断（宁可断在 fork）。
  左右两向由 `ExtendToRight(initial,0)` 与 `ExtendToRight(revcomp(initial),0)`
  分别调用（`:2448-2449`）。
- **`ExtendableSuccessor`**（`graphdigger.hpp:1657-1708`）：判定某后继是否"可扩展"，
  从它出发在图上至少走 `total_len=max(100,kmer_len)` 步（即**至少延伸到 100bp**）仍
  有路可走才算可扩展；沿途每步用 `FilterLowAbundanceNeighbors` 剔噪。这是"not
  extendable forks"过滤（`:1827`）的底层依据。
- **`jump`/`max_snp_len`**（`--max_snp_len` 默认 150）：开 SNP 时
  `max_extent=m_jump`（`:2276`），`DiscoverSNPCluster`（`:2170`）在 fork 处向前后
  两个方向探测一个 SNP 簇，验证 `step==step_back` 且长度一致后"跳过"多态区桥接
  （`:2324-2371`）；`--allow_snps` 开 `check_repeats` 额外做 SNP 发现（Rust
  `snp_discovery.rs`）。

## 5. 迭代组装编排（`assembler.hpp` / Rust `assembler.rs`）

`CDBGAssembler` 多轮迭代，k 从小到大逐步解重复：

1. **建首图 @ min_kmer**（默认 21）；算 read 平均长、genome size 估计
   （从 k-mer histogram 的 `CalculateGenomeSize`）。
2. **自动抬阈值**（`assembler.hpp:963-981`，Rust `assembler.rs:177-198`）：
   若 coverage 过高，`new_min_count = coverage/50`、`new_max_kmer_count =
   max(10, coverage/10)`（下限 10），并 `remove_low_count` 剪枝。
3. **GenerateNewSeeds → ImproveContigs**：`graph_digger` 保守组装出 seed contig
   （jump=0 的保守版，`ImproveContigs` 见 `assembler.hpp:713`）；有 `--seeds`
   则从种子扩展；`ConverToSContigAndMarkVisited`（`assembler.hpp:730`）把上一轮
   contig 用到的 k-mer 标为 visited，`GenerateNewSeeds`（`graphdigger.hpp:2871`，
   Rust 对应 `assemble_contigs_with_visited`，`graph_digger.rs:154`）在已组装区
   之外找新种子。
4. **max_kmer 估计**：`max_kmer = read_len+1 - (max_kmer_count/avg_count)×(read_len-min_kmer+1)`，
   clamp 到奇数（`assembler.rs:316-329`）。
5. **paired-end 连接**：
   - `estimate_insert_size`（抽样 10000 对，用首轮图估 insert N50，clamp 到
     `MAX_KMER`）；`paired_insert_limit = 3×N50`。
   - 若 `N50 > 1.5×max_kmer` 才启用**长 insert 双端迭代**（`use_long_paired_iterations`），
     否则直接 `connect_pairs` 在首图连 mate。
   - Rust `paired_reads.rs`：`connect_pairs` / `estimate_insert_size_full`。
6. **clean reads**：`clean_reads` 把完全落在已组装 contig 内的 read 剔除
   （`cleanup_min_contig_len = max(max_kmer, paired_insert_n50)`，Rust
   `assembler.rs:408`），防止陈旧 k-mer 污染下一轮 histogram。
7. **后续轮**：`max_kmer` 往上的每轮，用上一轮 contig 做"引导"，把未解重复区
   用更长 k 重连；`linked_contig.rs` 的 `ConnectFragments` 走连接链。
8. 输出 contigs（`--min-contig` 默认 200 过滤），可选 GFA。

> 关键工程点：**clean 的阈值语义**——min_contig_len 用的是 `max_kmer` 与
> `paired_insert_n50`（连接 mate 的 N50），**不是** insert 的 3×N50 上限；
> 用错会把 (N50, 3×N50) 区间已组装的 contig 排除在 kmer→contig 映射外，导致
> 陈旧 k-mer 残留（Rust `assembler.rs:400-407` 注释专门说明）。

### 5.1 关键参数与默认值（`skesa.cpp:318-357` 的 boost::program_options）

| 参数 | 默认 | 语义 |
|---|---|---|
| `--cores` | 0（=全部核） | 线程数；>硬件上限会 WARNING 并钳到硬件数（`skesa.cpp:445-456`） |
| `--memory` | 32 | 内存预算（GB，**仅排序计数**用，决定外部归并轮数） |
| `--kmer` | 21 | 最小 k（min_kmer，建首图的 k） |
| `--min_count` | 自动 | 保留 k-mer 的最小计数；高覆盖时自动抬到 `coverage/50`（`assembler.hpp:971`） |
| `--max_kmer` | 自动估计 | 最大 k；公式见 §5 第 4 步，clamp 到奇数 |
| `--max_kmer_count` | 自动 | 估 max_kmer 用的最低平均计数；自动抬到 `coverage/10`（`assembler.hpp:973`） |
| `--steps` | 11 | min→max k 的迭代轮数 |
| `--fraction` | 0.1 | 扩展允许的最大 noise/signal 比（`FilterLowAbundanceNeighbors` 的 `m_fraction`） |
| `--max_snp_len` | 150 | SNP 跳转的最大跨度（作为 digger 的 `jump`） |
| `--min_contig` | 200 | 输出 contig 最短长度 |
| `--vector_percent` | 0.05 | 含 19-mer 的 read 占比阈值，判定接头/载体（`1.` 关闭） |
| `--insert_size` | 自动估计 | 期望 insert size；`paired_insert_limit = 3×N50`，N50 由首图连接抽样 10000 对估得（`assembler.hpp:205,256`） |
| `--hash_count` | off | 用哈希计数（bloom + 分块表）取代排序计数；配套 `--estimated_kmers`（默认 100，百万）、`--skip_bloom_filter` |
| `--allow_snps` | off | 多一轮 SNP 感知遍历（`ImproveContigs(kmer,true)`，`assembler.hpp:301-306`） |
| `--use_paired_ends` | off | 单个 fasta/fastq 文件内即配对的 read（不逗号分隔） |

> 注：`estimate_min_count`（Rust `assembler.rs:177`）默认 true，对应 C++ 的
> `--min_count` 自动抬升逻辑；C++ 端该行为由 `GetGraph` 内的
> `total_seq>0 && genome_size>0` 触发（`assembler.hpp:963`）。

## 6. SAUTE / 引导组装（`guidedassembler.hpp`、`saute.cpp`）

- SAUTE 用**目标序列（参考）引导**：`guidedgraph`/`guidedpath_naa` 对目标区域
  做目标富集 de Bruijn 组装，输出 GFA + 两条 FASTA。
- Rust `guided_assembly.rs`/`guided_graph.rs`/`guided_path.rs` 仅为**简化版辅助**
  （README 明确 "full SAUTE parity is not yet implemented"）；`spider_graph.rs`
  对应 gfa-connector 的连接路径辅助，同样未完全对等。

## 7. 与 pgr 的关联 / 借鉴点

pgr 已有 `libs/kmer`（KmerTable：canonical 2-bit u128、精确计数、radix sort、
rayon 并行），是**精确计数路线**；SKESA/skesa-rs 提供两条互补路线 + 全套
de Bruijn 图遍历启发式：

1. **count 打包布局**（`DBGraph.hpp:180`）：total(32bit)+branch(8bit)+plus-fraction(16bit)
   一个 u64 同时承载计数/分支/链向——pgr 若扩展 k-mer 表，可参考这种打包省内存。
2. **canonical 存储 + Node 偶数/奇数编码**（`DBGraph.hpp:102-164`）：图只存
   canonical，用奇偶位表达链向、`Index()=m_node/2-1` 映射数组——pgr `KmerTable`
   已是 canonical key，若做 de Bruijn 图可直接套用该节点编码。
3. **分支位快速找后继**：8bit 分支掩码 + `(kmer<<2)&max` + 打包位跳过，避免
   逐碱基探测——对 pgr 未来 `asm` 类功能的图遍历是高性能范式。
4. **保守扩展 + fork 过滤启发式**（`graphdigger.hpp:1769-1887`）：低丰度 fork
   剔除、strand 特异 GGT 噪音剔除、strand 平衡检查、可逆性检查——四层过滤是
   "在重复区断裂"的实现核心，pgr 若实现 de Bruijn 组装应移植这层语义。
5. **迭代多 k + paired-end 连接 + read 清理**：从 min k 到 max k 渐进解重复，
   clean 阈值用 `max(max_kmer, paired_insert_n50)` 的细节值得照搬。
6. **确定性**：排序 + 稳定启发式保证输出确定（多线程不改变结果）——与 pgr
   "字节级一致"的硬约束一致，是排序计数优于哈希计数的理由之一。
7. **Rust 移植经验**（skesa-rs 对 pgr 的直接价值）：
   - `enum Kmer` + 宏展开 16 精度变体，替代 C++ boost::variant；
   - `Storage` enum 按精度选内联数组 vs boxed，兼顾缓存与正确性；
   - 排序阈值（>10000 并行 else 串行）保小规模确定性；
   - 用 `noodles`+`flate2` 替代 C++ NGS 库；把 CLI 做成 optional feature（库优先）。

### 7.1 OLC v1 借鉴映射（2026-08-12）

承接 `design/asm-olc.md` 的 v1 待决项，SKESA 提供三块直接素材：

1. **覆盖度/丰度驱动的 fork 过滤 → OLC repeat breaking 参数**：
   `FilterLowAbundanceNeighbors`（`graphdigger.hpp:1770`）的
   `abundance <= fraction × Σabundance`（`--fraction` 默认 0.1）与
   "`LowCount()==1` 且首后继丰度>5 时删丰度==1 尾巴"；`FilterNeighbors`
   的不可扩展 fork 剔除（`:1827`）与 strand 平衡检查（`:1863`，
   `min(plusf,minusf) < 0.1×fraction×max` 剔偏链）。这是"在重复区断裂"的
   成熟多层阈值语义——pgr layout 的 v0 repeat 检测只有 top2 近等边近似
   （`canu.md` §8.5 记录了 6× 低覆盖漏检案例），v1 应移植这层丰度阈值
   （pgr `asm unitig` 头部已带 `cov=`，可直接取用）。
2. **可逆性检查（`GetReversibleNodeSuccessors`，`graphdigger.hpp:1740`）→
   layout 延伸的"回得来"保证**：SKESA 扩展前验证后继的 revcomp 能回到原
   节点、否则断裂；与 pgr layout 的互惠 best edge（`canu.md` §8.5 连接端
   语义）是同一思想的两个实现——SKESA 在 k-mer 图、pgr 在 unitig 重叠图。
3. **迭代多 k + read 清理 → `asm olc` 的多 k 反馈**：pgr `asm olc` 目前
   各 k 独立出 unitig 再合并（无反馈）；SKESA 的"每轮用上一轮 contig 引导
   + `clean_reads`（阈值 `max(max_kmer, paired_insert_n50)`，
   `assembler.rs:407`）"与 metaMDBG 的 unitig 反馈（`metaMDBG.md` §4.1）
   同族，是 v2 候选。
4. **count 打包 + Node 奇偶编码（`DBGraph.hpp:180/102`）**：若 pgr 从
   bcalm 式哈希表切换到排序 DBG（kmer 表已走 radix 排序路线），
   total(32)+branch(8)+plus-fraction(16) 打包与 `Index()=m_node/2-1`
   是现成模式。
5. **"只沿唯一可回退路径"的扩展不变量（`ExtendToRight`，`graphdigger.hpp:2285-2294`）→
   layout 的唯一性约束**：SKESA 的简单扩展要求新节点前驱过滤后恰好 1 个、且其
   revcomp 能回到当前节点才步进；`ExtendableSuccessor` 再以"至少延伸 100bp"兜底
   （`:1659`）——这组"唯一 + 可回退 + 最小延伸长度"约束可直接映射 pgr `asm unitig`
   的延伸语义（`cov=` 已在头部），比单纯的 top2 近等边近似多一道方向性验证。

### 7.2 与 asm multik 的多 k 迭代对比（2026-08-14 重读补充）

本次重读 `assembler.hpp`（主流程）+ `graphdigger.hpp`（扩展）+ skesa-rs
对应文件（`assembler.rs` / `graph_digger.rs`），确认 SKESA 的迭代语义与
`asm multik`（metaMDBG 式）是**同一目标（大 k 解小 k 的重复）下的两个不同
实现族**。逐轮差异如下（行号均以本次阅读为准）：

| 维度 | SKESA / skesa-rs | `asm multik`（metaMDBG 式） |
|---|---|---|
| 每轮建图输入 | **清理后的未用 reads**（`CleanReads`，`assembler.hpp:685`；阈值 `max(max_kmer, paired_insert_n50)` + margin `max_kmer+50`）——read 集逐轮递减 | **全部 reads + 上一轮 unitig 序列**（metaMDBG §4.1.1）——数据只增不减 |
| 上一轮结果的反馈 | **contig 序列**：`ConverToSContigAndMarkVisited` 把旧 contig 的 k-mer 标记 visited（`assembler.hpp:730`），`GenerateNewSeeds` 只从**未 visited** 节点组装（`graphdigger.hpp:2871`），`ExtendContigsJob` 再从旧 contig 边缘 k-mer 扩展（`graphdigger.hpp:3226`） | **unitig 图 + 序列**：unitig 进计数表改 solid 集，图结构经跨接验证/回灌参与下一轮 |
| k 序列 | min→max 等距（`steps` 默认 11，`assembler.hpp:263-268`），长 insert 另加 `[1.25×max_kmer, 中点, insert_N50]` 三轮（`:286-289`） | metaMDBG 式逐 +1（k′-min-mer）或 `--kmer auto` 按读长 N50 生成 |
| bubble/fork 处理 | 默认 **断在 fork**（`ExtendToRight` 的非 SNP 分支即 break，`graphdigger.hpp:2308`）；开 `--allow_snps` 才用 `DiscoverSNPCluster` 前后双向探测、`step==step_back` 一致才跨过（`:2325-2371`） | `bridge_kmer` 跨接验证（跨 unitig 的 k-mer ≥2 才保留边）+ `progressive_filter` 中位 25% cutoff（直链保护）+ v6 **被删分支回灌**（低丰度分支继续参与） |
| 嵌合清理 | 无显式"整条删"，靠扩展不变量前置预防：简单扩展要求新节点前驱恰好 1 且可回退（`graphdigger.hpp:2285-2294`） | `remove_unsupported`：unitig 内部 k-mer 缺失率 >2% 则整条删 |
| 防 misassembly | 四层后继过滤：低丰度 fork（`≤fraction×Σ`）、GGT/ACC strand 噪音、不可扩展 fork（<100bp）、strand 平衡（`graphdigger.hpp:1769-1887`） | 最终阶段 `bridge_filter` + `split_by_bridge`（unitig 间/内 60-mer 探针 ≥2 验证） |
| 每轮产出 | contigs（SContig 链：connectors/extenders 通过 `m_left_link/m_right_link` 链回父 contig，`graphdigger.hpp:1287`） | unitig 图（unitigs + 跨接边） |

**两个关键共识**（也是 multik 设计里已吸收的）：

1. **宁断勿嵌合**：SKESA 默认断在 fork、扩展要求"唯一 + 可回退"；multik
   剪低丰度分支、删除无支撑跨接——都把正确性放在连通性前面，这是"无 N
   但不是无嵌合"的取舍（G37 misassemblies 8→0 的路线与此同源）。
2. **大 k 只解决小 k 未解决的部分**：SKESA 用 visited 标记让新种子只出现在
   旧 contig 之外、multik 的 unitig 反馈让已覆盖区天然有更高支持——都避免
   每轮从头重拼，计算量随轮次收敛（实测 multik 图结构逐轮递减：unitigs
   1345→396、edges 346→12，见 `benchmarks/multik-complexity.md`）。

**SKESA 细节的吸收情况（2026-08-14 实现验证）**：

- **"前驱恰好 1 + 可回退"不变量 → 已吸收（严格链唯一性）**：`merge_chains`
  / `recompact_graph` 原来用"先到先得"的端点占用检查（`right_of/left_of`
  `is_some()`），汇点（两个前驱）会被吞进先遍历到的链；改为**严格两端唯一**
  ——只有链段两端的定向度数都恰为 1 才合并（SKESA "predecessor == 1" 的
  unitig 图对应）。实现时发现一个易错点：`compute_links` 对同一 junction
  从两端各发一条 link，方向解析后是**同一条链段**，逐 link 度数统计会把
  对称 link 翻倍，必须用 HashSet 去重后再计度数。G37 回归：misassemblies
  保持 0、mismatches 27.7/100kbp（历史最佳）、N50 24.4K（-8%，宁断勿嵌合
  的正确性代价）、Genome fraction 95.86%（-0.13pp）。
- **read 清理（`CleanReads`）→ 确认不做**：与 multik 的 `remove_unsupported`
  机制冲突——multik 的 unitig 反馈以 1× 进计数表，剔除"完全落在 unitig 内"
  的 reads 后，已覆盖区内部 k-mer 计数降到 1×（< `min_count_extend`=2），
  `remove_unsupported` 会把真实 unitig 整条误删；而提高 unitig 反馈权重到
  ≥ threshold 又会让嵌合 unitig 自我支撑、嵌合清理失效。SKESA 能清理是因为
  它没有"内部 k-mer 必须 ≥ threshold"这一环（靠保守扩展不变量防嵌合）。
  且 multik 的耗时瓶颈是 `remove_unsupported` 的 O(总长×k)，不在计数输入量。
- **`allow_snps` 的双向 SNP 簇验证**（`DiscoverSNPCluster`）→ **暂缓**：
  multik v6 回灌保留低丰度分支但不合并成多态表示；SKESA 是"前后一致才合并
  为变体 chunk"。要表达菌株/单倍型多态需要引入 `CContigSequence` 的
  chunk+variants 表示，属输出表示层大改，等有真实菌株/宏基因组数据再定。

## 8. 局限

- C++ 依赖 Boost（variant）+ NGS 库（SRA），构建链较重；Rust 版已剥离。
- SAUTE / gfa-connector / SRA 在 Rust 版未对等移植。
- 哈希计数（bloom+分块表）有哈希碰撞/近似性，确定性弱于排序计数；pgr 的精确
  路线与哈希路线各有所长（内存 vs 精确）。
- SKESA 定位微生物短读组装；长读（HiFi/ONT）场景不在其路线（参照 metaMDBG）。
- skesa-rs 是 LLM 中介的"忠实翻译"，README 明示**不可完全信任**、复刻 bug、
  需自行验证——参考其工程手法而非当作权威实现。

## 9. 源码 quirks（异常/边界行为）

- **越界访问无检查**：`CDBGraph` 多处注释明示 "for all access with Node there is
  NO check that node is in range !!!!!!!!"（`DBGraph.hpp:173,425`）——调用方必须保证
  Node 有效，否则 `Index()=m_node/2-1` 可能下溢；`Node(0)` 是"无效/未命中"哨兵
  （`isValid()` 要求 `m_node>0`，`DBGraph.hpp:106`）。
- **Node 编码与奇偶**：偶数=正链、奇数=负链；`ReverseComplement` 就是 ±1（`:109-114`），
  `DropStrand` 归到偶数（`:115`）；`GetNode` 对图中不存在的 k-mer 返回 `Node(0)`
  （`:156-164`）。
- **`--cores` 边界**：负值直接 `exit(1)` 报错（`skesa.cpp:448-451`）；超过硬件上限只
  WARNING 并钳到硬件数（`:451-455`）；0 表示用全部核（`thread::hardware_concurrency`）。
- **输入去重**：`--reads`/`--fasta`/`--fastq`/`sra_run` 会 `sort+unique` 去重，有重复
  时打印 WARNING（`skesa.cpp:414-441`）——同输入会因**顺序无关**而稳定。
- **`--gz` 已废弃**：不再需要，自动识别 gzip，传了只打 WARNING（`skesa.cpp:384-385`）。
- **k 一律取奇**：`kmer_len -= 1-kmer_len%2`（`assembler.hpp:267`；Rust `assembler.rs:327-329`），
  保证 k 为奇数（回文 k-mer 无 canon 歧义）。
- **迭代 k 的取整**：主迭代 k 由 `min_kmer+step*alpha+0.5` 四舍五入再取奇，若 ≤上一
  张图的 k 则 `continue` 跳过（`assembler.hpp:264-269`）；长 insert 迭代固定用
  `[1.25×max_kmer, (1.25×max_kmer+max_kmer_paired)/2, max_kmer_paired]` 三个 k
  （`:286-289`）。
- **`ExtendableSuccessor` 的"至少 100bp"**：`total_len=max(100,kmer_len)`（`graphdigger.hpp:1659`），
  即使 k 很小也要求候选路径至少延伸 100bp 才算"可扩展"——对短读小 k 是保守阈值。
- **`ExtendToRight` 接受未拥有的节点**：签名注释 "initial_node may be not owned"
  （`graphdigger.hpp:2274`），且中途 `SetVisited` 失败即断（`:2297-2305`）。
- **GGT 噪声明辨**：strand 特异剔噪仅在 `GraphIsStranded()` 时生效；正链按末 3 碱基
  `== "GGT"` 找 target（`:1801`），负链按 `MostLikelySeq(suc,3)=="ACC"`（`:1845`，
  即 GGT 的 revcomp）——两处丰度阈值用 `1-PlusFraction`/`PlusFraction` 分别归一到
  对应链。
- **低丰度尾巴删除的触发条件**：`LowCount()==1 && 首后继丰度>5` 时才删丰度==1 的
  尾巴（`:1786-1789`），循环条件是 `j>0 && 当前末尾丰度==1`——避免过度删。
- **count 打包的"未用"字节随图类型漂移**：`CDBGraph` 中 `[40:47]` 是 8 bit 未用
  （`DBGraph.hpp:180`）；`CDBHashGraph` 把它改作 visited 控制位
  （`eVisited/eTemp/eMulti=1<<40/41/42`，`DBGraph.hpp:454`）——同一定义在不同图
  实现语义不同，移植/读码须留意。
- **insert 估计用 2000bp 上限**：`long_insert_size=2000` 作为连接步长上限
  （`assembler.hpp:240`），`m_insert_size=3×N50`（`:256`）。
- **自动抬 `min_count` 的封顶差异**：排序计数版无上限（`assembler.hpp:971`），哈希
  版封顶 `min(255,…)`（`:1010`，因哈希表计数是 8bit）；公式都带 `+0.5` 四舍五入。
- **空组装抛异常**：首轮 `ImproveContigs` 后 contig 为空即
  `throw runtime_error("Was not able to assemble anything")`（`assembler.hpp:182-183`）。

---

*参考来源: `SKESA-master/`（skesa.cpp、concurrenthash.hpp、counter.hpp、
DBGraph.hpp、graphdigger.hpp、assembler.hpp、guidedassembler.hpp、kmercounter.cpp、
gfa.hpp）+ `skesa-rs-main/`（src/kmer.rs、counter.rs、db_graph.rs、assembler.rs、
graph_digger.rs、cli.rs、Cargo.toml、README.md）*
