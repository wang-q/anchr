# `fq qc` 端到端基准：anchr vs FastQC vs Falco

> 2026-08-13。hyperfine 1.19.0，Lambda 数据（SRR5042715，108 bp）：
> `R1.2k.fq.gz`（2000 reads）与 `R1.fq.gz`（20000 reads）。
> 单线程对比（fastqc 默认 / falco `-t 1` / anchr `-p 1`，输出均为各自
> 完整 QC 报告）；另测多线程实际用法。

## 结果

### 2000 reads（R1.2k.fq.gz，warmup 2 / runs 5）

| 工具 | 时间（mean ± σ） | 相对 |
| :--- | ---: | ---: |
| fastqc 0.12.1（Java） | 1.119 s ± 0.003 s | 基线 |
| falco 2.0.1（C++）`-t 1` | 7.6 ms ± 0.5 ms | fastqc 的 147× |
| anchr（Rust）`-p 1` | 15.5 ms ± 1.1 ms | fastqc 的 72×、falco 的 0.49× |

### 20000 reads（R1.fq.gz，warmup 1 / runs 3）

| 工具 | 时间（mean ± σ） | 相对 |
| :--- | ---: | ---: |
| fastqc 0.12.1 | 1.129 s ± 0.008 s | 基线 |
| falco 2.0.1 `-t 1` | 25.1 ms ± 0.1 ms | fastqc 的 45× |
| anchr `-p 1` | 76.3 ms ± 2.6 ms | fastqc 的 15×、falco 的 0.33× |

### 多线程（20000 reads，warmup 1 / runs 3）

| 命令 | 时间 | 备注 |
| :--- | ---: | :--- |
| falco `-t 8` | 20.6 ms ± 0.7 ms | 相对单线程略快 |
| anchr `-p 8` | 72.7 ms ± 2.1 ms | 分块/合并开销抵消收益 |
| anchr `-p auto`（32 核） | 102.5 ms ± 14.3 ms | sys 703 ms，线程调度开销 |

## 分析

- **fastqc** 时间被 JVM 启动 + HTML/zip 报告主导（1.1s 基本与 reads 数无关），
  端到端最慢；
- **falco**（C++）单线程最快，anchr 慢约 3×——绝对值小（76 ms vs 25 ms），
  差距来自每 read 的 consume 开销（质量/GC/adapter indexOf/contaminant 标注）
  与 C++ 的 SIMD/内存布局差异；
- **多线程在小数据上是负优化**（20000 reads 分块粒度太细），真实规模
  （百万 reads）下 anchr 的 rayon 并行预期有收益，需大数据验证；
- 输出差异：fastqc 默认 html+zip，falco 三件套，anchr `-f txt` 只出文本
  （`-f both` 含 html，未在基准中测）。

## 备注

- 复现：`hyperfine --warmup 2 --runs 5 "fastqc -q <in> -o <d>" "falco -o <d> -t 1 <in>" "anchr fq qc <in> -o <d> -p 1 -f txt"`；
- fastqc 的 `-o` 目录需预先存在；
- 大数据（>1M reads）基准待做，用于评估并行扩展性与 3× 差距是否放大。

## G37（340k reads，2026-08-13 补充）

输入 `~/data/anchr/g37/2_illumina/R1.fq.gz`（E. coli G37，ENA ERR486835，
340322 reads × 150 bp），warmup 1 / runs 3：

| 命令 | 时间（mean ± σ） | 相对 anchr_p1 |
| :--- | ---: | ---: |
| fastqc | 3.151 s ± 0.002 s | 0.48× |
| falco `-t 1` | 390.6 ms ± 6.9 ms | **3.88×** |
| anchr `-p 1` | 1.515 s ± 0.018 s | 1× |
| falco `-t 8` | 247.2 ms ± 3.1 ms | 6.13× |
| anchr `-p 8` | 616.7 ms ± 5.3 ms | 2.46× |

**结论**：真实规模下单线程差距从 Lambda 的 3× **放大到 6.1×**（每 read
约 4.4 µs vs falco 1.1 µs）——anchr 的每 read consume 开销是主要瓶颈；
并行有效（`-p 8` 2.5×），但仍比 falco 单线程慢 1.58×。绝对值：340k reads
下 anchr 1.5 s vs falco 0.39 s；对几 M reads 的真实项目会到十几秒 vs 数秒，
**值得优化**（profile 定位热点，参照 pgr 的 SIMD 三步模式：
`notes/design/simd-optimization.md`，`pgr::libs::nt_simd::count_bases` 等
可复用；adapter 子串搜索可换 `memchr`）。

## 优化后（2026-08-13，perf 驱动）

针对 profile 热点（consume 内 4 项）：

1. `seq_counts` 的 `Vec<u8>` 堆分配 → 固定 `[u8; 50]` key（copy，无 alloc）；
2. 单线程流式处理（消除收集全部 reads 的 `SeqRecord::clone`）；
3. adapter `find_subseq` 首字节预筛 → 双字节预筛（memcmp 调用 ~16×↓，
   memcmp self 7.98% → 1.90%）；
4. per-read 三遍遍历（per-base + 平均质量 + GC 计数）合并为一遍；
5. seq_counts/kmer 的 HashMap 预分配容量 + per-base 双 match 合并。

优化后 G37 单线程 **1.515 s → 0.672 s（56%↓）**，8 线程 0.617 s → 0.546 s：

| 命令 | 优化前 | 优化后 |
| :--- | ---: | ---: |
| falco `-t 1` | 390.6 ms | 373.2 ms |
| anchr `-p 1` | 1.515 s | **0.672 s**（1.77× vs falco，原 6.13×） |
| anchr `-p 8` | 616.7 ms | 546.4 ms（1.44× vs falco） |

追加优化（第 2 轮）：`QualityCount` 128 槽 ASCII → 64 槽 Phred bin
（per-position 数组 150KB → 37KB，L1/L2 友好）；`seq_quality_hist`
BTreeMap → 扁平数组（消除每 read 一次树查找）。单线程再降
0.881 → 0.791 s（累计 48%）。

追加优化（第 3 轮）：adapter 匹配合并为**一次遍历查 6 个 adapter**
（前两字节互斥分支 + 早停，替代 6 次 `find_subseq` 扫描，memcmp
7.76% → 3.48%）；`seq_counts` 换 FNV-1a hasher + 更大容量；per-base
计数 u64 → u32（数组内存减半，L1 友好）。单线程 0.791 → 0.672 s
（累计 56%），差距 2.07× → 1.77×。

剩余差距来自 per-base 位置特异的质量/内容累积（本质工作，falco 同做，
C++ 内存布局/编译优化占优）；继续深挖需 SIMD 向量化位置统计
（pgr 三步模式），收益递减，暂缓。

## 第 4 轮：分支/解析器/数据结构深挖（2026-08-13 续）

用户要求继续挖（不接受"编译器优化差距"结论）。perf 揭示关键事实：
**指令数几乎相同（3.32B vs falco 3.39B）但 IPC 1.60 vs 2.97，分支未命中
23.8M vs 3.2M**——不是工作量问题，是分支预测/缓存效率问题。本轮改动：

1. **主循环 2 位置展开 + 静态查表替代 5 路跳表**：`BASE_TABLE`（256B，
   `(count_slot<<3)|pc_slot` 编码）替代间接跳转；A/C/G/T 计数合并为
   `base_counts: [u64;4]` 数组；每位置消掉 2 个边界检查（`qual[i]` 的
   索引证明丢失）+ 跳表间接分支（随机序列上误预测高）；
2. **trunc_gc 循环拆分**：按 fastqc 语义（GC 只数前 100/1000 bp）拆成
   两段循环，去掉每位置 `i < trunc_gc` 比较；
3. **len_hist BTreeMap → 扁平 Vec**：每 read 一次树查找（指针追逐）→
   O(1) 数组自增；
4. **kmer HashMap → 稠密网格**：`[u32; 16384 × L]` 行主序（key 主序，
   `kmer_stride` 增长时重排一次）；列主序版本因 64KB 步长缓存全 miss
   已回退；报告端遍历全部 16384 行求和；
5. **adapter 扫描**：2 字节前缀 → 64KB 查表 → 4-bit 打包 32KB 表（L1
   驻留）+ 每迭代 2/4 位置 + `[u8;65536]` 定长数组消边界检查 +
   `get_unchecked` 消 `sigs[a]`/`found[a]` 检查；12 bp 固定宽度比较
   （u16+u64 load）替代 libc memcmp；
6. **零拷贝 FASTQ 解析**：plain 文件 mmap + `next_record`（4 次 memchr
   + 切片，无逐行拷贝/分配），仿 falco；gz/stdin 保持 pgr 流式（整文件
   解压到 Vec 的匿名页错误开销 > 解析收益，实测倒退）；
7. **FNV 哈希 51B → 前 8 + 后 8 字节**（2 次乘法 vs 7 次；完整 key 仍
   存储比较，正确性不变；尾部折叠避免 adapter 前缀家族聚桶）；
8. **tile 解析去分配**：`split(':').collect()` → 直接数冒号 + 字节十进制
   解析，无 Vec 分配/utf8 转换。

回退实验（均实测更慢）：falco 位运算编码（-5%）、标量换行扫描替代
memchr（+47ms：memchr 的块级分支 miss 虽多但代价低，指令增量更贵）。

**G37 结果（12 runs，2026-08-13 晚）**：

| 命令 | plain（120MB） | gz（45MB） |
| :--- | ---: | ---: |
| falco 2.0.1 `-t 1` | 183.8 ms | 330.0 ms |
| anchr `-p 1` | **303.4 ms**（1.65×） | **376.2 ms**（1.14×） |

累计：初始 1.515 s → 376 ms（**4.0×**）。perf：consume 66.6%（主循环
~25%、adapter ~20%、seq_counts 5.5%），报告期 kmer_content+beta ~7%，
解析 ~4%；指令 3.23B / 1.67B cycles（IPC 1.93）；分支未命中 16.8M。
剩余差距主要是主循环每位置 3 次 RMW（质量槽 256B 步长、内容槽 20B
步长、计数）+ seq_counts 每 read 插入；falco 用独立简单循环换 IPC。

护栏：339 测试全绿；Lambda golden summary.txt 逐字节一致；fastqc_data
数值容差内一致（仅浮点显示格式差异，为既有状态）。

## Falco 实现方式分析（2026-08-13）

深入 `falco-2.0.1-Source/src/`（`results_collector.hpp` 的
`process_one_read`、`falco_utils.hpp` 的统计函数）：

1. **无 SIMD**：`count_nucs`/`count_gc`/`count_ns`/`count_quals` 全是标量
   `while` 循环（`++(*out_itr++)[encode(*seq_itr++)]`），无手写向量化；
2. **位运算编码**：`encode(c) = (c>>1)&3`、`is_gc(c) = (c>>1)&1`——
   但移植到 Rust 后**反而慢 5%**：LLVM 已把 anchr 的 5 路 ASCII `match`
   优化到等价的无分支查表，手动位运算 + N 守卫干扰了优化（已回退）；
3. **多次独立遍历**：每 read 4-5 遍（nucs/gc/ns/quals/dup/adapter）——
   anchr 合并为一遍（省遍历）；Falco 的独立简单循环在 C++ 下无分支
   也很快；
4. **质量数组 `[u64; 127]`**：和 anchr 优化前的 128 槽一样，Falco 未用
   紧凑存储；
5. **差距本质**：C++ `-O3`（无边界检查、无 trait 分发）对简单循环的
   激进编译优化。anchr 的算法层优化（合并遍历/紧凑存储/消除分配）已
   到位，剩余 ~1.7× 属语言/编译器级差距，算法层难以再压缩。
