# MEGAHIT（1.2.9）：succinct de Bruijn 图宏基因组组装器（源码分析）

> 2026-08-13 整理、2026-08-17 对照源码逐项复核修订，纯源码分析（`megahit-1.2.9/`，
> 版本 `v1.2.9`）。MEGAHIT 是面向
> NGS 短读（尤其宏基因组）的**超快、内存友好**组装器，论文 [MEGAHIT: An
> ultra-fast single-node solution for large and complex metagenomics assembly
> via succinct de Bruijn graph](https://doi.org/10.1093/bioinformatics/btv033)
> （Bioinformatics 2015，作者 Dinghua Li 等，香港大学 Lam 组）+ 后续 v1.0 论文
> （Methods 2016）。**核心创新**：把 de Bruijn 图压缩成 **succinct de Bruijn
> graph（SdBG）**——用位向量（rank/select）存图，边 = (k+1)-mer，内存占用远低于
> 传统哈希表/邻接表，从而能在单节点上组装超大宏基因组。**与 anchr 的关系**：它和
> anchr `asm contig/unitig`（tadpole/BCALM 移植）同属 k-mer 图组装路线，但走
> "外部排序计数 + 位向量存图 + 迭代多 k + 图上清洗"的磁盘友好范式，正好是 anchr
> 现在"内存中建 KmerTable"路线的对照物（见 §8）。
>
> 两层架构：**Python 驱动**（`src/megahit`）负责参数解析、checkpoint 断点续跑、
> 流程编排（反复调用 C++ `megahit_core` 子命令）；**C++ 核心**（`megahit_core`，
> `src/main.cpp` 分发）实现各底层步。GPLv3。

## 1. 概况

- **定位**：面向 Illumina 短读（单/双端，支持 gz/bz2），优化宏基因组，也适用于
  单基因组/单细胞。输出 `final.contigs.fa`（`--min-contig-len` 默认 200 过滤）。
- **双链语义**：图里存 canonical（k+1)-mer（k-mer 与其反向互补只占一条边），
  `w` 数组高位置 1 表达反向互补的字符（§3.1）。
- **迭代多 k**：k 从 `k_min`（默认 21）到 `k_max`（默认 141）递增，默认列表
  `[21,29,39,59,79,99,119,141]`（相邻差 ≤ 28，k 必须奇数）；每一轮用上一轮
  contig/bubble/local contig 引导重建更大 k 的图，逐步解重复。
- **磁盘友好**：k-mer 计数走"外部多趟排序"（§4.1，`kNumBuckets=65536` 桶按大小
  排序分批处理），不是全内存哈希表；SdBG 本身用位向量存。
- **构建**：C++14/17 + OpenMP，CMake；`git submodule update --init` 拉
  `kmlib`/`parallel_hashmap`/`xxhash` 子模块；运行时依赖 gzip/bzip2。
- **命令形态**：`megahit`（Python 驱动）+ `megahit_core <sub>`（C++ 底层步：
  `assemble/local/iterate/buildlib/count/read2sdbg/seq2sdbg/contig2fastg/
  readstat/filterbylen/checkcpu/checkpopcnt/checkbmi2/dumpversion/kmax`）。

### 1.1 核心算法与流程总览（先读这节）

**一句话**：reads 的 (k+1)-mer 经**外部多趟排序**计数（§4.1）滤出 solid 边 → 压成
**succinct de Bruijn 图**（位向量 + rank/select，§3.1）→ unitig 压缩 + 多轮**图清洗**
（tip/bubble/weak-link/low-depth，§5）→ **迭代多 k**：每轮用上一轮 contigs 引导建更大
k 的图、用 reads 补迭代边、对 contig 端点做本地组装（§4.3-4.6）→ 各轮产物合并输出。

```
buildlib → (count | read2sdbg[1pass]) → build_graph[k_min] → assemble[k_min]
  ┌─────────────────────────────┐
  while cur_k < k_max:          │
    local_assemble(cur_k)       │  contig 端点本地延伸（IDBA-UD 内核，k∈[11,next_k)）
    iterate(cur_k, step)        │  reads 补跨端点的 (k+step+1)-mer 迭代边
    build_graph(next_k, cur_k)  │  上一轮 contig/bubble/local 全长引导建新 k 图
    assemble(next_k)            │  unitig 化 + 四件套清洗
  merge_final(k_max)            │  cat *.final.contigs.fa + filterbylen
```

| 核心块 | 机制 | 详见 |
|---|---|---|
| 外部排序计数 | 65536 桶（前 8 碱基）→ Lv0 桶大小统计降序分批 → Lv1 4 字节差分偏移 → Lv2 排序计数滤 solid | §3.4、§4.1 |
| SdBG 位向量图 | 边 = canonical (k+1)-mer 排序数组，`w/last/tip` 位向量 + rank/select 实现 Forward/Backward | §3.1 |
| unitig 压缩 | `NextSimplePathEdge`：唯一出边且该边唯一入边才继续延伸 | §3.1 |
| 图清洗 | tip（2 起倍增阈值 + 8× 深度比）、naive/complex bubble（banded 相似度 0.95）、weak link（0.1× 断开 = 截短一格）、low depth（局部窗口 `min(min_depth, mean×ratio)`） | §5 |
| 多 k 引导 | 上一轮 contigs/bubbles/local **全长**喂 `seq2sdbg` 作新 k 图的种子；`iterate` 另从 reads 收集跨 contig 端点的迭代边 | §4.3、§4.6 |
| 本地组装 | reads 回帖 contig 端点（seed 31、sparsity 8、只信唯一最佳）→ insert size 定区间（封顶 650）→ IDBA 多 k 延伸 | §4.5 |

## 2. 仓库结构

```
megahit-1.2.9/
├── src/
│   ├── megahit                  # ★ Python 驱动：编排 + checkpoint（1100+ 行）
│   ├── main.cpp                 # C++ 入口，子命令分发（assemble/local/iterate/...）
│   ├── main_assemble.cpp        # ★ assemble：unitig 图构建 + 清洗 + 输出
│   ├── main_iterate.cpp         # ★ iterate：迭代边提取（多 k 桥接）
│   ├── main_local_assemble.cpp  # local：本地组装（contig 延伸）
│   ├── main_sdbg_build.cpp      # count / read2sdbg / seq2sdbg 入口
│   ├── main_buildlib.cpp        # buildlib：reads → 二进制库
│   ├── definitions.h            # 常量（kMaxK=255、word 布局、GenericKmer）
│   ├── sorting/                 # ★ 排序引擎 + 计数 + SdBG 构建
│   │   ├── base_engine.{h,cpp}  #   BaseSequenceSortingEngine（Lv0/Lv1/Lv2）
│   │   ├── kmer_counter.{h,cpp} #   count：(k+1)-mer 计数 → edges.* 文件
│   │   ├── read_to_sdbg_s1.cpp  #   1pass 阶段 1：solid 标记 + mercy 候选
│   │   ├── read_to_sdbg_s2.cpp  #   1pass 阶段 2：写 SdBG（含 $ 虚拟节点）
│   │   ├── seq_to_sdbg.{h,cpp}  #   seq2sdbg：从 contigs/bubbles/edges 建图
│   │   ├── edge_counter.h       #   EdgeMultiplicityRecorder
│   │   └── kmsort_selector.*    #   GPU/CPU 排序选择
│   ├── sdbg/                    # ★ SdBG 数据结构
│   │   ├── sdbg_def.h           #   常量（mul_t/small_mul_t/kMaxK/kAlphabetSize）
│   │   ├── sdbg.h               #   ★ SDBG 类：rank/select、Forward/Backward 等
│   │   ├── sdbg_raw_content.h   #   裸数据布局（w/last/tip/mul）
│   │   ├── sdbg_meta.{h,cpp}    #   bucket 元数据（串行化/反序列化）
│   │   └── sdbg_writer.{h,cpp}  #   SdBG 写盘
│   ├── assembly/                # ★ unitig 图 + 清洗
│   │   ├── unitig_graph.{h,cpp} #   UnitigGraph（vertices + 双适配器）
│   │   ├── unitig_graph_vertex.h#   UnitigGraphVertex（40 字节打包）
│   │   ├── tip_remover.cpp      #   尖端移除（2,4,8,... 倍增阈值）
│   │   ├── bubble_remover.{h,cpp}#   Naive/Complex 气泡解析
│   │   ├── weak_link_remover.cpp#   弱连接断开（depth 比例）
│   │   ├── low_depth_remover.cpp#   低深度修剪（局部/迭代）
│   │   ├── sdbg_pruning.{h,cpp} #   SdBG 层尖端移除 + InferMinDepth
│   │   └── contig_output.{h,cpp}#   ContigWriter（FASTA 输出）
│   ├── localasm/                # ★ 本地组装
│   │   ├── local_assemble.{h,cpp}#   contig 端点延伸（IDBA-UD 内核）
│   │   ├── hash_mapper.{h,cpp}  #   reads → contigs 哈希回帖
│   │   └── mapping_result_collector.h
│   ├── iterate/                 #   迭代边提取
│   │   ├── contig_flank_index.h #   contig 端点 k-mer 索引
│   │   └── kmer_collector.h     #   收集 reads 里的下一轮 k-mer
│   ├── sequence/                #   k-mer 编码、SeqPackage、IO
│   │   ├── kmer.h               #   Kmer<N,WordType> 2-bit 编码
│   │   └── io/                  #   fastx/binary/edge/contig 读写 + kseq
│   ├── kmlib/                   #   kmbitvector(rank/select)、kmsort、kmcompactvector
│   ├── parallel_hashmap/        #   phmap（flat_hash_map / parallel_flat_hash_map）
│   ├── idba/                    #   IDBA-UD 内核（本地组装复用）
│   └── tools/                   #   contig2fastg / read_stat / filter_by_len
└── test_data/                   # r1.il.fa.gz / r2.il.fa.bz2 / r3_1.fa / ...（--test）
```

## 3. 核心数据结构

### 3.1 SdBG（`sdbg.h` + `sdbg_raw_content.h` + `sdbg_def.h`）

图里 **每条边 = 一个 (k+1)-mer**（节点 = k-mer）。边按字典序（canonical 化后）
排序存储，用位向量记录结构信息：

```cpp
struct SdbgRawContent {
  SdbgMeta meta;                                    // bucket 元数据（§3.2）
  kmlib::CompactVector<kAlphabetSize, uint64_t> w;  // 每条边末字符（+RC 位）
  kmlib::CompactVector<1, uint64_t> last, tip;      // last/tip 位向量
  std::vector<small_mul_t> small_mul;               // 小多重度（uint8，254 封顶）
  std::vector<label_word_t> tip_lables;             // tip 边的完整标签
  phmap::parallel_flat_hash_map<uint64_t, mul_t> large_mul;  // 大多重度（稀疏）
  std::vector<mul_t> full_mul;                      // 全量多重度（可选）
};
```

- **`w`**：每条边的末字符，2-bit/碱基（A=0..T=3），`kAlphabetSize=4`；反向互补的
  边其末字符 = `原字符 + 4`（`kWAlphabetSize=9`，含哨兵），读取时 `a > 4 → a -= 4`
  （`sdbg.h` `Forward`/`GetLabel`）。
- **`last`**：某节点全部出边的最后一条置 1——`Forward`/`ComputeOutgoings` 靠它在
  排序序列里圈定"同一出边组"的范围。
- **`tip`**：无法沿 `Backward` 走到头的边（读取端点）置 1，其完整 k 标签单独存进
  `tip_lables`（`TipLabelStartPtr` + `CharAtTipLabel`），这样 `IndexBinarySearch`/
  `GetLabel` 在走到 tip 时改读标签而不用继续回退。
- **多重度**：`mul_t=uint16`（kMaxMul=65535）；常规用 `small_mul`（uint8，
  `kMaxSmallMul=254`，`kSmallMulSentinel=255` 表示"查 large_mul"）；超出走
  `large_mul`（phmap 稀疏哈希表）；`full_mul` 为可选的全量数组
  （`EdgeMultiplicity`，`sdbg.h:96-105`）。

**核心遍历原语**（都是 rank/select 位操作，`kmlib/kmbitvector.h`）：

- `Forward(edge_id)`：`rs_w_.rank(a, edge_id)` 数"本组内第几个同末字符边"，
  `rs_last_.select(...)` 定位该边指向的出边（`sdbg.h:107-114`）。
- `Backward(edge_id)`：用 `LastCharOf`（查累计表 `f_`）找该边所属组，再
  `rs_w_.select` 回第一条入边（`sdbg.h:116-121`）。
- `ComputeIncomings`/`ComputeOutgoings`：模板 flag（`kFlagMustEq0/1/WriteOut`）
  支持"收集入边/出边"、"判度是否为 0/1"三种模式（`sdbg.h:245-320`），供
  `EdgeIndegree/EdgeOutdegree/UniqueNextEdge/UniquePrevEdge/NextSimplePathEdge`
  等调用——`NextSimplePathEdge` 要求"唯一出边且其唯一入边"（`sdbg.h:419-427`），
  是 unitig 压缩"可继续走"的判据。
- `IndexBinarySearch`：把查询序列按前若干碱基（`prefix_look_up_` 前缀查找表）定位
  到 bucket 区间，二分（mid 处沿 `Backward` 逐碱基比对，遇 tip 改读标签），命中
  返回 `GetLastIndex(mid)`（`sdbg.h:141-207`）。
- `GetLabel`：从边 id 恢复 (k+1)-mer 全标签（`sdbg.h:214-231`）。
- `EdgeReverseComplement`：取标签 → 翻转 → 互补 → `IndexBinarySearch` 找到反向边
  （`sdbg.h:432-464`）。
- `SetInvalidEdge`/`invalid_`：`AtomicBitVector`，清洗阶段标记删除（`w==0` 的边
  加载时自动置无效，`sdbg.h:56-60`）。

> **前缀查找表**：`LoadFromFile` 建 `prefix_look_up_`（每 bucket 的 [first,last]
> 下标）与 `f_`（每字符累计边数），`rank_f_[a] = rs_last_.rank(f_[a]-1)` 供
> `Forward` 用（`sdbg.h:36-54`）。

### 3.2 bucket 元数据（`sdbg_meta.{h,cpp}`）

SdBG 按**边序列前若干碱基**分桶（`kBucketPrefixLength=8`，65536 桶），桶记录
`bucket_id/file_id/starting_offset/num_items/num_tips/num_large_mul`，以及全图
`item_count/tip_count/ones_in_last/w_count[]`。`FromBucketRecord` 排序后累加前缀
和（`accumulate_item_count`），按 file 顺序再排一次；`Serialize`/`Deserialize`
存成文本（`k`/`words_per_tip_label`/`num_buckets`/`num_files` + 每桶一行）。
`Read2SdbgS2` 用它把各线程写的 `*.sdbg.<tid>` 分片信息汇总。

### 3.3 k-mer 编码（`sequence/kmer.h` + `definitions.h`）

- 2-bit/碱基，`kCharsPerEdgeWord=16`（32-bit 字），`kBitsPerEdgeChar=2`，
  `kEdgeCharMask=0x3`（`definitions.h:33-36`）。
- `Kmer<N, WordType>` 模板：`N` 个字固定存，变长 k 用足够多的字；核心操作
  `ShiftAppend(c,k)`（向右移进）、`ShiftPreappend(c,k)`（向左移进，供 RC 同步）、
  `ReverseComplement(k)`（逐字位反转补 + 换序 + 移位校正）、`cmp`。
- `GenericKmer = Kmer<kUint32PerKmerMaxK, uint32_t>`，
  `kUint32PerKmerMaxK = (kMaxK+1+15)/16`，`kMaxK=255`（`definitions.h:45-49`、
  `sdbg_def.h:21`）——最大 k=255，`megahit_core kmax` 打印此值，驱动据此校验。

### 3.4 排序引擎（`base_engine.{h,cpp}`）

`BaseSequenceSortingEngine` 是 count / read2sdbg_s1/s2 / seq2sdbg 共用的
**外部多趟排序框架**，三段式（`base_engine.cpp:143-211` 的 `Run()`）：

1. **Lv0 准备**（`Initialize` + `Lv0CalcBucketSizeLaunchMt`）：reads 均分给线程，
   每线程统计各 bucket（65536 个 = 边序列前 8 个 2-bit 字符）的条数；
   `Lv0ReorderBuckets` 按条数**降序重排**（负载均衡，`base_engine.cpp:231-252`）；
   `AdjustMemory` 按 `host_mem`/`mem_flag` 折算 Lv1/Lv2 各装多少条
   （`base_engine.cpp:14-141`：`mem_flag=1` 适中、`0` 最小、其它用满）。
2. **Lv1 扫描**（`Lv1FindEndBuckets` → `Lv1FillOffsetsLaunchMt`）：每次取一批
   bucket，各线程把**所属 (k+1)-mer 的位置编码成 4 字节差分偏移**
   （`OffsetFiller::WriteNextOffset`，`base_engine.cpp:353-363`）；差分超过
   `kDifferentialLimit=(1<<31)-1` 的进 `special_offsets` 列表
   （`base_engine.h:181-204`）。按 bucket 大小自适应分批，保证 Lv1+Lv2 不超内存。
3. **Lv2 排序 + 后处理**（`Lv1FetchAndSortLaunchMt` → `Lv2Sort`）：每 bucket 用
   `OffsetFetcher` 逐条恢复完整偏移，`Lv2ExtractSubString` 拷出 (k+1)-mer 明文 +
   附加信息（前后碱基等），`substr_sort_`（`kmsort_selector`，GPU 或 CPU）排序，
   `Lv2Postprocess` 做计数/过滤/写边（见 §4）。

> **内存预算**：`AdjustItemNumbers` 反复迭代折减 Lv2 条数、按
> `bytes_per_lv2_item` 与 `kLv1BytePerItem=4` 反推 Lv1 条数，直到
> `num_lv2*4 <= num_lv1`（Lv1 至少 4 倍于 Lv2，保证扫描不吃紧，
> `base_engine.cpp:39-48`）。

## 4. 三段流水线（SdBG 构建 → 组装 → 迭代）

驱动主流程（`src/megahit` `main()`，`megahit:970-1038`）：

```
buildlib → (count | read2sdbg[kmin-1pass]) → build_graph[kmin] → assemble[kmin]
  ┌─────────────────────────────┐
  while cur_k < k_max:          │
    local_assemble(cur_k)       │  ← contig 端点本地延伸
    iterate(cur_k, step)        │  ← 提取迭代边（桥接多 k）
    build_graph(next_k, cur_k)  │  ← 用上一轮 contig/bubble/local 建更大 k 的图
    assemble(next_k)            │
  merge_final(k_max)            │  ← cat *.final.contigs.fa + filterbylen
```

### 4.1 count：k-mer 计数（`kmer_counter.{h,cpp}`）

`KmerCounter` 统计所有 **(k+1)-mer** 的频次，过滤出 solid（`>= solid_threshold`，
`--min-count` 默认 2），把边 + 多重度写进 `<prefix>.edges.<tid>`：

- `Initialize`：读 reads 进 `SeqPackage`，算 `words_per_substr_`（(k+1)-mer 字数）
  与 `words_per_edge_`（(k+1)-mer + 多重度字数），malloc
  `first_0_out_/last_0_in_`（每 read 首末 (k+1)-mer 在 edges 里的下标，供下一轮建
  图连边），初始化 `EdgeWriter`（`kmer_counter.cpp:60-112`）。
- `Lv0CalcBucketSize`/`Lv1FillOffsets`：对每个 (k+1)-mer 取其 canonical
  （`rev_edge.cmp(edge)<0` 取 RC），按前 8 字符分桶，记下 `EncodeOffset`
  （`read_id<<1 | strand`，`kmer_counter.cpp:114-206`）。
- `Lv2ExtractSubString`：恢复每个 (k+1)-mer 明文 + `read_info`（低位存前后碱基，
  RC 时互补，`kmer_counter.cpp:208-252`）。
- `Lv2Postprocess`：对排序后相同边的连续段 `[from_,to_)` 统计 `count`；用
  `count_prev/count_next` 统计前后碱基分布以决定该边是否"某个 read 的首/末
  (k+1)-mer"（更新 `first_0_out_/last_0_in_`）；`count >= solid_threshold` 的边
  `PackEdge` 写入 edges 文件（`kmer_counter.cpp:254-…`）。
- `PackEdge`：把 (k+1)-mer 拷到 `words_per_edge_` 个字的尾部，多余位清零，
  `dest[words_per_edge_-1] |= min(count, kMaxMul)`（`kmer_counter.cpp:32-52`）。

### 4.2 read2sdbg：1pass 建 SdBG（`read_to_sdbg_s1.cpp` / `s2.cpp`）

`--kmin-1pass`（或 `--min-count 1` 强制）时，不走 count+seq2sdbg，而是直接
`read2sdbg` 两阶段建 k_min 的图：

- **S1**（solid 标记）：同 count 一样数 (k+1)-mer，但结果用来给每个
  `<read_id,offset>` 打 `is_solid` 位（`SeqPkgWithSolidMarker.is_solid`），并把
  低丰度但能"桥接两个 solid 区"的 (k+1)-mer 记为 mercy 候选（`mercy_files_`，
  待 S2 加 mercy 边，`--need_mercy`）。
- **S2**（写图）：`words_per_substr_ = k*2 + kBWTCharNumBits(3) + 1` 位（k-mer +
  BWT 字符 + 方向位）；排序后按 (k-1)-mer 分组，写 `w/last/tip` 位向量；两 solid
  区之间用 **`$` 虚拟节点**（`kSentinelValue`）连接成 SdBG，必要时补 mercy 边
  （`Extract_a` 里 `non_dollar` 位区分真实碱基与哨兵，`read_to_sdbg_s2.cpp:69-89`）。

### 4.3 seq2sdbg：迭代建图（`seq_to_sdbg.cpp`）

`seq2sdbg -k kmer_size --kmer_from prev_k --contig ... --bubble ... --addi_contig
--local_contig --input_prefix prefix.edges.*`：把上一轮的
`contigs.fa / bubble_seq.fa / addi.fa / local.fa` 与 count 产出的 edges 合并，
用 `BaseSequenceSortingEngine` 排序去重，输出新 k 的 SdBG。它实现了
"多 k 迭代"的核心：**上一轮 contig 的碱基作为下一轮 (k+1)-mer 的来源**，从而把
低 k 已确定的骨架带进更高 k，重复区在更高 k 处逐渐可分。

### 4.4 assemble：unitig 图 + 清洗 + 输出（`main_assemble.cpp` + `assembly/`）

1. **加载 SdBG** → `sdbg_pruning::RemoveTips`（SdBG 层先删一波尖端）→
   `UnitigGraph graph(&dbg)` 建 unitig 图（§3.1 的 `NextSimplePathEdge` 压缩路径）。
2. **清洗轮**（`cleaning_rounds` 默认 5，`main_assemble.cpp:182-249`），每轮按序：
   - **tips**（round>1 才做）：`RemoveTips(graph, max_tip_len)`（§5 之 Tip）；
   - **naive bubbles**（`bubble_level>=1`）：`NaiveBubbleRemover::PopBubbles`；
   - **complex bubbles**（`bubble_level>=2`）：`ComplexBubbleRemover::PopBubbles`；
   - **weak links**：`DisconnectWeakLinks(graph, disconnect_ratio)`；
   - **excessive pruning**：`prune_level>=3` 时 `RemoveLowDepth(graph, min_depth)`
     （全图低于 `min_depth` 的 unitig 永久删）+ 再清一轮 bubble；
     `prune_level>=2` 时 `RemoveLocalLowDepth(graph, min_depth, max_tip_len,
     local_width, min(low_local_ratio,0.1))`（§5 之 LowDepth）。
   - 一轮内无任何变化（`!changed`）即提前终止。
3. **输出**：`OutputContigs` 写 `contigs.fa`（非最终轮）/ `addi.fa`（prune 后残余）
   / `final.contigs.fa`（`output_standalone` 或最终轮）；`careful_bubble` 时把被
   合并的气泡写 `bubble_seq.fa`（`main_assemble.cpp:251-301`，`contig_output`）。
   `careful_bubble` 由驱动**仅在非最终轮**附加、最终轮加 `--is_final_round`、
   `--no-local` 加 `--output_standalone`（`megahit:891-898`）。
   `min_standalone` 由驱动算：`max(min(k_max*3-1, min_contig_len*1.5),
   min_contig_len)`（`megahit:867`）。

### 4.5 local：本地组装（`main_local_assemble.cpp` + `localasm/`）

对上一轮 contig 的**端点**做局部延伸（C++ 默认 `kmin=11, kmax=41, step=6,
seed_kmer=31, sparsity=8, similarity=0.8, min_mapping_len=75`），完整流程：

1. `HashMapper` 建索引 + 回帖（`hash_mapper.cpp`）：
   - `LoadAndBuild`：contig 每 `sparsity=8` 取一个 seed k-mer（长度 `seed_kmer=31`）
     建 `index_`（canonical key → `EncodeContigOffset(contig_id, offset, strand)`，
     `hash_mapper.cpp:56-101`）；**同一 seed k-mer 出现多处时把 value 高位置 1 标记
     "多定位"**，回帖时跳过——只允许唯一回帖。
   - `TryMap`（`hash_mapper.cpp:135-268`）：`len<50` 或 `len<seed_kmer` 直接丢弃；
     逐碱基滑 seed（fwd/rc 同步），命中唯一索引后推出 query 在 contig 上的
     `[contig_from, contig_to]`，被剪贴（clip）的对齐要求长度 `>= min_mapped_len=75`；
     对所有候选做 `Match`（16-mer 批量比对 + POPCNT 数错配，相似度 `>= 0.8`，
     `hash_mapper.cpp:103-133`），取**唯一最佳**——并列最佳视为不可靠丢弃。
2. `EstimateInsertSize`（`local_assemble.cpp:83-138`）：抽样最多 2^18 对双端 read，
   只统计"配对双方回帖到**同一 contig 且不同链**"的对，且 `insert >= 两 read 长度`
   才入直方图；`Trim(0.01)` 去尾后取 mean/sd。
3. `LocalRange`（`local_assemble.cpp:140-153`）：基线 `max_read_len-1`；仅当 paired
   且 `mean >= max_read_len` 时取 `min(2*mean, mean+3*sd)`；一律封顶 `kMaxLocalRange=650`。
4. `MapToContigs` + `AssembleAndOutput`（`local_assemble.cpp:166-302`）：按 contig ×
   strand 收集落在端点局部区间的 read（`AddSingle`/`AddMate`）；`min_num_reads =
   local_range/max_read_len` 以下不值得组装；**同一映射位置最多取 3 条 read**
   （`pos_count<=3`，防覆盖度虚高）；`contig_end` = contig 端点 `local_range` 的序列。
   > 驱动只传 `--kmax`（= 下一轮的 k，`local_assemble(cur_k, next_k)`，
   > `megahit:908-913,1010`），其余用 C++ 默认——即 IDBA 本地组装的 k 区间实际是
   > `[11, next_k)`，专门补"下一轮 k 之前"的缺口证据。
5. `LaunchIDBA`（`local_assemble.cpp:28-81`）：对局部 read + `contig_end` 做多 k
   迭代组装（k 从 `kmin` 到 `min(kmax, max_read_len)`，`step` 递增）；每 k 用
   `HashGraph` 插入 k-mer，覆盖度直方图 `percentile(1 - local_range/num_vertices)`
   定阈值，`Assemble` 后经 `ContigGraph::RemoveDeadEnd(k*2)` → `RemoveBubble` →
   `IterateCoverage(k*2, 1, threshold)` 清洗；若只剩 1 条 contig 提前终止。
   产出延伸后的 `local.fa`。

### 4.6 iterate：迭代边提取（`main_iterate.cpp` + `iterate/`）

`iterate -c contigs.fa -b bubble_seq.fa -r reads.bin -k cur_k -s step -o prefix`：
用上一轮 contig/bubble 的端点 k-mer 建 `ContigFlankIndex`，扫 reads（`KmerCollector`
），把"能从 contig 端继续延伸的 (k+step+1)-mer"作为**迭代边**写
`<prefix>.edges.0`，供下一轮 `build_graph` 用（`main_iterate.cpp:117-222`）。
校验：`step` 为偶数且 `1<=step<=28`，`kmer_k+step` 小于
`max(Kmer<4>::max_size(), GenericKmer::max_size())`
（`main_iterate.cpp:77-96`）。这就是 MEGAHIT 多 k 桥接的另一种方式——**用 reads
证据补出高 k 才有的边**，避免丢失低 k 被过滤的序列。

### 4.7 输出与后处理

- `final.contigs.fa`：`merge_final` 把各 k 的 `*.final.contigs.fa` + 最终 k 的
  `contigs.fa` `cat` 后经 `filterbylen <min_contig_len>` 过滤（`megahit:918-937`）。
- `contig2fastg`（`tools/contigs_to_fastg.cpp`）：把某 k 的中间 contigs 转 FASTG。
- `read_stat`/`filter_by_len`（`tools/`）：reads 统计 / 按长度过滤 contigs。

## 5. 图清洗算法细节

### Tip（`tip_remover.cpp`）

`RemoveTips`：对每个长度 `< thre` 的 unitig（thre 从 2 起**倍增**直到
`max_tip_len`，`tip_remover.cpp:10-11`）：
- standalone（loop）直接删；
- 无入无出删；
- 出度 1 入度 0（或反之），且邻居平均深度 `> 8 ×` 自身深度时删——用深度比避免
  删掉低覆盖的真实短段。每轮 `graph.Refresh(false)` 后重扫。

> SdBG 层 `sdbg_pruning::RemoveTips`（`sdbg_pruning.cpp`）在 unitig 建图前先删
> 一波：遍历边，沿 `NextSimplePathEdge`/`PrevSimplePathEdge` 走到端点，长度 < 阈值
> 的标记 `SetInvalidEdge`。`InferMinDepth`（`sdbg_pruning.cpp`）用边多重度直方图
> 估最小深度，`--min_depth`（`--prune-depth`，默认 2）未显式给时用它。

### Bubble（`bubble_remover.{h,cpp}`）

`BaseBubbleRemover::SearchAndPopBubble` 识别"气泡"：当前节点的出度 ≥ 2，每个中间
unitig 入度 1、出度 1，且都汇到同一右节点（`right.b() == possible_right[0].b()`，
并校验 `right.canonical_id() >= adapter.canonical_id()` 防重复处理），长度 ≤
`max_len`（`bubble_remover.cpp:58-101`）。把中间按平均深度**降序**排，用
`checker(middle[0], middle[j])` 验证"最佳分支与其它分支可合并"，满足则删掉
`middle[1..]`（`SetToDelete`），仅保留深度最大的代表分支（`bubble_remover.cpp:103-133`）。

- **NaiveBubbleRemover**：checker 默认要求两分支长度（含 k-1 重叠）相等
  （`bubble_remover.h` 中 naive checker，`a.GetLength() == b.GetLength()`）。
- **ComplexBubbleRemover**：checker 用带通配编辑距离 `GetSimilarity`
  （banded DP，`bubble_remover.cpp:10-54`），要求
  `(b.len+k-1)*sim <= a.len+k-1`（长度相似 + 序列相似 ≥ `sim`），`max_len =
  lround(merge_level_*k/sim)`（默认 `--merge-level 20,0.95` → ≈21k）。`SetToDelete`
  的 unitig 若深度 ≥ 最佳分支的 `careful_threshold_`（0.2），连同左右端写入
  `bubble_seq.fa`（`careful_bubble`，`bubble_remover.cpp:109-132`）——下一轮
  `build_graph` 把这些"被合并的气泡"序列也喂进去，避免丢失多态序列。

### Weak link（`weak_link_remover.cpp`）

`DisconnectWeakLinks`：对出度 ≥ 2 的 unitig（跳过 standalone/回文），在**正反两
条链**上各看一遍，把"深度 ≤ `local_ratio`（0.1）× 邻居总深度"的邻居
`SetToDisconnect`（断开而非删除，`weak_link_remover.cpp:8-37`）。

断开后的实际处理在 `UnitigGraph::RefreshDisconnected`（`unitig_graph.cpp:140-208`）：
对每个 `to_disconnect` 的端点，`new_start = NextSimplePathEdge(old_start)`（沿出边走
**一格**）并把旧端点边 `SetInvalidEdge`；`new_length = old_length - 断开端数`，
`new_total_depth = lround(avg_depth × new_length)`；若 `length <= 断开端数` 则直接
`SetToDelete`——即断开是"**截短一格**"而非删整条，且维持平均深度不变。

### Low depth（`low_depth_remover.cpp`）

- `RemoveLocalLowDepth`：对长度 ≤ `max_len`、非 standalone 的 unitig，用
  `LocalDepth`（两端局部宽度 `local_width`（默认 1000）内**邻居**的深度加权均值，
  含反向互补，`low_depth_remover.cpp:10-35`）当参照。阈值取
  `min(min_depth, mean*local_ratio)`（代码写法：先 `threshold=min_depth`，若
  `min_depth < mean*local_ratio` 保持 `min_depth`、否则 `threshold=mean*local_ratio`，
  `low_depth_remover.cpp:63-68`）——即删 `depth < min(min_depth, mean*local_ratio)`
  的 unitig；且仅当 `indegree<=1 && outdegree<=1` 或入/出度为 0 时才考虑
  （`low_depth_remover.cpp:58`）。迭代时（`IterateLocalLowDepth`）`min_depth *= 1.1`
  递增，直到 `min_depth >= kMaxMul` 或一轮无变化（`low_depth_remover.cpp:88-102`）。
- `RemoveLowDepth`（prune_level≥3 的 excessive pruning）：全图删 `avg depth <
  min_depth` 的 unitig（`low_depth_remover.cpp:104-117`）。

## 6. 参数语义（`src/megahit` 驱动，`megahit:39-105` 与 `check_and_correct_option`）

| 参数 | 默认 | 语义 |
|---|---|---|
| `-1/-2/--12/-r` | — | 双端/交错/单端 reads（逗号分隔多库，gz/bz2 自动解压） |
| `-o/--out-dir` | `./megahit_out` | 输出目录；`--out-prefix` 改文件名 |
| `-m/--memory` | 0.9 | 内存（0-1 为总内存比例，否则字节数） |
| `-t/--num-cpu-threads` | 逻辑核数 | 线程；超硬件上限 WARNING 并钳制 |
| `--min-count` | 2 | solid (k+1)-mer 最小频次；==1 强制 `kmin-1pass`+`no-mercy` |
| `--k-list` | `21,29,39,59,79,99,119,141` | 显式 k 列表（全奇数，15..kmax，相邻差≤28） |
| `--k-min/--k-max/--k-step` | 21/141/10（帮助文本误写 `[12]`） | 替代 `--k-list`（step 偶数 ≤28） |
| `--kmin-1pass` | off | 1pass 建 k_min 图（低深度省内存） |
| `--no-mercy` | off | 不加 mercy 边 |
| `--max-tip-len` | -1(=2k) | 尖端长度上限 |
| `--bubble-level` | 2 | 气泡处理强度 0-2（驱动钳制） |
| `--merge-level` | 20,0.95 | 复杂气泡长度≤`l*k`、相似度≥s 才合并 |
| `--prune-level` | 2 | 低深度修剪强度 0-3 |
| `--prune-depth` | 2 | `min_depth`（unitig 平均 k-mer 深度下限） |
| `--disconnect-ratio` | 0.1 | 弱连接断开比例（[0,0.5]） |
| `--low-local-ratio` | 0.2 | 局部低深度比例 ((0,0.5]) |
| `--cleaning-rounds` | 5 | 清洗轮数（≥1） |
| `--no-local` | off | 禁用本地组装 |
| `--min-contig-len` | 200 | 最终 contig 最短长度 |
| `--presets` | — | `meta-sensitive`（min-count 1、k 21..141 步 10）/ `meta-large`（k 27..127 步 10） |
| `--continue` | off | 从 checkpoint 续跑（`checkpoints.txt` + `options.json`） |
| `--no-hw-accel` | off | 禁用 POPCNT/BMI2 硬件加速 |
| `--keep-tmp-files/--tmp-dir` | off | 保留中间文件 / 自定临时目录 |
| `--test` | off | 用 `test_data/` 玩具数据跑通 |

> 校验（`megahit:487-569`）：k 全奇数且在 [15, kmax]；相邻 k 差 ≤ 28；
> `k-step` 偶数；`prune-level∈[0,3]`；`merge-level` 相似度 ∈ [0,1]；
> `disconnect-ratio∈[0,0.5]`；`low-local-ratio∈(0,0.5]`；`cleaning-rounds≥1`；
> `min-count>0`。`set_max_k_by_lib` 会按 `max_read_len+20` 截掉过大的 k
> （`megahit:757-769`）。**k_max 实际封顶 255**（`kMaxK`，`megahit_core kmax`）。

## 7. 源码 quirks（异常/边界行为）

- **k=255 是编译期硬上限**：`kMaxK=255`（`sdbg_def.h:21`），`GenericKmer` 用
  `kUint32PerKmerMaxK=(255+1+15)/16=16` 个字（`definitions.h:45`）；驱动校验
  `k_list[-1] <= kmax`。
- **k 必须奇数**（驱动强制，`megahit:527-529`）；`iterate` 的 `step` 必须偶数且
  ≤28，`kmer_k+step < max(Kmer<4>::max_size(), GenericKmer::max_size())`
  （`main_iterate.cpp:92`）。
- **`--k-step` 帮助与代码默认不一致**：帮助文本 `[12]`（`megahit:61`），代码
  `k_step = 10`（`megahit:170`）；默认 `k_list` 显式给定，该默认仅在用户改用
  min/max/step 生成 k 列表时才生效。
- **驱动 `assemble_cmd` 重复传参**：`--cleaning_rounds` 在同一条命令里出现两次
  （`megahit:879,882`），值相同无影响——复制粘贴痕迹。
- **`seq2sdbg` 要求 k≥9**（`main_sdbg_build.cpp:203-205`）；`count`/`read2sdbg`
  要求显式 `--host_mem` 非 0（`main_sdbg_build.cpp:70-72,124-126`）。
- **`--min-count 1` 的隐式联动**：驱动设 `kmin_1pass=True` 且 `no_mercy=True`
  （`megahit:541-543`）——因为 1pass 模式下 S1 已按 1 阈值标记 solid，无需再走
  count+seq2sdbg 与 mercy。
- **`bubble-level` 语义漂移**：C++ `main_assemble` 帮助写 "0-3"
  （`main_assemble.cpp:79`），驱动只允许 0-2 并钳制（`megahit:564-569`）——
  实际运行以驱动为准。
- **清洗轮提前终止**：一轮内无任何 tips/bubbles/disconnect 变化就 `break`
  （`main_assemble.cpp:248`），所以 `cleaning_rounds` 是上限而非保证。
- **`--max-tip-len -1` 的三处默认不同**：SdBG 层 `RemoveTips`（
  `sdbg_pruning`）与 unitig 层（`main_assemble.cpp:143-145`）都是 `2k`；但驱动在
  `assemble()` 里若 `cur_k*3-1 > min_contig_len*1.5` 会把 `--max_tip_len` 设成
  `max(1, min_contig_len*1.5+1-cur_k)`（`megahit:886-887`）——当 k 较小时用
  `min_contig_len` 相关的阈值替代 2k。
- **`min_standalone` 计算**：`max(min(k_max*3-1, min_contig_len*1.5),
  min_contig_len)`（`megahit:867`）；`--max-tip-len>=0` 时改为
  `max(max_tip_len+k_max-1, min_contig_len)`（`megahit:868-869`）。
- **multplicity 饱和**：`mul_t=uint16`，计数超过 65535 封顶（`kMaxMul`，
  `PackEdge` 里 `min(count, kMaxMul)`）；`small_mul` 的 255 是"查 large_mul"哨兵
  （`sdbg_def.h:11-19`），所以常规小多重度不占 phmap 内存。
- **SdBG 加载自动置无效边**：`w==0` 的边在 `LoadFromFile` 里标记无效
  （`sdbg.h:56-60`）——0 是"哨兵字符"而非 A。
- **`FreeMultiplicity`**：可释放全部多重度数组换内存（调用后
  `EdgeMultiplicity` 失效，`sdbg.h:471-475`）——多轮迭代中已不需要旧图多重度。
- **checkpoint 续跑**：`Checkpoint` 装饰器每步完成写 `checkpoints.txt` 一行
  `<id>\tdone`，续跑时跳过已完成的步骤（`megahit:251-281,443-452`）；`--continue`
  忽略除 `-o` 外的所有选项（从 `options.json` 恢复）。
- **`EarlyTerminate`**：`build_graph` 发现 `file_size==0 && kmer_from!=0`
  （上一轮没有 contig/bubble 可喂）就提前结束迭代、直接 `merge_final`
  （`megahit:839-840,1019-1020`）——防止空转。
- **CPU 分发**：`checkcpu`（POPCNT+BMI2）/`checkpopcnt` 探测后选
  `megahit_core(_popcnt|_no_hw_accel)` 三个二进制之一（`megahit:613-631`）。
- **大差分偏移**：Lv1 差分超过 `2^31-1` 进 `special_offsets`，`xfatal` 于
  "Too many large difference items!"（`base_engine.h:184-185`）。

## 8. 与 anchr 的关联

anchr 现状：`asm contig/unitig` 用 pgr 的 `KmerTable`（canonical 2-bit u128、
精确计数、radix sort、rayon 并行）在**内存中**建 k-mer 图并做 unitig 压缩；
`asm olc` 走"多 k unitig → 精确 overlap → 贪心布局 → 一致序列"；`asm map`
是完美回帖。MEGAHIT 提供的是**同目标、不同工程路线**的完整参照：

1. **外部排序 + 分桶计数（§4.1，`base_engine`）**：MEGAHIT 用 65536 桶、Lv0 统计
   桶大小按降序分批、Lv1 差分偏移、Lv2 排序计数，全程受 `--memory` 约束——这是
   "低内存超大基因组"的教科书范式。anchr 的 `KmerTable` 全内存路线在基因组规模
   超内存时，可借鉴其"桶大小统计 → 降序分批 → 差分偏移"的落盘策略，而不必照搬
   GPU 排序。
2. **succinct 位向量图（§3.1）**：SdBG 用 `w/last/tip` 三个紧凑位向量 + rank/
   select 表达 de Bruijn 图，边 id 即下标，`Forward/Backward` 都是位操作；"tip 边
   的标签单独存"（避免从死端回退）与"前缀查找表二分"（`IndexBinarySearch`）是
   两个可独立移植的省内存技巧。anchr 若未来把 `asm unitig` 的图也落到磁盘格式，
   这套"排序边数组 + 位向量 + rank/select"是成熟模板（对比 SKESA 的
   total+branch+plus-fraction 打包，`skesa.md` §7.1）。
3. **unitig 压缩判据（§3.1 `NextSimplePathEdge`）**：MEGAHIT 的"唯一出边 + 该出边
   唯一入边"才继续压缩，与 anchr `asm unitig`（BCALM `graph3`）"唯一 solid 后继
   且其前驱唯一"语义一致（`docs/asm.md` unitig 一节）——可作交叉验证判据。
4. **清洗算法族（§5）**：MEGAHIT 把 tip（倍增阈值 + 8× 深度比）、naive/complex
   bubble（banded 编辑距离相似度）、weak link（0.1× 总深度断开）、low depth
   （局部窗口深度均值，`min_depth*=1.1` 迭代）做成一套可调参数——anchr `asm
   contig` 目前只做种子扩展 + 泡泡消除，若扩展清洗，这组默认值（tip=2k、
   disconnect-ratio=0.1、low-local-ratio=0.2、merge 20k/0.95）与"深度比防误删"
  的思路是现成参照。
5. **多 k 迭代 + 上一轮产物引导（§4.3/§4.6）**：MEGAHIT 用上一轮 contig/bubble/
   local contig 直接喂 `seq2sdbg`，另用 `iterate` 从 reads 补 (k+step+1)-mer 边；
   anchr `asm olc` 也是多 k，但目前各 k 独立出 unitig、无反馈（`notes/design/
   asm-olc.md`）——MEGAHIT 的"引导 + 迭代边"是 v2 反馈环的直接素材（与 SKESA 的
   `clean_reads` 反馈、metaMDBG 的 unitig 反馈同族）。**`asm multik`（2026-08-14）
   已实现反馈环**（metaMDBG 式图级反馈），与 megahit 的引导机制对比见 §8.6。

### 8.6 与 `asm multik` 的关联（2026-08-14 补充）

`asm multik` 已实现（借鉴 metaMDBG，`notes/design/asm-multik.md`），megahit 的
multi-k 迭代与之**同族但机制不同**：

| | megahit | `asm multik` |
|---|---|---|
| 迭代步 | k 列表 21→141（+8/+10/+20），每轮**重建**更大 k 的图 | auto 按 read N50 推导（`k_max=0.8×N50` clamp 31..256，150 bp→[50,70,90,110]；模板显式 KS=31..192），unitig 图结构保留 |
| 上一轮产物 | **序列级引导**：contigs/bubbles 喂 `seq2sdbg` 建新 k 图 | **图级反馈**：unitig 图 + compute_links 边保留，跨接验证选边 |
| reads 桥接 | `iterate`：contig 端点索引 + reads 回帖，提取跨端点的 (k+step+1)-mer 迭代边（**建图素材**） | `bridge_filter`：60-mer 探针验证 unitig 间连接（**验证边**） |
| unitig 生成 | 每轮 `assemble` 重新压缩（`NextSimplePathEdge`） | pass 0 一次，后续轮只验证/压实（不重新 unitig 化） |

**关键差异**：megahit 每轮"重建 + 引导"（旧 contigs 是新图的种子），multik 每轮
"验证 + 压实"（旧 unitig 图是验证对象）。两者都解决"更大 k 特异性的接入"：
megahit 用 iterate 边把跨 contig 的 reads 证据带进新图；multik 用 bridge_kmer/
探针验证跨 unitig 的连接。**megahit 的 iterate（contig 端点索引 + reads 回帖）
与 multik 的 bridge_filter（探针桥接）是同族机制**——megahit 在建图侧、multik
在验证侧。

**可借鉴**（multik 视角）：megahit 的引导把"已确定的骨架"（contigs）直接带进
下一轮，multik 的图级反馈保留了结构但**没利用 unitig 序列做更大 k 的计数引导**
（multik 的 count_at 把 unitigs 序列作为输入计数，但迭代边提取——跨 unitig 的
更大 k 连接——靠 bridge_kmer 验证而非"收集"）。若 multik 未来需要"把 contig
序列喂进下一轮建图"（megahit seq2sdbg 式），可用 `TadpoleTable::build_supermer`
对 unitigs 序列直接计数（已支持多序列输入）——这是 megahit 引导在 multik 的
最小映射。**已实现（2026-08-14，v6）**：multik 的渐进过滤删的低丰度分支作为
下一轮 unitigs 回灌（megahit bubble 回灌 + metaMDBG unitig 反馈结合），
G37 最长 contig +20%、misassemblies 保持 0（`notes/design/asm-multik.md`
§4.11）。

**可变 k 机制精读（2026-08-16 晚，针对 pgr K 上限 256 的对照）**：
* **k 序列**：megahit 默认 `[21,29,39,59,79,99,119,141]`（`src/megahit`
  `Options.__init__`），auto 模式只按读长裁剪——`set_max_k_by_lib` 把
  `k >= max_read_len + 20` 的项去掉，无更复杂的推导；相邻 k 差强制 ≤28
  （`--k-step` 校验，step 必须偶数、k 必须奇数）。
* **kMax = 255**（`sdbg/sdbg_def.h` `kMaxK`，`megahit_core kmax` 查询）——
  pgr 的 `Kmer::MAX_K = 256` 与其几乎一致（上游 128→256 是对标 megahit）。
* **迭代引导的两半**（缺一不可）：`seq2sdbg -k next --kmer_from cur` 的输入
  同时包含 **上一轮 contigs 全长**（`--contig`/`--bubble`/`--addi_contig`）
  和 **iterate 迭代边**（`--input_prefix`）：`iterate/contig_flank_index.h`
  `FeedBatchContigs` 存每个 contig 端点 (k+1)-mer + 延伸序列（`ext_seq`，
  最多 step-1 碱基），`FindNextKmersFromReads` 在 reads 里锚定端点后收集
  跨端点的 (k+step+1)-mer 作新 k 图的种子。
* **实验结论（G37 MRX40P000，450 bp reads）**：把上一轮 unitigs 的反馈从
  "全长"改成"仅两端各 k 片段"（简化版端点引导），或把 K31 端点片段重复喂给
  K160 主 K，unitig 输出均与基线完全一致（752/12,310；K160 100/72,998）
  ——高覆盖短读下 k-mer 本来就在 reads 里，碎片化是图结构分支而非缺种子，
  简单引导无增量。真正的 megahit 式引导（seq2sdbg 重建 + 迭代边）价值在
  低覆盖/长读场景（reads 支持不足，需低 k 结构引导），multik 暂无此类数据
  验证；`count_at` 保持全长反馈（与 megahit `--contig` 全长引导一致）。
* **引导实施与端到端（2026-08-16 深夜—08-17，已落地）**：`asm multik` 新增
  `--guide-contigs <fasta>`（cmd 层把上一主 K unitigs 全长作为伪 reads 写
  临时文件，重复到 solid 阈值后并入计数，对应 megahit `seq2sdbg --contig`）。
  单主验证：K192 用 K31 unitigs 引导，unitig N50 37.6K→**81.6K**（+117%）、
  K224 16.4K→37.9K、K256 6.0K→25.9K——**全长引导有效**（reads 高 k 覆盖
  不足时低 k 结构真实起引导作用）；500/1000 bp 端点片段引导无效（结构信息
  不足）。端到端首轮（G37 7 组，31..192 + 只 K192 引导）：GF 98.849→
  99.602% 但 Dup 1.210 + 1 mis。**三项修复后全链：N50 83.8K→318.1K、
  GF 98.849→99.642%、Dup 1.000、0 mis**：
  1. `consensus::coverage` 加 dominant 31-mer offset + banded identity 路径
     （容忍跨组 consensus 的大 indel——120 bp 插入让旧 3 种子锚点漏尾部，
     跨主 K 近似重复去重失效，Dup 1.210→1.001）；
  2. `multik` 单主模式也执行 `remove_unsupported`（内部 k-mer solidity +
     连续 unsupported 窗口判嵌合）——切断引导传播的嵌合 junction；
  3. `olc` 布局增加**目标端竞争检测**（`is_repeat` 不只查当前端，拼接前也
     查目标端——anchor_6 尾有 3 条 ~1.07K 竞争 overlap，mutual-best 只看
     best 指回漏掉，任选一条拼成 350473-389481+390412-429551 relocation）。
  模板：`6_unitigs` multik 分支 `KS` 含 192；K31 主 K 先跑（K_LIST 到
  192），其余主 K 并行 + `--guide-contigs unitigs_K31.fasta`。MG1655 验证
  因大 reads 计算成本过高暂缓（K192 引导 ~20 min/组 + 下游 olc 更慢），
  G37 端到端已支撑落地。
* **性能优化：验证轮稀疏化（2026-08-17）**：模板 multik 的 K_LIST 原来是
  "当前主 K + 所有更大 k"（11 个 k → 10 个验证轮，66 次全量计数/组——
  MG1655 单组 20+ min）。对照 SKESA（`steps=11` 但增量延伸）与 megahit
  （8 个稀疏 k）后实验：验证轮从 11 个减到 3 个（`VERIFY_KS="71 121 192"`），
  **单主 unitigs 完全一致**（G37 K31：751/12,310 相同）、时间 2m42s→
  **1m14s（2.2×）**；端到端（G37 7 组，稀疏验证 + K192 引导）质量保持
  （N50 317.5K、GF 99.563%、Dup 1.001、0 mis，vs 全量 318.1K/99.642/1.002）。
  配套修正 `coverage()` dominant-offset 的峰比门槛 0.7→0.6（120 bp indel
  会把 histogram 分成 66%/34% 两峰，主峰 ratio 0.664 被旧门槛跳过导致
  跨主 K 近重复去重失效、Dup 1.147；identity ≥99% 兜底）。
6. **本地组装（§4.5）**：MEGAHIT 用 reads 回帖 + IDBA-UD 内核做 contig 端点延伸，
   与 anchr `fq extend`（tadpole 沿图延伸）目标相近；其 `min_mapping_len=75`、
   `LocalRange = min(2*mean, mean+3*sd)` 封顶 650 等参数可对照。
7. **Python 驱动 + checkpoint（`src/megahit`）**：多步流水线用 checkpoint 文件
   断点续跑、`options.json` 固化参数——anchr 的模板流水线（`template` 命令）做
   大流程编排时可借鉴这种"每步可跳过、续跑只重做失败步"的工程形态。
8. **工程教训**：MEGAHIT 大量用 `phmap::parallel_flat_hash_map`（multiplicity、
   UnitigGraph 的 `id_map_`）、`AtomicWrapper`（位域原子 flag）、`#pragma omp
   parallel for reduction` 做无锁并行扫描——与 anchr 的 rayon 并行习惯可互证；
   但 MEGAHIT 的多重度是**饱和 uint16**（丢失 >65535 的精确计数），而 anchr
   `KmerTable` 是精确 u64 计数，这点上 anchr 更精确。

### 8.1 若要复刻类似效果，anchr 需要什么（差距清单）

把上面对照转成**落地清单**。前提：遵循 AGENTS.md——复杂逻辑进 `libs/`、`cmd/`
薄壳、**尽量不加新依赖**。anchr 现状：`asm contig/unitig` 走 `libs/asm/assemble.rs`
（tadpole contigMode，`TadpoleTable = HashMap<Kmer,u32>` 全内存），`asm olc` 走
`libs/olc/`（多 k unitig → overlap → layout → cns），`asm map` 是 `libs/map.rs`
完美回帖；pgr `KmerTable` 提供 packed bytes + radix sort + rayon 精确计数。
**anchr 的 `Kmer::MAX_K` 跟随 pgr（2026-08-16 起 256，对标 MEGAHIT 的 255）**
——k/step 校验动态取 `Kmer::MAX_K`，无需硬编码。

| 效果 | anchr 现状 | 需要补什么 | 放置 / 门槛 |
|---|---|---|---|
| 单 k 内存图组装 | ✅ 已具备（`assemble.rs`） | 无需补；MEGAHIT 的 unitig 判据（§3.1 `NextSimplePathEdge`）与 `assemble.rs` 的 `unique_solid_out/in` 语义一致，可加**交叉验证测试** | — / 低 |
| 图清洗（tip / weak link / low depth） | ⚠️ 只有 `pop_bubbles` | 新增 `libs/asm/clean.rs`：`RemoveTips`（倍增阈值 + `8×` 深度比，§5 Tip）、`DisconnectWeakLinks`（0.1× 总深度断开）、`RemoveLocalLowDepth`（端点局部宽度窗口深度参照）；深度直接复用 `Unitig` 的 `coverage/min_cov/max_cov`，邻接复用现有 `Link`/`EdgeRef` | `libs/asm/clean.rs` / 中 |
| 多 k 迭代 + 上一轮引导 | ⚠️ multik 图级反馈（unitigs 全长回灌计数 + 验证边），无 seq2sdbg 式"旧 contigs 建新 k 图" | 每轮 unitigs 序列回灌为下一轮 k 的输入（对应 `seq2sdbg`，§4.3，multik 已做全长回灌）；迭代边提取（对应 `iterate`，§4.6）用端点索引 + reads 锚定收集跨端点 k-mer——高覆盖短读验证无增量（§8.6 实验结论），价值在低覆盖/长读，待数据验证 | `libs/asm/iterate.rs` / 中 |
| 本地组装（contig 端点延伸） | ⚠️ `fq extend` 做 reads 延伸、`map.rs` 完美回帖，未组合 | 新增 `libs/asm/local.rs` 编排：`map.rs` 回帖（对应 `hash_mapper`）→ 按 insert 取端点局部区间（`min(2*mean, mean+3*sd)` 封顶 650，§4.5）→ 区间内跑 `assemble`（对应 IDBA-UD 内核）延伸端点 | `libs/asm/local.rs` / 低-中 |
| **低内存超大基因组（外部排序计数）** | ❌ `TadpoleTable`/`KmerTable` 全内存 | 新增 `libs/asm/externalsort.rs`：65536 桶（前 8 碱基）→ Lv0 统计桶大小降序分批 → Lv1 4 字节差分偏移（>2^31 进 `special_offsets`）→ Lv2 每桶取出明文排序计数、滤 `solid_threshold`（§4.1）；edges 按桶前缀分片落盘。排序用 std / 已有 radix，**无需新 crate** | `libs/asm/externalsort.rs` / **重（核心里程碑）** |
| succinct 位向量图（SdBG） | ❌ 无对应（图全内存） | 新增 `libs/asm/sdbg.rs`：自写 `Vec<u64>` + popcount 前缀表的 rank/select（约几百行，避免引入 `sucds` 等新依赖），`w/last/tip` 三数组 + RC 位 + `tip_lables` 旁存 + 前缀查找表二分（§3.1）。独立于外部排序，做磁盘图格式时才需要 | `libs/asm/sdbg.rs` / 重 |
| checkpoint 断点续跑 | ⚠️ `template` 命令生成脚本，无内置 checkpoint | 极简方案：`template` 产出的 shell 沿用"每步完成写一行 + 续跑跳过"模式（`megahit:251-281`），不进 Rust 核心 | `cmd/`（编排层）/ 低 |
| 多重度精确性 | ✅ anchr 精确 u32/u64 已胜过饱和 uint16 | 无需做；若复刻 SdBG，保留 `small_mul + large_mul` 两级思路但把上限提到 u32 | — / 低 |

> **落地顺序建议**：算法语义优先——②清洗 → ③多 k 反馈 → ④本地组装；低内存
> 里程碑（⑤外部排序 → ⑥SdBG）独立、按需启动，⑤是其中最有价值的一块（把 anchr
> 从"内存大小受限"里解放出来）。每一块的**成功判据**都建议先写测试：清洗对照
> MEGAHIT 默认值（tip=2k、disconnect-ratio=0.1、low-local-ratio=0.2，§6）在
> `tests/` 小数据集上对比；外部排序计数与 `KmerTable` 精确计数**逐 k-mer 一致**。

### 8.2 移植注意事项（本次精读确认，Rust 端）

- **"断开"是截短一格，不是删整条**：`DisconnectWeakLinks` 只标 `to_disconnect`，
  真正生效在 `RefreshDisconnected`（`unitig_graph.cpp:140-208`）——端点沿
  `NextSimplePathEdge` 走一格、旧端点边标 invalid、`length-1`、`total_depth =
  lround(avg × new_len)`（保持平均深度）。Rust 端若用可变图，需同时支持"删顶点"
  与"截短端点"两种操作，且截短时别把 total_depth 当新深度覆盖。
- **low depth 阈值是 `min` 不是 `max`**（易写反）：`depth < min(min_depth,
  mean*local_ratio)` 才删；且只在 `indeg+outdeg` 小的节点上判（`low_depth_remover
  .cpp:58`）。移植时对照 §5 Low depth 的实现写法。
- **SdBG 只读结构不必照搬**：MEGAHIT 因图是紧凑只读的，清洗只能靠位标记 + 每轮
  `Refresh` 全图重扫两遍（disconnect 一遍 + delete 一遍）。anchr 的 `assemble.rs`
  用可变 `HashMap`/`Vec`，可"原地标删 + 惰性清理"，比照搬 Refresh 更简单——只借
  清洗**判据与阈值**，不借重构机制。
- **local 回帖复用 `map.rs`**：`HashMapper` 丢弃多定位（value `>>63` 标记）与并列
  最佳（视为不可靠），本质是"只信唯一最佳"；anchr `map.rs` 完美回帖更严格，本地
  组装应复用 `map.rs` 而非照搬 seed-and-extend。可直接借的技巧：`sparsity=8` 抽样
  建索引、`pos_count<=3` 防覆盖度虚高、insert 只统计"同 contig 不同链"配对、
  `insert>=len` 入直方图、`Trim(0.01)` 去尾（§4.5）。
- **k 上限**：MEGAHIT 255 vs anchr `Kmer::MAX_K=256`（pgr，2026-08-16）——
  anchr 的 k/step 校验**动态取 `Kmer::MAX_K`**（multik 的 k 合法域
  `(1..=Kmer::MAX_K)`、`auto_ks` 的 `k_max` 亦 clamp 到 256），无硬编码数字。
- **多重度精度**：MEGAHIT 饱和 uint16（`kMaxMul`），anchr 是精确 u64 计数，更优；
  清洗阈值语义按 anchr 精确计数**重算默认值**，不照搬饱和行为。
- **unitig 判据交叉验证**：`NextSimplePathEdge`（唯一出边 + 该出边唯一入边，
  `sdbg.h`）与 `assemble.rs` 的 `unique_solid_out + unique_solid_in==1` 语义一致，
  可加"同一 k-mer 图两引擎产出相同 unitigs"的等价性测试锁定。

> 一句话：MEGAHIT 的**算法语义**（unitig 压缩判据、清洗阈值族、多 k 引导、本地
> 延伸）直接服务 anchr `asm` 的未来演进；它的**工程路线**（外部排序、位向量图、
> checkpoint）是 anchr 内存路线在"超大规模"场景下的对照与迁移蓝本，两者互补而非
> 竞争。

## 9. 局限

- **k 上限 255 且须奇数**：k>255 无法表示；中间迭代 k 用 `Kmer<N>` 模板实例
  逐个试选（`main_iterate.cpp:117-159`），k 越大字越多、位操作越慢。
- **多重度饱和**：边多重度 uint16 封顶，高覆盖（>65535×）时会丢失区分度（
  `kMaxMul`）。
- **两层架构复杂**：Python 驱动 + 3 个 C++ 二进制（hw-accel 分发），构建依赖
  git 子模块（kmlib/phmap/xxhash）；参数校验分散在驱动与各 main_* 里，存在
  `bubble-level` 帮助文本 0-3 与驱动 0-2 不一致这类漂移。
- **SdBG 是只读结构**：清洗靠位标记（`invalid_`）+ 每轮 `Refresh` 重扫，图本身不
  支持增量修改；`FreeMultiplicity` 释放多重度后无法恢复。
- **局部组装强依赖 IDBA-UD 内核**（`src/idba/` 一整块外部代码），占仓库体积大、
  与主数据结构（SdBG）不共享编码。
- **驱动对输入的校验是白名单式**：`getopt` 解析 + 若干 `raise Usage`，畸形输入
  走异常退出而非友好错误信息（对照 anchr 的 Zero-Panic 约定）。
- 面向短读组装：无长读（HiFi/ONT）路线（对照 metaMDBG）；纠错（mercy/read
  correction）仅做 k-mer 级别，不做 reads 级 polish。

---

*参考来源: 本项目源码 `megahit-1.2.9/`（src/megahit、main.cpp、main_assemble.cpp、
main_iterate.cpp、main_sdbg_build.cpp、sorting/{base_engine,kmer_counter,
read_to_sdbg_s1,s2,seq_to_sdbg}.{h,cpp}、sdbg/{sdbg.h,sdbg_raw_content.h,
sdbg_meta.cpp,sdbg_def.h}、assembly/{unitig_graph.{h,cpp},unitig_graph_vertex.h,
tip_remover.cpp,bubble_remover.{h,cpp},weak_link_remover.cpp,low_depth_remover.cpp,
sdbg_pruning.*}、localasm/{local_assemble.cpp,hash_mapper.*}、iterate/*.h、
sequence/kmer.h、definitions.h + README.md + CMakeLists.txt）*
