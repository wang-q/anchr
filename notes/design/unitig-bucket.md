# Unitig 计数分桶：已有基础设施调研与计划（设计）

> 2026-08-13。触发：`asm-assemble.md` §9.5 指出长 k 下 radix 排序成本随
> key_bytes 增长（small k=100 已 2.86 s），§9.3 第 4 项建议"minimizer
> 分桶"作为第二阶段。动手前先调研代码库——**分桶在 pgr/anchr 里已经有
> 成熟实现**，本文档盘点可复用面与真实缺口，再定计划。

## 1. 已有分桶基础设施（不需要重造）

### 1.1 `libs/fq/norm.rs`：k-mer 计数外部分桶（完整实现）

`fq norm` 的大内存输入走 `norm_buckets`（从 pgr 迁移，注释即引用
`pgr::kmer::count::count_keys` 的"memory-bounded bucket path"）：

* `bucket_of(key, buckets)`：**hash 分桶**（对打包字节做
  `h = h*131 + b`，再 `% buckets`），不是 minimizer；
* Pass A（`count_buckets`）：流式读 reads，`canonical_keys` 发射打包
  canonical key，按桶写 `bucket_{bi:05}.kmer` 文件（含重复）；
* Pass B：分**波**（wave）并行处理——每波取 `wave` 个桶，每个桶
  `count_keys`（radix 排序 + 分组）→ 序列化 `table_{bi:05}.tbl`；
  波大小由内存上限推导（`COUNT_BYTES_PER_KMER=21` +
  `TABLE_BYTES_PER_KMER=20`），`BUCKET_FRACTION=0.35`、`MAX_BUCKETS=4096`；
* 决策与输出与内存路径一致（逐桶打分 `score_in_chunks`，不需要全表
  同时在内存）。

### 1.2 `libs/fq/clump.rs`：pivot k-mer 外部分桶

`clump_buckets`（`SortMode::Bucket` / `--force-buckets`）：按 pivot k-mer
hash 分桶，桶内排序后按桶序输出；`MAX_BUCKETS` 与 norm 同源。同样是
hash 分桶 + 落盘 + 分波，无 minimizer。

## 2. 与 unitig 计数的关系（可复用面）

`TadpoleTable` 对齐 pgr 后（§9.5）的内存构成（G37 full，2.46 GB 峰值）：

| 组成部分 | 大小 |
| :--- | ---: |
| reads 双份拷贝（SeqRecord + `(seq, phred)`） | ~0.65 GB |
| 发射的中间打包 key（含重复，排序工作集） | ~1.0 GB |
| 最终 KmerTable（打包 key + counts，113 万条目） | ~13 MB |
| `sorted_entries` 惰性快照 `(Kmer, u32)` | ~54 MB |

中间 key 的 ~1 GB 是排序前的工作集。**分桶可把这一部分限到任意内存
上限内**：发射时按 `bucket_of` 落桶文件 → 分波 `count_keys` → 把各桶
的排序表 **k-way 合并**成完整 `KmerTable`。注意与 norm 不同：unitig 的
遍历（build_unitigs / scan_table）需要全表随机访问，不能像 norm 那样
逐桶打分后丢弃，所以合并是必须的（k-way merge 已排序桶表，代价低）。

可直接复用的部分：`bucket_count` 的内存推导、Pass A 落桶结构、
Pass B 分波并行、`bucket_of`、每个桶的 `count_keys`。需要新增的部分：
桶表 k-way 合并 + TadpoleTable 构建入口接入（质量门控发射保持）。

## 3. 真实缺口：minimizer 分桶（代码库里没有）

现有两个分桶都是 **hash 分桶**（目的只是内存有界化）。bcalm/FastK 的
**minimizer 分桶**是另一个维度：按 minimizer 聚类，桶内排序比较量随
minimizer 之后的部分变短——这正是长 k 下 radix 成本（§9.5：
k=31/64/100 = 1.89/2.42/2.86 s）的解法，**pgr/anchr 中不存在**，属于
新的算法工作。它可以叠加在 §2 的落盘/分波机制上（把 `bucket_of` 换成
minimizer 桶键），但 minimizer 的选择/窗口/桶键编码需要按 FastK/bcalm
思路新写。

### 3.1 FastK 实测与"MSD 排序"澄清（2026-08-14）

G37 full（`pe.cor.fa.gz`，k=31，同机）：

| 工具 | wall | RSS | 备注 |
| :--- | ---: | ---: | ---: |
| FastK `-t8` | **1.27 s** | **226 MB** | 只做计数 + profile，不做 unitig |
| anchr asm unitig `auto` | 2.18 s | 929 MB | 计数 + unitig |
| bcalm `-t8` | 2.38 s | 541 MB | unitig |
| cuttlefish `-t8` | 8.20 s | 1819 MB | unitig |

结论：

* **真正的差距在计数层，不在组装**——FastK 只花 1.27 s / 226 MB 完成
  计数，anchr 整套 2.18 s / 929 MB；
* **我们已经在使用 MSD radix 排序**（`pgr::libs::ds::radix_sort`，从
  Gene Myers FastGA `MSDsort.c` 移植、并行）；FastK 自己用的是 **LSD**
  （`LSDsort.c`）——排序方向不是差距；
* FastK 快的机制（`FastK.c`/`count.c` 注释）：
  1. **minimizer 分桶分发**（novel minimizer-based distribution），任意
     规模、磁盘友好；
  2. **两段式 "super-mer 然后 weighted k-mer" 排序**——低错误率
     （≤1%）数据下 super-mer 把待排序记录数压缩一个量级，这是关键；
  3. 位打包记录 + 多线程 LSD + 磁盘流式。
* 因此 §3 的"真实缺口"应明确为 **super-mer/minimizer 两段式计数**，
  而不仅是 minimizer 桶键。**pgr 已实现**（commit `769f82f`，
  `pgr::libs::kmer::supermer`，固定 m=12，输出与直接路径逐字节一致），
  anchr 已以 `asm unitig --supermer` 接入（`TadpoleTable::build_supermer`，
  无质量门控）；接入实测见 [asm-assemble.md](../design/asm-assemble.md)
  §12。

## 4. 计划（先定需求再实施）

1. **确认目标规模与内存上限**：当前 G37 full 中间 key ~1 GB，多大输入
   会撑爆？`fq norm` 已按 `--mem` 上限处理同类问题，unitig 应沿用同一
   种 cap 语义（或先只做内存路径 + 分桶兜底）；
2. **阶段 A（复用 norm 分桶）**：TadpoleTable 构建加内存有界路径——
   分桶落盘 + 分波 `count_keys` + k-way 合并，输出与现内存路径逐字节
   一致（golden 全绿为准）；
3. **阶段 B（super-mer/minimizer 两段式，长 k 治本，pgr 侧实现）**：
   在分桶机制上换 minimizer 桶键 + super-mer 预压缩，目标 k=100+ 的
   排序/比较量下降；先做小规模原型验证对 radix 成本的收益再全量接入；
4. 每一步先更新本文档/§9，再动代码。
