# cuttlefish (Cuttlefish 2)：紧凑 de Bruijn 图构建（源码分析）

> 2026-08-14 整理，纯源码分析（`cuttlefish-2.2.0/`）。Cuttlefish 是 Jamshed Khan &
> Rob Patro 等人的紧凑 de Bruijn 图（compacted de Bruijn graph, cdBG）构建工具：
> Cuttlefish 1 只接受**参考序列**（可输出 GFA），Cuttlefish 2 同时支持**测序 reads
> 与参考**（目前只输出 FASTA）。两篇论文：Cuttlefish (Bioinformatics 2021,
> btab309) 与 Cuttlefish 2 (Genome Biology 2022, s13059-022-02743-6)。
> 仓库内容构成：① `include/` 头文件 + `src/` 实现（`Application` 模板调度 →
> `CdBG`（1）/ `Read_CdBG`（2）→ 构造器/提取器/写入器）；② 内嵌第三方库
> （BBHash、compact_vector、KMC API、kseq、xxHash、spdlog、fmt、cxxopts、
> nlohmann/json、boost.preprocessor）；③ `data/` 示例数据；④ `scripts/` 合分支脚本。

## 1. 概况

- **定位**：从测序 reads（FASTQ）或参考序列（FASTA，可 gzip、多文件）构建**紧凑
  de Bruijn 图**——输出所有最大 unitig（最大非分支路径），并可选输出参考序列的
  unitig 拼贴（tiling）。以"极低内存 + 高并行 + 大规模可扩展"著称。
- **核心思想**（edge-centric）：把 dBG 的边当作 `(k+1)`-mer、顶点当作 `k`-mer
  来枚举。先由 **KMC3** 精确枚举/计数顶点与边（磁盘友好），再用 **BBHash 最小
  完美哈希（MPHF）** 把顶点集合映射到 `[0, n)`，在**紧凑位向量**里为每个顶点存
  一个 5/6 位的 **DFA 状态码**；随后并行扫描序列/边做**状态分类**（判定每个顶点
  是单入单出 / 多入单出 / 单入多出 / 多入多出），最后据此提取最大 unitig。
  全流程记忆体量极小：Cuttlefish 1 约 8.71 bits/vertex，Cuttlefish 2 约
  9.71 bits/vertex（`CdBG.hpp` 的 `bits_per_vertex`）。
- **双链语义**：顶点用 **canonical k-mer**（k-mer 与其反向互补视为同一顶点）；
  unitig 输出方向保证跨 run 稳定（canonical 形式）。
- **计数是精确的**（KMC3），不靠 Bloom filter；`cutoff` 只用于过滤低频
  (k+1)-mer（去测序错误）。
- **与 BCALM 2 的定位差异**：BCALM 2 用 minimizer 分桶 + 桶内压缩 + UF 全局拼接；
  Cuttlefish 不做 minimizer 分桶，改用"KMC 精确枚举 + MPHF + DFA 状态分类"，
  一次扫描即得每个顶点的度信息，unitig 提取时按状态即可确定起止。

## 2. 仓库结构

```
cuttlefish-2.2.0/
├── include/
│   ├── Application.hpp            # 顶层调度：按奇数 k 递归实例化 CdBG / Read_CdBG
│   ├── CdBG.hpp                   # ★ Cuttlefish 1 参考 dBG（构造 + 分类 + 输出）
│   ├── Read_CdBG.hpp              # ★ Cuttlefish 2 读/参考 dBG
│   ├── Read_CdBG_Constructor.hpp / Read_CdBG_Extractor.hpp   # 读模式构造/提取
│   ├── Kmer.hpp / Directed_Kmer.hpp / Annotated_Kmer.hpp / Kmer_Utility.hpp
│   ├── Directed_Vertex.hpp / Edge.hpp / Endpoint.hpp / Vertex.hpp
│   ├── Kmer_Hash_Table.hpp / Kmer_Hash_Entry_API.hpp / Kmer_Hasher.hpp
│   ├── State.hpp / State_Read_Space.hpp / Sparse_Lock.hpp / Spin_Lock.hpp
│   ├── kmer_Enumerator.hpp / kmer_Enumeration_Stats.hpp
│   ├── Kmer_Container.hpp / Kmer_Iterator.hpp / Kmer_SPMC_Iterator.hpp
│   ├── Thread_Pool.hpp / Task_Params.hpp / Job_Queue.hpp
│   ├── Maximal_Unitig_Scratch.hpp / Unitig_Scratch.hpp / Oriented_Unitig.hpp
│   ├── Unipaths_Meta_info.hpp / dBG_Info.hpp / dBG_Utilities.hpp / Data_Logistics.hpp
│   ├── Ref_Parser.hpp / Seq_Input.hpp / Validator.hpp / Build_Params.hpp
│   ├── globals.hpp                # MAX_K / INSTANCE_COUNT 宏与类型定义
│   ├── Input_Defaults.hpp / Output_Format.hpp / File_Extensions.hpp
│   ├── BBHash/ compact_vector/ kmc_api/ kseq/ xxHash/ spdlog/ fmt/ cxxopts/
│   │   nlohmann/ boost/preprocessor/        # 内嵌第三方库
│   └── DNA.hpp / DNA_Utility.hpp / utility.hpp / Progress_Tracker.hpp ...
├── src/
│   ├── main.cpp / commands.cpp     # CLI（cxxopts）
│   ├── Application.cpp             # Application<k,T_App>::execute 分派
│   ├── CdBG.cpp / CdBG_Builder.cpp / CdBG_Writer.cpp /
│   │   CdBG_GFA_Writer.cpp / CdBG_GFA_Reduced_Writer.cpp / CdBG_Plain_Writer.cpp
│   ├── Read_CdBG.cpp / Read_CdBG_Constructor.cpp / Read_CdBG_Extractor.cpp
│   ├── Kmer_Hash_Table.cpp / Kmer_Container.cpp / kmer_Enumerator.cpp ...
│   └── test.cpp                    # 开发者自测（kseq 读 FASTA、k-mer 去重检查等）
├── data/                           # 示例 reads.fq / refs1.fa / refs2.fa
├── scripts/                        # 合分支脚本（merge_into_master/develop）
├── CMakeLists.txt / cmake/         # INSTANCE_COUNT 控制最大 k
└── README.md
```

### 2.1 顶层调度（`Application.hpp` + `Application.cpp`）

`Application<k, T_App>` 是个**按奇数 k 递归**的模板类：每个实例持有一个
`app_next_level`（`Application<k-2>`）与一个 `T_App<k>* app`。构造时若传入的
`params.k() == k` 则实例化本层 app，否则把任务转给下一层 `k-2`；终止特化
`Application<1>`。`commands.cpp` 里按输入类型选 `T_App`：

```cpp
(params.is_read_graph() || params.is_ref_graph()) ?
    Application<MAX_K, Read_CdBG>(params).execute() :   // Cuttlefish 2
    Application<MAX_K, CdBG>(params).execute();          // Cuttlefish 1
```

`CdBG` 与 `Read_CdBG` 都通过 `ENUMERATE(INSTANCE_COUNT, INSTANTIATE, ...)`
宏（`globals.hpp`，基于 boost.preprocessor）对全部奇数 `k = 1, 3, ..., MAX_K`
做显式模板实例化。`INSTANCE_COUNT = (MAX_K+1)/2`，默认 32 → `MAX_K = 63`；
CMake `-DINSTANCE_COUNT=64` 可支持到 127。

## 3. 构建流水线总览

两套流水线都以 `construct()` 为入口（`CdBG.cpp` / `Read_CdBG.cpp`），开头都检查
`is_constructed()`——只要 `*.json` 元数据文件存在即视为已完成、直接返回（续跑
机制）。`dBG_Info` 收集"basic info / contigs info / short seqs / DCC info /
parameters info"五类信息，析构时 `dump_info()` 写 JSON。

### 3.1 Cuttlefish 1（`CdBG`，参考 dBG）

```
enumerate_vertices()   KMC3 数参考 FASTA 的顶点（k-mer），cutoff=1
construct_hash_table() 对顶点集构建 BBHash MPHF + 5 位紧凑位向量存状态
classify_vertices()    并行扫描参考序列，做 DFA 状态分类（§5）
output_maximal_unitigs() 按状态并行提取 unitig（plain / GFA / GFA-reduced）
```

- `enumerate_vertices`：`kmer_Enumerator<k>().enumerate(InputFileType::MULTILINE_FASTA, ...)`
  只数 k-mer 本身（`counter_max=1` 时 KMC 跳过计数，只去重）。
- `construct_hash_table`：`Kmer_Hash_Table<k, BITS_PER_REF_KMER>`（5 位/顶点），
  内存受限时用 `set_gamma` 调 BBHash 的 `gamma` 参数（§6.1）。
- 输出格式由 `-f` 决定（plain FASTA / GFA1 / GFA2 / GFA-reduced），三套输出
  写入器见 §7。

### 3.2 Cuttlefish 2（`Read_CdBG`，读/参考 dBG）

```
enumerate_edges()    KMC3 数 (k+1)-mer（FASTQ 或 FASTA），按 cutoff 过滤
enumerate_vertices() 以边库为输入（InputFileType::KMC）派生出 k-mer 顶点库
construct_hash_table() 对顶点集构建 BBHash MPHF + 6 位紧凑位向量存状态
compute_DFA_states()  Read_CdBG_Constructor：并行扫描边库做 DFA 状态计算（§5.2）
extract_maximal_unitigs() Read_CdBG_Extractor：并行扫描顶点库提取最大 unitig
```

- **边 → 顶点**：`enumerate_vertices` 用 `KMC::InputFileType::KMC` 把边库当输入，
  让 KMC 按新 k 重枚举，把 (k+1)-mer 边库投影成 k-mer 顶点库（`Read_CdBG.cpp:147-151`）。
- 计算完 DFA 状态后删边库，提取完 unitig 后删顶点库（除非 `--save-vertices`）。
- `--path-cover` 时提取"最大顶点不相交路径覆盖"而非最大 unitig（进度条与输出
  文案相应变化，`Read_CdBG_Extractor.cpp:39-40`）。

## 4. 核心数据结构

### 4.1 `Kmer<k>`（2-bit/碱基打包）

`include/Kmer.hpp`。`NUM_INTS = (k+31)/32` 个 64 位字，2-bit 编码，**后缀对齐**
（`kmer_data[0]` 存最末 32 个碱基）。定义了 `ODD_K` 宏声明"只用奇数 k"（保证
`Kmer<k>` 与 `Kmer<k+1>` 字长相同，顶点可从边直接 `from_prefix/from_suffix`
取前后缀，无需处理字边界）。关键操作：

- `reverse_complement()`：按**字节**查 `Kmer_Utility::REVERSE_COMPLEMENT_BYTE`
  表（4 碱基/字节），部分填充字节再位搬运，比逐碱基快。
- `roll_to_next_kmer(base, rev_compl)`：切掉首碱基、末尾补 `base`，同时增量更新
  reverse complement（`rev_compl.right_shift() + 高位补 complement(base)`），
  供序列滑窗时 O(1) 滚动。
- `roll_forward/roll_backward(edge)`：按 `Extended_Base`（边编码）向"右/左"滚动
  一个碱基——unitig 提取时沿边延展顶点用。
- `to_u64(seed)`：用 **xxHash3**（`XXH3_64bits_withSeed`）对打包数据哈希；配合
  `Kmer_Hasher`（沿用 BBHash 的种子 `0xAAAAAAAA55555555`）。
- `minimizer<l>()` / `count_lmers`：字典序 l-mer minimizer（2-bit 逐窗口比较），
  以及统计各 l-mer 频次——**注意**：Cuttlefish 2.2.0 里 minimizer 仅作 API 预留，
  **没有**被流水线使用（`Minimizer_Policy` 无引用，最小化分区未实现），与
  BCALM 2 的 minimizer 分桶不同。

### 4.2 `Directed_Kmer<k>` / `Directed_Vertex<k>` / `Edge<k>` / `Endpoint<k>`

- `Directed_Kmer<k>`（`include/Directed_Kmer.hpp`）：同时持有 `kmer_`、
  `rev_compl_`、`canonical_`（二者取小者）与 `dir_`（是否为正向）。滑窗
  `roll_to_next_kmer(next_base)` 一并对三者滚动。
- `Directed_Vertex<k>`（`include/Directed_Vertex.hpp`）：顶点实例。持有观测到的
  k-mer `kmer_`、其 RC `kmer_bar_`、指向 canonical 形式的指针 `kmer_hat_ptr`
  （指向 `kmer_` 或 `kmer_bar_` 之一），以及 canonical 的哈希 `h`。
  `in_canonical_form()` 即判断 `&kmer_ == kmer_hat_ptr`。`exit_side()/entrance_side()`
  依据"观测形式是否 canonical"决定该顶点作为边起点/终点时**入射边所在侧**
  （front/back）——这是双向边建模的关键。
- `Edge<k>`（`include/Edge.hpp`）：双向边实例 = `(k+1)-mer` 元组 `(u, s_u, v, s_v)`，
  持有 `e_`（边 (k+1)-mer）、`u_`/`v_`（`Endpoint`）。`configure()` 用 `e_`
  的前后缀配置两个端点；`is_loop()` 即 `u.canonical() == v.canonical()`。
- `Endpoint<k>`（`include/Endpoint.hpp`）：`{Directed_Vertex v, side_t s,
  edge_encoding_t e}`——端点顶点 + 边入射侧 + 边的 `Extended_Base` 编码。
  `neighbor_endpoint(e, hash)` 沿边滚动出邻居端点（`side_t::back` 时
  `roll_forward`，front 时 `roll_backward`）。

> `side_t`（front/back）与 `dir_t`（FWD/BWD）见 `globals.hpp`。`base_t` =
> `DNA::Base`（A/C/G/T/N，0-3 有效），`edge_encoding_t` = `DNA::Extended_Base`
> （E/A/C/G/T/N/OP_non_branch/OP_branching，见 `DNA.hpp`）。

### 4.3 `Kmer_Hash_Table<k, BITS_PER_KEY>`（MPHF + 紧凑位向量 + 稀疏锁）

`include/Kmer_Hash_Table.hpp`——Cuttlefish 内存效率的核心：

- **MPHF**：`boomphf::mphf<Kmer<k>, Kmer_Hasher<k>>`（BBHash）。`gamma` 参数控制
  bits/key 与构建速度的折衷（gamma 越大 MPHF 越大但构建/查询越快）；
  `bits_per_gamma[]` 表给出各 gamma 的经验 bits/key；`set_gamma(max_memory)` 在
  内存受限时自动选最大可用的 gamma。默认 gamma ∈ [2.0, 10.0]。
- **状态桶**：`compact::ts_vector<state_code_t, BITS_PER_KEY, uint64_t, ...>`
  （compact_vector 库的线程安全紧凑向量），以 `mph->lookup(kmer)` 为下标存每
  个顶点的状态码。`BITS_PER_REF_KMER=5`（参考）、`BITS_PER_READ_KMER=6`（读）。
- **并发**：`Sparse_Lock<Spin_Lock>`（65536 把自旋锁散列到桶下标）保护
  read-modify-write。`update(api)` 是 CAS 式更新：只有 `api.bv_entry` 仍等于读取
  时的旧值才写入新值，否则返回 false 让调用方重试（`CdBG_Builder` 里
  `while(!process_internal_kmer(...))` 重试）。`update_concurrent(api1, api2)`
  按桶下标排序加两把锁，避免死锁。
- **API 模式**：`operator[]`/`at(kmer)` 返回 `Kmer_Hash_Entry_API<BITS_PER_KEY>`
  ——包装"位向量条目 + 读取时的快照 + 可变副本"，读改写三段解耦。
- 支持 `save/load`（`--save-mph`/`--save-buckets`/`--save-vertices`）做续跑。

### 4.4 `State`（5 位，参考模式）与 `State_Read_Space`（6 位，读模式）

**`State`**（`include/State.hpp` + `src/State.cpp`）：一个 5 位状态码编码"未访问/
已访问 + 是否已输出 + 状态类 + （单入单出时）首/尾碱基"：

| code | 含义 |
|---|---|
| `0b00000` | 未访问 |
| `0b00011` | 多入多出（复杂节点，dead-end） |
| `0b001xx` | 多入单出（`xx` = back 碱基） |
| `0b010xx` | 单入多出（`xx` = front 碱基） |
| `0b01100-11` | 各状态类的**已输出**态（不存碱基） |
| `0b1xxxx` | 单入单出（`xx` = front/back 碱基） |

`0b00001/0b00010` 为非法码，任何构造/解码遇到都会 `exit(1)`。`decode()` 用
switch 把码还原成 `Vertex`（作者注释 TODO：换成查表可能更快）。
`is_dead_end()` = `is_visited() && state_class()==multi_in_multi_out`。
`outputted()` 把任意状态映射到对应"已输出"码。

**`State_Read_Space`**（`include/State_Read_Space.hpp`）：6 位 = front 3 位 +
back 3 位，各存一侧的 `Extended_Base` 边编码：`E`（无出边）、`A/C/G/T`（唯一
出边）、`N`（该侧分支/多出边）、`OP_non_branch`/`OP_branching`（已输出态）。
`mark_outputted()` 把两侧的 `N` 改成 `OP_branching`、其余改成 `OP_non_branch`。
读模式顶点无需存碱基——因为状态分类阶段已把"该侧是唯一边还是分支"编码进去了。

## 5. 状态分类（DFA 状态计算）

### 5.1 参考模式（`CdBG::classify_vertices`，`src/CdBG_Builder.cpp`）

`classify_vertices()` 用 `Ref_Parser`（kseq 封装）逐条读参考序列，按
`distribute_classification` 把 `[left_end, right_end]` 的 k-mer 区间切成
`thread_count` 份分发给 `Thread_Pool`（`Task_Type::classification`）。每条序列
的 k-mer 滑窗被拆成 **孤立 / 最左 / 最右 / 内部** 四种情形分别处理
（`process_isolated_kmer` / `process_leftmost_kmer` / `process_rightmost_kmer` /
`process_internal_kmer`）；对每个 k-mer 先求 canonical 与方向，再从哈希表取状态
做**状态转移**（CAS 失败则 `while` 重试）。核心转移逻辑（`process_internal_kmer`，
`CdBG_Builder.cpp:474-588`）：

1. 若 `state.is_dead_end()`（已判定多入多出）→ 直接返回（无需再转移）。
2. 若当前 k-mer 与下一 k-mer 的 canonical 相同（自环）→ `process_loop`。
3. 未访问（`!is_visited()`）→ 置为 `single_in_single_out(front=prev_base,
   back=next_base)`（按方向，BWD 时碱基取 complement）。
4. 已访问 → 解码成 `Vertex` 做**类间转移**：
   - `single_in_single_out`：`front==prev && back==next` → 不变；`front!=prev &&
     back!=next` → `multi_in_multi_out`；`front!=prev` → `multi_in_single_out`；
     否则（`back!=next`）→ `single_in_multi_out`。
   - `multi_in_single_out`：`back!=next` → 升格为 `multi_in_multi_out`。
   - `single_in_multi_out`：`front!=prev` → 升格为 `multi_in_multi_out`。
5. `state == old_state` 时跳过更新，否则 `hash_table->update(entry)`。

**自环处理**（`process_loop`，`CdBG_Builder.cpp:229-267`）：若环跨越顶点两侧
（直接重复导致的 crossing loop）或该 k-mer 是序列最左端，则直接判
`multi_in_multi_out`（复杂节点）；否则环被一侧完全包含（反向重复），该侧不能
再延展，等价于把 k-mer 当**最右 k-mer（哨兵）**处理——这正对应 bcalm 文档里
"单位化的顶点是单入单出；环侧阻塞延展"的语义。

**孤立 k-mer**（`process_isolated_kmer`）：序列中无任何邻居 → 直接判
`multi_in_multi_out`。

> 关键不变量：**一个顶点能作为 unitig 内部点，当且仅当它是单入单出**；多入多出
> 是 branch/junction（unitig 端点），单侧分支（多入单出 / 单入多出）是 unitig
> 的端点。这与 BCALM 的 unitig 定义（内部顶点 degree ≤ 2）等价，但 Cuttlefish
> 用"DFA 状态 + 单遍扫描"同时算出了方向和延展碱基，无需建显式图边。

### 5.2 读模式（`Read_CdBG_Constructor::compute_DFA_states`，
`src/Read_CdBG_Constructor.cpp`）

从边库（(k+1)-mer）出发，用 `Kmer_SPMC_Iterator`（单生产者多消费者迭代器）并行
扫每条边：取 `Directed_Vertex` 前缀/后缀作为两端点，对每端调
`add_incident_edge(endpoint)` 更新 `State_Read_Space`：

```cpp
cuttlefish::edge_encoding_t e_curr = state.edge_at(endpoint.side());
if(e_curr == N) return true;            // 该侧已标记为分支
cuttlefish::edge_encoding_t e_new = endpoint.edge();
if(e_curr != E) {                       // 该侧已有边
    if(e_new == e_curr) return true;    // 同一条边重复出现
    e_new = N;                          // 出现第二边 → 标记分支
}
state.update_edge_at(endpoint.side(), e_new);
return hash_table.update(bucket);       // CAS 更新
```

> 即：每侧首次出现的边记 `A/C/G/T`，第二次不同的边把该侧标成 `N`（分支）。
> 自环/回文等特例由 `Directed_Vertex::exit_side/entrance_side` 的侧判定自然
> 吸收。状态计算完成后边库即被删除。

## 6. k-mer 枚举与内存/磁盘管理

### 6.1 `kmer_Enumerator<k>`（`include/kmer_Enumerator.hpp`，KMC3 封装）

`enumerate(...)` 内部调用 KMC3（`kmc_runner.h`）两阶段：`Stage1` 近似统计
（`solid_kmer_count_approx` 预估 solid k-mer 数），`Stage2` 真正枚举去重。
常量：`min_memory=3`（GB，KMC3 要求）、`bin_count=2000`（临时文件数上限，README
建议 `ulimit -n 2048`）、`signature_len=11`、`counter_max=1`（`-cs 1` 跳过计数）。
`memory_limit(unique_kmer_count, bits_per_kmer)` 依据"预估唯一 k-mer 数 ×
bits/kmer"反推 KMC3 的内存预算，用于在 `--unrestrict-memory` 之外做严格内存控制。
`small_k_threshold=13` 是 KMC 的小 k 优化模式阈值。

### 6.2 内存/磁盘折衷

- `-m/--max-memory` 是**软上限**：Cuttlefish 会在"至少满足最小需求"的前提下
  尽量遵守；`--unrestrict-memory` 则放开（`strict_memory=false` 时给
  `Kmer_Hash_Table` 传 `gamma = DBL_MAX` 使 MPHF 最快、占内存最大）。
- `construct_hash_table` 里 `max_memory = max(process_peak_memory(),
  params.max_memory()) - parser_memory`（参考模式固定扣 256 MB 序列解析器，
  读模式扣 `Kmer_SPMC_Iterator::memory(thread_count)`）。
- 磁盘用量按阶段估计打印：`max_disk_usage(edge_stats, vertex_stats)`（读模式把
  边枚举期与顶点枚举期的临时文件 + 库文件相加取大，`Read_CdBG.cpp:209-216`）。

## 7. 输出

### 7.1 Cuttlefish 2（FASTA，`.fa`）

只支持 FASTA：`Read_CdBG_Extractor` 对每个**未输出**顶点从两侧延展拼出最大
unitig（`Maximal_Unitig_Scratch<k>`：把 unitig 在相遇顶点处切成 `u_b`/`u_f` 两段，
`\bar(u_f) glue_k u_b` 拼回），`finalize()` 时按端点哈希决定 canonical 方向并
赋唯一 ID；DCC（Detached Chordless Cycle，环状 unitig）用 `rotate_append_cycle`
旋转到最小顶点输出。统计由 `Unipaths_Meta_info` 汇总（条数/总长/最长最短/
DCC 数）。输出经 `Character_Buffer` 缓冲后写 `Output_Sink`。

### 7.2 Cuttlefish 1（GFA1 / GFA2 / GFA-reduced）

参考模式除 FASTA 外还可输出图与参考序列的 unitig 拼贴（`CdBG` 的 writer 家族，
`src/CdBG_Writer.cpp` 等，异步 spdlog 逐线程缓冲）：

- **GFA1（`.gfa1`）**：`S` Segment（unitig）+ `L` Link（连接）+ `P` Path（每条
  参考序列的 unitig 拼贴）。
- **GFA2（`.gfa2`）**：Segment + `E` Edge + `G` Gap + `O` Ordered Group。
- **GFA-reduced（`.cf_seg` + `.cf_seq`）**：`<id> <segment>` 段表 + `<id> <tiling>`
  拼贴表（unitig id 带 `+`/`-` 方向）。README 称 7 人基因组 k=31 时 GFA2 占
  112 GB、reduced 只占 29.3 GB。
- 跨线程的连接由 `first_unitig/second_unitig/last_unitig` 追踪后
  `write_inter_thread_connections` 补齐；路径拼贴由生产者线程把本序列各线程
  片段写临时文件、`write_sequence_tiling`（消费 `Job_Queue`）顺序归并。

**颜色语义**（GFA 的 `P`/`O`）：unitig 的"颜色集"= 包含它的输入参考集合，编码
在各序列的路径条目里；因此 unitig 天然"单色"（颜色集在 unitig 内部变化只可能
是分支/序列端点，都会截断 unitig）——README 明确这是 Cuttlefish 特定的颜色定义。

## 8. 参数语义

| 参数 | 含义 |
|---|---|
| `-s/--seq`, `-l/--list`, `-d/--dir` | 输入文件 / 文件列表 / 目录（可混用） |
| `-k/--kmer-len` | k-mer 长度（默认 27，**必须为奇数**，≤ MAX_K） |
| `-t/--threads` | 线程数（默认 = 硬件并发数 / 4） |
| `-o/--output` | 输出前缀（FASTA 追加 `.fa`，另产出 `.json`） |
| `-w/--work-dir` | 工作目录（**须已存在**，默认 `.`） |
| `-m/--max-memory` | 软内存上限 GB（默认 3，至少满足最小需求才遵守） |
| `--unrestrict-memory` | 放开内存限制（更快，MPHF gamma→max） |
| `-f/--format` | Cuttlefish 1：0 FASTA / 1 GFA1 / 2 GFA2 / 3 GFA-reduced |
| `--track-short-seqs` | 记录长度 < k 的序列到 JSON（Cuttlefish 1） |
| `--poly-N-stretch` | GFA / GFA-reduced 的拼贴（tiling）输出中记录 polyN 段（Cuttlefish 1） |
| `--read` / `--ref` | Cuttlefish 2 输入类型（FASTQ / FASTA）；两者都传 → `is_valid` 校验失败退出，都不传 → Cuttlefish 1 |
| `-c/--cutoff` | (k+1)-mer 频次阈值（默认 reads 2、refs 1） |
| `--path-cover` | 提取最大顶点不相交路径覆盖而非最大 unitig |
| `--save-mph` / `--save-buckets` / `--save-vertices` | 保存 MPHF / DFA 状态桶 / 顶点库（续跑） |
| `--vertex-set` / `--edge-set` | 调试：直接喂入 KMC 库前缀 |

校验（`Build_Params::is_valid`，`src/Build_Params.cpp`）：输入非空；k 须奇数且
≤ MAX_K；线程数 ≤ `hardware_concurrency`；输出目录与工作目录须存在。失败打印
原因并退出。

**更大 k**：源码安装默认 MAX_K=63；`cmake -DINSTANCE_COUNT=64` 支持到 127，
最多 255。只用 `(k+31)/32` 个字，故大 MAX_K 不影响小 k 的性能
（`Kmer<k>` 是编译期模板，只实例化用到的奇数 k）。

## 9. 与 anchr 的关联

- **anchr 现状**：`anchr asm unitig`（`src/cmd/asm/unitig.rs` +
  `src/libs/asm/assemble.rs::assemble_unitigs`）是 **BCALM 2 的移植**（`ograph.cpp`
  `graph3` 语义），不是 cuttlefish 的移植；cuttlefish 与 bcalm 是同一目标
  （cdBG 构建）的两条不同技术路线，见 `notes/design/asm-assemble.md` §10。
- **可借鉴点**：
  1. **边中心（edge-centric）建图**：以 `(k+1)-mer` 为边、`k-mer` 为顶点枚举，
     顶点可从边直接取前后缀（`Kmer<k>::from_prefix/from_suffix`，`ODD_K` 保证
     字长一致）——若 anchr 未来做"从精确计数直接求图"，可省一次独立建顶点集的
     扫描。
  2. **MPHF + 紧凑位向量存顶点状态**（§4.3）：BBHash 把顶点映射到 `[0,n)`，
     用 5/6 位/顶点存 DFA 状态，而非存邻接表/度数——这是 Cuttlefish 内存占用
     远小于 bcalm 的关键。anchr 的 `KmerTable` 若做内存内图，可用 MPHF 索引压缩
     顶点元数据。
  3. **状态分类代替显式建边**（§5）：单遍扫描即可确定每顶点的单入/多入/单出/
     多出，unitig 起止据此判定，无需真正物化边表。状态转移是纯函数式
     （读-改-CAS-写），天然可并行、可重试——比 bcalm 的桶内 `graph3` 双指针
     合并更简单直接，适合 anchr 的单机并行场景。
  4. **DFA 状态即"unitig 端点判据"**：单入单出 = 可穿越；多入多出 = 端点；
     单侧分支 = 端点。与 bcalm 形式化文档（`bidirected-graphs-in-bcalm2.md`）
     的 unitig 定义（内部顶点 degree ≤ 2）互相印证，可作为 `anchr asm unitig`
     验证逻辑的参考。
  5. **续跑机制**：以 `*.json` 元数据文件存在与否判定"是否已构建"，配合
     `--save-mph/--save-buckets/--save-vertices` 分阶段落盘复用。
  6. **路径覆盖变体**（`--path-cover`）：最大顶点不相交路径覆盖可作 unitig 的
     对照/替代输出。
  7. **验证器（`cuttlefish validate`，Cuttlefish 1/参考图专用）**：`Validator` 用另一份 BBHash 检查
     "unitig 的 k-mer 集合 == KMC 库 k-mer 集合"（概率性，MPHF 可能把额外 k-mer
     映射到合法哈希值）+ 逐参考行走验证完整覆盖——anchr 的 `asm unitig` 与
     bcalm 做逐字节对照时可借鉴这套"集合相等 + 序列覆盖"双层校验。
- **与 bcalm 的路线差异**（对照 `notes/references/bcalm.md`）：bcalm 走
  "DSK 计数 → minimizer 分桶 → 桶内 graph3 压缩 → UF 全局拼接 → 重算 L: 边"；
  cuttlefish 走 "KMC 精确枚举 → MPHF → DFA 状态分类 → 按状态提取"。前者把
  压缩并行化在"分桶"上、需要跨桶边落盘；后者把并行化在"每顶点状态 CAS 更新"
  上、无分桶无跨桶边。对 anchr 的单机内存场景，cuttlefish 的路线（省去外部
  排序与分桶落盘）比 bcalm 更适合参考；`anchr asm unitig` 现有实现对应 bcalm
  语义，若未来要加"极低内存大图"能力，可参考 cuttlefish 的 MPHF+状态码方案。

## 10. 局限

- **依赖 KMC3**：k-mer 枚举完全外包给 KMC（`kmc_runner.h`），编译与运行都依赖
  其二进制；KMC 阶段最少约 3 GB 内存、最多开 2000 个临时文件（需
  `ulimit -n 2048`，否则报 `Cannot open temporary file ./kmc_00000.bin`）。
- **k 必须为奇数且 ≤ MAX_K**（编译期 `INSTANCE_COUNT` 决定，需 cmake 重编）：
  偶数 k 直接拒绝（"even k not consistent with the theory"）。
- **Cuttlefish 2 尚无 GFA 输出**（README 列为 roadmap），只有 FASTA；GFA 输出
  是 Cuttlefish 1 专属（`-f`）。
- **内存模型依赖经验估计**：`bits_per_vertex`（8.71/9.71）与 `parser_memory`
  （256 MB）是经验常数；`--unrestrict-memory` 时 MPHF gamma 取 `DBL_MAX` 的
  行为依赖 BBHash 实现细节。
- 工作目录须预先存在；输出目录不存在即报错退出。
- `State` 的编解码用大量 switch-case（作者 TODO：换查表），对性能敏感路径是
  潜在优化点；`Kmer_Hash_Entry_API` 读取-修改-写回依赖 CAS 成功才更新，高
  竞争时会重试。
- 代码含 `#ifdef CF_DEVELOP_MODE` 调试分支（如直接喂 KMC 库、gamma 参数、
  `test.cpp` 自测），正式构建不启用。
- unitig 输出无丰度信息（FASTA 头仅 unitig id），读模式的 `cutoff` 过滤发生在
  KMC 枚举阶段，改阈值需重跑枚举（与 bcalm 把过滤推迟到消费侧不同）。

---

*参考来源: 本项目源码 `cuttlefish-2.2.0/`（include/Application.hpp、CdBG.hpp、
Read_CdBG.hpp、Kmer.hpp、Kmer_Hash_Table.hpp、Kmer_Hash_Entry_API.hpp、State.hpp、
State_Read_Space.hpp、Directed_Vertex.hpp、Edge.hpp、Endpoint.hpp、Thread_Pool.hpp、
Maximal_Unitig_Scratch.hpp、Unipaths_Meta_info.hpp、dBG_Info.hpp、kmer_Enumerator.hpp、
Build_Params.hpp、Input_Defaults.hpp、Validator.hpp + src/commands.cpp、Application.cpp、
CdBG.cpp、CdBG_Builder.cpp、Read_CdBG.cpp、Read_CdBG_Constructor.cpp、
Read_CdBG_Extractor.cpp、State.cpp、Build_Params.cpp、State_Read_Space.cpp +
README.md + CMakeLists.txt）*
