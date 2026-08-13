# `anchr asm contig`：tadpole contigMode 迁移（设计）

> 2026-08-11。目标：替代 anchr 中 tadpole 的**组装用途**（contigMode）：
> 2_insert_size 流程（硬依赖）与 unitigs 流程（`--unitigger tadpole` 可选
> 分支）。ecc/extend 已由 `fq ec-kmer`/`fq extend` 覆盖（见
> [fq-merge-replace.md](fq-merge-replace.md) §6-7）。
> 参考：BBTools-40.01 `assemble/Tadpole*.java`。

## 1. anchr 中的调用点（tadpole 组装）

| 流程 | 调用 | 参数 |
|---|---|---|
| `2_insert_size.tera.sh` | `tadpole.sh in=R1 in2=R2 out=PREFIX.tadpole.contig.fasta threads=N overwrite [prefilter]` | 默认 k=31 |
| `unitigs.tera.sh` | `tadpole.sh in=pe.cor.fa out=unitigs_K{k}.fasta threads=N k={k} overwrite` | k ∈ opt.kmer（如 "31 81"） |

两者都不传 `mode=` → **默认 contigMode**（`Tadpole.java`：无 ecc/extend/
toss 标志时 `processingMode=contigMode`）。下游：2_insert_size 接 bbmap +
reformat-ihist（另计）；unitigs 接 `anchr contained/orient/merge`（anchr
自有，不迁移）。**prefilter 默认关**（anchr `opt.prefilter=0`），本命令
不实现（同 ecc/extend）。

## 2. contigMode 语义（源码确认，逐条移植）

### 2.1 建表与输入

- 与 ecc/extend 相同的 canonical kmer 计数（`TadpoleTable`），minprob
  质量门控（FASTA 无质量 → 不过滤）。
- `kmerRangeMin=0`/`kmerRangeMax=MAX` 默认不做 kmer 范围过滤；
  `removeBubbles=false`/`removeDeadEnds=false`（shave/rinse 默认关，
  跳过）；`processContigs=false`（不建图/不 pop bubble）。

### 2.2 多轮种子（contigPasses=16，contigPassMult=1.7）

`BuildThread.run`（Tadpole2，k>31）：

```
for i = 15 .. 1:
    minCountSeedCurrent = max(3+i, floor(3 * 1.7^i * 0.92 - 0.25))
    扫描全表：count >= 阈值 且 未被 claim 的 kmer → processKmer
最终轮：minCountSeedCurrent = 3，再扫一遍
```

Tadpole1（k≤31）同构（`Tadpole1.BuildThread`）。**单线程下就是 16 轮
全表扫描**，每轮对未认领的种子 kmer 建 contig。

### 2.3 认领（ownership，contigMode 下 useOwnership=true）

- 每个 kmer 有 owner（-1 未认领 / 0..N-1 线程 id）。单线程 id=0：
  认领集合 = HashSet<Kmer>。
- `processCell`：count < 阈值跳过；已认领跳过；认领后 `processKmer`。
- 行走中每个新 kmer 也要认领；`owner==id`（本线程已认领）→ 环形检测，
  返回 `fbranch ? F_BRANCH : LOOP`；被其它线程认领 → BAD_OWNER（单线程
  不会发生）。
- `leftCounts` 在 BuildThread 里**非空**（区别于 extend 模式），因此
  行走时启用左 junction 与隐藏分支检查（`leftMaxPos != evicted` → 停）。

### 2.4 行走（extendToRight，contigMode 版）

入口：count(minCountSeed)=3、owner 检查 → 左/右计数（4 桶）→
`rightMax < minCountExtend(2)` → DEAD_END；`isJunction(rightMax,
rightSecond)` → `isJunction(leftMax,leftSecond) ? D_BRANCH : F_BRANCH`；
`isJunction(leftMax,leftSecond)` → B_BRANCH。循环：取 rightMaxPos 碱基
追加 → 新 kmer 认领/环形检测 → 若 `bbranch` 返回 `fbranch?D_BRANCH:
B_BRANCH`；`hbranch`（leftMaxPos!=evicted 且 branchMult1>0）同；追加后
`fbranch` → F_BRANCH；`rightMax<2` → DEAD_END。

isJunction（深度比阈值，与 ecc/extend 共用）：`second<1 || second*20<
max || (second<=3 && max>=max(2, second*3))` → false（非 junction）。

makeContig（Tadpole2）：种子 kmer → extendToRight → `reverseComplement`
再 extendToRight（即向左延伸）→ `doubleClaim`（单线程恒真）→
`trimEnds(0)` → 长度 `>= k+minExtension(2)` 且 `>= minContigLen` →
生成 Contig（`leftCode/rightCode/leftRatio/rightRatio`，canonical 方向）。
Tadpole1（k≤31）路径同构，junction 方向判定按
`kmer>rkmer`（Tadpole1）/`kmer<rkmer`（Tadpole2）翻转——与现有
`extend_to_right2` 的 `canonical_is_rc = k > 31` 一致。

### 2.5 覆盖度与输出

- `calcCoverage`：contig 每个 kmer（canonical count）均值
  `sum/(float)kmers`，并记 `minCov`/`maxCov`；`coverage<minCoverage(1)`
  或 `>maxCoverage` → 丢弃该 contig（默认不丢）。
- `minContigLen` 默认 `max(124, 2k)`（k=31→124，k=81→162）；`minExtension=2`。
- 输出 FASTA：`>contig_{id},len={len},cov={cov:.1},gc={gc:.3},min={min},
  max={max},hh={hh:.3},caga={caga:.3},left={code},right={code}` + 序列。
  gc/hh/caga 由 `calcScalarsFast` 计算（实现时逐行对照）。
- 排序：length 降序 → coverage 降序 → 序列字典序 → id。
- **顺序无关性**：行走受 leftCounts 隐藏分支约束，contig 集合由
  （图 + 阈值 + 认领）唯一确定，与种子扫描顺序无关 → anchr 的
  HashMap 迭代顺序不影响输出（黑盒验证确认后定案；排序含序列字典序
  兜底，id 只在完全重复时生效）。

## 3. 命令形状

```
anchr asm contig [OPTIONS] <infiles>...
  -k, --kmer <int>        默认 31（2_insert_size 用；unitigs 由模板按 k 循环调用）
  -o, --outfile <file>    输出 FASTA（默认 stdout）
  -p, --parallel <int|auto>  兼容参数，校验但不启用（确定性单线程，同 ecc/extend）
  --min-contig-len <int>  默认 auto = max(124, 2k)
```

- 单 k 每次调用；unitigs 模板按 k 循环（与现模板逐 k 调用 tadpole 一致）。
- 单线程确定性 = 与 `tadpole.sh threads=1` 黑盒对照的前提。

## 4. 验证

- Lambda `tests/bbtools/Lambda/pe.cor.fa.gz`（或 ecco 输出）k=31/81，
  `tadpole.sh threads=1` 生成 golden（contig FASTA），逐字节对照。
- 关键正确性锚点：contig 集合 + 序列 + 头字段（len/cov/gc/min/max/hh/
  caga/left/right）+ 排序。
- 已知风险：HashMap 迭代顺序与 Java 哈希表不同——若黑盒出现顺序相关
  差异，需复刻 Java 表布局或改确定性扫描顺序（先以实验定论）。

## 5. 不做

- prefilter（anchr 默认 0）；shave/rinse/pop/bubble；`mode=insert`；
  多线程；`trimCircular`/`trimEnds`（默认 0）；kmer 范围过滤。
- 2_insert_size 的 bbmap + reformat-ihist 仍属独立缺口（todo 挂账）。

## 6. 实现状态与已知偏差（2026-08-11 定案）

`anchr asm contig` 已实现：contig 构建（多轮种子/行走/认领）+ contig 图 +
BubblePopper + 排序重编号 + 输出，全部确定性与单线程等价。

**气泡开关（2026-08-11 定案）**：默认 `pop_bubbles=true`（tadpole
`popbubbles=t` 兼容，anchr 行为不变），新增 `--no-bubbles` 逃生门
（等价 tadpole `popbubbles=f`：跳过 process_contigs/pop_bubbles，直接对
pre-pop contig 集排序重编号输出）。Lambda 实测：pre-pop ≥124bp 77 条
（mincontig=1 全量 89 条），pop 后 66 条——合并路径让部分 <124bp 的
中间 contig 并入长 contig（输出总碱基反而 +105）。理由：气泡合并的
代表路径选择是任意启发式（依赖图布局），`--no-bubbles` 保留泡的两条
分支，语义更接近 cdBG unitig（见 `notes/references/bcalm.md` §4.3）。

### 已验证（逐字节）

- **pre-pop contig 集合 89/89 与 `tadpole.sh threads=1` 逐字节一致**
  （含短 contig，`popbubbles=f mincontiglen=1` 对照）。
- k=31 左端边行走在 rc 空间（Tadpole1 `processContigLeft` 交换
  kmer/rkmer 语义）——修复了"行走绕回自身生成自环边"的 bug。

### 已知偏差（bubble 解析顺序）

tadpole 的气泡消除**顺序相关**：expand 顺序决定重叠气泡中"谁吸收谁"，
进而决定链式合并。实验证据：

- 把 Java `popBubbles` 的 expand 迭代倒序，输出从 67 变 66 contig；
- Java 的哈希表 cell 顺序随 `-Xmx` 变化（Xmx1g/3g/8g → prime
  228983/213973/194057），但构建顺序 68/89、71/89 个 id 不同——即
  "逐字节一致"本身只对特定内存参数成立；
- anchr 用确定性扫描顺序（canonical kmer 排序）代替 Java 的哈希 cell
  顺序，输出确定且跨运行稳定，但 bubble 解析结果与 tadpole 有少量
  差异（Lambda 2000 对：67 vs 66 contig，总碱基差 ≤100，序列集合
  重合 ≥90%）。

**决策（用户确认，2026-08-11）**：不做哈希表布局复刻，接受确定性输出
+ 文档化偏差。理由：逐字节一致需复刻 BBTools 内存模型（`-Xmx` 相关
prime）+ 开放寻址插入顺序 + 溢出树，约几百行无生物学价值的"镜像
Java 内存布局"代码，且结果脆弱（换 `-Xmx` 即失效）。anchr 的 contig
集合/总碱基与 tadpole 一致，差异只是少数气泡的"走哪条路径"（两条
都是合法组装选择，序列质量等价），对 anchr 用途（insert-size 参考、
unitigs 组装）影响可忽略。

回归测试：`tests/cli_asm_contig.rs` + golden
`tests/bbtools/Lambda/golden/tadpole_contigs31.fasta.gz`（tadpole
默认输出 67 contig），断言确定性、总碱基差 ≤100、序列集合重合 ≥90%。

## 7. 性能优化（2026-08-11）

计数表与组装扫描按 `pgr::libs::kmer`（FastK/Myers 计数骨架）模式改造，输出
逐字节不变（golden 全绿）：

- **`TadpoleTable::sorted_entries`**：canonical k-mer 排序快照用
  `OnceLock` 缓存一次；`scan_table` 16 轮种子扫描改为线性迭代，去掉
  每轮 O(n log n) 的 collect+sort（原 HashMap 迭代需每轮排序保证确定性）。
- **并行构建**：`TadpoleTable::build` 按 4096 reads 分块 rayon 并行计数
  + 确定性合并；表内容与单线程一致（`pgr::libs::kmer::build_table` 同款模式）。
- **基准**（`benches/fq_assemble_benchmark.rs`，Lambda 20k reads，k=31，
  release）：assemble 全流程 576 ms →（sorted_entries）313 ms →
  （+并行 build）157 ms，~3.7×；build 247 ms → ~100 ms。
- **radix 化评估（2026-08-11，实测不做）**：曾尝试把 sorted_entries
  的排序换成 Myers radix（k≤64 用 u128 投影 + `radix_sort_u128`），
  Lambda 20k 实测反而更慢（比较排序 157 ms → radix 193 ms；用
  `mem::take` 消除占位分配后 164 ms，仍略慢）。结论：几十万唯一
  k-mer 规模下 `cmp_bases` 比较排序更优（缓存局部好、无投影/索引
  构建开销）；radix 的价值需数百万级 k-mer 才可能显现，届时再评估。
  k>64 的多 word radix 泛化同步搁置。

## 8. `anchr asm unitig` 命令（2026-08-11，借鉴 BCALM graph3）

新增独立命令 `anchr asm unitig`（**不从 assemble 加开关**）：不做种子
扩展/气泡，改为**最大 unitig 压缩**（`ograph.cpp` `graph3` 语义）。拆分
原因：`--no-bubbles`（tadpole 兼容参数）与 unitig 压缩语义不同但名字
相似，放在同一命令下造成困惑；独立命令让每个命令只有一种组装哲学
（assemble = tadpole 兼容 contig，unitigs = 严格图压缩）：

- solid 定义 = count ≥ `min_count_seed`（默认 3）；每个 solid k-mer 沿
  "唯一后继（out==1）且下一 k-mer 唯一前驱（in==1）" 双向延伸，分支/汇合/
  死端/环（`visited` 检测）处断开。
- **顺序无关**：unitig 集合由（k-mer 集 + solid 阈值）唯一确定，无认领/
  种子顺序依赖（对比 contig 模式的确定性排序只是"复刻"扫描顺序）。
- **无气泡**：平行路径各自成 unitig（测试
  `command_asm_unitig_keeps_branches` 验证 ≥4 条、不横跨）。
- 输出头沿用 contig 字段（len/cov/gc/min/max/hh/caga，无 left/right
  分支码），`>unitig_<id>`。
- 基准（Lambda 20k，release）：160 ms，与 contig 模式持平（计数占大头）；
  价值在语义不在速度。
- **适用**：高覆盖/已纠错输入（anchr unitigs 的 `pe.cor.fa`）；低覆盖
  原始 reads 上 unitig 会比 contig 碎（实测 2k Lambda：110 条/44823 bp
  vs contigs 77 条/48059 bp）。
- **环状处理（简化）**：纯环 k-mer 链由 `visited` 检测断开，输出近似环
  的路径而非闭合环（bcalm 用 expect_circular 兜底，未移植）。
- **待验证**：真实 `pe.cor.fa` 上 unitigs 与 bcalm 输出的对照（todo §5）。

### 8.1 `--links` / `--gfa`：unitig 间边输出（2026-08-11）

对齐 bcalm LinkTigs 语义：两条 unitig 共享端点 (k-1)-mer 即相连。
`--links` 在 FASTA 头追加 `L:<from±>:<to>:<to±>`（bcalm 格式）；
`--gfa` 输出 `H`/`S`/`L` 行（overlap `(k-1)M`）。方向规则（实际序列
匹配，简化自 LinkTigs 的 `beginInSameOrientation` 判定）：

- 源右端 `r` == 目标左端 `a` → `+`/`+`（3'→5' 正向出边）；
- `r` == rc(`a`)（目标左端）→ `+`/`-`；
- 3'-3' / 5'-5' 相遇 → 反链表示（`-`/`-` 或 `-`/`+`）。

边集合由（unitig 端点 (k-1)-mer + 阈值）唯一确定，输出排序去重后
确定性。单测 `links_directions_branch_and_rc` 锚定三种方向组合；
与 bcalm 真实输出对照（2026-08-12，MG1655 1M 纠错 reads，k=31）：

* **unitig 序列 100% 一致**（2403/2403，canonical 方向归一后逐条相同，
  聚合统计完全相同：总数/总长/N50/最长）——`asm unitig` 本体验证通过；
* **L: 边集对齐（2026-08-14，已解决）**：早前的 2026-08-12 对照（无向边
  3801 vs 3331、共同 2577）用的是简化语义与错误的对照口径。仓库里的
  bcalm v2.2.3 源码（`gatb-core/.../LinkTigs.cpp`）含完整 LinkTigs 实现，
  已按它精确重写 `compute_links`：**双侧 in/out（`L:-:` 从 begin、
  `L:+:` 从 end）各 4 种方向情形 + 偶数 `(k-1)` 回文特判 + 自连边**
  （bcalm 对首尾共享同一 (k-1)-mer 的 unitig 连到自己，如 poly-C）。
  G37 full 验证：**边集 0/1482 mismatch**；把 bcalm 的链向按各自 unitig
  存储方向（bcalm 为遍历方向、anchr 为 canonical 方向）规范化后，
  **链向编码也 0/1482 一致**；总 `L:` 条目 3974 = 3974。字节级差异仅剩
  unitig 编号/排序（bcalm 遍历序 vs anchr 长度序）与存储方向约定，
  图结构（unitigs + 边 + 链向）完全等价。`--links`/`--gfa` 现与 bcalm
  LinkTigs 语义一致。

### 8.2 GFA 版本支持分析（2026-08-14）

**结论：不需要主动支持其他版本，维持 GFA 1.0。** 当前输出恰好落在 1.0
的**核心子集**：`H`（`VN:Z:1.0` + `ks`）、`S`（id + seq）、`L`（overlap
`(k-1)M`），无 P/C/W/J 行。S 行带 `LN:i:`/`KC:i:`/`km:f:` 可选标签
（2026-08-14 起，与 bcalm `convertToGFA.py` 一致，见 §8.3）。

各版本差异与影响：

* **1.0 → 1.1**：只新增 `W` 行（walk，pangenome 用途）；S/L 行语法完全
  相同，可选字段（`LN`/`RC`/`FC`/`KC`/`SH`/`UR` 等）1.0 就存在，不需要
  升版本即可使用。因此声明 `VN:Z:1.1` 是**零成本切换**（S/L 一行不改），
  但目前没有下游工具要求它，切换无收益。
* **1.2**：新增 `J` 行（gap/距离连接）、P 行 `;` 分隔跳转与 `SC:i:1`
  快捷标签。anchr 没有 gap/scaffold 概念，`asm unitig` 是纯重叠图；
  只有未来 `asm cns/layout` 需要表达 gap 才需要。
* **2.0**：为 pangenome 语义重写的版本（W 行体系、行结构不同），unitig
  压实图不是它的目标场景；主流工具对 2.0 的适配也明显少于 1.x，不需要
  支持。

生态与对齐依据：

* GFA-spec 的 README 确认 Bandage、gfatools、vg、GraphAligner 等主流工具
  对 1.x 兼容；1.1/1.2 是向后兼容的扩展版 1，不会破坏 1.0 消费者。
* bcalm 的 `convertToGFA.py` 同样**写死 `VN:Z:1.0`**，并把 FASTA 头里的
  `LN:i:`/`KC:i:`/`km:f:` 原样带进 S 行。我们的 `--gfa` 现在输出相同的
  S 行标签（§8.3），保持与它的字段对齐。

触发升级的判据：未来引入 gap/scaffold 表达（需要 `J` 行）或明确要面向
pangenome 工作流时，再按下游工具实测选 1.1/1.2；在那之前保持 1.0。

### 8.3 bcalm 功能迁移清单（2026-08-14）

对照 bcalm 全部公开选项与 `scripts/` 的结论，需要迁移的三块已完成，
其余由 pgr/既有设计覆盖或暂不需要：

* **输入模式（已迁移）**：bcalm 的 `-in` 接受任意 FASTA/FASTQ（可 gz）
  且"不关心 paired"。`read_records` 改为**任意多文件顺序读取、不要求
  配对**，`asm unitig`/`asm contig`/`asm olc` 放开为 1 个以上文件，并加
  `--list-files`（每行一个序列文件路径，pgr `resolve_paths` 语义）。
  单文件奇数记录（如 bcalm `circular_unitigs_unittests/test1.fa`）不再
  报 `unpaired trailing read`。
* **环状 unitig（回归锚定）**：bcalm 有专门回归例
  `circular_unitigs_unittests`。实测 test1（单条 16 nt 环）输出
  `AAGTCCGCTAAGTCC` 与 bcalm **逐字节一致**；test3（环 + 一条短随机
  读）bcalm 输出 `GTCCGCTAAGTCCGC`、anchr 输出 `AAGTCCGCTAAGTCC`，
  长度相同、canonical k-mer 集相同（图等价），但旋转/切口不同——bcalm
  本身也不保证 orientation/切口稳定，测试以图等价为准。
* **丰度输出（已迁移）**：`--all-abundance-counts` 在 FASTA 头追加
  `ab:Z:<c1> <c2> ...`（序列顺序的逐 k-mer canonical 计数，bcalm
  `-all-abundance-counts` 格式）；`--gfa` 的 S 行输出
  `LN:i:`/`KC:i:`/`km:f:`（bcalm `convertToGFA.py` 字段）。
* **不需要搬**：`-histo`/`-histo2D`（pgr `kmer hist`/`plot spectra` 已
  覆盖）；`-nb-cores`/`-max-memory`/磁盘分桶（外部分桶见
  [unitig-bucket.md](unitig-bucket.md)，规模到了再上）；
  `split_unitigs.py`/`pufferize.py`（除非接 pufferfish/参考锚定图）；
  `unitigEvaluator.cpp`（bcalm 对照基准已承担验证职责）；
  minimizer/bloom/mphf 与 debug skip/redo 均为实现细节，不暴露。

## 9. 大规模性能基准与优化方向（2026-08-13）

基准（G37 full 144 MB / 656k reads，k=31）：`anchr asm unitig` 12.8 s /
5.3 GB，bcalm 2.2 s / 554 MB，Bifrost 2.7 s / 29 MB（详见
[unitig-bench.md](../benchmarks/unitig-bench.md)）。

### 9.1 时间分布（perf，full k=31）

* k-mer 计数循环（`TadpoleTable::build` per-chunk 计数）：**65.5%**
* 串行 map 合并：**15.9%**
* SipHash（RandomState 对 40 B key）：**~8.3%**
* rehash 1.4%、canonical 0.8%；真正压缩遍历（build_unitigs/排序）**<1%**

瓶颈是**计数 + 哈希**，不是组装算法本身。

### 9.2 内存 5.3 GB 根因（已确认）

`par_chunks(4096)` 把 656k reads 切成 160 块各建一个 `HashMap<Kmer,u32>`，
然后 `collect()` 全部 160 个 map 到 Vec 再串行 reduce——峰值就是 160 个
map 同时驻留（各 ~30 MB，合计 ~4.8 GB）。distinct k-mer 实测仅 113 万
（最终表 ~110 MB），reads 双份拷贝 ~0.65 GB。

诊断实验（临时把 chunk 调到 16384，40 个 map）：同一 full 输入
**RSS 5.28 → 2.86 GB（-46%），wall 12.8 → 6.1 s（-52%）**。

注意：§7 曾称 `TadpoleTable::build` 是"`pgr::libs::kmer::build_table`
同款模式"——实际并不同款：pgr 是 per-seq 打包 key + 并行基数排序 +
分组计数（FastK/Myers 骨架），anchr 是 HashMap + collect/reduce。
§7 的 radix 评估（20k Lambda 下 radix 更慢）是小规模结论，本节的
5.3 GB 根因说明**大规模下需要回到 pgr 模式**，两者不矛盾。

### 9.3 改进方向（按性价比）

> **前提**：pgr 的 k-mer 基础设施（`libs/kmer/key.rs`、`count.rs`、
> `mod.rs::canonical_keys`、`libs/ds/radix_sort.rs`）借鉴 Myers 的 FastK
> 项目，长期打磨，**不改 pgr 侧**；anchr 侧对齐/复用其模式。

1. **计数对齐 pgr（治本，首选）**：`TadpoleTable` 弃用
   `HashMap<Kmer,u32>`，改用打包 key 模式——per-seq 并行滚动窗口
   （`canonical_keys` 增量正/反互补 + `ceil(key_bytes/2)` 比较）+ 质量
   门控（tadpole 特性；pgr 的 `canonical_keys` 无质量过滤，门控留在
   anchr 的发射循环里）+ 全局并行基数排序 + 分组计数。一次性消除
   9.1 的三个热点与 5.3 GB 内存（打包 key 每条目 `ceil(k/4)` 字节 vs
   40 B，且无 collect 全部 map）；
2. **查询结构**：排序后的打包 key + 二分查找足够（遍历 <1%）；若将来
   遍历成瓶颈再加紧凑开放寻址索引；
3. **reads 双份拷贝减一份**（-0.65 GB）；
4. **算法级：minimizer 分桶**（bcalm 方案，长 k 大数据方向）：minimizer
   裁剪效果随 k 增强、桶可落盘——第二阶段，工作量最大；
5. 遍历并行化：收益小（<1%），不做。

> 分桶基础设施调研与计划见 [unitig-bucket.md](unitig-bucket.md)：anchr
> 已有 hash 分桶实现（`fq norm`/`fq clump` 外部分桶），minimizer 分桶
> 才是真实缺口。

过渡止血（可选，不改变方向）：`collect`+串行 reduce 换成 rayon 树状
reduce 即可立减 ~50% 内存/时间。

### 9.4 长 k 视角（规划是长 k，不止 31）

`key::Kmer` 固定 `k(8B)+[u8;32]` = 40 B（k≤128 恒定），条目内存不随
k 增长；滚动窗口 `push_right` 与哈希成本与 k 基本无关。k 缩放实测
（small 24 MB）：anchr unitig RSS 1245→1192→667 MB、wall 3.03→2.74→
2.45 s（k=31/64/100），bcalm 时间随 k 略升——**长 k 下差距收窄**
（anchr/bcalm 2.7× → 1.5×）。9.3 的对齐收益在长 k 下同样成立，且
内存问题相对更轻。

### 9.5 实施状态（2026-08-13）

9.3 第 1 项（计数对齐 pgr）已实施：`TadpoleTable` 存储改为
`pgr::libs::kmer::KmerTable`（打包 canonical key + u32 计数），构建改
为 per-chunk 质量门控发射打包字节 + rayon 树状拼接 + `count_keys`
（并行 MSD radix + 分组计数）。pgr 侧零改动；`get_count` 改为打包 key
二分，`sorted_entries` 惰性重建 `(Kmer, u32)` 快照（顺序与旧版一致：
`Kmer` 的 `Ord` 即 `to_bytes()` 字节序）。

验证与结果：

* 339 测试全绿；G37 full 的 unitig/contig 输出与旧实现**逐字节一致**
  （三档 N50/条数/总长全部相同）；
* full k=31：unitig 12.8 s / 5.3 GB → **4.4 s / 2.46 GB**（-65% /
  -54%）；small 3.03 s / 1.24 GB → 1.89 s / 547 MB；medium 5.1 s /
  2.9 GB → 2.3 s / 989 MB（详见 [unitig-bench.md](../benchmarks/unitig-bench.md)）；
* **长 k 新发现**：radix 排序成本随 key_bytes 增长（MSD 每字节一级），
  small 上 anchr unitig k=31/64/100 = 1.89/2.42/2.86 s，而旧 HashMap
  实现是随 k 下降（k-mer 变少）。k=100 比值回到 1.6×（vs bcalm）。
  长 k 的进一步优化走 minimizer 分桶（9.3 第 4 项，详见
  [unitig-bucket.md](unitig-bucket.md)），可同时降低有效排序长度与
  中间 key 总量。

### 9.6 发射与查询优化（2026-08-13 续）

perf 定位 k=100 的真实热点不在 radix（~8%）而在**发射循环与遍历查询**：

* 发射：`count_read_kmers_packed` 弃用每窗口 `canonical()`（重算 25 B
  rc，O(k)），改用 pgr `canonical_keys` 的增量模式——`win` 与 `win_rc`
  同步推进（`push_right(x)` + `push_left(3-x)`），canonical 判定只比
  前 `ceil(key_bytes/2)` 字节（FastK 镜像对称）；N 复位两者归零，空窗
  填充 k 步后自然对齐（已由 golden 验证）。发射热点 37.7% → 11.2%；
* 遍历查询：`get_count` 对 `k % 4 == 0`（整字节，如 k=100）用
  `REVCOMP_BYTE[256]` 表算打包 rc（kb 次查表 + 半字节比较），替代
  `key::Kmer::canonical` 的逐碱基 rc（O(k)）；`k % 4 != 0` 保持原路径
  （小 k 下遍历 <1%，无需优化）。遍历 canonical 热点 ~30% 消除。

结果（small，同批 runs 3）：k=31/64/100 = 1.85/1.86/1.98 s（原先
1.89/2.42/2.86 s），比值 vs bcalm 1.5×/1.2×/**1.1×**；内存 ~0.5 GB
持平。339 测试全绿，G37 输出逐字节一致。k 缩放基本拉平后，
minimizer 分桶（unitig-bucket.md）的优先级降低——内存有界化（阶段 A）
成为更相关的下一步。

### 9.7 中间 key 消除与流式读取（2026-08-13 续）

9.6 后 full 峰值仍 2.23 GB：reads 双份拷贝 ~0.6 GB + 含重复的全局
中间 key ~1 GB（132M 窗口 × 8 B）。两项改动：

* **流式读取**（`read_records` 直出 `(seq, phred)`）：只保留一个
  `SeqRecord` 缓冲，去掉全量 `Vec<SeqRecord>` 与二次拷贝（-0.25 GB）；
* **per-chunk 去重 + 树状合并**：每个 chunk 发射后立即 `count_keys`
  （排序去重成小表），rayon `reduce` 两两 `merge_tables`（合并相等
  key 的计数）——**含重复的全局 key 列表不再物化**；合并顺序与结果
  无关（确定性保持）。此即 unitig-bucket.md 阶段 A 的无盘版本，k-way
  合并正是其中标为"需要新增"的那块。

结果（full k=31）：2.23 GB / 4.1 s → **1.28 GB / 2.77 s**（内存 -43%，
时间反而更快——每 chunk 的排序工作总量小于全局排序）。339 测试全绿，
输出逐字节一致。最终对比：

| 指标 | 初始 | 现在 | bcalm |
| :--- | ---: | ---: | ---: |
| full wall | 12.8 s | **2.77 s** | 2.38 s（1.16×） |
| full RSS | 5.3 GB | **1.28 GB** | 555 MB（2.3×） |
| small k=64 比值 | — | **0.7×**（反超） | — |
| small k=100 比值 | — | **1.0×** | — |

剩余内存构成（1.28 GB）：reads ~0.29 GB、排序期 chunk 缓冲 ~0.4 GB、
最终表 + 快照 ~0.07 GB、其余为 contig 结构/分配器。进一步有界化走
unitig-bucket.md 阶段 A 的磁盘分桶（norm 机制），当前规模无需。

chunk 尺寸实测（full k=31，gz 输入）：2048/4096/8192/16384/32768/
65536 = 3.06/2.93/3.00/**2.87**/2.90/2.86 s、1393/1294/1160/**945**/
1119/1142 MB——**16384 为最优点**（排序总量与并发缓冲的折中），已设为
默认。perf 显示排序（partition+msd ≈ 48%）是当前最大项；若后续要再压
时间，方向是减少 per-chunk 排序总量（更大 chunk 或换排序结构），但
内存收益会回吐，需按目标规模权衡。

### 9.8 查询加速与 G37 序列级验证（2026-08-14）

遍历查询（`get_count` + claim 集合）在长 k 下占 ~32%（k=100），两项
改动：

* **前缀索引**：`TadpoleTable` 增加惰性 1-2 字节前缀桶偏移表（65537 /
  257 项，O(n) 一次扫描），`get_count` 在桶内二分——长 k 查询从
  log2(n) ≈ 20 次 25 B 比较降到 ~5 次；
* **FNV claim 集合**：`HashSet<Kmer>` 的 SipHash13 → 自带的 FNV-1a
  hasher（`KmerFnvHasher`，8 字节块乘法），visited/claimed 哈希
  6.6% → ~2%。

结果（small 同批 runs 3）：k=31/64/100 = 1.20/1.08/1.40 s，**全 k 反超
bcalm**（0.9×/0.4×/0.7×），k=64/100 与 Bifrost 打平；full k=31
2.05 s / 0.94 GB（vs bcalm 2.38 s / 555 MB，0.86×）。339 测试全绿，
输出逐字节一致。

**G37 全量序列级对照（2026-08-14）**：`anchr asm unitig` 默认输出
**1482 条 unitig，与 bcalm 全部输出规范化集合完全一致**（总长 622 758
相同）——无损等价（MG1655 的 2403/2403 对照之外的第二个数据集确认）。
此前默认 `min_contig_len = max(124, 2k)` 会把短 unitig 滤掉（只余
116 条），与 bcalm 无损压实语义冲突，已改为**默认不过滤**（见 §10.5）。

性能收官：初始 12.8 s / 5.3 GB → **2.05 s / 0.94 GB**（6.2× / 5.6×），
时间与内存均反超 bcalm。剩余：内存硬上限（阶段 A 磁盘分桶）仅在目标
规模更大时需要；遍历并行化因 BCALM claim 竞态需谨慎设计，收益有限，
暂缓。

## 10. bcalm 后处理语义分析（2026-08-14，`anchr asm unitig` 设计参考）

> 旧 `anchr contained` 与 `unitigs` 模板是**遗留物**（见 §10.4），本节
> 只分析 bcalm 语义本身，作为 `anchr asm unitig` 的对照；遗留组件的
> 实测仅作历史背景。

问题：清理被完全包含的 unitig 很简单，bcalm 为什么保留全部？

### 10.1 bcalm 的流程与输出契约（源码 + 形式化文档确认）

bcalm 2（仓库 `bcalm/`，GATB/bcalm）流程：`bcalm`（计数 + 按 minimizer
分区压实）→ `bglue`（分区 unitig 合并）→ `links`（输出 `L:` 边）。
**没有任何"清理被包含 unitig"的步骤**。

形式化定义（`bidirected-graphs-in-bcalm2.md`）：
* unitig = 图中不可扩展、不重复顶点的行走（最大非分支路径）；
* **"Maximal unitigs should be a vertex decomposition of the graph; in
  particular, two maximal unitigs should not share a vertex"**——两条
  unitig 不共享 k-mer 顶点；
* 因此**同链精确包含在语义上不可能**（若 A 的序列是 B 的子串，A 的
  k-mer 就是 B 路径的顶点，A 就不是最大 unitig）。

G37 实证：bcalm 1482 条中 ≥124 bp 的 116 条之间，精确包含对（同链或
rc）**为 0**。

### 10.2 bcalm 为什么留着"看起来被包含"的 unitig

1. **精确被包含的 unitig 根本不存在**（顶点分解性质，10.1）；
2. 对齐上"看起来被包含"的其实是两类，都不是图冗余：
   * **近相同重复拷贝**：k-mer 不同 → 图顶点不同 → 都是合法 unitig
     （bcalm 的顶点分解只约束精确 k-mer）；
   * **短 unitig**（如 <1000 bp）：不是被包含，只是短；
3. **输出契约 = 无损压实图**：每个 solid k-mer 恰好出现在一条 unitig。
   删除任何一条都会破坏顶点分解、丢失 k-mer 信息（有损）。删冗余是
   **组装层的决策**（组装要唯一长序列；pan-genome/图工具要全部节点），
   不是图压缩操作——bcalm 正确停留在无损层。

### 10.3 历史背景（旧 `anchr contained` 实测，不计入方向）

遗留模板在 bcalm 之后跑 `anchr contained`（`--len 1000 --idt 0.9999
--ratio 0.99999`），把 G37 的 1482 条砍到 84 条。实测其包含过滤
（idt/ratio）对 bcalm 输出**零删除**——1482 → 116 全由 `--len 124`、
116 → 84 全由 `--len 1000`。即旧步骤的实际作用只是长度过滤；该组件
属遗留物，不参与方向决策。

### 10.4 对方向的启示

* `anchr asm unitig` 与 bcalm 逐序列一致（§9.8）＝ 正确的无损压实语义，
  不内置任何有损清理（顶点分解契约）；
* 旧 `anchr contained`/`unitigs` 模板是遗留物，**不作为方向**——若要
  精简 unitig 供组装，是独立的新决策，语义（有损组装精简）与无损
  压实分开设计，不沿用旧命令；
* 方向以 `anchr asm unitig` / `asm contig` 为准（§8-§9），下游是否需要
  "精简输入"待定，不预设。

### 10.5 默认长度过滤的移除（2026-08-14）

随 §10 结论落地：`anchr asm unitig` 的默认 `min_contig_len` 从
`max(124, 2k)` 改为 **0（不过滤）**——与 bcalm 一致，输出完整的顶点
分解（G37 full：1482 条，与 bcalm 逐序列等价，§9.8）。

* `asm contig` 保留 tadpole 的 `mincontiglen` 自动默认 `max(124, 2k)`
  （tadpole 兼容语义，§2；由 contig 命令显式传入）；
* `asm olc`（多 k unitig OLC）随共享默认变为不过滤（unitig 无损输入）；
* `--min-contig-len` 选项保留：需要组装精简输出的用户可显式指定
  （`0` 以外的值即过滤）；
* 339 测试全绿（unitig/contig/olc 测试均显式传 `--min-contig-len`，
  不依赖旧默认）。

## 11. DFA 状态分类（借鉴 cuttlefish）评估（2026-08-14）

把 cuttlefish 的"先分类、后提取"路线按我们的语义移植成可选引擎
（**2026-08-14 起为 `asm unitig` 默认引擎**；`--no-dfa` 回退桶扫描，
实现见 `src/libs/asm/dfa.rs`）：

* **分类 pass**：对每个 solid canonical k-mer 预计算 4 个字段
  `VertexState{in_count, out_count, in_base, out_base}`（入/出度与唯一
  延展碱基，度 2+ 统一编码为 2）。分类只读计数表、按排序顶点表分块并行
  （rayon），**无 CAS/锁**——与 cuttlefish 的边扫描 CAS 不同，因为我们的
  计数已经物化在 `TadpoleTable` 里。
* **walk pass**：保持原 `build_unitigs` 的种子扫描顺序、visited 与
  环状处理；唯一变化是每次延展用 `out_base()`/`in_count()` 的 O(1) 状态
  查询替代 4 桶扫描。方向语义（forward/RC 互换 in/out、碱基取补）由
  canonical 判断推导，与旧遍历逐位等价。
* **阈值语义**：状态分类用的是 **k-mer 计数阈值**（`--min-count-seed`），
  不是 cuttlefish 的 (k+1)-mer 阈值，因此输出仍与 bcalm 逐序列一致。
* **全流水线绑定**（2026-08-14）：`--parallel N` 现在同时约束计数
  （`TadpoleTable::build_threaded`）与 DFA 分类两个可并行阶段；`auto`
  用满全部逻辑核（rayon 全局池）。**默认 `--parallel` 为 8**（2026-08-14
  改为 `min(逻辑核/2, 8)` 自适应，避免默认吃满逻辑核；`auto` 仍可显式
  指定）。walk 仍保持确定性单线程。

### 11.1 正确性

Lambda（k=31）与 G37 full（k=31/99）上，`--dfa --parallel 1/4/8` 的输出
与默认引擎 **逐字节一致**（G37 1482 / 568 条 unitig）。新增回归测试
`command_asm_unitig_dfa_matches_default`（线性 + 环状两组）。

### 11.2 性能与并行缩放（G37 full，hyperfine + `/usr/bin/time -v`）

| 引擎 | k=31 wall | k=31 RSS | k=99 wall | k=99 RSS |
| :--- | ---: | ---: | ---: | ---: |
| 默认 `auto`（全核） | 2.18 s | 929 MB | — | — |
| 默认 `-p 1` | 11.39 s | 373 MB | 12.74 s | 348 MB |
| 默认 `-p 4` | 4.15 s | 490 MB | 5.14 s | 747 MB |
| 默认 `-p 8` | 2.69 s | 680 MB | 3.82 s | 1019 MB |
| `--dfa -p 1` | 11.69 s | 474 MB | 12.68 s | 401 MB |
| `--dfa -p 4` | 3.77 s | 531 MB | 4.25 s | 744 MB |
| `--dfa -p 8` | **2.50 s** | 665 MB | **2.95 s** | 1017 MB |

分类阶段单独计时（`ANCHR_DFA_TIMING=1`，G37 full 单次）：

| k | t=1 | t=4 | t=8 | 8 线程加速 |
| :--- | ---: | ---: | ---: | ---: |
| 31 | 0.72 s | 0.25 s | 0.20 s | 3.7× |
| 99 | 1.34 s | 0.38 s | 0.23 s | 5.9× |

参考项目 cuttlefish 同机缩放（G37 full k=31，`-t 1/2/4/8/16`，runs 2）：
11.43 / 12.03 / 8.51 / 8.20 / 8.37 s——峰值仅 **1.39×@8**（t=2 反而
更慢，t=16 回落），总并行效率 ~17%。它的瓶颈阶段（KMC I/O、MPHF 构建、
提取）本身也不线性扩展。

分析：

* **`--parallel` 现在是全流水线语义**：1→4 线程加速约 2.7-3.0×（效率
  70-75%），4→8 约 1.5×（效率 ~40%），典型 Amdahl；`auto`（32 逻辑核）
  是默认的最快配置。不同线程数输出**逐字节一致**（计数合并与 walk 均
  与并行度无关）；
* **DFA vs 默认引擎同线程比较**：p=1 几乎持平（计数占绝对大头），
  p=4/8 快 8-23%（k=99 p8：2.95 vs 3.82 s）——分类把部分邻居查询变成
  O(1) 状态查找，长 k 收益更大；
* 分类阶段单段并行效率（t=4 73%/89%、t=8 46%/74%）**明显优于
  cuttlefish 的总缩放**；
* RSS：单线程时 `--dfa` 比默认多 ~15%（状态 4 B/顶点 + canonical→索引
  HashMap，G37 full ~113 万顶点 +20-70 MB）；线程数升高后计数缓冲占
  主导，差距消失（k=99 p8 两者都在 ~1.0 GB）；
* walk 仍单线程，遍历占比 <1%（§9.1），并行化收益有限且引入排序/确定性
  负担，暂不做。

> **并行 walk 实验已撤下（2026-08-14）**：曾实现"原子 seed 认领 +
> 最小 k-mer 拥有者保留"的并行 walk（`--parallel-walk`），但并行计数 +
> 分类 + walk 多层线程池叠加会把系统跑满（实测卡死），已整体移除；
> 默认 `--parallel` 改为自适应 `min(逻辑核/2, 8)`（本机为 8）。

### 11.3 后续方向

* **状态紧凑化**：`VertexState` 可压缩到 ~2-3 bit/顶点（度编码 + 唯一
  碱基），并用排序顶点表的二分/前缀索引替代 HashMap，为外存路径做准备；
* **覆盖度字段**：walk 后仍需 `calc_coverage` 从计数表取均值/min/max，
  若做"极低内存"模式需要压缩计数（u8/u16 饱和）或仅对 solid 顶点保留
  计数；
* **与分桶路径结合**：分类要求全局顶点表（类似 cuttlefish 的 KMC 全局
  枚举）；hash 分桶下 unitig 跨桶，仍需要 bcalm 式的 glue 阶段，状态
  分类应放在"全局排序/MPHF 之后"这一层。

## 12. supermer 两段计数接入（2026-08-14）

pgr 侧完成 FastK 式 super-mer/minimizer 两段计数
（`pgr::libs::kmer::supermer`，固定 m=12；stage 1 折叠 + 排序 +
stage 2 加权展开，输出与直接路径逐字节一致，测试覆盖 k=3..40、
重复/反向互补/N 处理）。anchr 以 `asm unitig --supermer` 接入
（`TadpoleTable::build_supermer`，无质量门控），Cargo.toml 的 pgr rev
升至 `769f82f`。

正确性：Lambda + G37 full（k=31/99）与默认引擎**逐字节一致**。

性能（G37 full，runs 2，`auto`）：

| 引擎 | k=31 wall | k=31 RSS | k=99 wall | k=99 RSS |
| :--- | ---: | ---: | ---: | ---: |
| 默认（直接排序） | 1.75 s | 976 MB | 2.61 s | 2010 MB |
| `--supermer` | 1.87 s | 999 MB | 3.62 s | 1442 MB |

结论：

* 输出一致；k=31 总耗时与默认**相当**（本批 +6%，噪声内），k=99 明显
  更慢（+39%）但 RSS 降 **28%**（2010 → 1442 MB）；
* pgr 自测（mg1655 663k reads）k=31 lib 计时 0.87 vs 1.74 s 的收益在
  anchr 端到端**没有兑现**：计数只是总耗时一部分，且本接入需克隆 reads
  序列、固定 m=12 对短读折叠有限；
* 与 pgr 结论一致：**super-mer 不是长 k 短读的解**（k=100 折叠仅
  ~1.1×）；k=31 的收益依赖数据冗余度；
* 后续：pgr 若支持质量门控与按数据自适应 minimizer，再重新评估；当前
  `--supermer` 作为实验开关保留，默认路径不变。

### 12.1 尝试背景与暂缓原因（记录）

* **动机**：FastK G37 实测 1.27 s / 226 MB（仅计数）vs anchr
  2.18 s / 929 MB（计数 + unitig）；[unitig-bucket.md](unitig-bucket.md)
  §3.1 确认差距在计数层，且我们已在用 MSD radix（FastK 用的是 LSD），
  缺口是 **super-mer/minimizer 两段式**，不是排序方向。
* **尝试**：pgr 实现 `supermer` 模块（commit `769f82f`，固定 m=12，
  stage 1 折叠 + 排序 + stage 2 加权展开，输出与直接路径逐字节一致，
  测试覆盖 k=3..40、重复/反向互补/N）→ anchr 以 `asm unitig --supermer`
  接入并验证（Lambda + G37 full k=31/99 逐字节一致），Cargo.toml pgr
  rev 同步升到 `769f82f`。
* **结果**：端到端无收益——k=31 与默认相当（±6%，噪声内），k=99
  +39% 但 RSS -28%；pgr 自测的 lib 级收益（mg1655 k=31：0.87 vs
  1.74 s）未兑现。原因：计数只是总耗时一部分、接入需克隆 reads 序列、
  固定 m=12 对短读折叠有限。
* **暂缓条件（重新评估的触发条件）**：
  1. pgr 在 supermer 路径支持**质量门控**——当前它无质量过滤，FASTQ
     语义与默认引擎不一致（FASTA 才等价）；有门控后才能无差别用于
     FASTQ；
  2. pgr 支持**按数据自适应 minimizer**——FastK 按输入训练 PAD_LEN，
     固定 m=12 未必匹配数据集；注意长 k + 短读时折叠本身只有 ~1.1×，
     自适应只能改善、不能根治。
  两者任一落地后重跑 §12 的基准表，再决定 `--supermer` 是否转正。

### 12.3 pgr 计数优化落地与 FASTA 默认切换（2026-08-14）

anchr 侧已收敛（DFA 默认 walk、流式 direct 计数、状态表内嵌计数、
单池复用、`--parallel` 默认 `min(逻辑核/2, 8)`）。剩余最大单一差距在
pgr 计数：G37 full k=31 约 **0.9 s** vs FastK **0.68 s**（同机、纯计数）。

给 pgr 的优化清单（按收益排序）：

1. **自适应 minimizer**：pgr supermer 固定 m=12；anchr 侧实测
   k=31 最优 m=8、k=63/99 最优 m=12、k=127 两者接近（启发式
   `min(12, max(5, k/4))` 已进 anchr CLI，pgr 可内置同样选择）；
2. **连续缓冲替代每序列 `Vec`**：supermer stage-1 目前 `per_seq: Vec<Vec<u8>>`
   逐个分配，可改为单块连续缓冲 + 偏移表（FastK 同款），减少分配与
   缓存未命中；
3. **避免克隆/移动开销**：`build_table` 接受 `&[Vec<u8>]`，调用方需先
   备好全量 seqs；可加"借用切片"API（`&[&[u8]]` 或迭代器）让 anchr
   流式路径直接喂；
4. **调 chunk / 排序阈值**：`radix_sort_bytes_par` 的 `PAR_SMALL`
   （1<<18）与 supermer stage-2 的展开粒度可按数据集调；
5. **（可选）并行读输入**：FastK `-T` 多线程读文件；pgr 读取单线程，
   gz 输入时差距更明显。

复现（同机、plain FASTA）：

```bash
# FastK（纯计数+profile）
/usr/bin/time -v FastK -T32 -t1 -k31 -N/tmp/fk31 -P/tmp /tmp/pe.cor.fa
# pgr 直接 / supermer
/usr/bin/time -v pgr kmer table /tmp/pe.cor.fa -k31 -o /tmp/d.pkt
/usr/bin/time -v pgr kmer table /tmp/pe.cor.fa -k31 --supermer -o /tmp/s.pkt
```

目标：pgr supermer k=31 端到端 ≤ 0.7 s、内存 ≤ 800 MB（当前 1.01 s /
854 MB，见 `notes/benchmarks/bench-supermer-vs-fastk.md`）。

**落地（pgr commit `b31af11`，2026-08-14）**：§12.3 清单全部完成——
自适应 minimizer（`min(12, max(5, ceil(k/4)))`）、stage-1 连续打包、
借用切片 API（`build_table_slices*`）、stage-2 group/expand 并行；
k=31 lib 0.87 → **0.595 s**（-32%）。anchr 接入新 API（bump rev +
slices 调用），并把 supermer 设为 **FASTA 默认计数**（FASTQ 自动回退
direct 保质量门控；`--no-supermer` 强制 direct）。G37 full half(8)：
k31 **1.358 s / 597 MB**（vs direct 2.103 s / 602 MB，-35%）、k99
**2.187 s / 1131 MB**（-13%）；`--parallel auto` k31 可到 1.056 s。

### 12.2 效率攻坚：walk 瓶颈与组合优化（2026-08-14）

分阶段计时（G37 full k=31，`--supermer --dfa -p8`）：

| 阶段 | 耗时 | 说明 |
| :--- | ---: | :--- |
| read_records | 0.12 s | 读 144 MB FASTA |
| supermer count | 0.85 s | pgr 两段计数（move+pack+sort+expand） |
| walk+classify+coverage | 0.47 s | dfa classify 0.10 s + walk/覆盖度/输出 |

关键发现：**walk+build 曾是隐藏瓶颈**（默认引擎 ~0.9 s，和计数同级）——
早期 §9.1 的"遍历 <1%"是在计数 12.8 s 时代测的，计数提速后 walk 占比
被放大。两个修复：

* **DFA 状态引擎**把每步 8 次邻居查询换成 O(1) 状态查询（classify 并行
  0.10 s @8 线程）；
* **`solid_entries`** 只物化 count ≥ threshold 的顶点，不再构造 55 万
  个低丰度 `Kmer`（walk+build 0.69 → 0.47 s）。
* **DFA 设为默认引擎（2026-08-14）**：full k31 **1.324 s / 948 MB** vs
  桶扫描 1.578 s / 1030 MB（-16% / -8%）；k99 **1.880 s / 1998 MB** vs
  2.622 s / 2095 MB（-28% / -5%），输出逐字节一致。`--no-dfa` 保留旧
  引擎作对照。
* **direct 计数流式化（2026-08-14）**：unitig 默认路径改为
  `TadpoleTable::build_streamed`——按 32768 记录/块经 bounded channel
  分发给 worker（std 线程），不保留全部 reads。G37 full：k31
  **1.43 s / 694 MB**（vs 内存路径 1.32 s / 948 MB：墙钟 +8%、内存
  -27%）；k99 **1.88 s / 1455 MB**（vs 1.88 s / 1998 MB：墙钟持平、
  内存 -27%）。`ANCHR_STREAM_CHUNK`/`ANCHR_STREAM_CAP` 可调
  （16384/cap2 更省内存但墙钟 +20%）。
* **supermer 成为 FASTA 默认计数（2026-08-14）**：pgr §12.3 落地后
  half(8) 下 k31 1.36 s vs direct 2.10 s（-35%）、k99 2.19 s vs
  2.52 s（-13%）；FASTQ 自动回退 direct 保 `min_prob` 语义，
  `--no-supermer` 强制 direct（见 §12.3）。
* **DFA 状态表内嵌计数（2026-08-14）**：分类时把 canonical 计数也存入
  `VertexStates.counts`，walk 的覆盖度/`ab:Z` 从 O(1) 索引取，不再每
  k-mer 重做 `get_count` 的 canonical + 二分。G37 full 默认组合进一步
  到 **k31 1.344 s / 670 MB**、**k99 1.69 s / ~1.6 GB**（walk+build
  0.62 → 0.36 s），输出逐字节一致。
* **分类条目复用 + visited 字节数组（2026-08-14）**：`VertexStates`
  保留排序后的 solid 条目（walk 不再二次构建 `solid_entries`），
  `visited` 从 FNV HashSet 换成按顶点索引的 `Vec<u8>`。G37 full 默认
  组合到 **k31 1.264 s / 699 MB**、**k99 1.625 s / 1348 MB**
  （walk+build 0.36 → 0.27 s），输出逐字节一致；`--no-dfa` 回退路径
  也已验证一致。
* **walk 覆盖度累积**：unitig 两条 walk 在扩展时直接累积 canonical
  计数（右循环 `right_counts` + 左循环 `left_counts`，按
  `reverse(left[1..]) + right` 拼回输出顺序），`calc_coverage` 与
  `calc_abundances` 的二次扫描已从 unitig 路径删除（contig 路径保留
  `calc_coverage`）。实测墙钟无显著变化——二次扫描本就不是大头，
  但逻辑上少了两遍滑窗 + canonical 重算，且 `ab:Z` 直接复用同一份
  计数。

minimizer 扫描（`--supermer --dfa -p8`，runs 2，输出均与默认逐字节一致）：

| k | m=5 | m=8 | m=10 | m=12 | 默认（直接+桶扫描） |
| :--- | ---: | ---: | ---: | ---: | ---: |
| 31 wall | 1.454 s | **1.442 s** | 1.475 s | 1.491 s | 1.574 s |
| 31 RSS | 779 MB | 791 MB | 831 MB | 841 MB | 945 MB |
| 99 wall | 2.861 s | 2.861 s | 2.868 s | **2.825 s** | 2.701 s |
| 99 RSS | 1338 MB | 1311 MB | 1301 MB | 1283 MB | 2044 MB |

当前最佳组合：**k=31 `--supermer --dfa -p8 -m8`**（1.44 s / 791 MB，
比默认快 8%、内存省 16%）；**k=99 `-m12`**（2.83 s / 1283 MB，比默认
慢 5%、内存省 37%）。剩余差距在计数（0.85 vs FastK 0.68 s）与
walk/覆盖度（0.37 s）；计数由 pgr 侧继续优化，覆盖度可在 anchr 侧用
walk 中累积计数替代二次扫描。
