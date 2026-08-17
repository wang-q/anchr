# OLC 组装设计：多 k unitig 层的 overlap → layout → consensus

> 状态：**已实现（2026-08-12，M1–M4 全部落地）**。命令：
> `anchr asm ovlp` / `layout` / `cns` / `olc`；库：`libs/olc/`
> （overlap / layout / consensus，anchr 业务库）。合成基因组端到端验证通过
> （30× 无错 reads → contigs 全部为基因组精确子串，最长覆盖 97.5%；
> 低覆盖 6× 出现重复区经典错装，符合预期）。
> **真实数据验证（Lambda，2026-08-12）**：见 §12——抓出并修复左向延伸
> 坐标下溢 bug；多 k 冗余 2.4×（v1 待消减）；参考菌株差异会误报
> "嵌合"，需 reads 侧验证。
> 需求来源：用户裁定"不对 reads 做 OLC，把不同 k 各自生成的 unitigs 当
> 伪 reads，在 unitig 层做 OLC 拼接"（`todo.md` §3、`references/canu.md` §8）。
> 参考源码：`canu-2.3/`（bogart + utgcns）、`wgs-8.3rc2/`（AS_BAT + AS_CNS），
> 分析笔记见 `references/canu.md` / `references/celera.md`。

## 1. 目标与范围

为 pgr 增加一条完整的 OLC 组装流程：

```
reads ──(多 k 各生成 unitigs)──> 伪 reads ──(overlap)──> PAF
  ──(layout)──> 路径 ──(consensus)──> contigs.fasta
```

**核心裁定（用户，2026-08-12）**：

* 不做 reads 级 OLC——reads 先用 DBG（`anchr asm unitig`）压缩成 unitigs；
* 多 k（默认 21/51/81）各生成一套 unitigs，合并后当伪 reads；
* unitig 语义 = 最大无分支路径（bcalm graph3 移植，无气泡），OLC 只处理
  unitig 间 overlap，不引入"平行路径选哪条"的启发式；
* 气泡/孤儿合并不处理（既有裁定）。

**成功标准**：合成基因组 reads → `anchr asm olc` 输出 contigs 覆盖完整基因组
（identity 100%，因 overlap 精确）；Lambda 真实 reads 冒烟测试输出合理
（contig 数 / N50 与 tadpole contig 同量级）；1701+ 测试全绿、fmt/clippy 干净。

## 2. 为什么 unitig 层 OLC 可行

* **数据量**：unitigs 数远小于 reads（宏基因组下少 1~2 个数量级），
  all-pairs overlap 成本从不可行变可行；
* **规避气泡**：unitig 无分支，重叠只发生在 unitig 间，天然无
  "平行路径投票"问题；
* **多 k 互补**：小 k 连通性好（低覆盖区、重复边界），大 k 特异性强
  （区分重复/菌株）；不同 k 的 unitigs 共享精确子串，overlap 证据天然存在；
* **精确性**：unitig 序列来自 DBG 固实 k-mer，无测序错误，overlap 可做
  全精确验证（不做 Myers/edlib 扩展，参考项目的高噪声机制不需要）。

## 3. 参考实现对照（详见 references/canu.md、celera.md）

| 环节 | Celera 8.3rc2 | Canu 2.3 | pgr 借鉴 |
|---|---|---|---|
| overlap | OlapFromSeedsOVL（k=9 seed + banded DP） | MHAP / overlapInCore（k-mer 哈希 + Myers） | seed→verify 骨架同源；pgr 是精确版（无错误模型） |
| layout | bogart（AS_BAT）：BestOverlapGraph + greedy 双向延伸 + 覆盖度证据 repeat split | bogart（同源）：互惠 best edge 种子 + 单向延伸 + markRepeatReads | greedy best-edge + 互惠种子；repeat 用"双定位"思想（v1） |
| consensus | AS_CNS：MA 列投票（BaseCallMajority） | utgcns：template stitch + edlib 重比对 + POA-DAG bestPath | 精确 overlap 下缝合即可；列投票 = 将来的鲁棒化方向 |

## 4. 管线设计（四阶段）

### S0 伪 reads 生成（复用，不新写）

每个 k 跑 `anchr asm unitig`（`libs/asm/assemble.rs::assemble_unitigs`，
bcalm graph3 压缩语义，默认 k=31 / solid ≥3）。产出物：unitig FASTA。

**命名**：`asm unitig` 的输出名恒为 `unitig_<id>`，多 k 合并必然撞名。
OLC 阶段统一重命名：`<tag>:<name>`，tag 默认取输入文件 stem
（仅保留 `[A-Za-z0-9_.-]`，空则用文件序号）——确定性且可回溯到 k。
（`anchr asm olc` 驱动器内部直接用 `k<k值>` 作 tag。）

### S1 overlap 检测（新 `libs/olc/overlap.rs`）

**算法**：seed → verify（与 `libs/map.rs` MapIndex 同构，精确版）：

1. 建索引：对全部 unitigs 做 canonical k-mer 索引
   （`MapIndex` 形态：`keys: Vec<u8>` 打包 FastK 字节 + `payloads:
   Vec<u64>` 存 `(cid<<32)|pos`，`radix_sort_bytes` 排序），seed k 默认 17
   （`--overlap-k`，≤ 目标最小 unitig 长度，越界则自动降到 min(k, len)）；
2. 候选：对每条 unitig q 查**边界 k-mer**——5' 端窗口 (0..k) 与
   3' 端窗口 (n-k..n)，各查正链与反互补（canonical 索引丢失方向，
   验证时双侧都试）；命中 → (cid, tpos)；
3. 验证：从 seed 对齐处向两端逐碱基扩展，得**最大精确 overlap**
   （含 seed 的 q∩t 最长精确段，长度 L ≥ k）；
4. 分类（按 q/t 覆盖关系）：
   * `dovetail`：q 5'/3' 端与 t 3'/5' 端重叠（两端各留出 >0 的非重叠段），
     或 q 完全包含于 t（contain，长度 ≥ L）；
   * `contain`：q ⊂ t 或 t ⊂ q——不参与延伸，留作共识覆盖证据；
5. 输出 PAF（复用 `pgr::libs::paf/record.rs` 12 列 + `ov:A:D|C` tag），
   去重（同一对多 seed 命中取最长）、排除自身（q==t 及回文 rc）。

**并行**：rayon 按 unitig 并行查询；索引构建与 `anchr asm map` 同路径。

### S2 layout（新 `libs/olc/layout.rs`）

**算法**：bogart 风格 greedy best-edge 路径延伸（简化为无 mate、无气泡）：

1. 只取 dovetail overlap 建**有向图**：node = unitig，edge = 一端 → 另一端
   （q 3' 端 → t 5' 端 / q 5' 端 → t 3' 端，方向由 overlap 坐标推出）；
2. 每 node 两端各选 best edge（最长 overlap L，平局按 (target, L) 字典序
   保证确定性）；**互惠**要求：seed 的 3' best 必须同时是对方相应端的
   best（Canu 互惠种子思想，防错装）；
3. greedy：按 unitig 长度降序取未放置 seed，沿 3' best 单向延伸，
   目标已放置 / 无 edge / 目标被标记 repeat 即停；反向延伸同理
   （通过 rc 复用同一逻辑）；
4. **repeat 标记（v0 简版）**：某端的 top2 edge 长度 ≥ 0.9×best 且指向
   不同 node → 该端标记 repeat，禁止从它延伸（Canu `markRepeatReads`
   的"双定位"思想的单元化近似；覆盖度证据版本留 v1）；
5. 输出 layout TSV（每 contig 一行一步）：
   `contig_id step unitig_name strand q_start q_end overlap_len`，
   `strand` 为 unitig 在 contig 中的方向，`q_start/q_end` 是其在
   contig 坐标系中的 0-based 区间。

### S3 consensus（新 `libs/olc/consensus.rs`）

overlap 全精确 ⇒ consensus = 沿 layout **精确缝合**：

1. 按路径顺序取每步 unitig 的方向片段，与上一步 overlap 部分对齐
   （坐标由 layout 记录），追加非重叠后缀；
2. 输出 FASTA：`>contig_<id>,len=...,cov=...`（cov = 路径上 unitigs
   平均覆盖深度，近似 = 参与步数），70 列换行，与 `asm contig` 输出风格一致；
3. `--min-contig-len` 过滤短 contig。

**列投票留 v1**：若未来引入错配 overlap 或真实数据暴露 junction 不一致，
再加 AS_CNS `BaseCallMajority` 式逐列投票 + min-coverage 修剪
（Canu `consensusNoSplit` 语义），复用 `anchr asm map` + `anchr sam to-rg` +
`pgr rg coverage` 的回放设施（`references/canu.md` §8.3 已论证）。

## 5. 命令设计

新增 `anchr asm` 三个叶子命令 + 一个驱动器（四层：`libs/olc/*` 管逻辑，
`cmd/asm/*` 薄壳）：

| 命令 | 输入 → 输出 | 逻辑 |
|---|---|---|
| `anchr asm ovlp` | unitig FASTA(s) → PAF | `libs/olc/overlap.rs` |
| `anchr asm layout` | PAF + unitig FASTA → layout TSV | `libs/olc/layout.rs` |
| `anchr asm cns` | layout TSV + unitig FASTA → contigs FASTA | `libs/olc/consensus.rs` |
| `anchr asm olc` | reads → contigs FASTA（驱动器） | 内部组合 S0–S3，阶段间走内存 |

`anchr asm olc` 参数：

```text
anchr asm olc <infiles>... -o contigs.fa \
    --kmer 21,51,81          # 逗号分隔，默认 21,51,81
    --min-count-seed 3       # 透传 asm unitig
    --overlap-k 17           # S1 seed k
    --min-overlap 34         # 最短接受的 overlap 长度
    --min-contig-len 500     # 输出过滤
    --keep-dir DIR           # 调试：落地中间文件（unitigs/ovlp/layout）
```

阶段命令天然可独立测试与组合（也支持用户自己跑
`asm unitig` → `asm ovlp` → `asm layout` → `asm cns` 的管道形态）。

## 6. 数据结构与格式

### overlap（PAF）

* 12 列标准 PAF（q 名 = `<tag>:<name>`，q 长 = unitig 长），
  `matches = block_length = L`，`mapq = 255`；
* tag：`ov:A:D`（dovetail）/ `ov:A:C`（contain）——不做 `cg:Z`（无错配无 CIGAR）。

### layout（TSV，无表头）

```text
contig_1<TAB>0<TAB>k21:unitig_5<TAB>+<TAB>0<TAB>2410<TAB>0
contig_1<TAB>1<TAB>k51:unitig_12<TAB>+<TAB>2410<TAB>4370<TAB>180
```

`q_start/q_end` 为该 unitig 在 contig 中的区间；`overlap_len` = 与上一步的
overlap（第 0 步恒 0）。同 contig 内区间连续（`q_end[i] == q_start[i+1]`）。

## 7. 现有基础设施复用映射

| 环节 | 复用 | 用途 |
|---|---|---|
| S0 | `libs/asm/assemble.rs::assemble_unitigs`（anchr 业务） | unitig 生成（命令层 `asm unitig` 已包装） |
| S1 | `libs/map.rs`（MapIndex 形态 + `canonical_keys` + radix，anchr 业务） | canonical k-mer 种子索引 |
| S1 | `pgr::libs::kmer/key.rs::Kmer` | 边界 k-mer 编解码 / rc / canonical |
| S1 | `pgr::libs::nt::rev_comp` | 方向验证 |
| S1 | `pgr::libs::paf/record.rs` | PAF 写出 |
| S1 | `pgr::libs::ds/radix_sort.rs::radix_sort_bytes` | 索引排序 |
| S2 | `pgr::libs::ds/dsu.rs`（仅若需要连通分量） | 布局分组（v0 可不用） |
| S3 | `pgr::libs::fmt/seq.rs::SeqReader` | unitig FASTA 读取 |
| 全部 | `pgr::libs::io.rs` reader/writer、`cmd/args.rs` 标准参数 | I/O 与 CLI 一致性 |
| 驱动 | `libs/asm/assemble.rs` + 上述各 libs | 内存组合，无中间文件 |

**不引入新依赖**（AGENTS.md 硬性要求）；k-mer 表示统一用 FastK 字节键
（`pgr::libs::kmer` 的唯一表示），与 `pgr kmer`/`pgi`/`anchr asm map` 同套。

## 8. 验证计划

### 单元测试（libs/olc/）

* overlap：构造两个已知精确 suffix/prefix overlap 的 unitigs → PAF 记录
  的坐标/L/方向正确；contain 与 dovetail 分类正确；rc 方向正确；
  重复 k-mer（polyA 区）不产生错误重叠；
* layout：线性链 / 分支（bubble）→ 只走一条；repeat unitig → 路径断开；
  确定性（相同输入两次运行逐字节一致）；
* consensus：缝合正确性（含跨多 k 的 contain 不引入重复碱基）。

### 集成测试（tests/cli_asm_olc.rs）

* 合成基因组（随机 ~2 kb 序列）→ 生成多份 reads（子串 + rc，覆盖 ~20×）
  → `anchr asm olc` → contigs 与基因组逐段精确一致（identity 100%）；
* 阶段管道形态（`asm unitig` ×3 k → `asm ovlp` → `asm layout` →
  `asm cns`）与驱动器输出一致；
* Lambda 真实 reads 冒烟（`tests/bbtools/Lambda/R1.2k.fq.gz` 等）：
  contig 数 / N50 合理、无 panic。

### 验收门

每里程碑 `cargo fmt` + `cargo clippy --all-targets -- -D warnings` +
`cargo test` 全绿；命令注册与 docs 由 `cli_consistency.rs` 约束。

## 9. 里程碑

| 里程碑 | 内容 | 验证 |
|---|---|---|
| M1 | `libs/olc/overlap.rs` + `anchr asm ovlp` | ✅ 5 单测 + 3 集成测试 |
| M2 | `libs/olc/layout.rs` + `anchr asm layout` | ✅ 5 单测（线性/反向/分支/互惠/contain） |
| M3 | `libs/olc/consensus.rs` + `anchr asm cns` | ✅ 4 单测（含不一致 overlap 友好报错） |
| M4 | `anchr asm olc` 驱动器 + 集成测试 + `docs/asm.md` + todo/笔记更新 | ✅ 合成基因组重建 + 确定性 + 阶段管道等价 |

### 实现说明（与设计定稿的偏差）

* layout 的互惠检查实现在**连接端**（target 的 junction end 的 best edge
  指回当前 unitig），而非自由端——线性链因此可连续延伸；
* 阶段命令 `layout` 需要 unitig FASTA（不止 PAF），用于长度与命名校验，
  `cns` 同理；命名逻辑抽到 `cmd/asm/common.rs` 三命令共用；
* 驱动器的 unitig 命名 `k<k>:unitig_<id>`（不是 `<stem>:`），阶段管道用
  文件 stem；两者互不冲突（驱动器 `--keep-dir` 产物可直接喂阶段命令）；
* `asm unitig` 新增内存版 `assemble_unitigs_buf`（`assemble.rs` 最小重构：
  核心逻辑抽出 `assemble_unitigs_core`，写盘版行为不变）。

## 10. 不做 / 待决

* **不做**：气泡/孤儿合并（用户裁定）；允许错配的 overlap（unitig 精确
  假设，DBG 错误在固实阈值下已排除；真实数据暴露再议）；scaffolding
  （无配对需求）；列投票 consensus（v1 待真实数据）。
* **待决（数据驱动）**：repeat breaking 的覆盖度证据阈值（Canu
  `SPURIOUS_COVERAGE_THRESHOLD=6` / `ISECT_NEEDED_TO_BREAK=15` 的
  单元化版本）；contain unitig 是否参与 consensus 投票；短 unitig
  （< seed k）的 overlap 缺失处理。
* **v1 素材来源（2026-08-12）**：repeat breaking 的覆盖度阈值与实现路径
  有两个成熟参考——`references/skesa.md` §7.1（`FilterLowAbundanceNeighbors`
  fraction=0.1 多层过滤 + 可逆性检查）与 `references/metaMDBG.md` §9
  （渐进丰度过滤 t=1.1/10% 步长 + RepeatRemover 的桥接 reads 证据，
  pgr 用 `anchr asm map` + `anchr sam to-rg` + `pgr rg coverage` 回放即等价设施）；
  多 k 反馈（SKESA clean_reads / metaMDBG unitig 反馈）为 v2 候选。
* **参考**：`canu-2.3/src/bogart/`、`wgs-8.3rc2/src/AS_BAT/` 源码随取随用；
  不引入其代码/依赖（Canu EOL）。

## 11. 相关文档

* `references/canu.md`（Canu OLC 源码分析 + §8 设计意图 + §8.5 实现后理解回写）
* `references/celera.md`（Celera 8.3rc2 源码分析 + 对照，§9 已按实现更新）
* `references/skesa.md` §7.1 + `references/metaMDBG.md` §9（v1 素材来源：
  fork 过滤/丰度过滤/桥接 reads repeat breaking）
* `references/metaMDBG.md` §4.1.1（metaMDBG 的**首要**借鉴点是 multi-k′
  迭代中的跨接验证——用当前 k 的 solid k′-min-mer 选择 unitig 连接/bubble
  分支，与 OLC 无关；OLC 只是其借鉴面之一）
* `design/asm-multik.md`（该机制在 anchr 的落地设计：`asm multik` 迭代
  驱动器，2026-08-14）
* `pgr: design/kmer.md` §11/§12（k 范围、FastK 字节键唯一表示，留 pgr）
* `design/asm-unitig.md` §8（`asm unitig` 语义与 L: 边）
* `todo.md` §3（多 k unitig OLC 挂账项，本项目承接）

## 12. 真实数据验证：Lambda（2026-08-12）

数据：`tests/bbtools/Lambda/R1.fq.gz` + `R2.fq.gz`（SRR5042715，108 bp ×
20k 对 = 40×，Illumina PE）；参考 `BBTools-40.01/resources/lambda.fa.gz`
（NC_001416.1，48,502 bp）；基线 golden `tadpole_contigs31.fasta.gz`
（BBTools 同 reads 组装，48,214 bp / N50 1199 / 最长 4258）。

### 12.1 抓出的 bug（已修复）

**左向延伸坐标下溢**：layout 坐标回填对 prepend 的首个 step 用占位
`q_end=0`，`prev_end − overlap_len` 下溢 panic（真实数据触发，合成数据
只测了右向延伸）。修复：首步坐标从自身长度算起；overlap > 前步末端改为
友好报错（零 panic 策略）。回归测试 `seed_extends_both_directions` /
`inconsistent_overlap_is_error`。

### 12.2 结果

| 实验 | 输入 | k | contigs | N50 | 最长 | 完美贴回参考 | 参考覆盖（正链） |
|---|---|---|---|---|---|---|---|
| A | 原始 reads 40× | 21,51,81 | 52 | 3409 | 19035 | 40/52 | 65.5% |
| B | 纠错 reads 9×（merge.ecco） | 21,51,81 | 202 | 459 | 2129 | 187/202 (92.6%) | 86.0% |
| C | 纠错 reads 9× | 21,31,41 | 248 | 469 | 2129 | 228/248 (91.9%) | 82.6% |
| D | 原始 reads 40× | 21,31,41 | 67 | 2233 | 8282 | 51/67 (76%) | 62.8% |

### 12.3 结论与教训

1. **长 contig ≠ 嵌合**：A 的最长 contig（19,035 bp）是**单个 k81 unitig**
   （cov=1.0），前 1708 bp 匹配参考（ref 29307–31015），随后跳出的序列
   **在 reads 中实锤存在**（100 bp 探针 fwd/rc 均命中）且两侧都是参考
   匹配区——是相对 NC_001416 的**菌株插入变异**（~1.3 kb at ref 31015），
   不是错装。教训：参考菌株与 reads 不同源时，"完美贴回"是错误判据，
   验证需 reads 侧证据（unitig 由 solid k-mer 建成本身即内部一致性证据）。
2. **多 k 冗余 2.4×**（A 总长 116 kb vs 基因组 48.5 kb）：不同 k 的 unitigs
   覆盖同一区域、contain 重叠被排除在延伸外 → 输出重复。v1 需消减
   （contain 去重，见 §13）。
3. **覆盖度 vs 纯度权衡**：40× 原始 reads 出长 contig 但变异区/重复区
   干扰贴回；9× 纠错 reads 纯度更高（92.6% 贴回）但碎片化（N50 459）。
   宏基因组真实数据的推荐路径待定（等数据）。
4. **k 选择**：108 bp reads 下 21/51/81 优于 21/31/41（原始 40×）——
   大 k 特异性在重复区更稳；k 应随读长自适应（设计默认 21/51/81 面向
   更长 reads）。
5. **合成数据验证的盲区**：合成测试全部是右向延伸链，漏了左向——真实
   数据验证的价值再次体现。

### 12.4 reads 回贴验证（2026-08-12，确认全长 contig 正确）

预过滤后 OLC 把 Lambda 拼成**单条 48,387 bp contig**（≈ 48,502 参考）。
reads 回贴验证（`anchr asm map`，完美匹配）：

* 40,000 reads 中 **34,697（86.7%）完美贴回** contig；对 NC_001416 参考
  只有 34,069（85.2%）——OLC contig 多捕获 628 条，正是参考缺失的
  变异区 reads；
* 覆盖剖面：平均深度 77.4、中位 78，**无 ≥50 bp 零深度缺口**、
  **无 ≥100 bp 低覆盖（<5）长段**——整条 contig 连续 reads 支持；
* 变异区（contig 1708–3000，相对 NC_001416 的插入）平均深度 64.5
  （min 21）——变异序列是 reads 实锤，不是错装。

结论：unitig 级预过滤 + 单条全长 contig 的正确性由 reads 侧证据确认；
后续 OLC 验证（宏基因组）沿用此口径（参考不匹配时看 reads 回贴而非
贴回率）。

## 13. v1 待办（真实数据驱动）

* **多 k 冗余消减：已完成（2026-08-12）**——两级：
  * 输出级：consensus 丢弃完全包含于更长 contig 的 contig（含 rc）；
  * unitig 级（`filter_contained`，布局前）：剔除被更长 unitig 完全包含
    的 unitig。Lambda 实测：unitigs 90→22（-76%），overlaps 386→50，
    layout 从 16 条碎片**合并为 1 条全长基因组 contig**（48,387 bp ≈
    48,502），16 条旧 contig 全部是它的子串（内容零丢失）。注意：
    过滤会改变 greedy 路径选择（这正是目的——多 k 冗余曾打断互惠链），
    "内容保留"而非"布局不变"。
* **repeat breaking 覆盖度证据**：桥接 reads 回放（`anchr asm map` +
  `anchr sam to-rg` + `pgr rg coverage`），阈值参考 SKESA fraction /
  metaMDBG 语义；
  需 reads 侧验证口径（参考菌株不匹配时不能只用贴回率）。
* 真实宏基因组数据验证 + 调参。

## 14. 设计评审与改造（2026-08-15，G37 multik 合并实验驱动）

### 14.1 输入约定的来源与"一代 vs 短读"辨析

**来源**：`asm olc` 驱动器"输入 reads"不是抄一代测序范式——2026-08-12
用户裁定明确"**不做 reads 级 OLC**：reads 先用 DBG（`asm unitig`）压缩成
unitigs"（§1），参考文档也写明"不是对 reads 做 OLC，把不同 k 的 unitigs
当伪 reads，在 unitig 层做 OLC"（`references/canu.md` §8.1）。
"输入 reads"只是 S0 的入口（内部生成 unitigs），**ovlp/layout/cns 操作
对象始终是 unitigs**。

**一代 vs 短读**：

* 直接对 reads 做全对 overlap 的 OLC 只有一代（Sanger：长 ~1kb、几万条）
  可行；短读（Illumina：百万级、150 bp、有错误）全对 overlap 不可行——
  这正是裁定"不做 reads 级 OLC"的原因（§2 数据量论证）；
* 一代也不是"reads 直接 layout"：Celera 原版（Sanger/short-read 时代）
  有 AS_BOG unitigger（Best Overlap Graph，先建 unitig 再布局，
  `references/celera.md` §6.1），只是它的 unitig 基于 overlap 图
  （mate 驱动）而非 DBG；
* 结论：unitig 层 OLC 是**短读/宏基因组导向**的设计，不是一代遗留。

### 14.2 G37 multik 合并暴露的三个问题

用 `anchr asm olc` 合并 23 组 multik 输出（G37，`benchmarks/
multik-allgroups.md`）后逐一定位：

1. **`asm olc` 驱动器对"已是 unitigs"的输入多做 S0 二次组装**：长 unitig
   在 k51 图里被多组序列"污染"（共享 51-mer 处分支）而切断（最长
   178,776 → 56,787）——用错命令（应走独立管道或新 `--unitigs` 模式）；
2. **撞名**：多组输出 cat 成单文件时 `unitig_<id>` 重名（`unitig_186`
   有 7 份），ovlp/layout/cns 按名字取序列会串——必须多文件分别传入
   （tag = 文件 stem，§5 设计如此，但 cat 是自然误区）；
3. **独立管道缺 `filter_contained`**：`filter_contained`（unitig 级去
   contain）只埋在驱动器（`olc.rs:164`），独立 `ovlp` 命令不调用 → 独立
   管道输出 87% contain（`ov:A:C`）的 unitigs 全部单步输出 → dup 2.659；
   `consensus::dedup_contained` 只丢"完全包含"，组间"部分重叠"（同区域
   不同边界）不去。

### 14.3 改造方案（适合现代流程：reads/multik 输出 → unitigs → OLC 合并）

1. **`asm olc --unitigs` 输入模式**：跳过 S0，输入直接是 unitigs（多文件，
   tag = stem 防撞名），内部走 filter_contained → ovlp → layout → cns——
   一条命令覆盖"多组 multik 输出合并"，避免误用驱动器和独立管道缺
   filter_contained 的问题；
2. **`asm ovlp --filter-contained`**（或 layout 前自动）：把 §13 的
   unitig 级冗余消减暴露给独立管道，使"用户自己跑 ovlp/layout/cns"也
   享受与驱动器一致的语义；
3. **cns 部分重叠去冗余（v1 列投票方向）**：`dedup_contained` 之外，
   对 overlap 链上近似包含（边界/内容微差）的 contig 做对齐投票合并，
   AS_CNS `BaseCallMajority` 语义（§S3 预留）——去冗余同时保留覆盖证据
   （G37：输入 contained 去重 dup 1.201 但 GF -0.55pp，两全需此步）；
4. 帮助文本/文档强化：`asm olc` 输入 reads、`--unitigs` 输入 unitigs、
   独立管道供高级用法。

**改造后的目标流程**：

```text
reads ──multik/unitig──> 多组 unitigs ──asm olc --unitigs──>
  filter_contained → ovlp → layout → cns（含部分重叠去冗余）
  ──> contigs（0 mis / 0 N / dup≈1 / 覆盖互补）
```

实施顺序：先 14.3-1（--unitigs 模式，复用现有 libs，改动最小、收益最大：
G37 合并 dup 2.659 → 期待 ~1.2 且不损失 GF），再 14.3-2，最后 14.3-3。

### 14.4 改造实施状态（2026-08-15）

**14.3-1 已实现**：`anchr asm olc --unitigs`（跳过 S0，输入直接是
unitigs/contigs，多文件 tag = stem 防撞名；内部 filter_contained → ovlp →
layout → cns）。集成测试 `command_asm_olc_unitigs_merges_files`（dovetail
合并 + contain 过滤 + 多文件 tag）。

**14.3-2 已实现**：`anchr asm layout --filter-contained`（独立管道在布局前
调用与驱动器相同的 `filter_contained`，unitig 级去精确包含）。集成测试
`command_asm_layout_filter_contained`。G37 验证：独立管道
ovlp → layout --filter-contained → cns 与 `--unitigs` 结果完全一致
（N50 39,496 / dup 2.050 / GF 96.54%）。

**14.3-3 已实现（`--dedup-ratio`）**：cns 输出后做 **contig 级近似包含
去冗余**——`consensus_with_ratio` + `coverage`（seed 锚定 + 允许 ~1%
错配的扩展，对齐 `anchr contained --idt 0.99` 语义）。`asm cns
--dedup-ratio <f>`（默认 1.0 = 精确子串，向后兼容）；`asm olc --unitigs`
默认 0.99（现代流程自动去重）。单元测试
`dedups_approximate_contained_contigs`（边界差异的近似重复：ratio 1.0
保留两条、0.95 合并为长条）。

**G37 最终结果**（23 组 multik 输出，`--unitigs` 一条命令，dedup 0.99）：

| 指标 | 改造前（独立管道） | `--unitigs` v1 | `--unitigs` 最终 |
|---|---:|---:|---:|
| # contigs | 65 | 48 | **21** |
| N50 | 39,104 | 39,496 | **54,964** |
| Largest | 179,610 | 179,610 | 179,610 |
| Duplication ratio | 2.659 | 2.050 | **1.094** |
| Genome fraction (%) | 96.44 | 96.54 | **96.54** |
| # misassemblies | 0 | 1 | 1* |
| # N's / 100 kbp | 0 | 0 | 0 |

\* 唯一 mis 是 560 bp、cov 1.4 的低覆盖碎片（filter_contained 后 layout
新拼接的 relocation，`--min-contig-len 1000` 可滤除）——不是合并引入的
系统性问题。

**改造完成度**：14.3-1/2/3 全部落地，`reads → multik（多组）→ asm olc
--unitigs → contigs` 的现代流程成立：0 N、0 大 mis、GF 96.54%、N50 54.9K、
dup 1.094（接近单组 MR 的 1.000）。剩余：560 bp 碎片 mis 的覆盖度门槛
（可选）、真实宏基因组/长读验证（todo）。

### 14.5 `anchr asm anchor`（2026-08-15，用户建议：reads mapping → anchors → OLC）

用户建议："像老流程一样，用真实的 reads mapping 一遍得到 anchors，再对
anchors 做 OLC"——老流程 anchors = bbwrap 回贴 + [lower, upper] 覆盖区间
过滤（`references/anchr-legacy-pipeline.md` §2.4）。已实现为
**`anchr asm anchor`**：`asm map`（完美匹配，bbmap perfectmode 等价）→
逐碱基覆盖度（差分数组 sweep line）→ 老流程公式
`lower = max(mincov, (median − mscale×MAD)/lscale)`、
`upper = (median + mscale×MAD)×uscale` → 连续覆盖区 = anchors。

* 库：`libs/olc/anchor.rs`（AnchorOptions / coverage_from_alignments /
  anchor_thresholds / anchor_regions / extract_anchors），单元测试 2 个；
* 命令：`anchr asm anchor <unitigs.fa> <reads...> -o anchors.fa`，参数
  `--mincov 5 --mscale 3 --lscale 3 --uscale 2 --min-anchor-len 500
  --kmer 31 --parallel`；集成测试 `tests/cli_asm_anchor.rs`；
* **注意**：reads 必须用**产生该组 unitigs 的同一覆盖子集**（不是全量
  reads——全量会得到 ~180× 覆盖，median 失真）。

**G37 验证**（与 python 概念实验逐位一致）：

| 方案 | N50 | dup | GF (%) | mis | mm/100k |
|---|---:|---:|---:|---:|---:|
| MRX40P001 单组 unitigs | 54,888 | 1.000 | 96.11 | 0 | 30.7 |
| MRX40P001 单组 anchors | 54,841 | 1.000 | 95.66 | 0 | **27.21** |
| 7 组 MR `--unitigs`（无 anchors） | 54,964 | 1.618 | 96.37 | 0 | 32.95 |
| **7 组 MR anchors → `asm olc --unitigs`** | 54,858 | **1.002** | 96.04 | **0** | **28.12** |

**现代流程最终形态**：

```text
reads（每覆盖度子集）
  → asm multik → 每组 unitigs
  → asm anchor（同子集 reads 回贴 + 覆盖过滤）→ 每组可靠 anchors
  → 所有组 anchors（cat，名字唯一）→ asm olc --unitigs → contigs
```

覆盖过滤同时解决 `--unitigs` 的三大遗留：dup（upper 排除重复区多版本）、
低覆盖碎片 mis（lower 排除）、mm（过滤掉高/低覆盖区错配）。代价：GF
96.37 → 96.04（过滤掉的问题区本身是低质量区）。

### 14.6 时间分析（2026-08-15，用户确认可接受）

`anchr asm anchor` 实测（G37，MR 组，release，`-p 16`）：

| 阶段 | 单组耗时 | 占比 |
|---|---:|---:|
| 总时间 | 0.24s | — |
| reads 完美回贴（`asm map`） | ~0.14s | 58% |
| SAM 写盘+读回+解析 | ~0.05s | 21% |
| 覆盖度+阈值+区间 | ~0.03s | 12% |

7 组 MR 串行 1.6s（每组 0.24–0.28s；80× 组 230K reads 只比 40× 组慢
17%——mapping 的并行 verify 扩展性良好）。**现代流程的时间大头是
`asm multik`**（单组 MR 40× 2.35s，占 ~65%），anchor 占 ~7%、OLC 合并
~28%；23 组全跑（8 路并行）multik ~30s、anchor ~2s、OLC ~3s。用户
2026-08-15 确认"时间还可以"（此前觉得偏长是错觉）。

**保留的优化点**（非当前瓶颈，宏基因组数据时再评估）：SAM 中间文件
内存化（map libs `map_read` 已算对齐、只写盘读盘，宏基因组 SAM GB 级时
收益显著）；multik 性能（计数复用 / `remove_unsupported` 查表化 /
轮数裁剪，见 `benchmarks/multik-complexity.md` 待办）。

**G37 实测**（23 组 multik 输出，`--unitigs` 一条命令）：

| 指标 | `--unitigs` | 独立管道（多文件） | 输入 contained + 管道 |
|---|---:|---:|---:|
| N50 | 39,496 | 39,104 | 39,189 |
| Largest | 179,610 | 179,610 | 179,610 |
| # misassemblies | **1** | 0 | 0 |
| Genome fraction (%) | **96.54** | 96.44 | 95.89 |
| Duplication ratio | 2.050 | 2.659 | **1.201** |
| # mismatches / 100 kbp | 33.54 | 31.46 | 33.69 |

**已知问题（1 mis）**：`contig_48`（560 bp、cov 1.4）是 filter_contained
后 layout 新拼接的低覆盖碎片（relocation）——filter_contained 改变 greedy
路径选择（§13 已注明是设计意图），低覆盖碎片被链错。独立管道无此问题
（不含该序列）。影响小（560 bp），彻底解决需 layout 的覆盖度证据（14.3-3
方向）。

**dup 2.050 仍未到 1.2**：filter_contained 只去**精确包含**（ov:A:C），
多组 unitigs 的冗余主要是**近似包含/部分重叠**（边界/内容微差）——
`anchr contained --idt/--ratio`（近似判据）能到 1.201 但损失覆盖
（GF -0.65pp）。**两全方案 = 14.3-3（cns 部分重叠去冗余/列投票）**，
待实施。

## 15. 现代组装流程总结（2026-08-15 会话收官）

> 2026-08-17 并入自 `asm-olc-modern-flow.md`。原文是长会话
> （SKESA/multik/OLC 改造/anchors/替代盘点）的完整总结，供后续会话
> 快速恢复上下文。快速交接看 §15.1 现代流程与 §15.2 的用户裁定
> （老流程理解修正）；机制细节见 §14、
> `references/anchr-legacy-pipeline.md` §5、
> `benchmarks/multik-allgroups.md`。

### 15.1 现代流程最终形态

```text
reads（每覆盖度子集）
  → fq 家族：clean（trim）/ clump（dedupe）/ merge --ecphase（merge+纠错）
             / s-filter（quorum 替代：丢弃坏 reads、保留原样）
  → asm multik（每组 unitigs，metaMDBG 式 multi-k 迭代）
  → asm anchor（同子集 reads 完美回贴 + [lower, upper] 覆盖过滤 → 可靠 anchors）
  → asm olc --unitigs（所有组 anchors 合并 → contigs）
  → quast（外部，最终质量确认）
```

目标：**无 N 染色体**（无 gap、无嵌合、完整覆盖）。正确性硬指标由 multik
（0 mis / 0 N）+ anchor 覆盖过滤 + OLC 合并共同保证。

### 15.2 会话完成的工作

#### SKESA 借鉴（multik）

* 读 SKESA-master（C++）+ skesa-rs-main（Rust 移植），核对笔记并补充
  §7.2（与 multik 的多 k 迭代对比）；
* 落地"严格链唯一性"：`merge_chains`/`recompact_graph` 从先到先得占用
  改为**严格两端唯一**（SKESA "predecessor == 1"），含对称 link 去重修正
  （`compute_links` 双向 link 导致度数翻倍）。
* G37：misassemblies 0、mm 27.7/100kbp 历史最佳、N50 24.4K（-8% 宁断勿嵌合）。

#### 老流程理解修正（用户多次纠正，重要！）

* 主干 = **reads 拆多覆盖度部分 → 各部分组装 → 取可靠 anchors → OLC 式合并**，
  不是一条线；
* **合并（contained/orient/merge）才是经典 OLC**（`anchr merge` = overlap →
  有向图 → 拓扑排序 → 拼接）；**fill/glue 不是主 OLC**（用 2GS 长 contigs
  补 gap 的辅助）；
* G37 文档记录 **40×/80×**（`results/model.md` statQuorum 表）；30×/60× 是
  用户口头经验（初验：单跑都达标）；
* **pe.cor.fa 未纠错**：quorum 丢弃被修正 reads（`:sub:`/`trunc` 标记），
  pe.cor 是未修正原始序列；
* **quorum 替代 = `fq s-filter`**（检查的正是 quorum 的信号：anchor/
  truncation/substitution + Poisson 碰撞，保留原样或丢弃）；
* **fastk/GeneScope 替代 = `pgr kmer gsize --model`**（GenomeScope 拟合，
  `libs/kmer/genomescope.rs`）；
* **spades/megahit = 可选参考**（保留用于对照，不实现不替代）；
* 唯一保留外部组件：**quast**（用户要求最终质量确认）；bcalm 作为
  对照 unitigger 随 `anchr template --unitigger` 恢复（2026-08-15，见
  `todo.md`），用于与现代 multik 做 unitig/最终 N50 的 A/B 对照。

#### multik 验证与基准

* 全分组复核：老流程 23 组全部 0 mis / 0 N / dup≤1.001，MR 组全面更优
  （N50 34-55K）、X40 优于 X80、全量 reads 反而最差（`multik-allgroups.md`）；
* 覆盖度实验：30×/60× 单跑质量达标（60× mismatch 最优 25.9）
  （`multik-cov.md`）。

#### OLC 改造（现代流程的命令支撑）

* **`asm olc --unitigs`**：跳过 S0（不二次组装切断长 unitig）、多文件 tag
  防撞名、内部 filter_contained；
* **`asm layout --filter-contained`**：独立管道与驱动器语义对齐；
* **`--dedup-ratio`**（cns）：contig 级近似包含去冗余（允许 ~1% 错配，
  对齐 `anchr contained --idt 0.99`）；
* G37 23 组合并：dup 2.659 → 1.094、GF 96.54%、N50 54.9K、0 N、0 大 mis。

#### `anchr asm anchor`（用户建议：reads mapping → anchors → OLC）

* 实现：`asm map`（完美回贴）+ 逐碱基覆盖度 + 老流程公式
  （`lower = max(mincov, (median−mscale·MAD)/lscale)`、
  `upper = (median+mscale·MAD)·uscale`）→ 连续覆盖区 = anchors；
* 库 `libs/olc/anchor.rs`、命令 `cmd/asm/anchor.rs`、测试
  `tests/cli_asm_anchor.rs`；
* G37 7 组 MR anchors → OLC：**dup 1.002、mm 28.12、0 mis、GF 96.04%**
  ——覆盖过滤一次解决 dup（upper 排重复区）、560bp 碎片 mis（lower 排
  低覆盖）、mm（滤错配区）；
* **注意**：reads 必须用产生该组 unitigs 的同一覆盖子集（全量 reads 会
  得到 ~180× 覆盖，median 失真）。

#### 时间分析（用户确认可接受）

* anchor 单组 0.24s（map 58%、SAM I/O 21%、覆盖 12%），7 组串行 1.6s；
* 现代流程大头是 **multik**（单组 2.35s，占 ~65%）；
* 优化点保留：SAM 内存化（宏基因组时）、multik 性能（multik-complexity.md）。

### 15.3 当前状态

* 测试：**29 套件全绿、0 失败**；fmt/clippy 干净；
* 工作区：**17 个文件未提交**（14 修改 + 3 新增：`anchor.rs`×2 + 测试），
  `.git` 只读需用户本机 commit；
* 文档：asm-olc.md §14（改造+时间分析）、anchr-legacy-pipeline.md §5
  （替代对照）、multik 基准、todo 全部同步。

* 2026-08-15 追加：`--unitigger` 恢复 bcalm（multik 默认；bcalm 每 k
  unitigs + `asm olc --unitigs` 跨 k 合并，G37 N50 31199 与 legacy 一致）；
  MG1655 对照显示 multik unitigs 比 bcalm 短 2-2.5×（MRX40：21.2K/1455
  条 vs 53.6K/158 条），最终 N50 28K vs legacy 63-97K —— 差距来源待
  bcalm 链端到端确认。

### 15.4 下一步（详见 todo）

1. 真实宏基因组/长读数据验证（决定拆分合并路线的最终定位）；
2. anchor 补洞逻辑（老流程 fill）；
3. 560bp 碎片 mis 的覆盖度门槛（`--min-contig-len 1000` 可滤，可选）；
4. SAM 内存化（宏基因组时再做）；
5. multik 性能优化（计数复用 / remove_unsupported 查表化 / 轮数裁剪）。

### 15.5 `--cross-validate` 跨组嵌合投票（2026-08-18）

背景（G37 MRX40P002 归因，详见 `results/model_org.md` 2026-08-17/18
两节）：单组内 k≥101 的图会连上某些"低丰度菌株 A-B 相邻"结构
（junction 桥接 reads 9 条，真实存在），quast 按主参考报 relocation；
而跨组 `olc --unitigs` 的 `filter_contained` 会把被嵌合 contig
前缀整条包含的其他组正确 contigs 删掉，嵌合反而存活到最后。

机制（`libs/olc/overlap.rs::drop_cross_chimeras`，CLI
`olc --unitigs --cross-validate`，默认关；跨组模板
`7_merge_anchors.tera.sh` 已开）：

* contig 两端（flank = min_overlap）各自被 ≥2 个其他**文件**的
  contigs 覆盖（按文件 stem 计票），且中部 junction 窗口
  （span = min_overlap/2）无任何其他文件 contig 横跨 → 删除该
  contig；序列由其他组的分开 contigs 完整提供，GF 不降（G37 实测
  反升 98.774→98.797）；
* 横跨判定先**按来源 contig 合并对齐区间**：精确重叠链在每个错配
  处断裂，组0 的 236 kb 等价 contig 表现为两段
  [0,48851)+[48757,236423)，单段检查会误判（MG1655 首版实测误删，
  N50 118731→113444，合并后恢复）；
* 必须在 `filter_contained` **之前**跑（否则正确 contigs 已被
  前缀包含删除）；
* 单文件输入时无跨文件 covers，自然无操作。

语义边界：这是跨样本多数投票（其余组分开 ↔ 本组连上），只能压制
"单组私有的连接"；各组一致的真连接（有横跨者）保留。精确重叠的
固有盲区（两条等价 contig 端部错位、对齐不含任何端点时检测不到）
由 min_groups≥2 门槛兜底，实测 MG1655/G37 双门禁 0 mis。
