# anchr 近期待办

> 依据 `project-understanding.md`、`notes/design/fq-asm-migrate.md` 与
> 各审计/基准记录整理。已完成条目只留一行结论，细节见链接文档。
> 按类型组织（已完成 / 待实现 / 挂账待决 / 待验证 / 低风险审计 /
> 技术债 / 明确不做），不按会话轮次。

## 0. 会话交接（2026-08-15，长会话收官）

**主题**：现代 OLC 组装流程定型——`fq 家族 → asm multik → asm anchor →
asm olc --unitigs → quast`。完整总结见
[design/asm-olc-modern-flow.md](design/asm-olc-modern-flow.md)。

**本次会话完成**：
* multik 严格链唯一性（SKESA 借鉴）+ G37 全分组复核（23/23 组 0 mis/0 N）；
* OLC 改造：`--unitigs`（跳过 S0、防撞名）、`layout --filter-contained`、
  `--dedup-ratio`（近似去冗余，dup 2.659 → 1.094）；
* **`anchr asm anchor`**（reads 回贴 + 覆盖过滤 → anchors；7 组 MR →
  OLC = dup 1.002、mm 28.12、0 mis）；
* 老流程替代盘点 + 时间分析（用户确认可接受）。

**用户关键裁定（勿再犯）**：
* 老流程 = 多覆盖度拆分 → 各部分组装 → **取可靠 anchors → OLC 合并**；
  fill/glue 不是主 OLC（合并才是）；
* G37 文档记录 **40×/80×**（30×/60× 是口头经验）；pe.cor 未纠错；
* quorum 替代 = **`fq s-filter`**；fastk/GeneScope 替代 =
  **`pgr kmer gsize --model`**；spades/megahit = **可选参考**；
* 唯一保留外部组件：**quast**。

**工作区状态**：17 个文件未提交（14 修改 + 3 新增，`.git` 只读需本机
commit）；测试 29 套件全绿、fmt/clippy 干净。

**下一步优先级**：真实宏基因组/长读验证 > anchor 补洞 > 560bp 碎片门槛
（可选）> SAM 内存化（宏基因组时）> multik 性能。

## 1. 已完成（一行结论，细节见链接）

- fq/asm 迁移阶段 1-4：25 命令 + 业务 libs + 21 测试 + golden 数据
  （[fq-asm-migrate.md](design/fq-asm-migrate.md)）；
- pgr 依赖规范化：rev 锁定 + 本地 patch（`Cargo.toml` +
  `.cargo/config.toml`，后者 gitignore）；
- `bio`/`intspan` 移除，covered 改 `pgr::libs::runlist`
  （[project-understanding.md](project-understanding.md) §2.3）；
- 基础设施：AGENTS.md / rust-toolchain / CI / docs 目录统一
  （[AGENTS.md](../../AGENTS.md)）；
- 项目理解与索引：`notes/project-understanding.md`、`notes/todo.md`；
- 基准：`fq_assemble`/`asm_map`/`covered` benchmarks +
  `bbtools-vs-anchr.md`（[benchmarks](benchmarks/bbtools-vs-anchr.md)）；
- QC M1-M4 完成（`anchr fq qc`，FastQC/Falco 底层数值一致，见
  [qc.md](design/qc.md)）；
- cuttlefish 2.2.0 调研 + 基准：不整体借鉴，DFA 状态分类思路已吸收
  （[references/cuttlefish.md](references/cuttlefish.md)、
  [asm-assemble.md](design/asm-assemble.md) §11/§12.2）；
- `asm unitig` 效率攻坚：DFA 默认 walk、流式计数、单池复用、状态表
  内嵌计数、`--parallel` 默认 `min(逻辑核/2, 8)`；输出与旧引擎逐字节
  一致（[asm-assemble.md](design/asm-assemble.md) §12.2）；
- supermer 转正为 **FASTA 默认计数**（pgr `b31af11`：自适应 m + slices
  API；FASTQ 自动回退 direct 保 `min_prob`；`--no-supermer` 强制），
  新默认 half(8) 基线 k31 1.36 s/597 MB、k99 2.19 s/1131 MB
  （[asm-assemble.md](design/asm-assemble.md) §12.3、
  [unitig-bench.md](benchmarks/unitig-bench.md)）；
- **`asm multik` v4（2026-08-14）**：multi-k 迭代组装落地——unitig 图
  跨轮验证（桥接 k-mer ≥2 + 内部 k-mer 嵌合清理）+ **渐进丰度过滤 +
  recompact**（主路径保护：只删分支/孤立节点，直链永远保留）+ **k 序列
  自适应**（`--kmer auto`，读长 N50 驱动）。
  **Lambda 20k：38 contigs / 最长 46467 / N50 46467 / 总长 49359**——最长
  ≈ 参考 48502 的 95.8%，reads 100% 覆盖、零缺口（短读数据下第一条近整条
  染色体级无 N contig；OLC 对照最长 19035）。合成长读：**完整覆盖（环状）
  → 单条 100,000 bp = 基因组 100%（k-mer 多重集一致）**；0.1% 错误 +
  2 kb 重复 ×2 → k-mer 覆盖 97.4%（重复区完整）。
- **`asm multik` v5-v7（2026-08-14）**：真实数据验证 + 防错连 +
  megahit 借鉴——G37 端到端（misassemblies 8→**0**、N50 26.6K、
  最长 52.8K）、被删分支回灌（megahit bubble 回灌 + metaMDBG unitig
  反馈）、tip/weak_link 清洗（megahit 算法族）、Lambda/20k 不回归、
  374 测试绿。设计 [asm-multik.md](design/asm-multik.md) §4.10-§4.12、
  [asm-multik.md](design/asm-multik.md) §9/§10、
  [megahit.md](references/megahit.md) §8.6；
- 双轨核对 22/22（历史）：`scripts/verify-migrate.sh` 只对"删除
  fq/asm/sam 之前"的 pgr 有效，新版 pgr 下全 FAIL 属预期（脚本头已
  注明），勿当 bug。

## 2. 待实现

- **`fq range` 的 BGZF `.gzi` 索引 CLI 化**：目前 BGZF 输入需要外部生成
  `.gzi`（pgr `libs/bgzf::build_gzi_index`），plain 文本输入自动建 `.loc`；
  若要在 anchr 侧支持 BGZF range，需在 pgr 补 CLI 或封装
  （[fq-range.md](design/fq-range.md)）；
- **模板链（`templates/*.tera.sh`）端到端替换验证**：trim.tera.sh 等模板
  仍引用 BBTools/dazzler 工具，逐步切到 anchr fq/asm 命令并按
  `fq-trim-replace.md`/`fq-merge-replace.md` golden 核对；
- **`fq norm` 精确 vs 近似定稿**（pgr 移交，`fq-trim-replace.md` §4.8
  未定）：anchr 走精确表 + 外部桶；bbnorm `bits=16` 近似表结果依赖
  `-Xmx`。差异 = 定义差异不是 bug，需在文档中定稿并记录边界差异；
- **`dep`/`ena`/`template` 命令的外部工具版本核对**：依赖
  dazzler/hnsm 系统工具，CI/容器环境预装清单待整理；
- **pgr 并行读 gz（可选）**：asm-assemble.md §12.3 第 5 条，pgr 侧
  决定是否做（风险>收益暂缓）。

## 3. 挂账 / 待决

- **pgr rev 更新流程**：pgr 基础 API 变更后需要手动 bump rev + 重跑
  `verify-migrate.sh`；本地 patch 可能掩盖发布版差异，发布前需确认；
- **双轨遗留**：pgr 侧删除后，`notes/design/fq-asm-migrate.md` 的 anchr
  副本需随 pgr 更新（阶段 4 完成标注），`project-understanding.md` §6.3
  的"待补全"条目逐项销账；
- **audit 文档增量**：`notes/audit/audit-fq.md`/`audit-asm.md` 是迁移时
  的审计快照；后续 anchr 侧对 fq/asm 的修改应更新审计记录而非 pgr 侧；
- **pgr supermer 质量门控**：当前 FASTQ 自动回退 direct；若 pgr 给
  supermer 补质量门控，可去掉回退并统一计数路径
  （[asm-assemble.md](design/asm-assemble.md) §12.1）。

## 4. 待验证 / 等数据或场景到位

- **`asm multik` 端到端**：Lambda 短读（最长 46467 ≈ 参考 95.8%，零缺口）、
  合成长读（完整环状覆盖 → **单条 100% 基因组**；1 Mb 基准 9.4 s / 816 MB）、
  **G37 真实基因组（2026-08-14）**：misassemblies **0**、N50 26.6K、最长
  52.8K（回灌后）、Genome fraction 95.33%（设计 §4.10-§4.12）已过；
  **tip/weak_link 的宏基因组价值**（多菌株弱连接）与真实宏基因组/长读
  无 N 判据（覆盖完整 + 无 gap + 无嵌合）待数据验证；
  `--parallel` 扩展性复测与 `--min-count-extend` 阈值调参（设计 §7）；
- **`asm multik` 迭代性能**（2026-08-14 基准，`benchmarks/multik-complexity.md`）：
  图结构递减 ✓（unitigs 1345→396、edges 346→12）但耗时递增
  （remove_unsupported O(总长×k) 为瓶颈，迭代 3.25 s > 直接大 k 1.07 s）；
  优化：计数复用（minimizer 流）、remove_unsupported 查表化、轮数裁剪；
- **`asm anchor` SAM 内存化（可选，2026-08-15 时间分析）**：map libs
  `map_read` 已算出对齐（cid/pos/rc）只写 SAM 再读回（单组 G37 占 21%、
  宏基因组 SAM GB 级时收益显著）；当前时间可接受（用户确认），等
  宏基因组数据再评估；
- **`asm multik` 防 misassembly（层次 3 落地，2026-08-14）**：
  `bridge_filter`（unitig 间探针）+ `split_by_bridge`（unitig 内部 60-mer
  窗口切分）——**G37 misassemblies 8 → 0**、N50 24.5K → 26.6K、Lambda/
  20k 不回归（设计 `asm-multik.md` §9）。待办：探针长度/阈值
  调参、<500 bp dropped 碎片输出策略；
- **multik 多覆盖度拆分（2026-08-15 初验，`benchmarks/multik-cov.md`）**：
  老流程"reads 拆成多个覆盖度部分（G37 文档记录 40×/80×；用户口头提过
  30×/60× 合适）→ 各部分组装 → OLC 拼装"（`references/anchr-legacy-pipeline.md`
  §2.2/§2.7）。初验：30×/60× 单跑 quality 达标（0 mis/~96%，60× mismatch
  最优 25.9）；**30×+60× contained 合并 N50=40× 水平但 duplication 1.124**
  ——多部分合并的去冗余机制（anchors upper / orient+merge / glue OLC）
  待设计；**全分组复核已完成**（`benchmarks/multik-allgroups.md`：老流程
  全部 23 组 40×/80× × Q/L × P + MR，23/23 组 0 mis、0 N、dup≤1.001，
  MR 组 N50 34-55K 全面优于非 MR 17-24K）；**合并实验**：contained 合并
  23 组输出 N50 37.5K 但 dup 1.827；1× 重新压实引入 1 mis（无 reads 证据）；
  全量 reads 反而最差（N50 14.5K）；单组 MR（54.9K）已接近老流程合并
  （55.0K）；**正确合并方式（2026-08-15 定位）**：`asm olc` 驱动器的
  设计输入是 reads（内部 S0 做 unitig 重组装），喂 multik 输出会切断长
  序列 + cat 单文件会撞名；应直接用**独立管道 `asm ovlp → layout → cns`
  （多文件分别传入）**——23 组合并：最长 179,610 保持、0 mis、GF 96.44%
  全组最优、mm 31.46 全组最优（dup 2.659 未去冗余，contained 后 1.899）
  （`benchmarks/multik-allgroups.md` "为什么之前变短"节）——外部 unitigs
  反馈（--unitigs）降级为可选优化；
- **megahit 借鉴待续**（`megahit.md` §8 差距清单）：low depth 局部窗口
  判定（宏基因组覆盖不均）、本地组装（contig 端点 reads 延伸）——
  有真实宏基因组数据后再定优先级；
- **SKESA 借鉴落地**（`skesa.md` §7.2，2026-08-14 重读本体 + Rust
  移植）："前驱恰好 1"不变量已吸收为 **merge/recompact 严格两端唯一**
  （`multik.rs`，含对称 link 去重修正）；read 清理确认不做（与
  `remove_unsupported` 机制冲突，详 §7.2）；SNP 簇变体表示暂缓（等菌株/
  宏基因组数据）。G37 回归：misassemblies 0、N50 24.4K（-8%）、
  mismatches 27.7/100kbp 历史最佳；
- **大规模真实数据**：Lambda 20k/40k reads 之外，用真实染色体数据跑
  `fq → asm → map → template` 全链，核对统计（覆盖量/unitig 数/PSL 行数）；
- **多线程与内存**：`fq norm` 外部 hash-bucket 路径（`--mem`）、
  `fq trim-adapter --parallel` 在 50 万-pair 合成数据上的扩展性复测；
- **`asm olc` 参数验证**（pgr 移交）：overlap 是否允许少量错配、不同 k
  unitig 冗余去重、repeat breaking 覆盖度证据阈值（Canu 6/15 单元化
  版本）、列投票 consensus（v1）——见 `notes/references/canu.md` §8；
- **`asm cns` contain 吸收（v1，2026-08-15 multik 合并实验暴露）**：
  多组 unitigs 合并时 ovlp 87% 是 contain（ov:A:C），v0 cns 只缝合 layout
  的 dovetail、contain 序列各自单步输出 → dup 2.659；输入侧 contained 可
  降到 1.201 但损失覆盖（GF -0.55pp）。v1 应实现 contain 序列对齐主
  contig 的列投票/吸收（AS_CNS BaseCallMajority 语义，`design/asm-olc.md`
  §S3 已预留），既去冗余又保留覆盖证据；
- **`asm olc` 现代流程改造（2026-08-15 完成，`design/asm-olc.md` §14）**：
  `--unitigs` 输入模式（跳过 S0、多文件 tag 防撞名、filter_contained）、
  `asm layout --filter-contained`（独立管道对齐）、`--dedup-ratio`（contig
  级近似包含去冗余，允许 ~1% 错配）。G37 23 组 multik 合并最终：
  **dup 2.659 → 1.094、GF 96.44 → 96.54%、N50 39.1K → 54.9K、0 N、0 大
  mis**（唯一 mis 为 560 bp 低覆盖碎片，`--min-contig-len 1000` 可滤）。
  剩余：560 bp 碎片 mis 的覆盖度门槛（可选）、真实宏基因组/长读验证；
- **`anchr asm anchor`（2026-08-15 实现，`design/asm-olc.md` §14.5）**：
  reads 完美回贴（`asm map`）+ 覆盖度区间 [lower, upper] 过滤 → 可靠
  anchors（老流程 anchors 的现代命令）。G37 7 组 MR：
  anchors → `asm olc --unitigs` = **dup 1.002、mm 28.12、0 mis、GF 96.04%**
  （覆盖过滤同时解决 dup / 560bp 碎片 mis / mm 三大遗留）。现代流程：
  `reads(每子集) → multik → anchor → 各组 anchors → olc --unitigs`。
  待办：补洞逻辑（老流程 fill）、others 输出、真实宏基因组/长读验证；
- **OLC 宏基因组/长读真实数据验证**（pgr 移交）：`asm olc` 四命令在
  宏基因组/长读数据上的端到端验证；
- **fq range BGZF**：真实 BGZF 输入（含 `.gzi`）的端到端验证
  （测试已覆盖 plain + BGZF 小数据，见 `tests/cli_fq_range.rs`）；
- **gz/大输入回归**：默认 supermer 路径在 gz 输入下的全链回归 +
  更大数据集的峰值内存（当前基准为 plain 144 MB 级别）。

## 5. 低风险审计记录项（可顺手修）

- `scripts/verify-migrate.sh` 的 `asm_olc` 用例用 Lambda 数据（约 6 s），
  可考虑缩小输入加速日常核对；
- ~~`notes/benchmarks/` 目录索引~~：已补 `README.md`（2026-08-14，含
  `multik.md`）。

## 6. 技术债（有空再议）

- **`.cargo/config.toml` 的 patch 与发布形态差异**：本地走 `../pgr`，
  发布走 GitHub rev；若 pgr 未 push 的本地改动被 anchr 依赖，构建会
  在外部环境失败——需要约定"依赖 pgr 新 API 必须先 push 再 bump rev"；
- **外部工具链依赖**：流程命令强依赖 dazzler/hnsm（PATH），测试缺失时
  跳过可能掩盖回归（`tests/cli_pl.rs` 的 LAshow 模式）；
- **golden 数据体积**：`tests/bbtools/Lambda/` 约 13 MB，若仓库膨胀
  可考虑外部数据源或生成脚本（当前随测试提交）。

## 7. 明确不做（避免重复立项）

- 不重复基础层：格式 I/O、Phred、k-mer 表、PAF 解析一律走 pgr
  （`pgr::libs::*`）；
- 不搬 pgr 专属基准：`benches/` 33 个 .rs 与 `notes/benchmarks/` 12 篇
  均服务 pgr 剩余命令/基础层，留在 pgr；
- 不内置 pgr 的比对/索引/遮蔽命令（`align`/`pgi`/`pbit`/`rept` 等）；
- 不重新引入 `bio`/`intspan` 外部 crate（pgr 已内嵌等价实现）；
- 不做并行 walk（`--parallel-walk` 已撤：多阶段线程池叠加会卡死系统，
  教训已记入 [asm-assemble.md](design/asm-assemble.md) §12.2）。
