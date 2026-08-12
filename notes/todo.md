# anchr 近期待办

> 依据 `project-understanding.md`、`notes/design/fq-asm-migrate.md` 与
> 各审计/基准记录整理。已完成条目只留一行结论，细节见链接文档。
> 按类型组织（已完成 / 待实现 / 挂账待决 / 待验证 / 低风险审计 /
> 技术债 / 明确不做），不按会话轮次。

## 0. 会话交接（2026-08-13，fq/asm 迁移收尾）

> 会话交接材料，供下一次会话恢复上下文；读取后按用户指示清理。

**当前状态**：289 测试通过；`scripts/verify-migrate.sh` 22/22 双轨核对通过；
fq/asm 迁移阶段 1-4 完成（pgr 侧 `ace4ee8` 已删除对应代码/文档，零代码差异）。
工作树 2 处未提交：`Cargo.toml`（pgr rev → `ace4ee8`）、
`notes/design/fq-asm-migrate.md`（同步 pgr 的基准迁移说明）。

**最近提交**：`88fac59`（project-understanding + 基准引用修正）、
`a8edc18`（pgr dep 更新 + fq_assemble/asm_map benchmarks 迁入）、
`0abbd3b`（fq/asm 全套工具迁入）、`9060e73`（covered 改用 pgr runlist）、
`9f05c71`（bio/intspan 依赖移除）。

**本会话成果**：
- fq/asm 业务 + 命令 + 测试 + golden 全部迁入 anchr，依赖 pgr 基础层
  （git rev 锁定 + 本地 `.cargo/config.toml` patch）；
- `bio`/`intspan` 外部 crate 移除：序列 I/O、IntSpan、分层覆盖改用 pgr
  （covered 经基准验证从 vendored Coverage 迁到 `pgr::libs::runlist`，
  sweep 路径快 2-30×，`benches/covered_benchmark.rs`）；
- 基础设施对齐 pgr：AGENTS.md、rust-toolchain、CI（zigbuild）、
  `docs/`（原 `doc/`）、notes 目录、`notes/project-understanding.md`、
  `notes/todo.md`（本文）、`notes/benchmarks/bbtools-vs-anchr.md`；
- 双轨 golden 核对脚本 `scripts/verify-migrate.sh`（22 命令，含 stderr
  规范化）；核对中发现并修复：dispatch `unwrap` → `?`（Zero Panic）、
  env_logger 缺失（日志不可见）、核对脚本只比 stdout 的漏洞。

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
- 双轨核对 22/22：`scripts/verify-migrate.sh`（[verify-migrate.sh](../../scripts/verify-migrate.sh)）。

## 2. 待实现

- **`fq range` 的 BGZF `.gzi` 索引 CLI 化**：目前 BGZF 输入需要外部生成
  `.gzi`（pgr `libs/bgzf::build_gzi_index`），plain 文本输入自动建 `.loc`；
  若要在 anchr 侧支持 BGZF range，需在 pgr 补 CLI 或封装
  （[fq-index.md](design/fq-index.md)）；
- **模板链（`templates/*.tera.sh`）端到端替换验证**：trim.tera.sh 等模板
  仍引用 BBTools/dazzler 工具，逐步切到 anchr fq/asm 命令并按
  `anchr-trim-replace.md`/`anchr-merge-replace.md` golden 核对；
- **`fq norm` 精确 vs 近似定稿**（pgr 移交，`anchr-trim-replace.md` §4.8
  未定）：anchr 走精确表 + 外部桶；bbnorm `bits=16` 近似表结果依赖
  `-Xmx`。差异 = 定义差异不是 bug，需在文档中定稿并记录边界差异；
- **`dep`/`ena`/`template` 命令的外部工具版本核对**：依赖
  dazzler/hnsm 系统工具，CI/容器环境预装清单待整理。

## 3. 挂账 / 待决

- **pgr rev 更新流程**：pgr 基础 API 变更后需要手动 bump rev + 重跑
  `verify-migrate.sh`；本地 patch 可能掩盖发布版差异，发布前需确认；
- **双轨遗留**：pgr 侧删除后，`notes/design/fq-asm-migrate.md` 的 anchr
  副本需随 pgr 更新（阶段 4 完成标注），`project-understanding.md` §6.3
  的"待补全"条目逐项销账；
- **audit 文档增量**：`notes/audit/audit-fq.md`/`audit-asm.md` 是迁移时
  的审计快照；后续 anchr 侧对 fq/asm 的修改应更新审计记录而非 pgr 侧。

## 4. 待验证 / 等数据或场景到位

- **大规模真实数据**：Lambda 20k/40k reads 之外，用真实染色体数据跑
  `fq → asm → map → template` 全链，核对统计（覆盖量/unitig 数/PSL 行数）；
- **多线程与内存**：`fq norm` 外部 hash-bucket 路径（`--mem`）、
  `fq trim-adapter --parallel` 在 50 万-pair 合成数据上的扩展性复测；
- **`asm olc` 参数验证**（pgr 移交）：overlap 是否允许少量错配、不同 k
  unitig 冗余去重、repeat breaking 覆盖度证据阈值（Canu 6/15 单元化
  版本）、列投票 consensus（v1）——见 `notes/references/canu.md` §8；
- **OLC 宏基因组/长读真实数据验证**（pgr 移交）：`asm olc` 四命令在
  宏基因组/长读数据上的端到端验证；
- **fq range BGZF**：真实 BGZF 输入（含 `.gzi`）的端到端验证
  （测试已覆盖 plain + BGZF 小数据，见 `tests/cli_fq_range.rs`）。

## 5. 低风险审计记录项（可顺手修）

- 既有 8 处 warning：`Overlap` dead_code（libs/overlap.rs）、
  `tadpole` 未读字段（`error_extension_pincer`/`error_extension_tail`）、
  `HashMap` unused import（overlap2.rs）、lib.rs `#[macro_use]` 等；
- `scripts/verify-migrate.sh` 的 `asm_olc` 用例用 Lambda 数据（约 6 s），
  可考虑缩小输入加速日常核对；
- `notes/benchmarks/` 目录索引：迁移 `bbtools-vs-anchr.md` 后补一篇
  README 或索引（当前只有一篇）。

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
- 不重新引入 `bio`/`intspan` 外部 crate（pgr 已内嵌等价实现）。
