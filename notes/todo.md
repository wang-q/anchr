# anchr 待办

> 仅保留 actionable 待办；历史结论、已完成项与细节见各 design/benchmark/
> audit 文档（索引见 `project-understanding.md`；2026-08-15 会话交接见
> `design/asm-olc-modern-flow.md`；明确不做清单见
> `project-understanding.md` §6.4）。

## 待实现

- ~~G37 megahit/spades 策略对比分析~~（2026-08-16 完成，见
  `benchmarks/g37-megahit-spades.md`）：结论——GF 缺口 76% 是低复杂度/重复区
  （anchor 上界过滤排除）；外部组装器的 mis 全部在重复区边界；最大单项改进
  是提高主 K 上限（MR 链 31..81→31..121，7 组全链实验 N50 55.0K→83.6K、
  GF 96.96→98.63%、0 mis 保持）。
- ~~模板主 K 上限~~（2026-08-16 完成：`{4,6}_unitigs.tera.sh` 三 unitigger
  分支统一 `KS`，4_ 链 31..91、6_ MR 链 31..121；G37 MR 链全链复测 N50
  55.0K→83.6K、GF +1.7 pp、0 mis，见 `g37-megahit-spades.md` §7/§8）。
- ~~multik 气泡合并~~（2026-08-16 完成：`multik.rs` `bubble_merge` 按
  megahit `SearchAndPopBubble` 语义（两端唯一 + banded 相似度 ≥0.95 +
  ≤20k，branch 节点不参与，备选路径独立输出），CLI `--merge-similar`/
  `--merge-len`；G37 MR 链实测 contigs 13→12、Dup 1.005→1.001、0 mis、
  N50 持平、GF -0.05 pp；5 个单元测试 + CLI 冒烟测试，见 §8）。
- ~~unitig 端点本地组装/缺口闭合~~（2026-08-16 完成：新命令
  `anchr asm extend CONTIGS READS...`，base-by-base reads k-mer 图游走 +
  严格多数（`--min-support` 且 `>=2x` 次高）防错连；`{4,6}_unitigs`
  三 unitigger 分支在 olc 合并后自动调用。G37 MR 链 GF 98.585→98.809%、
  未覆盖 8,210→6,907 bp、0 mis，见 `g37-megahit-spades.md` §8/§10）；
- 重复区解析（P2，剩余 7,928 bp：高 k 单主取用 + 双端连接证据裁决，歧义
  保持断开）——§8；
- QUAST 口径：mm/indel 按覆盖 bp 归一化 + reads-vs-consensus 一致率列
  （P2，G37 高 k 链 mm 上升是 reads 与 NC_000908 在低复杂度区的真实差异，
  非 consensus 劣化）——§8。
- ~~`fq range` 剩余~~（2026-08-15 全部完成：双端感知 S2 已核对 +
  BGZF `.gzi` 自动生成已实现，见 `design/fq-range.md` §7）；
- ~~`dep`/`ena`/`template` 的外部工具版本核对~~（2026-08-15 完成：
  `check_dep.sh` 已更新现代依赖清单；`ena` 已 Rust 化（`ena meta`/
  `ena manifest`，ureq 客户端）；`template` 无外部工具依赖）。
- ~~`anchr contained` 命令删除~~（2026-08-15 完成：`src/cmd/contained.rs`
  与 `src/libs/overlap.rs` 一并删除；8_* 模板直接用 `contigs.fasta`/
  `final.contigs.fa` 产出 anchors，`9_quast`/`9_busco`/`9_stat_final`
  同步去掉 non-contained 引用；G37 全流程验证通过）。
- ~~G37 全流程端到端验证 + `results/model_org.md`~~（2026-08-15 完成：
  现代模板链（multik → anchor → olc --unitigs）在 G37 全量数据跑通，
  文档已生成；过程中修掉：2_quorum 对 Q0L0 未压缩 `.fq` 的兼容、
  4/6_anchors 的 reads 路径、8_spades repair 模式（fint→rp）、
  8_megahit 奇数交错 reads（`--12` 需偶数）、7_merge_anchors 第三次
  合并的 find 遗漏 `anchor.merge.fasta`、`asm anchor` 新增 `--stats`
  输出（Mapped/median/MAD/lower/upper/SumOthers，替代旧 env.json 统计）。
- ~~`anchr template --unitigger` 恢复 bcalm + 加入自研 unitig~~（2026-08-15
  完成：`--unitigger` 支持 multik（默认）/ unitig / bcalm 空格分隔列表；
  4/6_unitigs、0_master/0_bsub、0_cleanup、9_quast 全部按 unitigger 参数化；
  unitig/bcalm 每 k 独立 unitigs（31..81）后用 `asm olc --unitigs` 跨 k
  合并——G37 bcalm 链 N50 31199 与 legacy 一致；MG1655 全链验证 `asm
  unitig` 与 bcalm 等价（见下）。

## 挂账 / 待决

- `fq merge` 嵌合 reads（暂不修；机制、失败方案与 merge/no-merge 取舍见
  `fq-merge-replace.md` §9）：短读 merge 无法区分重复序列重叠与真实 insert
  重叠，reads 层无信号，由 multik 图防御兜底；merge 保留（N50 +54%）；
  真实宏基因组/长读数据到位后重审。
- pgr rev 更新流程 + `.cargo/config.toml` patch 发布差异
  （`project-understanding.md` §8.1）；
- pgr supermer 质量门控：落地后 FASTQ 可去掉 direct 回退
  （`asm-assemble.md` §12.1）；
- pgr 并行读 gz（`asm-assemble.md` §12.3，风险>收益暂缓）；

## 待验证 / 等数据到位

- 真实宏基因组/长读（决定拆分合并路线定位）：multik tip/weak_link 价值、
  `asm olc` 四命令端到端、anchor 补洞（老流程 fill）+ others 输出、
  megahit low depth/本地组装与 SKESA SNP 簇变体优先级
  （`megahit.md` §8、`skesa.md` §7.2）、SAM 内存化评估（`asm-olc.md` §14.6）；
- multik `--parallel` 扩展性复测、`--min-count-extend` 调参
  （`asm-multik.md` §7）；
- multik 性能：计数复用 / remove_unsupported 查表化 / 轮数裁剪
  （`benchmarks/multik-complexity.md`）；
- multik 防 misassembly：bridge 探针长度/阈值调参、<500 bp 碎片输出策略
  （`asm-multik.md` §9）；
- `asm cns` contain 吸收（`asm-olc.md` §14.3-3/§S3：列投票去冗余同时保留
  覆盖证据）；
- `asm olc` 参数验证：overlap 少量错配、不同 k unitig 冗余去重、repeat
  breaking 覆盖度阈值（`references/canu.md` §8）；
- ~~大规模真实数据全链~~（2026-08-15 G37 全量：`fq → asm → template`
  跑通，结果见 `results/model_org.md`；统计核对完成）；
- ~~MG1655 multik vs bcalm 初步对照~~（2026-08-15 初步数据，详见下条
  端到端对照）
- ~~MG1655 bcalm 链端到端对照~~（2026-08-15 完成，5 组同输入
  MRX40P000/P001/P002 + MRX80P000/P001，唯一变量 unitigger）：unitig N50
  bcalm 42-58K、`asm unitig` 51-61K vs multik 19-21K；最终 merge N50
  **unitig 链 105.7K / bcalm 链 95.5K vs multik 23.4K**——N50 差距基本
  全部来自 multik unitig 碎片化，`asm unitig` 与外部 bcalm 等价且略优
  （自研可替代外部依赖）；mis 4→1（legacy 0）。
- ~~merge 近似重叠去重（dup 1.07-1.08）~~（2026-08-15 完成：
  `asm olc` consensus 后新增 `merge_overlapping_contigs`：跨组 anchors
  边界不一致导致 exact overlap 连不上、残留同区域不同边界的 contig 对；
  现用 31-mer 定位主导 offset + 头部锚定 + 重叠区 ≥99% 一致才合并，嵌合
  多块对齐拒绝——unitig 链 Dup 1.079→1.000、bcalm 链 1.068→1.001，
  GF 97.77% 不变，unitig 链 N50 提到 110.2K；mis 仍 1（contig_26
  relocation，属 multik 防嵌合任务）。
- ~~multik 碎片化修复（auto k 序列 k0 太低）~~（2026-08-15 完成：
  multik 的 unitig 骨架冻结在 pass 0 的 k0，k0=21/31 时 MG1655 N50 仅
  21K（迭代轮并不切碎，1661→1661 条稳定；单跑 asm unitig k=81 就有
  53.5K/705 条）。`auto_ks` k0 从 `N50/10` 改为 `N50/3`（clamp 31..51）：
  MG1655 auto → 50/70/90/110，unitig N50 **21K → 58K**（687 条），
  5 组全链 merge N50 **23.4K → 65.8K**（128 contigs），GF 97.36%；
  395 测试全绿。k0 越高越低覆盖丢失风险，宏基因组需再验证。
- ~~multik 防嵌合（mis 4→0）~~（2026-08-16 完成：机制、修复、验证与
  unitig/bcalm/multik 三链对比见 `asm-multik.md` §9.8 +
  `benchmarks/mg1655-unitig-bcalm-multik.md`；397 测试全绿）
- multik 多主 K 架构（2026-08-15 落地）：用户裁定 Rust 内多骨架并行太慢，
  改为**主 K/从 K**：multik 一次只跑一个主 K（ks[0] 骨架，更大 k 验证），
  模板（4/6_unitigs）用 bash 并行跑每个主 K（31..81），`asm olc --unitigs`
  跨主 K 合并。G37：N50 39K/427 条 → 55.4K/**20 条**（~1 分钟）；
  MG1655 单组：21K/1455 → **95.5K/108 条**（~8 分钟，6 主 K 并行 -p 4）。
- **关键发现**（2026-08-15）：multik 的 mis 随主 K 增大而减少——K=31 主
  3 mis、K=51 主 1 mis、**K=81 主 0 mis**（N50 58.8K、GF 97.88% 最高、
  Sum 4.62M≈基因组 99.4%）；多主 K 合并（olc）引入 mis（全 K 6 个、仅
  高 K 61/71/81 也有 3 个——跨主 K 的 unitigs 经重复序列 exact-overlap
  错连）。**高 K 单主（无迭代）已是正确性最优**，迭代验证的价值需在
  低覆盖/长读数据再确认；方向 2（重复区断开）待做。
- ~~mis 根因：现代 olc `--min-overlap 34` 太松~~（2026-08-15 修正，用户
  指出 legacy E. coli 0 mis 后复查）：legacy `orient/merge` 要求
  **overlap ≥ 1000 bp + idt ≥ 99.9%**（daligner 全长比对），现代
  `asm olc --unitigs` 默认**端点 17-mer 种子 + 34 bp exact overlap**——
  短重复序列（46bp 反向重复）被当 overlap → inversion mis（contig_73
  实测）。**min-overlap 提到 1000 后**：K=81 单主链 6→**0 mis**、unitig
  链 6 k 合并 1→**0 mis**（N50 53.6K/146 条、GF 97.66%）；N50 略降
  （短 overlap 不合并，宁断勿错，与 legacy 一致）。模板 7_merge 与
  4/6_unitigs 合并已统一 `--min-overlap 1000`。低主 K（31/41）unitigs
  自身仍有 mis（multik 迭代错连，非 olc 合并），方向 2 仍待做。
- gz/大输入回归：默认 supermer 路径全链回归 + 峰值内存；
- 560 bp 碎片 mis 覆盖度门槛（`asm-olc.md` §14.3，`--min-contig-len 1000`
  可滤，可选）。
- `tsv-sample` 无种子（2_quorum/2_merge 的 shuffle）：down-sampling 随机
  分桶会让 anchor/merge 有运行间波动。已修掉由此暴露的 dup 问题：
  `consensus::coverage()` 从"最长单段覆盖率"改为"多块区间累计覆盖率"
  （合并输出中整条被另一 contig 分两段覆盖的冗余 contig 现在会被
  `dedup_contained_ratio` 0.99 丢弃），G37 merge_anchors dup 稳定回到
  **1.000**（此前随机差的一轮为 1.205），0 mis、0 N、GF ~95.9% 不变；
  新增测试 `dedups_multi_block_contained_contigs`；另加 `asm anchor`
  边界校正（contig 两端 readl/2 窗口线性缩放，GF 95.79→95.99%）。
  完全可复现仍可给 `tsv-sample` 加固定 `--seed`（待决）。

## 低风险审计（可顺手修）

- ~~`scripts/verify-migrate.sh` 的 `asm_olc` 用例加速~~（2026-08-15 评估：
  输入已是 `R1.2k`（2k reads 缩小版），约 6 s 是双 k OLC 合理开销；脚本
  已历史化（仅旧 pgr 二进制可用），加速价值低，不做）。
- warning 清理（2026-08-15 完成）：`tadpole` 的
  `error_extension_pincer/tail` 未读字段已删；`Overlap` 的 dead_code 已
  随 `contained` 命令删除（`src/libs/overlap.rs` 移除）一并解决。

## 技术债

- golden 数据体积：`tests/bbtools/Lambda/` 约 13 MB，可考虑外部数据源或
  生成脚本；
- 外部工具链依赖：流程命令强依赖 dazzler/hnsm（PATH），测试缺失时跳过
  可能掩盖回归（`project-understanding.md` §8.3）。
