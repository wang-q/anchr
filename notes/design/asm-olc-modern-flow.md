# 现代组装流程总结（2026-08-15 会话收官）

> 本文是长会话（SKESA/multik/OLC 改造/anchors/替代盘点）的完整总结，
> 供后续会话快速恢复上下文。快速交接看 §1 现代流程与 §2.2 用户裁定；
> 机制细节见 `asm-olc.md` §14、`references/anchr-legacy-pipeline.md` §5、
> `benchmarks/multik-allgroups.md`。

## 1. 现代流程最终形态

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

## 2. 会话完成的工作

### 2.1 SKESA 借鉴（multik）

* 读 SKESA-master（C++）+ skesa-rs-main（Rust 移植），核对笔记并补充
  §7.2（与 multik 的多 k 迭代对比）；
* 落地"严格链唯一性"：`merge_chains`/`recompact_graph` 从先到先得占用
  改为**严格两端唯一**（SKESA "predecessor == 1"），含对称 link 去重修正
  （`compute_links` 双向 link 导致度数翻倍）。
* G37：misassemblies 0、mm 27.7/100kbp 历史最佳、N50 24.4K（-8% 宁断勿嵌合）。

### 2.2 老流程理解修正（用户多次纠正，重要！）

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
* 唯一保留外部组件：**quast**（用户要求最终质量确认）。

### 2.3 multik 验证与基准

* 全分组复核：老流程 23 组全部 0 mis / 0 N / dup≤1.001，MR 组全面更优
  （N50 34-55K）、X40 优于 X80、全量 reads 反而最差（`multik-allgroups.md`）；
* 覆盖度实验：30×/60× 单跑质量达标（60× mismatch 最优 25.9）
  （`multik-cov.md`）。

### 2.4 OLC 改造（现代流程的命令支撑）

* **`asm olc --unitigs`**：跳过 S0（不二次组装切断长 unitig）、多文件 tag
  防撞名、内部 filter_contained；
* **`asm layout --filter-contained`**：独立管道与驱动器语义对齐；
* **`--dedup-ratio`**（cns）：contig 级近似包含去冗余（允许 ~1% 错配，
  对齐 `anchr contained --idt 0.99`）；
* G37 23 组合并：dup 2.659 → 1.094、GF 96.54%、N50 54.9K、0 N、0 大 mis。

### 2.5 `anchr asm anchor`（用户建议：reads mapping → anchors → OLC）

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

### 2.6 时间分析（用户确认可接受）

* anchor 单组 0.24s（map 58%、SAM I/O 21%、覆盖 12%），7 组串行 1.6s；
* 现代流程大头是 **multik**（单组 2.35s，占 ~65%）；
* 优化点保留：SAM 内存化（宏基因组时）、multik 性能（multik-complexity.md）。

## 3. 当前状态

* 测试：**29 套件全绿、0 失败**；fmt/clippy 干净；
* 工作区：**17 个文件未提交**（14 修改 + 3 新增：`anchor.rs`×2 + 测试），
  `.git` 只读需用户本机 commit；
* 文档：asm-olc.md §14（改造+时间分析）、anchr-legacy-pipeline.md §5
  （替代对照）、multik 基准、todo 全部同步。

## 4. 下一步（详见 todo）

1. 真实宏基因组/长读数据验证（决定拆分合并路线的最终定位）；
2. anchor 补洞逻辑（老流程 fill）；
3. 560bp 碎片 mis 的覆盖度门槛（`--min-contig-len 1000` 可滤，可选）；
4. SAM 内存化（宏基因组时再做）；
5. multik 性能优化（计数复用 / remove_unsupported 查表化 / 轮数裁剪）。
