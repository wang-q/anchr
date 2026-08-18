# DH5alpha relocation 嵌合调查记录

> 维护约定：本文件是 DH5alpha 多组多 k 流程中 relocation 嵌合的完整调查快照。
> **任何后续尝试本问题的方案前，必须先读本文件**，避免重复已完成/已排除的手段。
> 命令级门禁数据以 `results/model_org.md` 为唯一权威；本文件只记录调查过程与结论。

## 1. 问题现象与目标

**现象**：DH5alpha 13 组下采样链（模型门禁链）最终 `quast -m500` 报 **2 条 relocation
misassemblies**（最初为 contig_3 / contig_18，之后演变为 contig_4 / contig_22）。

**目标**：消除这 2 条 relocation 嵌合，使 misassemblies 降为 0（与 G37 / MG1655 基线一致）。

**现状（2026-08-18）**：**misassemblies 已归 0**。`contig_4` 由 extend 低覆盖缝检查
消除；`contig_22` 由 extend 跨 contig 所有权护栏消除（见 §3.5）。

## 2. 关键文件与当前改动（工作区未提交）

| 文件 | 改动 | 作用 |
|------|------|------|
| `src/libs/asm/multik/graph.rs` | `internal_repeat_bridge_split`：从高覆盖尖峰(RATIO=4)改为**低覆盖缝**(LOW_RATIO=0.3, MIN_RUN=5, MAX_RUN=200, GAP=3, MIN_MEDIAN=8, SPAN=2000)；`is_repeat_bridge` 按 junction probe count / unitig 中位数覆盖率比(1.5x)阻断跨重复合并；`recompact_graph` 重复桥检测 | multik 内拆/阻断嵌合 |
| `src/libs/asm/multik/master.rs` / `schedule.rs` | 31-mer 短 k 重复表（探测 SNP 型串联重复） | 捕获 k-长 k 表漏掉的重复桥 |
| `src/libs/asm/extend.rs` | walk 覆盖率检查：`LOW_RATIO=0.3`, `MIN_LOW_RUN=5`, `contig_median_count()` | 连续≥5 步 k-mer 计数 < 0.3×中位数即回滚停止，阻断跨低覆盖缝延伸 |
| `results/model_org.md` | DH5alpha 模板 `--unitigger` 去掉 `bcalm`、去掉非法 `--redo` | 模板修正 |

## 3. 已尝试且最后被采纳的方案（有效）

### 3.1 multik 内部低覆盖缝切分（`internal_repeat_bridge_split`）
- **v1 高覆盖尖峰**（4×median, FLANK 隔离检查）→ 过度切割，28 条 unitig 被切，
  输出碎片化（最长 2.5 kb）→ 废弃。
- **v2 低覆盖缝**（0.3×median）→ 每条仅切 4-5 处，保留 258 kb 大 unitig → 采纳。

### 3.2 重复桥阻断（`is_repeat_bridge` / `recompact_graph`）
- junction probe count > 1.5× 该 unitig 中位数覆盖 → 视为重复桥，不合并。
- 成功消除 MRX40P001 的 95 kb contig_18。

### 3.3 31-mer 短 k 重复表
- 长 k 重复表漏掉 SNP 型串联重复；31-mer 表补上，阻断 MRX80P001 的 u4→u314 重复桥。

### 3.4 extend 低覆盖缝检查（采纳，消除 contig_4）
- walk 连续 5 步计数 < 0.3×中位数即回滚，消除 `contig_4`。

### 3.5 extend 跨 contig 所有权护栏（2026-08-18 采纳，消除 contig_22 → 0 mis）
- **根因**：contig_23 的 3' 端 walk 越过强覆盖保守重复区，接进远程 3922.2 kb 位点
  （该序列已被同组 contig_19 装配占有）。接合处覆盖不低（74-104×），低覆盖检查
  无效——这是真正的高覆盖保守重复区 relocation。
- **实现**（`src/libs/asm/extend.rs`）：`cross_contig_kmers` 预建所有输入 contig 的
  规范 k-mer 归属索引（`sole` = 仅属于单条 contig 的 k-mer → 其索引；`multi` = 出现
  于 ≥2 条不同 contig 的 k-mer）。walk 每步检查新窗口：若 `multi` 含它、或其 `sole`
  属主 ≠ 当前 contig（即 foreign），连续 ≥ `MIN_FOREIGN_RUN(5)` 步即回滚停止，去掉
  进入他人领地的延伸（回滚到最后一个自属窗口）。
- **验证**（`/tmp/dh5alpha_gate6`，仅重跑 MRX80P001 extend+anchor，其余 12 组复用
  gate5 anchor）：**misassemblies 1 → 0**，N50 82.9 k 持平、contigs 119→122（无碎片化）、
  Genome fraction 98.84 → 98.34（-0.5 pp，宁断勿错的覆盖代价，约 ~23 kb 唯一序列）。
  local misassembly 2 者均为 1（既有，未新增）。
- **回归门禁**（`/tmp/regress_extend.sh g37|mg1655`，仅改 extend.rs，其余步骤/数据不变，
  复用各 pre-extend unitig，比较新旧二进制的 quast）：

  | 数据集 | misassemblies | Genome fraction | N50 |
  |--------|---------------|-----------------|-----|
  | G37  基线 (gate) → 新 (gate2) | 0 → **0** | 97.979 → 97.743 (-0.24 pp) | 55170 → 48850 |
  | MG1655 基线 (fix) → 新 (fix2) | 0 → **0** | 98.856 → 98.418 (-0.44 pp) | 96447 → 95499 |

  两数据集 mis 均保持 **0**（不高于基线），未引入新错装；GF 代价 ~0.2-0.4 pp，
  与 DH5alpha 的 -0.5 pp 同量级，属宁断勿错的覆盖代价，回归门禁通过。
- **单测**：`stops_when_crossing_another_contig`（contig A 的 overhang 接进仅被 contig
  B 占有的区域 → 截断）；现有 5 个单 contig 测试不受影响（sole/multi 为空）。

## 4. 已排除的手段（不要重复）

| 手段 | 结果 | 原因 |
|------|------|------|
| 调 multik `min_count_extend` | 无效 | 嵌合不源自 k-mer 阈值 |
| 调 `probe_half` | 无效 | 同上 |
| 重复区支持阈值 RATIO=1.5 + 130-mer probe | 失败 | 未消除嵌合 |
| 高覆盖尖峰切分 v1 (RATIO=4, FLANK) | 碎片化严重 | 误切真实区域 |
| `--cross-validate` 跨组 OLC `drop_cross_chimeras` | 无法消除 contig_22 | 嵌合嵌在**单个** anchor 内部，跨文件级校验无效 |

## 5. 残差 contig_22 的根因（2026-08-18 定位）

**起源链（MRX80P001）**，用 `/tmp/probe_real.fa`（4498.5 kb 位点 200bp probe）与
`/tmp/probe_chi.fa`（3922.2 kb 位点 200bp probe）扫描各阶段：

| 阶段 | 文件 | 结果 |
|------|------|------|
| multik | `MRX80P001/unitigs_all.fasta` | 仅 real，**干净**（unitig_125-130/216/217） |
| per-group OLC | `MRX80P001/unitigs.fasta` contig_23 | 仅 real，干净 |
| **per-group extend** | `MRX80P001/unitigs.ext.fasta` contig_23 | **同时含 real + chi + 嵌合在此形成** |
| anchor | `MRX80P001/anchor.fasta` anchor_24_contig_23_173-64837 | 携带嵌合进后续 |
| cross OLC+extend | `merge.fasta`→`merge.ext.fasta` | 在嵌合 anchor 上再 +1000 bp（两端各 500） |

**结论**：嵌合由**单组 `asm extend`** 造成。contig_23 的 3' 端 k-mer walk 从真实位点
延伸越过一个**强覆盖保守重复区**，接到远端 3922.2 kb 位点。

**为什么低覆盖 seam 检查拦不住**：MRX80P001 是高覆盖样本（~80x），3922.2 kb 重复位点
在它自己的 reads 中覆盖极高，接合处 k-mer 计数不低、有单一明确路径，**并非低覆盖 seam**。
这与已消除的 contig_4 本质不同（contig_4 的接合位点在 MRX40P000 reads 中近零覆盖，故
低覆盖检查生效）。这是真正的保守重复区强覆盖 relocation，单一覆盖率度量无法与正确延伸
区分。

**关联**：这正是 `results/model_org.md` 记载的 DH5alpha 跨组保守重复区 relocation 已知
边界。G37/MG1655 基线为 0 不受影响。

## 6. 后续方向（已实现跨 contig 所有权，见 §3.5）

延伸 walk 进入已被其他 contig 占有的 k-mer 区判定为嵌合这一方向**已落地**（§3.5），
跨 contig relocation 在 DH5alpha 上归零。以下为仍可考虑的保守化余量（勿在未 A/B 前
直接改公式）：

1. **重复上下文抑制**：extend 时检测 seed/延伸 k-mer 的多重性，在重复上下文中进一步
   缩减延伸。当前跨 contig 所有权已等价覆盖大部分重复区场景，可作为更激进路线。
2. **高覆盖样本保守化**：对高覆盖样本限制 `--max-extend` / 提升 `--min-support`。

## 7. 手动分析脚本（/tmp，仅参考）

- `/tmp/find_seq.py`：在 fasta 中精确检索序列子串，定位嵌合来源。
- `/tmp/prof_terminal.py`：对 unitig 做 31-mer 覆盖剖面，找低覆盖缝。
- `/tmp/map_tail.py`：将 unitig 3' 端按 24-mer seed 比对参考，定位嵌合接合。
- `/tmp/scan_source.py`：确定嵌合 contig 的源组/unitig 归属。
- 门禁链：`/tmp/dh5alpha_gate4/run.sh`、`/tmp/dh5alpha_gate5/run.sh`
  （gate4=低覆盖缝 split 代码，gate5=含 extend 保守检查）、`/tmp/dh5alpha_gate6`
  （含 extend 跨 contig 所有权；仅重跑 MRX80P001，其余组复用 gate5 anchor）。
- quast 结果：gate5 `quast/report.tsv`（1 mis：contig_22）、gate6 `quast/report.tsv`
  （**0 mis**）。

## 8. 需更新的下游记录

- 若最终接受 1 mis 或设法降 0，须把结论追加进 `results/model_org.md`（quast 参数保持
  既有表格一致，禁止变更导致纵向对比失效）。
- multik 低覆盖缝 / extend 保守检查的相关参数若继续调整，需按实验纪律在 G37 与 MG1655
  上过质量门禁（misassemblies 不得高于基线 0）。