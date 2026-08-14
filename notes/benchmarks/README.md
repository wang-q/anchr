# Benchmarks 索引

| 文档 | 内容 |
| :--- | :--- |
| `bbtools-vs-anchr.md` | fq 命令 vs BBTools 39.38 CLI 基准（随迁自 pgr，`kmer hist` 行属 pgr） |
| `qc-bench.md` | `fq qc` 端到端基准（anchr vs FastQC vs Falco） |
| `unitig-bench.md` | `asm unitig` 效率基准（supermer 默认计数、DFA walk） |
| `multik.md` | `asm multik` 吞吐 sanity check（1 Mb 合成基因组 → 单条 100%，9.4 s / 816 MB） |
| `multik-complexity.md` | 复杂度验证：前期小 k 压实 vs 大 k 图小（G37 每轮分解；图结构递减 ✓、耗时未递减，remove_unsupported 为瓶颈） |
| `multik-g37-quast.md` | G37 Quast 质检：multik 序列质量最优（mismatch 32/indel 2.5/0 N），misassemblies 8 → **0**（bridge_filter + split_by_bridge），N50 26.6K |
| `multik-cov.md` | 覆盖度实验（2026-08-15）：30×/60× 单跑质量达标（0 mis、~96% 覆盖，60× mismatch 最优 25.9），40× N50 略优；30×+60× contained 合并 N50=40× 水平但 duplication 1.124（多部分合并需完整 anchors+OLC 机制） |
| `multik-allgroups.md` | 全分组复核（2026-08-15）：按老流程 model.md 的 Q/L × X40/X80 × P + MR 全部 23 组跑 multik——**23/23 组 0 mis、0 N、dup≤1.001**；MR 组全面更优（N50 34-55K、GF 96.1%），X40 N50 优于 X80，Q/L 档位影响小；**`anchr asm olc` 合并多组输出（须先 contained 去重）N50 54.9K ≈ 老流程 7_merge_anchors 55.0K、0 mis；单组 MR 已近天花板，合并且只补覆盖（GF +0.3-0.5pp）、mm 反升** |
