# Benchmarks 索引

| 文档 | 内容 |
| :--- | :--- |
| `bbtools-vs-anchr.md` | fq 命令 vs BBTools 39.38 CLI 基准（随迁自 pgr，`kmer hist` 行属 pgr） |
| `qc-bench.md` | `fq qc` 端到端基准（anchr vs FastQC vs Falco） |
| `unitig-bench.md` | `asm unitig` 效率基准（supermer 默认计数、DFA walk） |
| `multik.md` | `asm multik` 吞吐 sanity check（1 Mb 合成基因组 → 单条 100%，9.4 s / 816 MB） |
| `multik-complexity.md` | 复杂度验证：前期小 k 压实 vs 大 k 图小（G37 每轮分解；图结构递减 ✓、耗时未递减，remove_unsupported 为瓶颈） |
| `multik-g37-quast.md` | G37 Quast 质检：multik 序列质量最优（mismatch 32/indel 2.5/0 N），misassemblies 8 → **0**（bridge_filter + split_by_bridge），N50 26.6K |
