# MG1655 unitig / bcalm / multik 三链对比（2026-08-16）

> 同一批 5 组 MG1655 reads（`6_down_sampling/MRX40P000/P001/P002`、
> `MRX80P000/P001`）上，用当前代码对比三条组装链的端到端质量。

## 命令（统一口径）

* 每条链每 k / 每主 K 独立产出 unitigs：
  - multik：`anchr asm multik -k <K,后续 k 列表> -p 8`，主 K = 31..81
    （模板 6_unitigs multik 分支）；
  - unitig：`anchr asm unitig -k <K>`，K = 31..81；
  - bcalm：外部 bcalm，K = 31..81；
* 跨 k 合并：`anchr asm olc --unitigs unitigs_K*.fasta --min-overlap 1000
  --min-contig-len 1000` → `unitigs.fasta`；
* anchors：`anchr asm anchor unitigs.fasta reads.fa --mincov 5 --mscale 3
  --lscale 3 --uscale 2 -p 8`；
* 最终合并：`anchr asm olc --unitigs --list-files anchors.list
  --min-overlap 1000 --min-contig-len 1000`；
* 质检：quast，参考 `ref/mg1655/genome.fa`。

## 结果（quast，contigs ≥ 1 kb）

| 指标 | multik | unitig | bcalm |
| :--- | ---: | ---: | ---: |
| # contigs (≥ 1 kb) | **107** | 121 | 121 |
| # contigs (≥ 10 kb) | **78** | 88 | 88 |
| N50 | **95,478** | 67,358 | 67,358 |
| # misassemblies | 0 | 0 | 0 |
| Genome fraction (%) | **97.61** | 97.54 | 97.54 |
| mismatches / 100 kbp | 0.00 | 0.00 | 0.00 |

## 结论

* **multik 全面占优**：N50 比 unitig/bcalm 高 42%（95.5K vs 67.4K），
  contigs 更少、GF 略高，三链均 0 mis——multik 迭代验证（更大 k 桥接
  验证 + 跨主 K 合并）的价值得到端到端确认；
* **unitig 与 bcalm 等价**：N50/GF/contigs 完全一致（`asm unitig` 是
  bcalm 的 Rust 复刻，自研可替代外部依赖）；
* **多主 K > 单主 K**：multik 跨主 K（31..81）N50 95.5K，优于单主 51
  （60.3K，见 `asm-multik.md` §9.8）；
* **代价**：multik 迭代最耗时（每组 6 主 K 并行 ~1-2 分钟 vs unitig/bcalm
  每 k 独立几秒）。

## 备注

* 三链均 0 mis：multik 的 4 mis 由 §9.8 修复归零；unitig/bcalm 旧报告
  的 1 mis 是 2026-08-15 旧 OLC 生成方式遗留（80 bp 重复 overlap 被旧版
  接受），当前 `olc --min-overlap 1000` 重跑已归零（见 `asm-multik.md`
  §9.8）；
* merge vs no-merge 对照与 `fq merge` 嵌合分析见 `fq-merge-replace.md`
  §9。
