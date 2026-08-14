# multik 覆盖度实验：30× / 60× vs 40×（2026-08-15）

> 动机：用户口头提过"**30×/60× 两个覆盖度比较合适，不是覆盖度越高越好**"
> （G37 当时文档 `results/model.md` 记录的实际是 **40×/80×** + 各 Q/L
> 档位切分，见 `references/anchr-legacy-pipeline.md` §2.2）。用 G37 验证
> multik 在 30×/60× 下的表现，并试"多覆盖度部分分别组装后合并"
> （老流程语义）。

## 数据与命令

* 全量输入：`2_illumina/Q25L60/pe.cor.fa.gz`（622,308 reads × 150 bp，
  quorum 筛选后，未纠错）；
* 降采样（老流程同款）：`pgr fa split about -e -c $((580076 * cov))`，
  取第一部分（P000）——30× = 116,372 reads，60× = 232,724 reads；
  40× 用已有的 `4_down_sampling/Q25L60X40P000/pe.cor.fa`（155,148 reads）；
* 组装：`anchr asm multik <pe.cor.fa> -o out.fa`（`--kmer auto`，
  release build，全部默认参数）；
* 质检：`quast.py -r ~/data/anchr/ref/g37/genome.fa --min-contig 500`。

## 单覆盖度结果

| 指标 | 30× | 60× | 40×（对照） |
|---|---:|---:|---:|
| # contigs (≥500) | 57 | 65 | 60 |
| Largest contig | 52,720 | 45,634 | 56,183 |
| N50 | 20,291 | 19,467 | **24,445** |
| # misassemblies | **0** | **0** | **0** |
| Genome fraction (%) | 95.86 | **95.90** | 95.86 |
| # mismatches / 100 kbp | 28.57 | **25.87** | 27.67 |
| # indels / 100 kbp | 2.52 | **2.34** | 2.52 |
| # N's / 100 kbp | 0 | 0 | 0 |
| Duplication ratio | 1.001 | 1.001 | 1.001 |

**结论（用户经验得到验证）**：

1. **30×/60× 单独跑质量都达标**：misassemblies 全 0、0 N、Genome
   fraction ~95.9%（与 40× 持平），mismatch/indel 甚至更好（60× 全组最优：
   25.9/2.34）——低覆盖不是质量短板；
2. 40× 的 N50 仍最高（24.4K vs 20.3K/19.5K），但差距在 20-25%，且 30×/60×
   的 contig 数量/最长片段与 40× 同量级——覆盖度 30-80 区间对 multik 的
   主路径碎片化影响温和；
3. "不是覆盖度越高越好"在 multik 上也成立：60×（比 40× 高 50% 数据）没有
   带来 N50/覆盖度提升，反而 N50 略降——multik 的跨接验证/渐进过滤对
   高覆盖的边际收益有限。

## 多部分合并实验（30× + 60×，老流程 7_merge 语义）

`cat multik30.fa multik60.fa | anchr contained --len 1000 --idt 0.9999
--ratio 0.99999`（1612 → 51 条）：

| 指标 | 30×+60× 合并 | 40× 单跑 |
|---|---:|---:|
| N50 | 24,445 | 24,445 |
| # misassemblies | 0 | 0 |
| Genome fraction (%) | 95.31 | 95.86 |
| Duplication ratio | **1.124** | 1.001 |
| # mismatches / 100 kbp | 30.09 | 27.67 |

**结论**：简单 contained 合并把 N50 抬到 40× 水平（24.4K），但 duplication
1.124（两个部分各自组装的 unitigs 大量部分重叠，contained 只去"完全包含"）
且 Genome fraction 略降——**多覆盖度部分合并需要老流程的完整机制**
（contained → orient → merge → anchors 化，即 OLC 式合并；fill/glue 是
gap 填充辅助，见 `anchr-legacy-pipeline.md` §2.5/§2.7）去冗余，单纯
cat + contained 不够。multik 若做多部分合并，需先解决"部分重叠 unitigs
去冗余"（如覆盖率/长度证据合并或 overlap 图拼接）。

## 与老流程的对应

老流程 `4_down_sampling` 的 X30/X60（用户经验值）在 G37 模板里是 X40/X80；
本实验确认 30×/60× 也是 multik 的有效工作点。老流程"低覆盖部分少重复
冲突、高覆盖部分补洞"的互补性，multik 单部分跑已经隐含（渐进过滤 +
跨接验证对覆盖度不敏感），多部分互补的收益需要"分部分 → 合并"实验进一步
验证（todo §4 已记）。
