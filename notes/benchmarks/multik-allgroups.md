# multik 全分组复核：按老流程 40×/80× × Q/L × P 每组（2026-08-15）

> 背景：老流程（`results/model.md`）按 **Q/L 档位（Q0L0/Q25L60/Q30L60）×
> 覆盖度（X40/X80）× P 随机副本**切分 reads，另加 merged reads 的 MR
> 版本（6_down_sampling），每组独立组装。用 multik 把 **全部 23 组** 都
> 复核一遍（对齐老流程 statQuorum/statUnitigs 的逐行记录）。

## 命令

* 输入：`4_down_sampling/<Q L X P>/pe.cor.fa`（16 组）+ `6_down_sampling/
  <MRX P>/pe.cor.fa`（7 组），全部为 quorum 筛选后（未纠错）reads；
* 组装：`anchr asm multik <pe.cor.fa> --parallel 4 -o multik.fa`
  （`--kmer auto`，release build，默认参数；8 路并行 × 4 线程）；
* 质检：`quast.py <23×multik.fa> -r ref/g37/genome.fa --min-contig 500`。

## 结果（quast，contigs ≥ 500 bp）

| 组 | # contigs | N50 | Largest | mis | GF (%) | mm/100k | indel/100k | dup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| MRX40P000 | 28 | 39,148 | 79,567 | 0 | 96.20 | 29.56 | 2.51 | 1.000 |
| MRX40P001 | 26 | 54,888 | 178,776 | 0 | 96.11 | 30.67 | 2.33 | 1.000 |
| MRX40P002 | 32 | 34,086 | 73,192 | 0 | 96.20 | 30.46 | 3.76 | 1.000 |
| MRX40P003 | 28 | 39,104 | 68,099 | 0 | 96.14 | 29.76 | 3.05 | 1.000 |
| MRX40P004 | 24 | 37,552 | 179,610 | 0 | 96.09 | 29.78 | 2.87 | 1.000 |
| MRX80P000 | 35 | 34,012 | 113,997 | 0 | 96.04 | 30.15 | 2.33 | 1.000 |
| MRX80P001 | 35 | 54,964 | 178,776 | 0 | 96.08 | 30.86 | 3.05 | 1.000 |
| Q0L0X40P000 | 56 | 19,467 | 70,464 | 0 | 95.68 | 26.65 | 2.70 | 1.001 |
| Q0L0X40P001 | 52 | 24,527 | 68,315 | 0 | 95.75 | 27.71 | 2.52 | 1.001 |
| Q0L0X40P002 | 55 | 19,467 | 39,920 | 0 | 95.60 | 27.57 | 2.52 | 1.001 |
| Q0L0X40P003 | 60 | 19,615 | 79,461 | 0 | 95.84 | 26.60 | 2.16 | 1.001 |
| Q0L0X80P000 | 70 | 21,494 | 67,515 | 0 | 95.49 | 25.80 | 2.16 | 1.001 |
| Q0L0X80P001 | 74 | 16,719 | 54,460 | 0 | 95.52 | 23.62 | 1.98 | 1.001 |
| Q25L60X40P000 | 60 | 24,445 | 56,183 | 0 | 95.86 | 27.67 | 2.52 | 1.001 |
| Q25L60X40P001 | 51 | 19,467 | 56,941 | 0 | 95.36 | 26.92 | 2.53 | 1.001 |
| Q25L60X40P002 | 60 | 24,445 | 79,544 | 0 | 95.85 | 27.86 | 2.34 | 1.001 |
| Q25L60X40P003 | 60 | 18,347 | 52,023 | 0 | 95.86 | 26.78 | 2.34 | 1.001 |
| Q25L60X80P000 | 67 | 20,147 | 61,072 | 0 | 95.53 | 26.68 | 2.16 | 1.001 |
| Q25L60X80P001 | 77 | 17,700 | 52,315 | 0 | 95.88 | 24.61 | 1.98 | 1.001 |
| Q30L60X40P000 | 62 | 20,291 | 47,698 | 0 | 95.96 | 26.92 | 2.69 | 1.001 |
| Q30L60X40P001 | 55 | 20,281 | 43,730 | 0 | 95.79 | 26.98 | 2.34 | 1.001 |
| Q30L60X40P002 | 58 | 18,361 | 68,315 | 0 | 95.88 | 27.31 | 2.52 | 1.001 |
| Q30L60X80P000 | 76 | 17,277 | 49,592 | 0 | 95.64 | 23.95 | 2.34 | 1.001 |

（mis = # misassemblies，GF = Genome fraction，mm = mismatches/100kbp，
dup = Duplication ratio；全部 0 N / 100 kbp。）

## 结论

1. **23/23 组正确性全达标**：misassemblies 全 0、0 N、Duplication ≤ 1.001
   ——multik 在所有 Q/L × 覆盖度 × 副本组合下都稳定，不挑数据子集；
2. **MR（merged reads）组全面优于非 MR**：N50 34-55K vs 17-24K、GF
   96.0-96.2% vs 95.4-95.9%、最长 179K vs 79K——merged reads 的长读长
   让 multik 的 auto k 更大、跨接验证更有效（老流程 statUnitigs 也显示
   MR 的 N50Anchor 更高：31-34K vs 7-25K，方向一致）；
3. **非 MR 组：X40 的 N50 优于 X80**（19-24K vs 17-21K），X80 的
   mismatch/indel 略低（23.6-26.7 vs 26.6-27.9）——覆盖度不是越高越好，
   与老流程 statQuorum（X80 的 anchors 更多、N50Anchor 更低）一致；
4. **Q/L 档位影响小**：Q0L0/Q25L60/Q30L60 三档的 N50/GF/mismatch 同量级
   （trim 档位主要影响 reads 数 677K→638K，multik 对质量过滤不敏感）；
5. **P 副本有波动**：同组不同 P 的 N50 差 ~6K（如 Q25L60X40 18.3-24.4K），
   但正确性不变——单组结果解读时要看副本范围；
6. **对照老流程单组 anchors**（statQuorum：N50Anchor 7.3-25.4K、anchors
   39-103 条、Sum 546-562K）：multik 每组 ≥500 bp contigs（24-77 条）
   的 N50 相当或更高、条数更少，Sum 对应 GF 95-96% 的参考覆盖——multik
   单组输出已不弱于老流程单组 anchors，老流程的优势在跨组合并
   （7_merge 的 OLC 式合并，fill/glue 是 gap 填充辅助）。

## 与 multik-cov.md 的关系

`multik-cov.md` 的 30×/60× 是用户口头建议的补充工作点（G37 文档记录是
40×/80×）；本表覆盖 model.md 记录的全部 40×/80× × Q/L × P + MR 组。

## 合并实验（2026-08-15 补：老流程"合并后效果好"的 multik 对应）

老流程合并效果（`results/model.md` statMergeAnchors）：单组 anchors N50
7.3-25.4K → **7_merge_anchors 合并后 N50 55,041 / 17 条**。multik 的
对应实验：

| 方案 | N50 | Largest | mis | GF (%) | dup | mm/100k |
|---|---:|---:|---:|---:|---:|---:|
| 单组非 MR（Q25L60X40P000） | 24,445 | 56,183 | 0 | 95.86 | 1.001 | 27.7 |
| 单组 MR（MRX40P001，全组最优） | **54,888** | 178,776 | 0 | 96.11 | 1.000 | 30.7 |
| 全量 Q25L60（~150×，622K reads） | 14,536 | 50,525 | 0 | 95.57 | 1.001 | 24.5 |
| 23 组合并 `cat + contained`（老流程 7_merge 第一步） | 37,552 | 179,610 | 0 | 95.89 | **1.827** | 32.8 |
| 合并后 1× 重新压实（unitigs 当 reads，`--min-count-seed 1`） | 19,615 | 94,906 | **1** | 96.05 | 1.001 | 29.1 |
| **`anchr asm olc` 合并 23 组 multik 输出** | **54,964** | 148,437 | 0 | 95.86 | **1.000** | 33.4 |

**结论**：

1. **multik 单组 MR 已接近老流程合并水平**（N50 54.9K vs 7_merge_anchors
   55.0K）——merged reads 是 multik 在 G37 上的最佳输入，单组即达老流程
   合并后的长度水平（0 mis、GF 96.1% 还更高）；
2. **全量 reads 反而最差**（N50 14.5K、GF 95.57%）——高覆盖（~150×）下
   渐进过滤 cutoff 升高/噪音增多，"不是覆盖度越高越好"在 multik 上强成立，
   与老流程降采样动机一致；
3. **contained 合并**把非 MR 组 N50 从 ~20K 抬到 37.5K（0 mis），但
   **duplication 1.827**——23 组输出大量部分重叠，contained 只去完全包含；
4. **1× 重新压实**去掉了冗余（dup 1.001、GF 96.05% 全组最高）但 N50 降且
   引入 **1 个 misassembly**——无 reads 覆盖证据的 unitig 压实会错连；
5. **multik 复现"合并效果好"的正确路径**：多组 unitigs 反馈 + reads 支撑
   的重新压实（对应老流程"合并后重新 anchors 化用全量 reads"），需要
   multik 支持外部 unitigs 输入（每轮 `count_at` 加入外部 unitigs 反馈），
   当前 CLI 无此入口——记入 todo。

## `anchr asm olc` 合并（2026-08-15 补：用户建议，已验证）

用户建议直接用现有的 `anchr asm olc`（unitig → overlap → layout → cns）
合并多组 multik 输出，对应老流程"合并 = 经典 OLC"。命令：

```bash
cat 4_down_sampling/*P0*/multik.fa 6_down_sampling/*P0*/multik.fa \
    > multik_all_cat.fa
anchr asm olc multik_all_cat.fa \
    --min-count-seed 1 --min-overlap 34 --min-contig-len 500 \
    -o multik_olc_merged.fa
```

（`--min-count-seed 1`：multik 输出是 1× 无错 unitigs，k-mer 计数 1 即
solid。**输入预处理很关键**：contained 去重后（43 条）与全 cat
（16,584 条含短碎片）结果**不同**——短碎片会引入 misassembly 和 N50 崩，
必须先去重（见下表），对应老流程"先取可靠 anchors 再合并"。）

**结果**（quast，contigs ≥ 500 bp）：

| 指标 | asm olc（43 条输入） | asm olc（全 cat 16,584） | asm olc（MR-only 7 组） | 单组 MR P001 | 老流程 7_merge |
|---|---:|---:|---:|---:|---:|
| N50 | **54,964** | 31,730 | 48,850 | 54,888 | 55,041 |
| Largest | 148,437 | 56,787 | 130,840 | **178,776** | 179,685 |
| # misassemblies | 0 | **1** | 0 | 0 | — |
| Genome fraction (%) | 95.86 | **96.47** | 96.37 | 96.11 | — |
| Duplication ratio | 1.000 | 1.001 | 1.000 | 1.000 | — |
| # mismatches / 100 kbp | 33.4 | 36.1 | 35.6 | **30.7** | — |

**结论**：`anchr asm olc` 合并多组 multik 输出，N50 54.9K **与老流程
7_merge_anchors（55.0K）基本持平**，0 mis、dup 1.000、0 N——"多组 multik
输出 → asm olc 合并"是 multik 侧复现老流程"取可靠 anchors → OLC 合并"
的现成路径，无需新增 --unitigs 反馈机制（todo 该项降级为可选优化）。
**但注意（用户质疑点，2026-08-15）**：单组 MR（P001）在 mm（30.7）、
Largest（178.8K）、GF（96.11%）上仍优于 asm olc 合并（33.4/148K/95.86），
只有 N50 打平——**合并没有变得更好，为什么？** 分析见下节。

### 为什么单组 MR 比整套流程（23 组 + asm olc 合并）还好

1. **单组 MR 已接近 G37 组装的"天花板"**：G37 是单菌株、基因组 580K、
   重复区仅 ~3.6K（`1_genome/repetitive/`），组装简单；merged reads 读长
   长 + 40× 适中覆盖 → multik auto k 大 → 单组直接产出最长 178.8K、
   N50 55K 的长 unitigs。**asm olc 合并的 N50 54,964 与 MRX80P001 单组
   的 N50 54,964 相同不是巧合**：已比对合并输出的 contig_3（54,964）=
   MRX80P001 单组 unitig_3 的 reverse complement——**合并输出继承了 MR
   组的长 unitigs**，N50 打平是因为合并没有破坏它们，而不是合并变好了。
2. **合并引入的代价**：
   * **mm 升高（30.7 → 33.4/36.1）**：合并把多组 unitigs 的重叠区做
     consensus 拼接，不同组（不同 reads 子集）在同一区域有微小序列差异，
     拼接处产生错配；单组内部自洽，无跨组差异；
   * **Largest 变短（178.8K → 148K/57K/131K）**：多组 unitigs 引入更多
     ambiguous junction（重复区边界不一致、短 unitigs 的分支），layout 在
     ambiguous 处停止；全 cat 输入短碎片最多 → 最长只有 56.8K 且 1 mis；
   * **GF 是唯一稳定变好的指标**（全 cat 96.47%、MR-only 96.37% > 单组
     96.11%）——覆盖互补确实有效，但以 mm/N50/Largest 为代价。
3. **为什么老流程合并价值大、multik 合并价值小**：老流程单组产出是
   **可靠 anchors（覆盖筛选片段）**，N50 仅 7-25K，合并（OLC 拼装）把
   碎片拼成整条 → 55K，提升 2-7 倍；multik 单组产出已是**完整 unitigs**
   （MR 55K），再合并是"把多条整条再拼"——边际收益小、风险高（错配、
   ambiguous）。**老流程"先取可靠 anchors"的精髓 = 合并输入必须去重后
   的一致片段**：全 cat 直接 OLC（16,584 条）会因短碎片引入 1 mis 和
   N50 崩（31.7K），contained 去重后（43 条）才 0 mis、N50 54.9K。
4. **对 multik 的启示**：G37 这类单菌株数据，**multik 单组 MR 就是最佳
   工作点**；多组合并只做"覆盖互补"（GF +0.3-0.5pp），且必须先 contained
   去重再 OLC；需要更长 contig 时应从 reads 端入手（更长合并 reads、
   更高 k），而不是堆多组输出。

### 为什么之前"最长的 contig 变短"——用错了命令 + 撞名（2026-08-15 定位）

**现象**：单组 MR 输出最长 178,776（如 MRX40P001 的 unitig_1）；多组合并
（MR-only 7 组 / 全 cat 23 组）经 `anchr asm olc` 后最长掉到 130,840 /
56,787。**已排除"单组就被切"**：4 个单组分别过 asm olc，最长全部保持
（79567→79567、73192→73192、113997→113997、178776→178776）——切断只
发生在**多组混合**时。

**根因（用户追问后核对原始设计文档 `design/asm-olc.md` §5）**：

* `anchr asm olc` 驱动器的设计输入是 **reads**：内部 S0 对每个 k 跑
  `asm unitig` 生成伪 reads → S1 ovlp → S2 layout → S3 cns。独立命令
  `asm ovlp` / `layout` / `cns` 才是给"已经是 unitigs 的输入"用的管道
  （文档原文："支持用户自己跑 `asm unitig` → `asm ovlp` → `asm layout`
  → `asm cns` 的管道形态"）——**我把 multik 输出喂给 `asm olc` 驱动器，
  等于多做了一遍 S0 unitig 重组装**；
* 且把 23 个文件 cat 成单文件导致 **unitig_<id> 撞名**：设计文档 §S0
  的 `<tag>:<name>` 机制要求 tag = 输入文件 stem、多文件分别传入；
  `multik_mr_cat.fa` 里 `>unitig_186` 有 7 份（7 组同名不同序列），
  cat 单文件管道 cns 报 "overlapping bases disagree"。

**切断机制（`asm olc` 驱动器 S0 阶段，基于 `--keep-dir` 中间产物）**：

1. S0 把输入（1× unitigs）当 reads 在 k51 的 de Bruijn 图上重新走唯一
   路径；**单组输入**内部 51-mer 唯一（同一 reads 集自洽）→ 保持完整；
2. **多组输入**：其他组的 unitigs 与长 unitig 存在组间重叠 / 51-mer
   共享（同一基因组区域在不同组的组装边界不同，或重复序列），k51 图里
   长 unitig 内部出现第二个后继/前驱 → 分支 → 唯一路径**在分支处断开**；
3. 证据：MR-only 7 组时 MRX40P004 的 179,610 被切成 130,840
   （`k51:unitig_0`）+ 其他片段，130,840 末端与一个 93 bp 短 unitig
   （其他组）重叠 80 bp——断点引入分支的就是它；
4. 切断片段不会拼回（断点是分叉不是 overlap），组越多切越狠（MR-only
   7 组 → 130,840；全 cat 23 组 → 56,787）。

**正确用法**（多文件分别传入，不 cat）：

```bash
FILES=$(ls 4_down_sampling/*P0*/multik.fa 6_down_sampling/*P0*/multik.fa)
anchr asm ovlp $FILES -o all.ovlp.paf
anchr asm layout all.ovlp.paf $FILES -o all.layout.tsv
anchr asm cns all.layout.tsv $FILES -o all.contigs.fa
```

**正确管道结果**（quast ≥500 bp）：

| 方案 | N50 | Largest | mis | GF (%) | dup | mm/100k |
|---|---:|---:|---:|---:|---:|---:|
| 单组 MR P001（multik 直接输出） | 54,888 | 178,776 | 0 | 96.11 | 1.000 | 30.7 |
| 管道 MR-only 7 组 | 54,964 | **179,610** | 0 | 96.37 | 1.618 | 32.95 |
| 管道全部 23 组 | 39,104 | **179,610** | 0 | **96.44** | 2.659 | **31.46** |
| 管道 23 组 + 输出 contained | 39,104 | 179,610 | 0 | 95.90 | 1.899 | 31.99 |
| **`asm olc --unitigs`（改造后，dedup 0.99）** | **54,964** | 179,610 | 1* | **96.54** | **1.094** | 36.91 |
| （对比）asm olc 驱动器 23 组 | 31,730 | 56,787 | 1 | 96.47 | 1.001 | 36.05 |

\* 唯一 mis 是 560 bp、cov 1.4 的低覆盖碎片（`--min-contig-len 1000` 可滤）。
改造（`design/asm-olc.md` §14.4）：`--unitigs` 输入模式 + `--dedup-ratio`
近似包含去冗余，dup 2.659 → 1.094、GF 96.54%（全组最优）、N50 54.9K。

**结论**：

1. **正确管道不再切长序列**：最长 179,610 保持（单组里最长的 178,776 /
   179,610 都保留）；之前 asm olc 驱动器的切断 = S0 重组装被多组序列
   "污染"（共享 51-mer 处分支）——是**用错命令**造成，不是 OLC 本身；
2. **管道合并确实带来覆盖互补**：23 组 GF 96.44%（全组最优，比单组 MR
   96.11 高 0.33pp）、mm 31.46（全组最优）；MR-only 7 组 N50 54.96K 保持
   MR 单组水平且 GF 96.37%——多组合并的收益真实存在；
3. **dup 偏高（1.6-2.7）是未去冗余**：管道不做 S0 重组装，组间部分重叠
   的 unitigs 都保留；contained 只能去"完全包含"（dup 1.899）——彻底
   去冗余需要输入侧 contained 或老流程式重新 anchors 化；
4. **撞名是独立管道的必踩坑**：多组 unitigs 文件必须分别传入（tag =
   stem），cat 单文件会串序列（cns 报错是防御机制在起作用）。

### reads mapping → anchors → OLC（2026-08-15 用户建议，已验证）

用户建议："像老流程一样，用真实的 reads mapping 一遍得到 anchors，再对
anchors 做 OLC"（老流程 anchors = bbwrap 回贴 + [lower, upper] 覆盖区间
过滤，`references/anchr-legacy-pipeline.md` §2.4）。现代等价用
`anchr asm map`（完美匹配）+ `sam to-rg` + 自算每碱基覆盖度。

**单组验证（MRX40P001，40× 组自己的 reads）**：median 39、lower 7、
upper 114（40× 合理）→ 25 个 anchors（≥500 bp）、554.9K bp（95.7%）。
anchors 质量 vs 原 unitigs：

| 指标 | 原 unitigs | anchors |
|---|---:|---:|
| N50 | 54,888 | 54,841 |
| # misassemblies | 0 | 0 |
| Genome fraction (%) | 96.11 | 95.66 |
| # mismatches / 100 kbp | 30.7 | **27.21** |

**多组验证（7 组 MR → 各自 anchors → `asm olc --unitigs` 合并，201 条
anchors → 20 条 contigs）**：

| 方案 | N50 | dup | GF (%) | mis | mm/100k |
|---|---:|---:|---:|---:|---:|
| MR-only `--unitigs`（无 anchors） | 54,964 | 1.618 | 96.37 | 0 | 32.95 |
| MR-only **anchors → OLC** | 54,858 | **1.002** | 96.04 | **0** | **28.12** |
| 23 组 `--unitigs`（dedup 0.99） | 54,964 | 1.094 | 96.54 | 1* | 36.91 |

**结论（用户建议的价值）**：reads mapping + 覆盖度过滤（[lower, upper]）
同时解决了 `--unitigs` 的三大遗留问题——**dup 1.618 → 1.002**（upper 排除
重复区多版本）、**mm 32.95 → 28.12**（过滤掉高/低覆盖区的错配）、**560 bp
碎片 mis 消失**（lower 排除低覆盖碎片）。代价：GF 96.37 → 96.04（过滤掉
的问题区本身就是低质量区）。N50 保持 54.9K。

**集成建议**：把"reads mapping + 覆盖过滤 + 区间提取"封装进命令（如
`asm olc --anchor-reads <reads>`：内部先出 anchors 再 OLC），现代流程变为
`reads → multik（多组）→ anchors（覆盖筛选）→ asm olc --unitigs`——与
老流程"取可靠 anchors → OLC 合并"对齐。

### 为什么 cns 之后还有冗余（用户追问，2026-08-15）

**现象**：直接管道（23 组多文件）输出 dup 2.659。用户预期 OLC 的
consensus（C）应该把重叠的 unitigs 合并去冗余——为什么还有？

**数据证据**（all23.ovlp.paf / all23.layout.tsv）：

* ovlp 分类：138,627 条 **contain（ov:A:C，87%）** vs 20,016 条 dovetail
  （ov:A:D）——多组 multik 输出对同一基因组区域高度重叠，短 unitigs 被
  长 unitigs 包含是主流关系；
* layout：16,584 条输入 unitigs 里 **16,123 条是单步 layout**（各自成
  contig），只有 227 条进了多步链——contain 不参与 layout 延伸
  （`design/asm-olc.md` §S1："contain：q ⊂ t 或 t ⊂ q——不参与延伸，
  留作共识覆盖证据"）；
* 因此 v0 `cns`（精确缝合）只缝合 layout 路径上的 dovetail unitigs，
  **contain 的 unitigs 没有被吸收也没有被丢弃**——各自作为单步 contig
  原样输出 → 同一区域多条 contig → dup 2.659。

**这是 OLC 的 v0 设计缺口**：`design/asm-olc.md` §S3 明说"列投票留 v1：
若未来引入错配 overlap 或真实数据暴露 junction 不一致，再加 AS_CNS
`BaseCallMajority` 式逐列投票 + min-coverage 修剪"——v1 的 consensus
才应该把 contain 序列对齐到主 contig 参与投票并**丢弃冗余**；v0 精确
缝合没有这一步。

**缓解实验**（输入侧先 contained，多文件）：

| 方案 | N50 | Largest | mis | GF (%) | dup | mm/100k |
|---|---:|---:|---:|---:|---:|---:|
| 直接管道 23 组 | 39,104 | 179,610 | 0 | **96.44** | 2.659 | **31.46** |
| 输入 contained 后管道 | 39,189 | 179,610 | 0 | 95.89 | **1.201** | 33.69 |

输入侧 contained 把 dup 2.659 → 1.201，但 GF 96.44 → 95.89（-0.55pp）、
mm 31.46 → 33.69——被删的"被包含 unitigs"里有组间差异序列（近似包含、
边界/内容微差），删掉损失覆盖证据。**两全方案 = cns v1 的 contain 吸收**
（短序列对齐主 contig 投票，既去冗余又保留覆盖证据），记入 todo。
