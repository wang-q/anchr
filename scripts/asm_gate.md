# asm-gate baselines

> 门禁基准记录（`scripts/asm-gate.sh` 的唯一事实来源）。与
> `results/model_org.md` 的错装基线一致；此处只为分层门禁存机器可读基线。
>
> 复核：mis 以 quast 全比对为准（AGENTS.md 错装判定纪律）。

## smoke (L1)

G37 `MRX40P000` 整组 `asm multik --all-masters auto`（参考 580 kb、40×，
几秒 / <1 GB）。multik 确定性：同参数输出逐字节一致。

- golden-md5 `1dc52b1f4d1c989eb6497be07c0eacb0`
- count = 2745
- N50 = 37655
- total = 5936843

复捕：`bash scripts/asm-gate.sh smoke --write`（先把本文件 golden-md5 更新为现产物的 md5）。

## single (L2)

每数据集取 `MRX40P000` 组 multik→olc→extend，只看趋势（N50/count），
mis 不在此层判定（单组 mis 属预期，由多组 anchor 投票消解）。

计数口径：`count N50 total`（count=unitig 数、N50、total bp）。
暖色线：big N50/Total drop 或 unitig count 异常爆炸/塌缩 → 提示检查。

| dataset | stage | count | N50 | total(bp) | 说明 |
| ------- | ----- | ----: | --: | --------: | ---- |
| G37    | multik(unitigs_all) | 2745 | 37655 | 5936843 | 与 smoke 同源 |
| G37    | olc+extend (final) | 34 | 81699 | 587581 | 近整环 |
| MG1655 | multik(unitigs_all) | 8562 | 40990 | 46221276 | 参考 4.6 Mb |
| MG1655 | olc+extend (final) | 275 | 58128 | 4635367 | ~参考大小 |

## full (L3)

来自 `results/model_org.md` §end-multiplicity 门控全链门禁（quast
`--min-contig 10`，错装权威判定）。

| dataset | 组数 | #contigs | Largest | Total | N50 | #mis | GF% | Dup |
| ------- | --- | ------: | ------: | ----: | ---: | ---: | ---: | ---: |
| G37     | 7 | 10 | 187498 | 579008 | 121382 | 0 | 98.548 | 1.000 |
| MG1655  | 5 | 91 | 268281 | 4617679 | 112557 | 0 | 98.197 | 1.013 |
| DH5alpha| 13 | 105 | 258601 | 4496026 | 99473 | 0 | 97.800 | 1.003 |

> 用 `scripts/asm-gate.sh full <g37|mg1655|dh5alpha>` 复跑的实测值（quast
> report.tsv，`--min-contig 10`）会自动打印 #contigs/Total/N50/mis/GF，
> mis 为硬红线。工作目录默认 `/tmp/asm-gate-full/<ds>`（可 `FULL_WORK=`
> 覆盖），已完成的各组 anchor 会被复用，避免每次重跑 multik。
>
> 2026-08-20 L3 全链重跑（三数据集）：三者 #mis 均 0，MG1655 与
> DH5alpha 逐指标与上表一致；G37 现为 **10 contigs**（merge 用
> `--min-contig-len 1000` 丢弃早前未过滤纪录中包含的 3 条 <1000bp 尾部
> 子区 contig_11/12/13，共 2331 bp；N50 不变 121382）。该差异为测量口径
> 差异（gate 统一 merge 过滤 <1000bp），非代码回归，且合并基准与模板
> `9_quast.sh` 一致。

---

## gate history（由 `results/model_org.md` 迁入）

> 历史门禁记录，原为各数据集的"门禁"章节（multik 单次调用 / 性能优化 /
> end-multiplicity 门控全链）。`results/model_org.md` 只保留数据集定义与
> run（质量基线表格），此处保存算法改动验证与复现方法。

### g37: multik 单次调用（--all-masters）全流程门禁（2026-08-17 追加）

`asm multik --all-masters`（auto 阶梯：N50 408 → 31..192）单次调用，
7 组 MR anchor 合并口径，guide 与否结果完全相同：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 单次调用 auto（guide 与无 guide 一致） | 58 | 179800 | 580756 | 55170 | 0 | 97.943 | 1.003 | 0.00 | 108.81 | 28.78 |
| 08-16 旧链基线（merge_multik） | 15 | 179712 | 563707 | 55098 | 0 | 96.997 | 1.000 | 0.00 | 76.44 | 19.20 |

要点：
* **0 mis 保持**，N50 与旧 multik 链持平（55.2K vs 55.1K），GF +0.95；
* guide 与无 guide 输出完全一致，与 MG1655 结论互相印证；
* 本流程为快速门禁复现（手写 anchor 合并脚本），contigs 数多于模板
  7_merge 链属预期；单组验证 N50 121K（auto 31..192）。

### g37: multik 性能优化门禁（2026-08-17 晚追加）

walk 滚动窗口 + succ 索引（删除 400 MB HashMap index）、
`remove_unsupported` 按 unitig 并行、pass0 与 rounds 并发（详见
`notes/design/asm-multik.md` §性能）。multik 输出与优化前字节级一致，
全链复跑：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | -------: | ---------: |
| 性能优化后（auto 31..192） | 58 | 179800 | 581098 | 55170 | 0 | 97.979 | 107.51 | 28.59 |
| 08-17 门禁记录（同口径） | 58 | 179800 | 580756 | 55170 | 0 | 97.943 | 108.81 | 28.78 |

* **0 mis 保持**，逐指标与 08-17 门禁一致（GF +0.036 源于其间合入的
  `remove_unsupported` run>=2 过剪修复，非本次性能改动）；
* 单组 multik 4.2–6.5 s（-p 8、2 组并发），峰值 ~2.0 GB/进程。

### g37: end-multiplicity 门控全链门禁（2026-08-20 追加）

针对 Q25L60X80P001 unitig_411 的 end-multiplicity 图级门控（共享
(k-1)-mer ≥3 末端阻断压缩/融合，`graph.rs`；详见
`notes/design/asm-multik.md` §2026-08-20）合入后的 G37 全链门禁：7 组 MR
multik(--all-masters auto 31..160) → olc → extend → anchor → 跨组
`asm olc --unitigs` merge，QUAST `--min-contig 10`（同模板 9_quast.sh）。

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 08-16 merge_mr_multik 基线 | 17 | 187458 | 586244 | 121369 | 0 | 98.666 | 1.005 | 0.00 | 241.36 | 70.08 |
| **本次 end-multiplicity（修复后）** | 13 | 187498 | 581339 | 121382 | 0 | 98.674 | 1.001 | 0.00 | 238.19 | 70.15 |

* **0 mis 保持**；N50 +13 bp（+0.01%）、GF +0.008pp，无回归；
* Total 581.3K vs 基线 586.2K（−0.8%）、Dup 1.001 vs 1.005：门控切开
  少量 3-末端重复区错误融合，冗余更低、更贴近参考大小；
* 单级门禁（Q25L60X80P001 同组口径）4 mis 与修复前持平——单级 unitig
  层 mis 由下游多组 anchor 投票消解，与全链 0 mis 不矛盾。

### mg1655: multik 单次调用（--all-masters）门禁（2026-08-17 追加）

`asm multik --all-masters` 单次调用取代模板 per-master 循环（k-major
顺序共享每 k 计数表、rounds 跨 master 并行、计数表即建即弃）。同 5 组
anchor 合并口径，三个变体与 08-16 基线对比：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 08-16 基线（per-master 链） |  90 | 268333 | 4584494 | 123988 | 0 | 98.523 | 1.002 | 0.00 | 0.22 | 0.09 |
| 单次调用 guide 31..192 | 181 | 174306 | 4603907 |  79617 | 0 | 98.880 | 1.003 | 0.00 | 0.83 | 0.07 |
| 单次调用 无guide 31..192 | 182 | 174306 | 4604139 |  79625 | 0 | 98.877 | 1.003 | 0.00 | 0.83 | 0.07 |
| **单次调用 auto（模板默认）** | 178 | 174306 | 4603590 |  79617 | 0 | 98.852 | 1.003 | 0.00 | 1.80 | 0.07 |

要点：
* **三变体全部 0 mis**；guide 与否逐指标几乎相同（GF 差 0.003）→ 5 组
  anchor 投票下 guide 无贡献，模板移除 `--use-guide`（速度 ~2.1×）；
* **单组对照无回归**：新旧同 k（31..128）单组质量完全一致（均 2 mis、
  GF 99.40、N50 112,514），重构本身质量中性；
* **N50 79.6K vs 08-16 链 124.0K**：~~归因为验证密度~~（08-17 下午证伪，
  实为跨组 olc 缺 `--unitigs` 的测量错误，见下节）；
* 计时（release、-p 8、单组）：auto 140–172 s（无 guide）vs guide
  ~360 s vs 旧 per-master 串行 7:23（-p 4、9 次全量计数）；
* 内存：计数表即建即弃后峰值 6.6 GB（40×组）/ 10.8 GB（80×组）；
  曾因缓存全部 K 张表达 26.6 GB 并在 5 组并发时 OOM，已改；
* auto 阶梯改为固定梯 `31..192` 截断于 `clamp(N50/2, 81, 192)`
  （本数据 N50 339 → 31..160）。旧公式 0.8×N50 给出 51..251，高 k
  master 被残余错误打碎（N50 9.4K、5 mis），不可用。

### mg1655: multik 性能优化门禁（2026-08-17 晚追加）

`--all-masters` 单次调用的三个性能改动（walk 滚动 fw/rc 窗口 + 分类期
预联唯一后继下标、`remove_unsupported` 按 unitig 并行、pass0(k) 与
earlier-masters rounds(k) rayon::join 并发），multik 输出与优化前
**字节级一致**，全链复跑（跨组 olc 带 `--unitigs`，同上节口径）：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | -------: | ---------: |
| 性能优化后（auto 31..160） | 90 | 268842 | 4624071 | 118731 | 0 | 99.072 | 2.06 | 0.06 |
| 上节基线（优化前同口径） | 90 | 268842 | 4624071 | 118731 | 0 | 99.072 | 2.06 | 0.06 |

* 逐指标完全相同（multik 字节级一致 → 下游全同），**0 mis 保持**；
* 单组计时（release、-p 8、单进程）：105.8 s → **48.9 s**（2.2×）；
  5 组链 2 并发下每组 59–65 s；峰值内存不变（~11.7 GB/进程）；
* 分解：DFA walk 7.3→1.9 s（k=160 隔离）、classify 1.5→1.1 s（删
  HashMap 构建）、46 轮 `remove_unsupported` CPU 104.9 s 转入并行。

### mg1655: end-multiplicity 门控全链门禁（2026-08-20 追加）

针对 Q25L60X80P001 unitig_411 的 end-multiplicity 图级门控（共享
(k-1)-mer ≥3 末端阻断压缩/融合，`graph.rs`；详见
`notes/design/asm-multik.md` §2026-08-20）合入后的 MG1655 全链门禁：
标准 5 组 MR（MRX40P000/P001/P002 + MRX80P000/P001）multik(--all-masters
auto 31..160) → olc → extend → anchor → 跨组 `asm olc --unitigs` merge，
QUAST `--min-contig 10`（同模板 9_quast.sh 口径）。

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 08-17 性能优化门禁基线（auto 31..160） | 90 | 268842 | 4624071 | 118731 | 0 | 99.072 | 1.005 | 0.00 | 2.06 | 0.06 |
| 08-18 现行复现（无 cv，c9da0ce 合入后） | 107 | 268272 | 4632159 |  95706 | 0 | 98.325 | 1.001 | — | — | — |
| **本次 end-multiplicity（5 组，修复后）** | 91 | 268281 | 4617679 | 112557 | 0 | 98.197 | 1.013 | 0.00 | 0.00 | 0.04 |

要点：
* **0 mis 保持**；N50 112.6K、GF 98.197 与 08-18 现行复现（c9da0ce
  合入后 5 组无 cv：N50 110.5K、GF 98.346）同代吻合，低于 08-17 基线
  （118.7K、GF 99.072）——该差归因 c9da0ce（DH5alpha relocation 修复）
  合入后 extend 跨 contig 护栏，非本次端门控；门控在 MG1655 上不新增
  回归；
* 单级门禁（Q25L60X80P001 组）4 mis 与修复前持平——单级 unitig 层
  mis 由多组 anchor 投票消解，与全链 0 mis 不矛盾；
* 13 组全跑口径 am：N50 117.8K、GF 98.206、0 mis（contigs 数更多），
  与 5 组口径结论一致。

### dh5alpha: multik 单次调用（--all-masters）全流程门禁（2026-08-18 追加）

承接 `g37/mg1655: olc --cross-validate 跨组嵌合投票`后的第三条基准
数据集验证。`asm multik --all-masters`（auto 阶梯按 reads N50 生成，
本数据 N50 250 → 31..160）单次调用，13 组 MR anchor 合并口径，
跨组 `olc --unitigs --cross-validate`：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 13 组 anchor + cross-validate + extend 跨 contig 护栏（现行） | 122 | 259019 | 4581642 |  82991 | 0 | 98.342 | 1.016 | 0.00 | 1.27 | 0.20 |
| 13 组 anchor + cross-validate（消除前，2 mis） | 105 | 259450 | 4619600 | 112966 | 2 | 98.848 | 1.019 | 0.00 | 2.71 | 0.24 |
| spades（旧运行） | 78 | 132337 | 4490000 | 112448 | — | — | — | — | — | — |
| mr_spades（旧运行） | 59 | 178373 | 4510000 | 132590 | — | — | — | — | — | — |
| mr_megahit（旧运行） | 70 | 133730 | 4520000 | 132754 | — | — | — | — | — | — |

要点（单组计时与资源）：
* 13 组 MR（MRX40P000-008 + MRX80P000-003）每组合计 multik 20–33 s
  （-p 8、2 组并发）、olc ~75–95 s、extend ~2 s、anchor ~2–3 s；
  峰值内存 8.4–14.9 GB/进程（40× 组 ~8.6–10.5 GB、80× 组
  ~10.8–14.9 GB），2 组并发峰值 ~30 GB（< 机器 88 GB 的 1/2）；
* 跨组 `olc --unitigs --cross-validate` 96.4 s / 820 MB，extend 1.4 s；
  QUAST 用 `quast.py -m 500 -r genome.fa --min-contig 200`。

### dh5alpha: end-multiplicity 门控全链门禁（2026-08-20 追加）

针对 Q25L60X80P001 unitig_411 的 end-multiplicity 图级门控（共享
(k-1)-mer ≥3 末端阻断压缩/融合，`graph.rs`；详见
`notes/design/asm-multik.md` §2026-08-20）合入后的 DH5alpha 全链门禁：
13 组 MR（MRX40P000-008 + MRX80P000-003）multik(--all-masters auto
31..160) → olc → extend → anchor → 跨组 `asm olc --unitigs` merge
（无 cv，模板默认），QUAST `--min-contig 10`（同模板 9_quast.sh）。

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 08-18 门禁 cv + 护栏（cv 口径） | 122 | 259019 | 4581642 |  82991 | 0 | 98.342 | 1.016 | 0.00 | 1.27 | 0.20 |
| 08-18 门禁 13 组 merge（无 cv） | 117 | 258601 | 4505173 |  99473 | 0 | 97.942 | 1.004 | — | — | — |
| **本次 end-multiplicity（13 组无 cv，修复后）** | 105 | 258601 | 4496026 |  99473 | 0 | 97.800 | 1.003 | 0.00 | 0.18 | 0.16 |

要点：
* **0 mis 保持**；N50 99.5K 与 08-18 门禁无 cv 记录逐位一致（99,473），
  Largest 一致（258,601），GF -0.14 pp、Total -9K（contigs 105 vs 117，
  分组投票消元更干净）——端门控在 DH5alpha 上无新增回归；
* 对位取**无 cv** 口径：cv（`--cross-validate`）在现行 anchor 上会假拆
  真连接（N50 82.9K，-17%，见 `notes/design/asm-multik.md` §跨组嵌合），
  且 AGENTS.md 错装判定以 quast 为准；无 cv 0 mis 即通过门禁；
* 单级门禁（Q25L60X80P001 组）4 mis 与修复前持平——单级 unitig 层
  mis 由下游多组 anchor 投票消解，与全链 0 mis 不矛盾。

### bcer: 现代流程门禁（2026-08-20 已跑）

* 数据已就位（本地 R1/R2 + ref），用 `/tmp/bcer_full` 完整模板链初跑
  （`--unitigger "multik unitig"`，非弃用的 `bcalm` 写法）；
* **首个革兰氏阳性 / 低 GC（~38%）/ 多复制子（染色体 + pBc10987 质粒）**
  视角：multik 相关组装（merge_multik / merge_mr_multik / merge_unitig /
  merge_mr_unitig / merge_anchors）**全部 0 mis**；`merge_anchors`
  N50 62,043、GF 98.517、Dup 1.018；`merge_mr_multik` N50 44,305、
  GF 98.488。对比组装器 megahit 7 mis / mr_megahit 1 mis；
* 对照老流程基线（见 `results/model_org.md` §bcer statQuast）：多复制子上
  Dup>1 由质粒/重复区正常放大；无 cv（单样本，非共组装场景）；
* 工作目录 `/tmp/bcer_full`（quast `--min-contig 100`，同模板口径）。

### rsph: 现代流程门禁（待跑）

* 需先下载 reads（ccb.jhu.edu），再按 template/run 节命令执行；
* 关注点：7 复制子（2 染色体 + 5 质粒）下 GF/dup 统计会被质粒放大，
  QUAST 用 `--no-check`（老流程即如此），mis 判定看染色体级骨架；
  `--cross-validate` 跨组投票对高重复（12.8%）菌株是核心验证。