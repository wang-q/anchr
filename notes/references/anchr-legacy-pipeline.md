# 老手工流程理解（G37 例子，2026-08-14/15，基于完整代码阅读）

> 本文记录我对 `results/model.md` 的 G37 手动模拟流程的理解。**每个结论
> 都来自完整阅读对应代码**：流程模板 `templates/*.tera.sh`（生成实际脚本
> 的源）、命令实现 `src/cmd/*.rs`、以及 G37 产物脚本
> `~/data/anchr/g37/`（`0_script/*.sh` 和各目录的 `unitigs.sh`/
> `anchors.sh`/`quorum.sh`）。
>
> **2026-08-15 修正（用户纠正）**：之前我把流程理解成"一部分 reads 从前往后
> 组装得到一个东西"（单条线），漏掉了老流程的主干——**把 reads 拆分成多个
> 覆盖度部分，每个部分分别组装，最后用经典 OLC（overlap → layout）拼装
> 起来**。本文 §2.2/§2.6 已重写；§3 关键结论已更新。

## 1. 流程总览（`0_master.sh` 调用顺序 + `model.md` template 参数）

`model.md` 的 `anchr template` 参数（决定哪些脚本生成）：

```text
--genome 580076 --parallel 8 --xmx 12g --repetitive
--fastqc --insertsize --fastk
--trim "--dedupe --cutoff 30 --cutk 31" --qual "25 30" --len "60"
--filter "adapter artifact" --quorum --merge --ecphase "1 2 3"
--cov "40 80" --unitigger "bcalm bifrost superreads tadpole"
--statp 2 --readl 125 --uscale 2 --lscale 3 --redo
```

**老流程的主干结构**（用户 2026-08-15 纠正，详见 §2.2/§2.6）：

```
reads（全量）
  → 拆分成多个覆盖度部分（4_down_sampling X40/X80 × P 副本；6_down_sampling MR 版）
  → 每个部分分别组装（unitigs：4 工具 × 6 k → contained）
  → 每个部分出 anchors（bbwrap 回贴 + 覆盖度区间过滤）
  → 所有部分的 anchors 合并（7_merge_anchors：contained + 重新 anchors 化）
  → 经典 OLC 拼装（7_extend/7_glue/7_fill：overlap2 → group → layout）
  → 对照路线 8_spades/8_megahit（全量 reads 直接组装）
  → 9_stat_final / 9_quast 比较全部产物
```

`0_master.sh` 按序调用（存在才跑）：

```
2_fastqc / 2_insert_size / 2_fastk / 2_trim / 9_stat_reads
2_merge（mergeread + ecphase） / 2_quorum
4_down_sampling → 4_unitigs_{bcalm,bifrost,superreads,tadpole}
    → 4_anchors → 9_stat_anchors
6_down_sampling → 6_unitigs_* → 6_anchors → 9_stat_mr_anchors
7_merge_anchors（4/6 × 4 unitigger + 最终 7_merge）
8_spades / 8_mr_spades / 8_megahit / 8_mr_megahit
9_stat_other_anchors / 9_stat_final / 9_quast
```

> 注：`0_master.sh` **不包含** 7_extend/7_glue/7_fill——它们是
> `templates/0_bsub.tera.sh`（集群版模板）里 `--extend 1` 时自动衔接的
> 手动阶段（G37 本地流程里手动执行，产物存在，9_stat_final/9_quast 引用）。

## 2. 各阶段（依据：模板 + 产物脚本）

### 2.1 reads 准备

* **trim**（`templates/2_trim.tera.sh`）：
  `anchr trim --dedupe --cutoff 30 --cutk 31 --qual "25 30" --len "60"
  --filter "adapter artifact"`；**输入前缀 R/S/T**（模板 `for PREFIX in
  R S T`；G37 只有 R = 原始 PE，S/T 未用）。`2_illumina/Q{1}L{2}/` 是
  `trim/Q{1}L{2}/` 的符号链接（Q0L0/Q25L60/Q30L60 = 质量/长度档位）。
* **merge**（`templates/2_merge.tera.sh`）：`anchr mergeread <R1 R2 Rs>
  --ecphase "1 2 3"`，前缀映射 R→M/U、S→N/V、T→O/W；产物在
  `2_illumina/merge/`（statMergeReads 的 clumped/ecco/eccc/ecct/extended/
  merged.raw/unmerged.raw/M1/U1/U2/Us/M.cor）。**ecphase 在 mergeread 里**。
* **quorum**（`templates/quorum.tera.sh`）：`anchr quorum <R1 R2 [Rs]>`
  对每个 trim 档位跑；**关键语义（用户纠正）**：
  `quorum_error_correct_reads` 给修正过的 read 头追加 `:sub:`/`trunc`
  标记，脚本把带标记的 reads 全部丢弃（`fa some -i discard.lst`）——所以
  **pe.cor.fa 的内容 = 未被修正的原始序列**（被 quorum 改过的都丢了），
  "cor" 只是 quorum 输出文件后缀（`{{prefix}}.cor.fa`），不是 corrected。
  `trim/pe.cor.fa.gz` 是全量（未降采样）的这份输出。

### 2.2 降采样（`templates/4_down_sampling.tera.sh` / `6_down_sampling.tera.sh`）

**这是"把 reads 拆分成多个部分"的一步**（不是一条线）：

* **4_**：输入 `2_illumina/Q{1}L{2}/pe.cor.fa.gz`，`pgr fa split about -c
  genome×cov` 按目标覆盖切分 → `4_down_sampling/Q{1}L{2}X{40/80}P{000...}/
  pe.cor.fa`。组合 = **3 个 trim 档位（Q0L0/Q25L60/Q30L60）× 2 个覆盖度
  （X40/X80）× P 随机副本**（P 数量由 statp=2 约束）——每个部分是一份
  独立 reads 子集。
* **6_**：输入 `2_illumina/merge/pe.cor.fa.gz`（merged reads 的 quorum
  输出）→ `6_down_sampling/MRX{40/80}P{...}/pe.cor.fa`（MR 版本，X40/X80
  × P 副本）——**merged reads 是另一条平行的 reads 拆分线**。
* **覆盖度选择（model.md 记录为准，2026-08-15 核对）**：G37 当时文档
  （`results/model.md`）记录的就是 **40×/80×**（template `--cov "40 80"`，
  statQuorum/statUnitigs 表格全部为 X40/X80）——**不是 30×/60×**。
  用户口头提过"30×/60× 比较合适"是另一处的经验（2026-08-14 讨论，已用
  30×/60× 做了 multik 初验，见 `benchmarks/multik-cov.md`），但 G37
  model.md 的事实是 40×/80×。
* **Q/L 档位也切分**：`4_down_sampling` 对 **每个 trim 档位**（Q0L0 /
  Q25L60 / Q30L60）都做 X40/X80 降采样（+ P 随机副本），`results/model.md`
  statQuorum 表记录了全部 14 行（如 Q25L60X40P000/P001/P002、
  Q30L60X80P000）；MR 版本（merged reads）同样 MRX40/MRX80 × P。

每个部分**独立走完 unitigs → anchors**（§2.3/§2.4），最后才在
§2.5/§2.6 合并。

### 2.3 unitigs（`templates/4_unitigs.tera.sh` + `templates/unitigs.tera.sh`）

* `anchr unitigs <pe.cor.fa> <env.json> -u {bcalm|bifrost|superreads|tadpole}
  --kmer '31 41 51 61 71 81' -o unitigs.sh`，然后 `bash unitigs.sh`；
* **unitigs.sh 内容**（`templates/unitigs.tera.sh`）：逐 k 生成 unitigs
  （bcalm/tadpole/bifrost/superreads 各用自己的命令）→ `anchr contained
  unitigs_K*.fasta --len min --idt 0.9999 --ratio 0.99999` 合并 6 个 k；
  若 `--merge` 选项开（G37 未开）：再 `orient --idt 0.999` + `merge
  --idt 0.9999` → unitigs.fasta；否则 contained 后直接输出。

### 2.4 anchors（`templates/4_anchors.tera.sh` + `templates/anchors.tera.sh`）

`anchr anchors <unitigs.fasta> <pe.cor.fa> --readl 125 --uscale 2 --lscale 3
-p 8 -o anchors.sh && bash anchors.sh`。**anchors.sh 内容**
（`templates/anchors.tera.sh`，`anchr anchors` 的 about 是 "Select anchors
(proper covered regions) from contigs"）：

1. unitigs → `UT.fasta`；**bbwrap perfectmode**（`maxindel=0
   strictmaxindel perfectmode ambiguous=all`）把 reads 映射回 UT.fasta →
   `basecov.txt`（每碱基覆盖度）；
2. 覆盖统计 → 阈值 `lower = max(mincov=5, (median − mscale×MAD)/lscale)`、
   `upper = (median + mscale×MAD)×uscale`（mscale=3/lscale=3/uscale=2）；
3. 保留覆盖度 ∈ [lower, upper] 的碱基（低覆盖=错误区、高覆盖=重复区都
   排除）→ 连续覆盖区 `contig.covered.txt`；
4. **补洞**：覆盖 ≥ratio 且 ≥min 的 contig 填全部洞（`fill×10` →
   excise min），其余填小洞（fill）→ union → `contig.proper.json`；
5. `spanr convert` → regions；`pgr fa range UT.fasta -r ...` **从 unitigs
   按 region 提取** → `pe.anchor.fa`；others = 补集 → `pe.others.fa`；
6. `anchr contained pe.anchor.fa --idt 0.9999 --ratio 0.99999` → `orient
   --idt 0.999` → `merge --idt 0.9999` → `contained --idt 0.98 --ratio 0.99`
   → `anchor.fasta`。

**结论**：anchors = unitigs 中被 reads 良好覆盖（[lower, upper] 区间、
滤低覆盖错误区和重复区）的连续片段，从 unitigs 提取后 contained/
orient/merge。**不是 reads 延伸**。

### 2.5 merge anchors（`templates/7_merge_anchors.tera.sh`）

* **`--redo` 分支**（G37 model.md 带 `--redo`，走此路径）：
  1. `anchr contained $(find ... -name anchor.fasta | sort -r)` 合并所有
     unitigger × 降采样的 anchor.fasta（Q30L60X80 优先）→
     `anchor.non-contained.fasta`；
  2. **重新 anchors 化**：`anchr anchors anchor.non-contained.fasta
     ${BASH_DIR}/../2_illumina/trim/pe.cor.fa.gz --keepedge --ratio 0.98`
     ——输入是 trim 目录的全量 quorum 筛选后 reads（脚本 log_info 写
     "anchors with Q0L0 reads" 但实际路径是 `trim/pe.cor.fa.gz`）→
     `anchor.merge.fasta`；
  3. others（pe.others.fa）也 contained 合并（--len 500）。
* **`--redo 0` 分支**：直接 `orient --idt 0.999` + `merge --idt 0.9999` +
  `contained --idt 0.98 --ratio 0.99`（不重新 anchors 化）。

### 2.6 最终组装（`templates/8_spades.tera.sh` 等）

* `pgr fa filter --min-len 60 ${DIR_READS}/pe.cor.fa.gz | repair.sh` →
  re-pair（R1/R2/Rs）；
* `spades.py -t 8 --only-assembler -k 21,33,55,77 -1 re-pair/R1.fa
  -2 re-pair/R2.fa`（**直接用 reads 组装，不是 anchors 输入**）；megahit
  同理；
* `contigs.fasta` → `anchr contained --len 1000 --idt 0.98 --ratio 0.99999
  | faops filter -a 1000` → `spades.non-contained.fasta`；
* 再 **anchors 化**：`anchr anchors spades.non-contained.fasta
  pe.cor.fa.gz --ratio 0.98` → `8_spades/anchor/anchor.fasta`。

### 2.7 经典 OLC 拼装（`7_extend` / `7_glue` / `7_fill`，2026-08-15 补）

这是用户纠正的核心："每一个部分拼接完之后，**再用经典的 OLC 把它们拼装
起来**"。调用链在 `templates/0_bsub.tera.sh`（`--extend 1` 分支），本地
G37 手动执行（产物在 `7_glue_anchors/`、`7_fill_anchors/`，
9_stat_final/9_quast 都引用）：

```bash
# 1. 长序列集（FILE_LONG）：8_* 各组装器的 contigs 合并（≥1000 bp）
cat 8_spades/spades.non-contained.fasta \
    8_megahit/megahit.non-contained.fasta \
    8_mr_spades/spades.non-contained.fasta \
    8_mr_megahit/megahit.non-contained.fasta \
    | anchr dazzname --no-replace stdin \
    | pgr fa filter --min-len 1000 stdin -o 7_extend_anchors/contigs.2GS.fasta

# 2. glue：anchors × long contigs 做 OLC，输出 contig.fasta
bash 0_script/7_glue_anchors.sh 7_merge_anchors/anchor.merge.fasta \
    7_extend_anchors/contigs.2GS.fasta 3

# 3. fill：glue 的 contig × long contigs 再 OLC（填 gap）
bash 0_script/7_fill_anchors.sh 7_glue_anchors/contig.fasta \
    7_extend_anchors/contigs.2GS.fasta 3
```

**fill/glue 脚本内部 = 标准 OLC 三段**（`templates/7_glue_anchors.tera.sh`
与 `7_fill_anchors.tera.sh` 相同骨架）：

1. **Overlap**：`dazz overlap2 ${FILE_ANCHOR} ${FILE_LONG} -b 50 --len 1000
   --idt 0.999 --all`（anchor × long 两两比对）→ `anchorLong.ovlp.tsv`；
   glue 额外做 anchor 内部 overlap（`anchr overlap --serial --len 30
   --idt 0.9999`，只保留每端度数 ≤2 的）→ `anchor.ovlp.tsv`；
2. **Group**：`dazz group anchorLong.db anchorLong.ovlp.tsv
   [--oa anchor.ovlp.tsv] --range 1-N --len 1000 --idt 0.999
   [--max 100 | --max -30] -c GAP_COV`——按 anchor 分组建图（-c = gap
   覆盖度阈值）；
3. **Layout（每组内）**：`anchr orient`（方向一致化）→ `anchr overlap`
   （组内精确 overlap，`--len 1000 --idt 0.9999`）→ `anchr restrict`
   （按 restrict.tsv 约束）→ **`dazz layout`** → `group/*.contig.fasta`；
   汇总 `non_grouped.fasta + *.contig.fasta` → `contig.fasta`
   （`faops filter -a 1000`）。

**fill vs glue 的差异**：glue 用 anchor 内部 overlap 辅助分组（`--oa`）、
`--max -30`（允许每个 anchor 进多个组）；fill 用 `--keep --max 100`，
`-c 2`（bsub 版统一传 3）。fill 的输入 anchor = glue 的 contig——即
**先 glue 把 anchors 和长 contigs 拼起来，再 fill 用长 contigs 填 gap**。

**G37 本地产物与 bsub 版参数不同**：本地 `7_fill_anchors/anchor.fasta`
（16 条，md5 ≠ `7_merge_anchors/anchor.merge.fasta` 的 17 条）说明手动
实验调整过输入，但流程结构一致。

## 3. 关键结论（修正我此前的理解）

0. **老流程是"多覆盖度拆分 → 各部分组装 → OLC 拼装"**（用户纠正，
   2026-08-15）：reads 按覆盖度拆成多个部分（G37：3 trim 档 × X40/X80 ×
   P 副本 + MR 版；用户经验 30×/60× 合适，不是越高越好），各部分独立
   unitigs → anchors，合并后再用经典 OLC（overlap2 → group → layout）
   拼装（§2.2/§2.7）——**不是"一部分 reads 从前往后组装得到一个东西"**；
1. **anchors ≠ reads 延伸**（用户纠正）：anchors = unitigs 的覆盖质量筛选
   片段（bbwrap 回贴 + [lower, upper] 覆盖过滤 + 补洞 + pgr fa range 提取）；
2. **pe.cor.fa 序列未纠错**（用户纠正）：quorum 修正过的 reads 带
   `:sub:`/`trunc` 标记且被丢弃，pe.cor 保留未被修正的原始序列；"cor" 是
   quorum 输出后缀，不代表 corrected；
3. **8_* 不是"以 anchors 为输入"**：spades/megahit 直接用 quorum 筛选后的
   全量 reads（trim/pe.cor.fa.gz 或 merge/pe.cor.fa.gz）组装，contigs 再
   anchors 化；
4. **7_merge 的重新 anchors 化**：`--redo` 分支用 trim/pe.cor.fa.gz +
   `--keepedge --ratio 0.98`（不是 Q0L0——log_info 消息有误导）；
5. **statFinal / quast 是多条路径的对比**：unitigs→anchors→merge→
   glue→fill（7_*）vs reads→spades/megahit→contigs→anchors（8_*），
   mr_* = merged reads 版本；**7_extend 的 2GS contigs（8_* 合并）同时是
   glue/fill 的 FILE_LONG**——7_* 与 8_* 不是孤立的两条线，glue/fill 把
   它们合流；
6. **multik 的正确对比对象**是老流程的最终产物集合（quast 里的
   merge_anchors / glue_anchors / fill_anchors / spades / mr_spades /
   megahit / mr_megahit），不是中间 unitigs；
7. **G37 文档记录的是 40×/80×**（2026-08-15 核对 `results/model.md`）：
   覆盖度不是越高越好，高低两个覆盖度各司其职（低覆盖 = 少重复冲突、
   高覆盖 = 补低覆盖区的洞），这是"拆分多个部分"的动机；30×/60× 是用户
   口头提过的另一经验值（multik 初验见 `benchmarks/multik-cov.md`），
   不要与 model.md 的 40×/80× 混记。

## 4. 待讨论（我的初步看法，非结论）

* multik 的渐进丰度过滤只做"低覆盖排除"，anchors 的 `upper`（高覆盖 =
  重复区排除）没有对应——老流程的"多覆盖度部分"（30×/60×）本质也是这个
  思想的另一形态：低覆盖部分天然避开重复区冲突；multik 若按覆盖度拆多个
  部分跑再合并，或引入覆盖度区间过滤，等价于老流程的 anchors 语义；
* 老流程把"每个部分"的 unitigs → anchors → **OLC 拼装（overlap2 → group
  → layout）**；multik 的 unitig 图跨接验证 + 链压实是"OLC 拼装"的
  图论替代（边有 k-mer 计数证据，而非仅 overlap 长度）——两套语义的
  对应关系值得展开（§2.7 的 dazz layout vs multik merge_chains）；
* merged reads（MR）路径对 multik 的意义（multik 直接吃 merged reads 的
  可行性）；
* 7_extend 的 2GS contigs（8_* 合并）作为 glue/fill 的 long 输入：multik
  输出作为 FILE_LONG 与老流程 anchors 做 OLC 拼装的混搭可行性。
