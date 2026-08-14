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
  → 每个部分取可靠的 anchors（bbwrap 回贴 + 覆盖度区间过滤，排除错误区/重复区）
  → 所有部分的可靠 anchors 合并（7_merge_anchors：contained + orient/merge 或
    重新 anchors 化）——**合并本身是经典 OLC**（overlap → layout，见 §2.5）
  → fill/glue 辅助：用 8_* 的长 contigs（2GS）填充/粘合 gap（非主 OLC，见 §2.7）
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

### 2.5 merge anchors（`templates/7_merge_anchors.tera.sh`）——**经典 OLC 的合并**

**用户纠正（2026-08-15）**："先取其中可靠的 anchors，然后再合并，这个合并
其实非常像经典的 OLC"——**合并（7_merge_anchors）才是老流程的 OLC 主步骤**，
fill/glue 不是（见 §2.7）。OLC 语义落在 `anchr contained/orient/merge` 三个
命令上（都基于 `anchr overlap`）：

* `anchr merge`（"Merge overlapped unitigs"）实现 = **overlap → 有向图 →
  拓扑排序去环 → 按序拼接**（`src/cmd/merge.rs:125-260`：`anchr overlap`
  产出所有重叠 → 重叠关系建 DiGraphMap（边权 = append length）→ 环删除、
  分支保留 → 拓扑序按序把后序序列的非重叠部分接上）——这就是 OLC 的
  **overlap + layout**（consensus 简化：高 idt 0.9999 下直接按 overlap 拼接）；
* `anchr orient` = 重叠序列方向一致化（layout 前置）；
* `anchr contained` = 丢弃被包含的片段（overlap 关系去冗余）。

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
* **`--redo 0` 分支**：`contained → orient --idt 0.999 → merge --idt 0.9999
  → contained --idt 0.98 --ratio 0.99`（**不重新 anchors 化**）——这一串就是
  完整的 OLC 式合并：overlap（contained/merge 内部）→ layout（orient +
  拓扑拼接）→ 去冗余。

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

### 2.7 gap 填充/粘合辅助（`7_extend` / `7_glue` / `7_fill`，非主 OLC）

**用户纠正（2026-08-15）**：fill/glue **不是**老流程的"经典 OLC"——主 OLC
在 §2.5 的合并（contained/orient/merge）。fill/glue 是合并之后用
**8_* 的长 contigs（2GS）**对 anchors 做 gap 填充/粘合的辅助步骤：它的
`dazz overlap2` 是 **anchor × long**（不是片段间两两 overlap），`dazz
group/layout` 把 long 序列挂到 anchor 分组上补 gap，而不是从零做
overlap-layout。调用链在 `templates/0_bsub.tera.sh`（`--extend 1` 分支），
本地 G37 手动执行（产物在 `7_glue_anchors/`、`7_fill_anchors/`，
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

**fill/glue 脚本内部**（`templates/7_glue_anchors.tera.sh` 与
`7_fill_anchors.tera.sh` 相同骨架）：

1. **Overlap**：`dazz overlap2 ${FILE_ANCHOR} ${FILE_LONG} -b 50 --len 1000
   --idt 0.999 --all`（anchor × long 两两比对）→ `anchorLong.ovlp.tsv`；
   glue 额外做 anchor 内部 overlap（`anchr overlap --serial --len 30
   --idt 0.9999`，只保留每端度数 ≤2 的）→ `anchor.ovlp.tsv`；
2. **Group**：`dazz group anchorLong.db anchorLong.ovlp.tsv
   [--oa anchor.ovlp.tsv] --range 1-N --len 1000 --idt 0.999
   [--max 100 | --max -30] -c GAP_COV`——按 anchor 分组建图（-c = gap
   覆盖度阈值）；
3. **Layout（每组内，拼 long 补 gap）**：`anchr orient`（方向一致化）→ `anchr overlap`
   （组内精确 overlap，`--len 1000 --idt 0.9999`）→ `anchr restrict`
   （按 restrict.tsv 约束，只允许 anchor-anchor / anchor-long 的合法连接）
   → **`dazz layout`** → `group/*.contig.fasta`；
   汇总 `non_grouped.fasta + *.contig.fasta` → `contig.fasta`
   （`faops filter -a 1000`）。

**fill vs glue 的差异**：glue 用 anchor 内部 overlap 辅助分组（`--oa`）、
`--max -30`（允许每个 anchor 进多个组）；fill 用 `--keep --max 100`、
`-c 2`（bsub 版统一传 3）。fill 的输入 anchor = glue 的 contig——即
**先 glue 把 anchors 和长 contigs 粘合，再 fill 用长 contigs 填 gap**；
两者都是合并（§2.5）之后的精修，不改变"合并 = 主 OLC"的定位。

**G37 本地产物与 bsub 版参数不同**：本地 `7_fill_anchors/anchor.fasta`
（16 条，md5 ≠ `7_merge_anchors/anchor.merge.fasta` 的 17 条）说明手动
实验调整过输入，但流程结构一致。

## 3. 关键结论（修正我此前的理解）

0. **老流程是"多覆盖度拆分 → 各部分组装 → 取可靠 anchors → OLC 式合并"**
   （用户纠正，2026-08-15）：
   老流程**先取可靠的 anchors**（覆盖度区间过滤，排除错误区/重复区），再
   **合并**（7_merge_anchors：contained/orient/merge 或重新 anchors 化）——
   **合并本身是经典 OLC**（`anchr merge` = overlap → 有向图 → 拓扑排序 →
   拼接，见 §2.5）；**fill/glue 不是主 OLC**（是用 2GS 长 contigs 做 gap
   填充/粘合的辅助步骤，见 §2.7）；reads 按覆盖度拆成多个部分（G37：
   3 trim 档 × X40/X80 × P 副本 + MR 版；用户经验 30×/60× 合适，不是
   越高越好），各部分独立 unitigs → anchors——**不是"一部分 reads 从前往后
   组装得到一个东西"**；
   - **0a. "取可靠 anchors" = 覆盖度区间筛选**（用户强调）：anchors 步骤
     （bbwrap 回贴 + [lower, upper] 过滤）先选出可靠片段再合并——低覆盖
     （错误区）和高覆盖（重复区）都排除，合并的输入是"可靠 anchors"而非
     全部 unitigs；
   - **0b. 合并 = 经典 OLC**：`anchr merge`（`src/cmd/merge.rs`）= overlap
     → 有向图（边权 = 追加长度）→ 拓扑排序去环 → 按序拼接；`--redo 0`
     分支 `contained → orient → merge → contained` 是完整 OLC（overlap +
     layout + 去冗余）；`--redo`（G37）用 contained 合并 + 重新 anchors 化，
     同一骨架；
   - **0c. fill/glue = gap 填充/粘合辅助**（非主 OLC）：输入 = 合并的
     anchors × 2GS 长 contigs，`dazz overlap2` 是 anchor×long 比对、
     group/layout 把 long 挂到 anchor 分组补 gap——不是从零的片段 OLC。
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
* 老流程把"每个部分"的 unitigs → **可靠 anchors（覆盖筛选）→ OLC 式合并**
  （`anchr merge` = overlap → 拓扑拼接，§2.5）；multik 的 unitig 图跨接
  验证 + 链压实是"合并"的图论替代（边有 k-mer 计数证据，而非仅 overlap
  长度）——两套语义的对应关系值得展开（`merge.rs` 的拓扑拼接 vs multik
  `merge_chains`）；fill/glue（§2.7）的 dazz layout 是 gap 填充辅助，不是
  主 OLC，对比时优先级低；
* merged reads（MR）路径对 multik 的意义（multik 直接吃 merged reads 的
  可行性）；
* 7_extend 的 2GS contigs（8_* 合并）作为 glue/fill 的 long 输入：multik
  输出作为 FILE_LONG 与老流程 anchors 做 OLC 拼装的混搭可行性。

## 5. 现代替代对照（2026-08-15 盘点：老流程每个组件是否都有自实现替代）

| 老流程组件 | 旧工具/命令 | 现代替代（anchr/pgr） | 状态 |
|---|---|---|---|
| trim（dedupe/质量/长度/adapter） | bbtools bbduk + clumpify | `fq clean`（bbduk 等价：adapter k-mer/质量/组成过滤）+ `fq clump`（clumpify-compatible） | **自实现 ✓** |
| PE 合并 + **纠错（ecphase）** | bbtools bbmerge（`anchr mergeread --ecphase "1 2 3"`：ecco/ecct/eccc） | `fq merge`（bbmerge-compatible，含 ecco 等 phase；`--ihist` 写 insert 直方图）/ `fq ec-overlap`（= merge phase 1，ecco） | **自实现 ✓** |
| quorum 步骤（**2026-08-15 核对修正**） | 外部 quorum（k-mer 计数纠错器）；老流程用法 = **丢弃带 `:sub:`/`trunc` 标记的被修正 reads**，pe.cor 保留**未被修正的原始序列**——实际效果是 **reads 筛选**，不是纠错输出 | **`fq s-filter`**（自实现 ✓）：帮助/实现明确"用输入 reads 自身当参考，检查 **quorum 的信号**（无高质量 anchor / truncation / 会 substitution 的碱基 + Poisson 碰撞测试），**不产出修正序列——保留原样或丢弃**"——正是老流程 quorum 的语义（`pgr::libs::kmer::qcheck`） | **自实现 ✓**（用户 2026-08-15 确认；multik 的 solid 阈值是另一层冗余过滤） |
| QC | FastQC | `fq qc`（FastQC-compatible） | **自实现 ✓** |
| insert size 统计 | bbtools/picard statInsertSize | `fq merge --ihist`（bbmerge ihist 格式） | **自实现 ✓** |
| fastk（k-mer 谱/基因组特征，**2026-08-15 核对修正**） | 外部 fastk（-NTable/-Histex）+ GeneScope R | **`pgr kmer table/hist`**（k-mer 表/频率直方图，替代 fastk 计数与 Histex）+ **`pgr kmer gsize --model --plot`**（GenomeScope 模型拟合：kmercov/het/genome size + 谱图，替代 GeneScope） | **自实现 ✓**（`libs/kmer/genomescope.rs`、`cmd_pgr/kmer/gsize.rs`；另有 `profile`/`qhist`） |
| 降采样 | `pgr fa split about` | `pgr fa split`（pgr 库） | **自实现 ✓** |
| unitigs | bcalm / bifrost / superreads / tadpole | `asm unitig`（bcalm graph3 移植 + supermer/DFA 优化） | **自实现+优化 ✓** |
| 迭代组装（现代新增） | — | `asm multik`（metaMDBG 跨接验证 + SKESA 严格链 + megahit 清洗） | **自实现 ✓** |
| anchors（覆盖筛选） | bbwrap perfectmode + basecov + R 阈值 | `asm map`（完美匹配，替代 bbwrap）+ `asm anchor`（覆盖区间过滤） | **自实现+优化 ✓**（补洞逻辑待办） |
| merge anchors | `anchr contained/orient/merge`（自实现） | 同上 + `asm olc --unitigs`（现代合并） | **自实现+优化 ✓** |
| glue/fill | dazz overlap2 / group / layout | `asm ovlp / layout / cns`（精确 overlap 替代） | **自实现 ✓** |
| spades / megahit（**可选参考**） | 外部组装器 | 现代主路线 = multik + OLC；**spades/megahit 保留为可选参考**（用户 2026-08-15：用于对照，如 `multik-g37-quast.md` 的 mr_spades/mr_megahit 对比列） | 可选参考 ✓（不实现、不替代） |
| quast（最终质检） | 外部 quast | 仍用外部 quast（用户明确要求） | 有意保留 ✗ |
| repetitive（重复区提取） | `pgr fa` | `pgr fa range/filter` | **自实现 ✓** |

**结论**：老流程的每个**功能性组件**（trim/merge/纠错/unitigs/anchors/合并/
OLC）现在都有自实现的现代替代（fq 家族 + asm 家族 + pgr fa），且带优化
（supermer/DFA、精确 overlap、图论合并）。**quorum 步骤的替代是
`fq s-filter`**（用户 2026-08-15 确认：s-filter 检查的正是 quorum 的信号、
且"保留原样或丢弃"语义一致），老流程的纠错在 merge 的 ecphase
（`fq merge`/`fq ec-overlap` 已替代），fastk/GeneScope 的替代是
**`pgr kmer gsize --model`**（GenomeScope 拟合，用户 2026-08-15 核对确认）。
**保留外部的组件：quast（最终质量确认）+ spades/megahit（可选参考，
用于对照，用户 2026-08-15 确认）**。
