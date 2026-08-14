# 老手工流程理解（G37 例子，2026-08-14，基于完整代码阅读）

> 本文记录我对 `results/model.md` 的 G37 手动模拟流程的理解。**每个结论
> 都来自完整阅读对应代码**：流程模板 `templates/*.tera.sh`（生成实际脚本
> 的源）、命令实现 `src/cmd/*.rs`、以及 G37 产物脚本
> `~/data/anchr/g37/`（`0_script/*.sh` 和各目录的 `unitigs.sh`/
> `anchors.sh`/`quorum.sh`）。

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

* **4_**：输入 `2_illumina/Q{1}L{2}/pe.cor.fa.gz`，`pgr fa split about -c
  genome×cov` 按目标覆盖切分 → `4_down_sampling/Q{1}L{2}X{40/80}P{000...}/
  pe.cor.fa`（P = 随机降采样副本，数量由 statp=2 约束）。
* **6_**：输入 `2_illumina/merge/pe.cor.fa.gz`（merged reads 的 quorum
  输出）→ `6_down_sampling/MRX{40/80}P{...}/pe.cor.fa`（MR 版本）。

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

## 3. 关键结论（修正我此前的理解）

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
5. **statFinal 是两条并行路径的对比**：unitigs→anchors→merge（7_*）vs
   reads→spades/megahit→contigs→anchors（8_*），mr_* = merged reads 版本；
6. **multik 的正确对比对象**是老流程的 `7_merge_anchors` / `8_*` 最终产物
   （statFinal），不是中间 unitigs。

## 4. 待讨论（我的初步看法，非结论）

* multik 的长 unitigs 是否天然包含 anchors 的"覆盖筛选"语义（bbwrap 回贴 +
  覆盖区间过滤）？metaMDBG 的渐进丰度过滤覆盖"低覆盖排除"，但 anchors 的
  `upper`（重复区排除）在 multik 没有对应；
* 老流程的 `anchr contained/orient/merge`（anchors 后处理）与 multik 的
  压实/去冗余的对应关系；
* merged reads（MR）路径对 multik 的意义（multik 直接吃 merged reads 的
  可行性）。
