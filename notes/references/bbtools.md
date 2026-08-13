# BBTools: reads 处理与组装的主参考工具包

> 整理于 2026-08-13，源自对仓库根 `BBTools-40.01/`（2026-02-11）源码的分析。
> BBTools 是 fq/asm 迁移的主参考：`fq` 家族以 BBTools 39.38 逐字节核对
> （golden 见 [audit-fq.md](../audit/audit-fq.md)），`asm` 家族对照 BBTools
> tadpole（contig）与 bbmap `perfectmode`（map，见
> [audit-asm.md](../audit/audit-asm.md)）。本文汇总与 anchr 相关的模块与
> 语义要点，作为已随迁设计文档的索引与补全。

## 1. 简介

`BBTools`（Brian Bushnell, JGI/DOE）是 Java 实现的快速序列处理工具包，
以单机大内存、多线程为设计目标。本仓库内的 `BBTools-40.01/` 是其完整源码
发行（脚本壳 + Java 源码 + 资源），目录自带 `.gitignore`（`*`），仅作参考，
不被 git 跟踪。

- **版本基线**：本地实际安装 39.38（`~/.cbp/libexec/bbtools/`，`bbduk.sh`
  2024-11-18，入口 `jgi.BBDuk`）；仓库根 40.01 是更新版（入口
  `bbduk.BBDukS` 家族）。运行对照以 39.38 为准，40.01 只作行为参考。
- **运行方式**：`*.sh` 脚本壳调用 `java -cp BBTools.jar <入口类>`；
  `current/` 下为按包组织的 Java 源码，`jni/` 下为原生扩展（两版本均禁用）。

## 2. 目录结构

| 路径 | 内容 |
| :--- | :--- |
| `*.sh`（顶层 309 个） | 命令入口脚本，每个工具一个 `.sh` |
| `current/` | Java 源码，按包组织（`clump/`、`jgi/`、`assemble/`、`align2/`、`bbduk/` 等） |
| `resources/` | 参考数据（如 `lambda.fa.gz`，OLC 冒烟用） |
| `jni/` | JNI 原生扩展（`BBMergeOverlapper.c`、`BandedAlignerJNI.c` 等，已禁用） |

## 3. 与 anchr 相关的模块映射

BBTools 入口 → anchr/pgr 目标（映射出处 [fq-trim-replace.md](../design/fq-trim-replace.md) §4.2、
[fq-merge-replace.md](../design/fq-merge-replace.md)、
[asm-assemble.md](../design/asm-assemble.md)）：

| BBTools 工具 | 入口类 | 源码位置 | anchr 目标 |
| :--- | :--- | :--- | :--- |
| clumpify | `clump.Clumpify` | `current/clump/Clumpify.java` | `fq clump` |
| bbnorm | `jgi.KmerNormalize` | `current/jgi/KmerNormalize.java` | `fq norm` |
| reformat 降采样 | `jgi.ReformatReads` | `current/jgi/ReformatReads.java` | `fq sample` |
| bbduk | `jgi.BBDuk`（39.38）/ `bbduk.BBDukS`（40.01） | `current/jgi/BBDuk.java`、`current/bbduk/` | `fq trim-adapter` / `fq trim` |
| bbmerge / tbo/tpe | `jgi.BBMerge`（委托 `BBMergeOverlapper.mateByOverlapRatio`） | `current/jgi/BBMergeOverlapper.java` | `fq merge` / `fq overlap` |
| repair | `jgi.SplitPairsAndSingles rp` | `current/jgi/SplitPairsAndSingles.java` | `fq split` |
| kmercountexact | `jgi.KmerCountExact` | `current/jgi/KmerCountExact.java` | k-mer 基础（留 pgr） |
| tadpole | `assemble.Tadpole` | `current/assemble/Tadpole.java` | `asm contig`/`asm unitig`、`fq extend`/`fq ec-kmer` |
| bbmap | `align2.BBMap`（核心映射线程 `align2.AbstractMapThread`） | `current/align2/AbstractMapThread.java` | `asm map`（`perfectmode`） |

## 4. 已确认的关键语义要点

以下结论来自已随迁的 audit/design 文档，实现与 golden 核对时直接引用：

1. **质量值存储**：`Read.quality` 存 phred 值（输入时减 ASCII_OFFSET），输出再加回
   （[fq-trim-replace.md](../design/fq-trim-replace.md) §4.1）。
2. **clumpify 桶序**：外部 bucket 路径输出为"按桶拼接"序，与 BBTools 大数据行为
   一致，但与内存路径不同（[audit-fq.md](../audit/audit-fq.md)）；`--dedupe
   --dupesubs 0` 整对去重，N 通配精确匹配，保留期望错误更少的一对。
3. **tadpole 内存模型**：多字 `Kmer`（`Vec<u64>` 共 2k 位，镜像 BBTools long
   array），k 无上限；逐字节一致需复刻其 `-Xmx` 相关内存语义，已文档化偏差
   （[asm-assemble.md](../design/asm-assemble.md)）。
4. **bbmap perfectmode**：种子-验证、完美匹配（无错配无缺口）、`ambiguous=all`
   语义在 `current/align2/AbstractMapThread.java:1371` 确认（`maxMismatches=
   (PERFECTMODE||SEMIPERFECTMODE)?0:...`，另见 :810 门控"imperfect 不可完美映射"）
   （[asm-map.md](../design/asm-map.md)）。
5. **golden 数据**：`tests/bbtools/Lambda/`（`R1.2k.fq.gz`、`golden/` 等）是
   fq 逐字节对照的基准，随 fq 测试迁移批次从 pgr 迁入。

## 4.5 关键工具源码级要点（默认值 + 算法）

> 本节直接读 `BBTools-40.01/` 的 `*.sh` 脚本（usage help 即各工具默认值文档）与
> `current/` Java 源码，给 anchr 各目标命令提供**调参与语义对照**。`*.sh` 的默认值
> 是作者在脚本里写死的，比 README 更可信；Java 源码再确认算法语义。
> 已覆盖：bbnorm、bbmap、tadpole、bbduk、clumpify、reformat、bbmerge（§4.5.1–4.5.7）。

### 4.5.1 bbnorm（`fq norm` 血统）——`jgi/KmerNormalize`

脚本默认值（`bbnorm.sh` usage）：

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `k` | 31 | kmer 长度（<32 最省内存，可任意高） |
| `bits` | 32（脚本硬编码） | bloom filter 每 cell 位数（2/4/8/16/32）；最大可记录深度 = 2^bits；2-pass 自动降为 16 |
| `hashes` | 3 | 每个 kmer 哈希并存储的次数 |
| `minq` / `minprob` | 6 / 0.5 | 忽略含质量 <minq 碱基的 kmer / 整体正确概率低于此的 kmer |
| `target` | 100 | 目标归一化深度（**kmer 深度**，非 read 深度） |
| `maxdepth` / `mindepth` | -1 / 5 | 低于 maxdepth 不再降采样；低于 mindepth 的 kmer 不参与 read 深度计算 |
| `minkmers`（mgkpr） | 15 | 每条 read 至少需这么多 >mindepth 的 kmer 才保留 |
| `percentile`（dp） | 54.0 | read 深度取 kmer 深度数组的**第 54 百分位**（1-100 可调） |
| `uselowerdepth`（uld） | t | 双端时取两端较低者作深度代理 |
| `passes` | 2 | 2-pass 更高精度/纠错/深度控制 |
| `deterministic` | t | 确定性随机数，保证多次运行输出一致 |

算法语义：
- **归一化判定**：先建 kmer 计数 bloom filter，再逐 read 收集其 kmer 深度数组，
  取 `percentile`（默认 54%）作为该 read 的深度代理；深度 > target 的 read 降采样/丢弃。
  深度换算 `Dr = Dk*(R/(R-K+1))`（`bbnorm.sh:56`），源码里
  `readDepth = median_all*(avgReadLen/(avgReadLen-k+1))`（`KmerNormalize.java:2149`）。
- `percentile` 解析：>1 且 ≤100 时除以 100（`KmerNormalize.java:331`），故给 54 等价 0.54。
- **错误检测**（可选，默认 `tossbadreads=f`）：按 kmer 深度高/低百分位（`hdp=90`/`ldp=25`）
  与 `errordetectratio=125` 判错 read。**纠错**（`ecc=f`，作者明示 "Tadpole is now preferred"）
  有 `ecclimit=3`、`errorcorrectratio=140` 等。
  - 硬编码常量：`HIGH_PERCENTILE=0.90f`、`LOW_PERCENTILE=0.25f`（`KmerNormalize.java:3798-3800`）、
    `ERROR_DETECT_RATIO=125`、`HTHRESH=12`、`LTHRESH=3`（:3807-3811）、
    `ERROR_CORRECT_RATIO=140`（:3814）。
  - 判错规则（`:3135`）：read 判为错误 read 当且仅当 `high<=LTHRESH`
    （高深度 kmer ≤3）**或**（`high>=HTHRESH` **且** `low<=LTHRESH`）**或**
    `high >= low*ERROR_DETECT_RATIO`（高/低深度比 ≥125）；其中 high/low 是该 read 的
    kmer 深度数组在 `hdp`/`ldp` 百分位处的取值（HTHRESH 默认 12，见 :3135-3139 判定计数）。
  - `percentile` 解析 `>1 且 ≤100` 时除以 100（`KmerNormalize.java:331`），hdp/ldp 同理（:335/:339）。
- 硬编码：`bbnorm.sh:152` 启动命令固定追加 `bits=32`——这是脚本层写死的默认，
  与 usage 里的 `bits=32` 一致。

### 4.5.2 bbmap（`asm map` 的 `perfectmode` 血统）——`align2/BBMap` + `AbstractMapThread`

脚本默认值（`bbmap.sh` usage）+ 硬编码：

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `k` | 13 | 索引/映射 kmer（范围 8-15） |
| `fastareadlen` | 500（脚本硬编码，`bbmap.sh:341`） | 长 FASTA read 在此长度处切开 |
| `minid` | 0.76 | 最低比对 identity（越高越快越不敏感） |
| `minhits` | 1 | 候选位点所需最小 seed 命中数 |
| `perfectmode` / `semiperfectmode` | f / f | 仅允许完美映射（见下） |
| `ambiguous` | best | 多最高分位点时的行为（best/toss/random/all） |
| `tipsearch` / `maxindel` | 100 / 16000 | 末端缺失暴力搜索范围 / indel 上限 |
| 硬编码追加（`bbmap.sh:341`） | `build=1 overwrite=true fastareadlen=500` | 启动时固定参数 |

`perfectmode` 语义（`AbstractMapThread.java`）：
- `maxMismatches = (PERFECTMODE || SEMIPERFECTMODE) ? 0 : ...`（`:1371`）——完美/半完美
  模式强制 0 错配，即 **seed-verify 精确匹配**，对应 `pgr asm map` 的完美匹配路线。
- `PERFECTMODE||SEMIPERFECTMODE` 时，含未定义碱基（N）的 read 直接判不可映射
  （`r.containsUndefined()` → return -1，`:810`）——`semiperfectmode` 的"完美"允许参考
  侧 N，但不允许 read 侧 N。
- `ambiguous=all` 保留所有最高分位点（对应 `pgr asm map` 的多 hit 语义）。

### 4.5.3 tadpole（`asm contig`/`unitig`、`fq extend`/`ec-kmer` 血统）——`assemble/Tadpole`

脚本默认值（`tadpole.sh` usage）：

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `k` | 31 | kmer 长度（1 到无穷，内存随 k 增） |
| `minprob` | 0.5 | 忽略整体正确概率低于此的 kmer |
| `rcomp` | t | kmer 与其反向互补合并计数 |
| `mincountseed`（mcs） | 3 | 播种/开始延伸的最小 kmer 计数 |
| `mincountextend`（mce） | 2 | 继续延伸的最小 kmer 计数（低覆盖宏基因组建议 1） |
| `mincountretain` / `maxcountretain` | 0 / inf | 低于/高于此计数的 kmer 丢弃 |
| `branchmult1` / `branchmult2` / `branchlower` | 20 / 3 / 3 | 高/低深度下的分支比阈值 |
| `minextension` / `mincontig` / `mincoverage` | 2 / auto / 1 | 延伸/写 contig 的最小长度与覆盖度 |
| `contigpasses` / `contigpassmult` | 16 / 1.7 | 递减 seed 深度建 contig 的轮次 / 每轮比值 |
| `popbubbles` | t | 弹泡（需 `processcontigs=t`） |
| `mode` | contig | contig/extend/correct 等 |

> 内存模型：多字 `Kmer`（`Vec<u64>` 共 2k 位）镜像 BBTools long array，k 无上限——
> 逐字节一致需复刻其 `-Xmx` 内存语义（详见 `asm-assemble.md`）。

**延伸/分支（junction）判定**（`assemble/Tadpole.java`）：
- 参数硬编码默认：`minExtension=2`、`minCoverage=1`、`minCountSeed=3`、`minCountExtend=2`、
  `branchMult1=20`、`branchMult2=3`、`branchLowerConst=3`（`Tadpole.java:2689-2703`）；
  `minContigLen=-1` 时自动取 `max(124, 2*k)`（:582-583）。
- 分支判定 `isJunction(max, second)`（`Tadpole.java:2580`）：**不认为分支**（继续直行）当
  `second<1` **或** `second*branchMult1 < max` **或**（`second<=branchLowerConst` 且
  `max>=max(minCountExtend, second*branchMult2)`）；否则判为分支/断点。两条方向任一判分支
  即算 junction（:2573-2576）。
- `errorPath/errorMult1/errorMult2` 等纠错延伸阈值见 `Tadpole.java:2705-2709`，供 ec-kmer
  扩展参考。

### 4.5.4 bbduk（`fq trim-adapter`/`trim` 血统）——`jgi/BBDuk`（39.38）/ `bbduk.BBDukS`（40.01）

脚本默认值（`bbduk.sh` usage）：

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `k` | 31 | 找污染物的 kmer 长度（比 `bbnorm.sh` 的 k 语义不同：这是参考 kmer） |
| `ways` | 8 | 参考 kmer 索引分片数（7 或 2 的幂） |
| `rcomp` | t | 除正向外也查反向互补 |
| `maskmiddle`（mm） | t | 把 kmer 中间碱基当通配符提敏感度（=mm=1 奇长 / mm=2 偶长） |
| `hammingdistance`（hdist） | 0 | 参考 kmer 的最大汉明距（仅替换；内存 ∝ (3K)^hdist） |
| `editdistance`（edist） | 0 | 最大编辑距（替换+indel；内存 ∝ (8K)^edist） |
| `minkmerhits` / `minkmerfraction` | 1 / 0.0 | 判匹配的最小命中 kmer 数 / 分数 |
| `mincovfraction`（mcf） | 0.0 | 覆盖分数（给则覆盖 mkh/mkf） |
| `copyundefined`（cu） | f | 对参考的非 AGCT IUPAC 碱基展开所有可能拷贝 |
| `removeifeitherbad`（rieb） | t | 双端任一端匹配即丢 |
| `ecco` | f | 交叠双端先用 BBMerge 纠错 |

> bbduk 的 `--qtrim`/`trimq` 是 **BBDuk 质量修剪**模式（`r/l/rl/w/f`），**不是** cutadapt
> 的 Mott（见 `cutadapt.md` §5.4 澄清）。

**qtrim 质量修剪算法**（`bbduk.sh:208-217`，qtrim 取值 `rl/f/r/l/w`，默认 `f`；`trimq=6` 默认）：
- 实际修剪委托 `shared/TrimRead.trimFast`（`bbduk/BBDukProcessorS.java:1297` 调用）。
  `TrimRead.java` 有三种修剪器，由静态开关选择：
  - **`optimalMode=true`（默认，`TrimRead.java:954`）**：走 `testOptimal`
    （`TrimRead.java:348`）——这是一个 **Kadane 最大子数组**算法（源码注释 :354-357 明示）。
    对每个碱基 `delta = avgErrorRate - probError`（`avgErrorRate = phredToProbError(trimq)`，
    `probError` 查 `QualityTools.PROB_ERROR[q]`，N 或 q<1 用 `nprob`，见 :375）；累计和
    非负则并入当前段，取累计和最大的连续区间为"保留段"，其左/右侧全部修剪
    （`:398-404`，全坏 read 时 `maxScore=0 → left=0,right=len` 全剪）。
  - **`windowMode=false`（`TrimRead.java:956`）**：`testRightWindow`（:438）滑动窗口，
    `windowLength=4`（:961），窗口质量均值 <trimq 即剪。
  - **朴素模式**：`testLeft`/`testRight`（:416/:457）要求 **`minGoodInterval=2`**（:949）个
    连续 `q>trimq` 的好碱基，否则从最后一个坏碱基处修剪；质量数组为 null 时退化为按 N
    修剪（`testLeftN`/`testRightN`，:477/:491）。
- `anchr fq trim` 的质量修剪血统即此 `trimFast`，逐字节 golden 时注意默认走 **optimal
  (Kadane)** 分支，而非朴素的"连续 2 好碱基"——两者结果不同。

### 4.5.5 clumpify（`fq clump` 血统）——`clump/Clumpify`

脚本默认值（`clumpify.sh` usage）：

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `k` | 31 | clump/纠错 kmer（1-31） |
| `hashes` | 4 | 哈希掩码数（0=raw kmer） |
| `border` | 1 | 距 read 端此 bp 内的 kmer 不用 |
| `dedupe` | **f** | 去重默认**关闭**（双端需两端都匹配才去重） |
| `subs`（s） | 2 | 两条重复序列间允许的最大替换数 |
| `subrate`（dsr） | 0.0 | 若设，则允许替换数 = max(subs, subrate*min(len1,len2)) |
| `allowns` | t | no-called（N）碱基不算替换 |
| `scanlimit` | 5 | 遇非重复后继续扫描这么多 read（提高近似去重检测） |
| `optical` / `dupedist` | f / 40 | 光学重复检测 / 距离 |
| `groups` | auto | 中间文件数（按文件大小估，防内存溢出） |

> 注意：`dedupe` 脚本默认是 `f`（不开启去重）——需显式 `dedupe=t`。anchr `fq clump`
> 的 `--dedupe --dupesubs 0` 整对去重语义即对应此处 `subs=0` 的精确去重。

**去重判定算法**（`clump/Clump.java`）：
- 默认 `maxSubstitutions=2`（`Clump.java:1398`）、`dupeSubRate=0`（:1400）；
  `dupesubs` 若设则 `maxSubstitutions` 取 `max(subs, subrate*min(len1,len2))`。
- 序列相等判定 `equals(a,b,maxSubs)`（`Clump.java:607`）：长度不等即非重复；逐位比较，
  累计替换数 `subs>maxSubs` 即判非重复（:611-626）。**`allowNs=t` 时 N 视为通配**
  （与任何碱基相等，:612-618）；`allowns=f` 则 N 也计为替换。
- `scanlimit=5`（`Clump.java:1396`）遇非重复后仍继续扫描 N 个 read 提高近似检测，
  与脚本 usage 一致。
- 纠错阈值硬编码见 `Clump.java:1415-1429`（`minCountCorrect=4`、`minRatio=30.0`、
  `minRatioMult=1.80`、`minRatioOffset=1.9` 等），供 clump 纠错扩展参考。

### 4.5.6 reformat（`fq sample` 血统）——`jgi/ReformatReads`

脚本默认值（`reformat.sh` usage）：

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `samplerate` | 1 | 随机输出该比例的 read（1=不采样） |
| `samplereadstarget`（srt） | 0 | 精确输出 read（或 pair）数 |
| `samplebasestarget`（sbt） | 0 | 精确输出碱基数 |
| `sampleseed` | -1 | 采样 PRNG 种子（正数=确定性采样） |
| `upsample` | f | 目标大于输入时允许复制 read |
| `verifypaired`（vpair） | f | 校验 read 名是否配对 |

> `pgr fq sample` 做确定性降采样时，`sampleseed` 正数即确定性语义的来源。

### 4.5.7 bbmerge（`fq merge`/`fq overlap` 血统）——`jgi/BBMergeOverlapper`

> §3 中 bbmerge 映射到 `fq merge`/`fq overlap`，此前 §4.5 未覆盖，此处补齐。

脚本默认值（`bbmerge.sh` usage）：

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `mininsert` / `mininsert0` | 15 / 12 | 可合并的最小 insert size / 第一轮筛选下限（`bbmerge.sh:84-86`） |
| `minoverlap` / `minoverlap0` | 12 / 8 | 可合并的最小重叠 bp / 预筛选下限（:87-89） |
| `maxratio` | 0.09 | 最大错误率（越高合并越多，`:128`） |
| `ratiomargin` | 5.5 | 比值判定裕度（越低合并越多，min=1，:129） |
| `minsecondratio` | 0.1 | 次优重叠的比值下限（判是否歧义，:135） |
| `maxexpectederrors`（mee） | 0 | 若 >0，期望错误超限的 read 不合并（:59） |
| `ecco` | f | 只纠错重叠部分不合并（:79） |
| `minprob` | 0.5 | 忽略正确概率低于此的 kmer（:191） |

算法语义（`jgi/BBMergeOverlapper.java`）：
- 主入口 `mateByOverlapRatio`（:98）→ 纯 Java 分支 `mateByOverlapRatioJava`（:411），
  JNI 分支（`Shared.USE_JNI`，:114）两版本均禁用（anchr 用纯 Rust 移植）。
- `mateByOverlapRatioJava`（:411）两阶段：
  1. `findBestRatio`（:425）先在各 insert 位置打分，`x>maxRatio` 则直接判无合并（:429-433）；
  2. 主循环（:460 起）按 insert 从大到小遍历，对每个重叠区累计 good（`gIncr`）与 bad（`bIncr`），
     坏碱基数超 `badlimit = extraMult*(min(bestRatio,maxRatio)*margin*overlapLength)+1+extraBadlimit`
     （:468）即提前终止；记录最优与次优，`minSecondRatio` 用于歧义判断。
- 返回 `-1` 表示无有效合并；`rvector[2]/[4]` 记录失败原因（:430-432）。
- `expectedMismatches(a,b,overlap)`（:1117）供 `maxexpectederrors` 判定用。
- 无质量数组时走 `mateByOverlapRatioJava`（:125），有质量走
  `mateByOverlapRatioJava_WithQualities`（:158），把质量转概率加权匹配/错配。

## 5. 版本差异与注意事项

- 39.38 → 40.01 的主要变化：`bbduk` 入口从 `jgi.BBDuk` 改为 `bbduk.BBDukS`；
  其余 7 个 trim 流水线工具的入口类不变（[fq-trim-replace.md](../design/fq-trim-replace.md) §4.2）。
- 各 `*.sh` 在 launch 时除类名外还会硬编码默认参数，比对语义时应留意，如
  `bbnorm.sh` 追加 `bits=32`（`bbnorm.sh:152`）、`bbmap.sh` 追加
  `build=1 overwrite=true fastareadlen=500`（`bbmap.sh:341`）。脚本内的 usage
  help 文本还有大量工具默认值可直接参考，如 bbnorm 默认 `k=31 target=100
  percentile=54.0 passes=2 minq=6 mindepth=5`（`bbnorm.sh`），bbmap 默认
  `k=13 perfectmode=f ambiguous=best minid=0.76`（`bbmap.sh`）——`pgr fq norm`
  （bbnorm 血统）与 `asm map`（bbmap perfectmode）调参时可直接对照这些默认值。
- `BBMergeOverlapper` 的 Java fallback 有 JNI 版本，anchr 侧用纯 Rust 移植，
  两版 JNI 均不采用。
- 性能对照见 [bbtools-vs-anchr.md](../benchmarks/bbtools-vs-anchr.md)
  （hyperfine 实测 BBTools 39.38 8 线程 vs anchr release 单线程）。
