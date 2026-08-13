# Unitig 生成基准：anchr（unitig / contig / contig --no-bubbles）vs bcalm vs Bifrost vs cuttlefish

> 2026-08-13。hyperfine 1.19.0、`/usr/bin/time -v`（peak RSS）。
> 输入为 G37（E. coli 模拟，ENA ERR486835）的 corrected reads：
> `4_down_sampling/` 的 40×（24 MB）与 80×（48 MB）抽样，
> 以及全量 `2_illumina/merge/pe.cor.fa.gz`（解压后 144 MB，656k reads）。
> k=31，bcalm/Bifrost 用 8 线程；每次运行在全新临时目录。

对比对象：`anchr asm unitig`（anchr 的 k-mer 图压缩实现，BCALM 2
`ograph.cpp graph3` 语义）、`anchr asm contig`（种子化贪心遍历 + 弹泡，
tadpole 兼容）、`anchr asm contig --no-bubbles`（同上但不弹泡）、bcalm 2、
Bifrost 1.3.5、cuttlefish 2.2.0（KMC3 + MPHF + DFA 状态路线）。都从
solid k-mer 生成 unitig/contig，参数按各自常规用法（存在不对称，见备注）：

| 工具 | 命令要点 | 固实阈值 | 长度过滤 | 线程 |
| :--- | :--- | :--- | :--- | :--- |
| `anchr asm unitig` | `--kmer 31 --min-count-seed 3` | ≥3 | 默认 0（不过滤） | 计数并行（默认 8 线程；auto=全核）/ walk 单线程、确定性 |
| `anchr asm contig` | `--kmer 31`（种子化遍历 + 弹泡） | ≥3 | 同上 | 同上 |
| `anchr asm contig --no-bubbles` | `--kmer 31`（种子化遍历，不弹泡） | ≥3 | 同上 | 同上 |
| bcalm | `-kmer-size 31 -abundance-min 3` | ≥3 | 无 | 8 |
| Bifrost | `--kmer-length 31 --clip-tips --del-isolated` | 默认 1 | 无（只剪 tip） | 8 |
| cuttlefish | `-k 31 -t 8 --ref -c 3` | (k+1)-mer ≥3 | 无 | 8 |

## 结果（warmup 1 / runs 3）

### small：24 MB（Q0L0X40P000，40×）

| 工具 | wall（mean ± σ） | peak RSS | 序列数 | N50 | 总长 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| anchr asm unitig | 1.127 s ± 0.045 s | 260 MB | 1048 | 9032 | 606 566 |
| anchr asm contig | 1.393 s ± 0.024 s | 260 MB | 77 | 14 203 | 561 887 |
| anchr asm contig --no-bubbles | 1.340 s ± 0.029 s | 261 MB | 75 | 14 203 | 561 490 |
| bcalm | 1.332 s ± 0.004 s | 231 MB | 1048 | 9032 | 606 566 |
| Bifrost | 899 ms ± 93 ms | 28 MB | 1191 | 7745 | 612 634 |
| cuttlefish | 6.950 s ± 0.103 s | 1211 MB | 1018 | 9032 | 605 663 |

### medium：48 MB（Q0L0X80P000，80×）

| 工具 | wall | peak RSS | 序列数 | N50 | 总长 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| anchr asm unitig | 1.361 s ± 0.041 s | 401 MB | 1554 | 4677 | 626 022 |
| anchr asm contig | 1.583 s ± 0.019 s | 419 MB | 149 | 7633 | 561 886 |
| anchr asm contig --no-bubbles | 1.456 s ± 0.083 s | 470 MB | 145 | 7633 | 561 098 |
| bcalm | 1.556 s ± 0.025 s | 264 MB | 1554 | 4677 | 626 022 |
| Bifrost | 1.090 s ± 0.013 s | 28 MB | 1679 | 4656 | 631 026 |
| cuttlefish | 7.103 s ± 0.037 s | 1267 MB | 1526 | 4677 | 625 181 |

### full：144 MB（2_illumina/merge，全量 656k reads）

| 工具 | wall | peak RSS | 序列数 | N50 | 总长 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| anchr asm unitig | 2.088 s ± 0.075 s | 909 MB | 1482 | 8790 | 622 758 |
| anchr asm contig | 2.466 s ± 0.031 s | 916 MB | 97 | 11 657 | 561 352 |
| anchr asm contig --no-bubbles | 2.397 s ± 0.030 s | 910 MB | 96 | 11 657 | 561 163 |
| bcalm | 2.382 s ± 0.051 s | 541 MB | 1482 | 8790 | 622 758 |
| Bifrost | 2.803 s ± 0.105 s | 28 MB | 1505 | 8350 | 623 546 |
| cuttlefish | 7.685 s ± 0.082 s | 1697 MB | 1459 | 8790 | 622 067 |

## 分析

- **速度**：full k=31 上 anchr unitig 2.09 s vs bcalm 2.38 s（0.88×），
  但这个优势来自**核数差**：anchr 计数用 rayon 全局池（32 逻辑核），
  bcalm 只用了 8 线程。**同核数（1 线程）下 bcalm 反而快 ~10%**
  （10.40 vs 11.46 s，full k=31）。cuttlefish 2.2.0（KMC3 + MPHF +
  DFA 状态）7.69 s，约 anchr unitig 的 3.7×、bcalm 的 3.2×；其总并行
  效率也很低（§11.2）。anchr 的计数实现（per-chunk `count_keys` + 树状
  合并 + 增量 canonical + 前缀索引查询）真正的优势在**内存**而非速度；
- **内存是最大短板**：anchr 峰值 RSS 随输入近似线性且基数很大——
  24 MB 输入 260 MB、48 MB 473 MB、144 MB **0.91 GB**（分别为 bcalm 的
  1.1×/1.7×/1.7×，Bifrost 的 9×/16×/33×）；cuttlefish full 为
  **1.70 GB**（anchr 的 1.9×、bcalm 的 3.1×），KMC3 的 3 GB 内存预算
  与 BBHash MPHF 是主要来源。per-chunk 去重 + 树状合并
  消除了含重复的全局 key 中间体后相对旧实现（1.24/2.9/5.3 GB）大幅
  改善，chunk 尺寸 4096 → 16384 再降 ~25%（根因与后续方向见
  [asm-assemble.md](../design/asm-assemble.md)
  §9）；
- **contig vs unitig**：三种 anchr 模式共享同一 k-mer 表（RSS 几乎
  相同），种子化贪心遍历比严格压缩合并得更多——full 上 contig 97 条、
  unitig 1482 条（无损压实，不过滤），N50 11 657 vs 8790；
- **弹泡影响很小**：contig 与 --no-bubbles 的输出几乎一致（full 97 vs
  96 条、N50 相同 11 657，各档差异仅 1-4 条；耗时上 --no-bubbles 反而
  略慢（14.0 vs 13.6 s，保留路径后排序/处理略多）。对该数据集
  （纠错 reads，气泡少）弹泡与否的影响可忽略；
- **产出一致性**：`anchr asm unitig` 默认不过滤（§10.5），输出与 bcalm
  **逐序列完全相同**（三档 1048/1554/1482 条，N50/总长全部一致，规范化
  集合双向相等）。contig 模式合并更多（full 97 条、N50 11 657）。
  cuttlefish 输出与 bcalm **接近但不逐序列相同**：full 1459 条
  （bcalm 1482）、总长 622 067（bcalm 622 758）、N50 相同 8790——差异
  来自它按 (k+1)-mer 过滤（`-c 3`）且走边中心/DFA 语义，低覆盖端点和
  环状 unitig 的表示与 bcalm 的 k-mer 过滤/压缩不完全一致；canonical
  k-mer 集合核对几乎相等（578 297 vs 578 298，仅缺 1 个）。
  bcalm 保留全部 unitig 的语义分析（顶点分解、无损压实，含旧
  `contained` 的历史实测）见 [asm-assemble.md](../design/asm-assemble.md)
  §10；
- **确定性与功能**：anchr 的 walk 单线程、扫描顺序无关，输出确定性；
  计数阶段并行但结果与并行度无关。原生支持 gz 输入、FASTQ/FASTA、
  GFA 输出、`--links` 边信息——这些是 bcalm/Bifrost 命令行不具备的
  便利性，但速度上并没有真正的单核优势。

## 备注

- 复现：`bash scripts/unitig-bench.sh [small|medium|full] [runs]`
  （脚本自动用 `mktemp` 隔离每次运行，输出 wall/RSS/统计表）；
- 本表 small/full 同批、medium 同批（各批内机器状态一致）；绝对时间随
  机器负载漂移，**相对比值**才是稳定指标；
- 参数不对称：anchr/bcalm 固实阈值 3，Bifrost 默认 1（无 `--min-abundance`
  等价参数，保留低丰度 k-mer），后者 unitig 数略多；anchr 默认不过滤
  （与 bcalm 无损压实一致）；cuttlefish 用 `--ref -c 3`，
  但它的阈值作用在 (k+1)-mer 上（bcalm/anchr 作用在 k-mer 上），因此
  数量只可近似比较；
- 单线程对拍（full k=31，runs 2）：anchr（`RAYON_NUM_THREADS=1`）
  11.46 s / 331 MB，bcalm（`-nb-cores 1`）10.40 s / 654 MB；
  anchr 的 32 核并行版本 2.09 s / 909 MB；
- RSS 用 `/usr/bin/time -v` 单次测量，不是 hyperfine 均值；输入统一为
  明文 FASTA（脚本对 gz 先解压）。
- cuttlefish 由源码构建（`cuttlefish-2.2.0/build`，依赖外部拉取的
  KMC 3.2.1 与 jemalloc 5.2.1，GCC 15 下需补 `<stdexcept>`/`<cstdint>`
  等头文件）；长 k 版用 `-DINSTANCE_COUNT=64` 重编（k 上限 127，cuttlefish
  只支持奇数 k）；脚本通过 `CUTTLEFISH` 环境变量指定二进制路径。

## K 缩放（small 24 MB，runs 2；`K=64` / `K=100 bash scripts/unitig-bench.sh small 2`）

| 工具 | k=31 | k=64 | k=100 |
| :--- | ---: | ---: | ---: |
| anchr unitig wall | 1.20 s | 1.09 s | 1.29 s |
| anchr unitig RSS | 267 MB | 342 MB | 330 MB |
| bcalm wall / RSS | 1.27 s / 237 MB | 2.31 s / 254 MB | 2.10 s / 242 MB |
| Bifrost wall / RSS | 0.88 s / 29 MB | 1.06 s / 29 MB | 1.15 s / 29 MB |
| anchr/bcalm 比值 | **0.9×** | **0.5×** | **0.6×** |

前缀索引查询（`get_count` 二分从 ~20 次比较降到 ~5 次）+ FNV claim 集合
后，k 缩放拉平且**全 k 反超 bcalm**（k=31 0.9×、k=64 0.5×、k=100
0.6×），各 k 输出与 bcalm 逐序列相同（1048/486/251 条）；k=64/100 与
Bifrost 打平。内存各 k 均 ~0.3 GB。contig 模式在长 k 下合并优势更明显
（k=100 N50：contig 36 028 vs unitig 14 304 vs Bifrost 37 787）。

> 注：本节为较早批次的偶数 k 数据；下节是 2026-08-14 的奇数 k 重测
> （含 cuttlefish），同批内 bcalm 在 k=99/127 更快，绝对数值以同批
> 内相对比值为准。

## 长 K（small 24 MB，runs 2，奇数 k 对齐 cuttlefish；`K=63|99|127 bash scripts/unitig-bench.sh small 2`）

cuttlefish 只接受奇数 k，因此本组统一用 63/99/127 对比。anchr unitig 与
bcalm 输出逐序列一致（各 k 计数相同）；cuttlefish 用 `--ref -c 3`
（(k+1)-mer 阈值）。

| 工具 | k=63 | k=99 | k=127 |
| :--- | ---: | ---: | ---: |
| anchr unitig wall | 1.70 s | 2.32 s | 2.49 s |
| anchr unitig RSS | 341 MB | 345 MB | 242 MB |
| anchr unitig 输出（条/N50/总长） | 495 / 11 496 / 610 751 | 249 / 14 057 / 604 449 | 5156 / 245 / 1 194 913 |
| bcalm wall / RSS | 1.70 s / 217 MB | 2.16 s / 238 MB | 1.84 s / 235 MB |
| bcalm 输出（条/N50/总长） | 同 anchr | 同 anchr | 同 anchr |
| Bifrost wall / RSS | 1.08 s / 28 MB | 1.16 s / 28 MB | 0.43 s / 28 MB |
| Bifrost 输出（条/N50/总长） | 449 / 18 957 / 608 642 | 105 / 37 785 / 589 056 | 1072 / 733 / 665 590 |
| cuttlefish wall | 7.29 s | 7.80 s | 7.79 s |
| cuttlefish RSS | 1204 MB | 1200 MB | 1194 MB |
| cuttlefish 输出（条/N50/总长） | 486 / 11 496 / 610 188 | 250 / 14 303 / 604 543 | 6040 / 223 / 1 305 829 |

分析：cuttlefish 耗时与 RSS 对 k 几乎**不敏感**（7.3-7.8 s、~1.2 GB，
KMC 枚举 + MPHF 构建占大头），k=63/99 输出与 bcalm 接近（条数差 ≤1，
N50 相同或近同）；k=127 时明显分化（6040 vs 5156 条、总长 1 305 829 vs
1 194 913）——长 k 下 reads 相对 k 变短，(k+1)-mer 阈值与 k-mer 阈值、
环状/端点语义的差异被放大。anchr 长 k 表现：k=63 与 bcalm 打平
（1.70 vs 1.70 s），k=99/127 反而慢（2.32 vs 2.16 s、2.49 vs 1.84 s）；
cuttlefish 全程约 anchr 的 3-4.5×。

## supermer + DFA 组合（2026-08-14，`--supermer --dfa -p8 -m<k/4>`）

计数（pgr supermer）+ 状态分类（DFA）组合 vs 默认引擎，输出全部逐字节
一致。m 扫描：k=31 最优 8，k=63/99 最优 12，k=127 两者接近；启发式
`m = min(12, max(5, k/4))` 覆盖各档。

> 本节"默认"列是 2026-08-14 前的桶扫描 + 内存计数引擎；此后默认已是
> DFA walk + 流式计数 + 状态表内嵌计数（默认线程 = `min(逻辑核/2, 8)`，
> 最新数字见下节"当前默认引擎实测"）。
> small/medium 行仍为旧引擎对照。

| 规模 | k | 默认 wall / RSS | 组合 wall / RSS | 变化 |
| :--- | ---: | ---: | ---: | ---: |
| small | 31 | 911 ms / 261 MB | **746 ms** / 270 MB | -18% / +3% |
| medium | 31 | 1000 ms / 480 MB | **949 ms** / 360 MB | -5% / -25% |
| full | 31 | 1.26 s / 699 MB | 1.442 s / 791 MB | +14% / +13% |
| small | 63 | 1.268 s / 332 MB | **1.034 s** / 313 MB | -18% / -6% |
| small | 127 | 1.877 s / 240 MB | **1.265 s** / 297 MB | -33% / +24% |
| full | 99 | 1.63 s / 1348 MB | 2.825 s / 1283 MB | +73% / -5% |

结论：**当前默认（流式 direct 计数 + DFA walk）是最佳配置**；supermer
组合在旧引擎对照下多数档位更快，但在新默认下 full k31 无优势、k99
墙钟 +50%（内存仅省 12%）——本节为 pgr §12.3 优化前的旧结论；
**pgr 优化后（2026-08-14）supermer 已反超 direct 并成为 FASTA 默认**，
见下节。

## 当前默认引擎实测（2026-08-14 优化后）

默认 = 流式 direct 计数（chunk 32768）+ DFA walk + 状态表内嵌计数 +
visited 字节数组；**FASTA 默认计数 = pgr supermer（§12.3 优化后）**，
FASTQ 自动回退 direct 保质量门控语义。`--parallel` 默认
`min(逻辑核/2, 8)`（本机 32 核 → 8 线程），单池复用、线程数稳定。
runs 2，输出与旧引擎逐字节一致。

| 规模 | k | wall | RSS | vs 旧默认（桶扫描+内存计数） |
| :--- | ---: | ---: | ---: | ---: |
| small | 31 | 590 ms | 171 MB | -35% / -34% |
| medium | 31 | 809 ms | 242 MB | -19% / -50% |
| full | 31 | 1.358 s | 597 MB | -14% / -37% |
| small | 63 | 768 ms | 269 MB | -39% / -19% |
| small | 127 | 899 ms | 239 MB | -52% / 0% |
| full | 99 | 2.187 s | 1131 MB | -19% / -45% |

> 旧默认引擎用满 32 逻辑核，墙钟更快（full k31 1.21 s）；默认降为
> `half(≤8)` 后大输入墙钟回升，但内存大幅下降、线程占用稳定（不
> 吃满机器）。需要全速时显式 `--parallel auto`（k31 可到 ~1.06 s）。

至此 anchr 侧已把能压的压完：walk 隐藏瓶颈（DFA）、计数流式化
（内存）、状态表内嵌计数 + visited 字节数组（walk O(1)）。剩余最大
单一差距是 pgr 计数本身（~0.9 s vs FastK 0.68 s，pgr 侧推进）。
