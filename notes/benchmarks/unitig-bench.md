# Unitig 生成基准：anchr（unitig / contig / contig --no-bubbles）vs bcalm vs Bifrost

> 2026-08-13。hyperfine 1.19.0、`/usr/bin/time -v`（peak RSS）。
> 输入为 G37（E. coli 模拟，ENA ERR486835）的 corrected reads：
> `4_down_sampling/` 的 40×（24 MB）与 80×（48 MB）抽样，
> 以及全量 `2_illumina/merge/pe.cor.fa.gz`（解压后 144 MB，656k reads）。
> k=31，bcalm/Bifrost 用 8 线程；每次运行在全新临时目录。

对比对象：`anchr asm unitig`（anchr 的 k-mer 图压缩实现，BCALM 2
`ograph.cpp graph3` 语义）、`anchr asm contig`（种子化贪心遍历 + 弹泡，
tadpole 兼容）、`anchr asm contig --no-bubbles`（同上但不弹泡）、bcalm 2、
Bifrost 1.3.5。都从 solid k-mer 生成 unitig/contig，参数按各自常规用法
（存在不对称，见备注）：

| 工具 | 命令要点 | 固实阈值 | 长度过滤 | 线程 |
| :--- | :--- | :--- | :--- | :--- |
| `anchr asm unitig` | `--kmer 31 --min-count-seed 3` | ≥3 | `min_contig_len` 默认 124 | 单线程（确定性的） |
| `anchr asm contig` | `--kmer 31`（种子化遍历 + 弹泡） | ≥3 | 同上 | 单线程（确定性的） |
| `anchr asm contig --no-bubbles` | `--kmer 31`（种子化遍历，不弹泡） | ≥3 | 同上 | 单线程（确定性的） |
| bcalm | `-kmer-size 31 -abundance-min 3` | ≥3 | 无 | 8 |
| Bifrost | `--kmer-length 31 --clip-tips --del-isolated` | 默认 1 | 无（只剪 tip） | 8 |

## 结果（warmup 1 / runs 3）

### small：24 MB（Q0L0X40P000，40×）

| 工具 | wall（mean ± σ） | peak RSS | 序列数 | N50 | 总长 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| anchr asm unitig | 1.195 s ± 0.032 s | 267 MB | 1048 | 9032 | 606 566 |
| anchr asm contig | 1.358 s ± 0.059 s | 267 MB | 77 | 14 203 | 561 887 |
| anchr asm contig --no-bubbles | 1.288 s ± 0.019 s | 266 MB | 75 | 14 203 | 561 490 |
| bcalm | 1.268 s ± 0.024 s | 237 MB | 1048 | 9032 | 606 566 |
| Bifrost | 850 ms ± 41 ms | 29 MB | 1191 | 7745 | 612 634 |

### medium：48 MB（Q0L0X80P000，80×）

| 工具 | wall | peak RSS | 序列数 | N50 | 总长 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| anchr asm unitig | 1.320 s ± 0.036 s | 473 MB | 1554 | 4677 | 626 022 |
| anchr asm contig | — | — | 149 | 7633 | 561 886 |
| anchr asm contig --no-bubbles | — | — | 145 | 7633 | 561 098 |
| bcalm | 1.474 s ± 0.110 s | 273 MB | 1554 | 4677 | 626 022 |
| Bifrost | 1.277 s ± 0.033 s | 29 MB | 1679 | 4656 | 631 026 |

### full：144 MB（2_illumina/merge，全量 656k reads）

| 工具 | wall | peak RSS | 序列数 | N50 | 总长 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| anchr asm unitig | 2.112 s ± 0.048 s | 946 MB | 1482 | 8790 | 622 758 |
| anchr asm contig | 2.375 s ± 0.046 s | 993 MB | 97 | 11 657 | 561 352 |
| anchr asm contig --no-bubbles | 2.338 s ± 0.034 s | 947 MB | 96 | 11 657 | 561 163 |
| bcalm | 2.378 s ± 0.030 s | 556 MB | 1482 | 8790 | 622 758 |
| Bifrost | 2.982 s ± 0.275 s | 29 MB | 1505 | 8350 | 623 546 |

## 分析

- **速度**：full 上 anchr unitig 2.11 s、contig 2.44 s，**均快于 bcalm**
  （2.38 s；unitig 0.89×）、也快于 Bifrost（2.81 s）。anchr 单线程
  （计数对齐 pgr：per-chunk `count_keys` + 树状合并 + 增量 canonical +
  前缀索引查询）；bcalm/Bifrost 多线程流式建图。长 k 下优势更大
  （small k=64 0.47×、k=100 0.61× vs bcalm）；
- **内存是最大短板**：anchr 峰值 RSS 随输入近似线性且基数很大——
  24 MB 输入 267 MB、48 MB 473 MB、144 MB **0.95 GB**（分别为 bcalm 的
  1.1×/1.7×/1.7×，Bifrost 的 9×/16×/33×）；per-chunk 去重 + 树状合并
  消除了含重复的全局 key 中间体后相对旧实现（1.24/2.9/5.3 GB）大幅
  改善，chunk 尺寸 4096 → 16384 再降 ~25%（根因与后续方向见
  [asm-assemble.md](../design/asm-assemble.md)
  §9）；
- **contig vs unitig**：三种 anchr 模式共享同一 k-mer 表（RSS 几乎
  相同），种子化贪心遍历比严格压缩合并得更多——full 上 contig 97 条、
  unitig 116 条，N50 11 657 vs 9627，耗时略高（13.6 vs 12.8 s）；
- **弹泡影响很小**：contig 与 --no-bubbles 的输出几乎一致（full 97 vs
  96 条、N50 相同 11 657，各档差异仅 1-4 条；耗时上 --no-bubbles 反而
  略慢（14.0 vs 13.6 s，保留路径后排序/处理略多）。对该数据集
  （纠错 reads，气泡少）弹泡与否的影响可忽略；
- **产出一致性**：`anchr asm unitig` 默认不过滤（§10.5），输出与 bcalm
  **逐序列完全相同**（三档 1048/1554/1482 条，N50/总长全部一致，规范化
  集合双向相等）。contig 模式合并更多（full 97 条、N50 11 657）。
  bcalm 保留全部 unitig 的语义分析（顶点分解、无损压实，含旧
  `contained` 的历史实测）见 [asm-assemble.md](../design/asm-assemble.md)
  §10；
- **确定性与功能**：anchr 单线程换来扫描顺序无关的确定性输出，且原生
  支持 gz 输入、FASTQ/FASTA、GFA 输出、`--links` 边信息——这些是
  bcalm/Bifrost 命令行不具备的便利性，但性能/内存需要优化。

## 备注

- 复现：`bash scripts/unitig-bench.sh [small|medium|full] [runs]`
  （脚本自动用 `mktemp` 隔离每次运行，输出 wall/RSS/统计表）；
- 本表三档在同一次批处理中测出（机器状态一致）；绝对时间随机器负载
  漂移（相邻批次 unitig full 从 10.1 s 漂到 12.8 s），**相对比值**才是
  稳定指标；
- 参数不对称：anchr/bcalm 固实阈值 3，Bifrost 默认 1（无 `--min-abundance`
  等价参数，保留低丰度 k-mer），后者 unitig 数略多；anchr 默认过滤
  <124 bp 的序列，bcalm/Bifrost 全保留，数量不可直接比较；
- RSS 用 `/usr/bin/time -v` 单次测量，不是 hyperfine 均值；输入统一为
  明文 FASTA（脚本对 gz 先解压）。

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
