# `asm multik` 基准（2026-08-14）

> `anchr asm multik`（multi-k 迭代组装，设计见 `design/asm-multik.md`）的
> 吞吐 sanity check。release build（`cargo build --release`），单线程计数
> 路径（`--parallel 0`），合成长读数据。

## 数据与命令

* 1 Mb 随机基因组（无重复，固定种子）；
* 15 kb HiFi 式长 reads：40× 随机 + 150 条环状覆盖（跨起点/从起点），
  0.1% 替换错误；
* 命令：`anchr asm multik m1mb.fa -o m1mb_out.fa`（`--kmer auto` → 31/61/91/121）。

## 结果

| 指标 | 值 |
|---|---|
| 输出 | **单条 1,000,000 bp = 基因组 100%**（N50 1,000,000，总长 1,000,000） |
| 墙钟 | 9.41 s（release，`/usr/bin/time -v`） |
| CPU | 15.35 s user + 0.91 s system |
| 峰值内存 | 816 MB（RSS） |
| 吞吐 | ≈ 106 kb 基因组 / s |

## 解读

* **正确性**：1 Mb 单条零缺口（此前已验证 k-mer 多重集与基因组一致 =
  环状旋转等价）——multik 在 1 Mb 规模保持 100% 完整无 N；
* **扩展性**：内存 ≈ 0.8 GB / 1 Mb（计数表 + 图 + reads 缓冲）；真实细菌
  基因组（~5 Mb）预计 ~45 s / ~4 GB，可接受；宏基因组（Gb 级）需外部
  分桶（`design/asm-assemble.md` §13，todo §2）；
* **对比**：`asm olc`（并行多 k + 启发式）Lambda 20k 最长 19 kb vs
  multik 46.5 kb（短读）/ 100%（长读）——迭代路线在连通性与完整性上
  全面占优（无冗余 + 计数证据选边）。

## 待补（真实数据）

* 真实宏基因组/长读数据的端到端验证（覆盖完整 + 无 gap + 无嵌合）；
* `--parallel` 计数并行在更大数据上的扩展性复测（todo §4）。
