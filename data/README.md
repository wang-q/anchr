# Reference sequence data

统一存放 anchr 命令使用的参考序列（adapter / contaminant / artifact）与
对应许可声明。所有命令通过 `include_str!` 从本目录内嵌，不散落在
`templates/`、`src/` 或测试数据中。

## 文件索引

| 文件 | 用途 | 命令 | 来源 | 许可 |
| :--- | :--- | :--- | :--- | :--- |
| `illumina_adapters.fa.gz`（158 条） | 接头去除/过滤 | `anchr trim`、`fq clean`、`fq filter`（`--ref`） | BBTools `resources/adapters.fa` | BBTools free |
| `sequencing_artifacts.fa.gz`（164 条 + PhiX） | 污染/假象过滤 | `anchr trim`（`--artifact`） | BBTools `sequencing_artifacts.fa.gz` + NCBI PhiX174 | BBTools free + 公共领域 |
| `adapter_list.txt`（6 条） | 接头检测（Adapter Content） | `fq qc` | FastQC `Configuration/adapter_list.txt` | GPL v3 |
| `contaminant_list.txt`（151 条） | 过度代表序列来源标注 | `fq qc`（Overrepresented） | FastQC `Configuration/contaminant_list.txt` | GPL v3 |
| `FastQC_DATA_LICENSE` | GPL v3 全文 | — | GNU（随 FastQC 数据分发） | GPL v3 |

## 许可说明

- anchr 本身是 MIT；本目录中 **FastQC 来源的序列数据（adapter_list.txt、
  contaminant_list.txt）按 GPL v3 单独声明**（对齐 Falco 的做法：MIT 代码 +
  `data/Configuration/FastQC_DATA_LICENSE`）。
- BBTools 数据（`illumina_adapters.fa`、`sequencing_artifacts.fa` 的 BBTools
  部分）按 BBTools "free for unlimited use" 使用；PhiX174 序列来自 NCBI
  （公共领域）。
- 修改本目录数据时保持来源与许可声明。
- 大序列（两个 `.fa.gz`）压缩存储：`trim.rs` 用 `include_bytes!` 内嵌 +
  flate2 解压，运行时写出平文件（`--adapter`/`--artifact` 默认行为不变）。

## 添加新参考序列

按用途放本目录（命名保持 `*.fa` / `*.txt` 后缀），在 README 表格补一行，
并注明来源与许可。代码侧用 `include_str!("../../data/<file>")` 内嵌。
