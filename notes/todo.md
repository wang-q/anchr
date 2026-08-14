# anchr 待办

> 仅保留 actionable 待办；历史结论、已完成项与细节见各 design/benchmark/
> audit 文档（索引见 `project-understanding.md`；2026-08-15 会话交接见
> `design/asm-olc-modern-flow.md`；明确不做清单见
> `project-understanding.md` §6.4）。

## 待实现

- `fq range` 剩余：BGZF `.gzi` 免预生成（封装
  `pgr::libs::bgzf::build_gzi_index` 或 pgr 补 CLI）；
  双端感知 S2 已完成（2026-08-15 核对，`design/fq-range.md` §7）；
- `dep`/`ena`/`template` 的外部工具版本核对（dazzler/hnsm 预装清单）。

## 挂账 / 待决

- pgr rev 更新流程 + `.cargo/config.toml` patch 发布差异
  （`project-understanding.md` §8.1）；
- pgr supermer 质量门控：落地后 FASTQ 可去掉 direct 回退
  （`asm-assemble.md` §12.1）；
- pgr 并行读 gz（`asm-assemble.md` §12.3，风险>收益暂缓）；

## 待验证 / 等数据到位

- 真实宏基因组/长读（决定拆分合并路线定位）：multik tip/weak_link 价值、
  `asm olc` 四命令端到端、anchor 补洞（老流程 fill）+ others 输出、
  megahit low depth/本地组装与 SKESA SNP 簇变体优先级
  （`megahit.md` §8、`skesa.md` §7.2）、SAM 内存化评估（`asm-olc.md` §14.6）；
- multik `--parallel` 扩展性复测、`--min-count-extend` 调参
  （`asm-multik.md` §7）；
- multik 性能：计数复用 / remove_unsupported 查表化 / 轮数裁剪
  （`benchmarks/multik-complexity.md`）；
- multik 防 misassembly：bridge 探针长度/阈值调参、<500 bp 碎片输出策略
  （`asm-multik.md` §9）；
- `asm cns` contain 吸收（`asm-olc.md` §14.3-3/§S3：列投票去冗余同时保留
  覆盖证据）；
- `asm olc` 参数验证：overlap 少量错配、不同 k unitig 冗余去重、repeat
  breaking 覆盖度阈值（`references/canu.md` §8）；
- 大规模真实数据全链：`fq → asm → map → template` 统计核对（覆盖量/
  unitig 数/PSL 行数）；
- gz/大输入回归：默认 supermer 路径全链回归 + 峰值内存；
- 560 bp 碎片 mis 覆盖度门槛（`asm-olc.md` §14.3，`--min-contig-len 1000`
  可滤，可选）。

## 低风险审计（可顺手修）

- ~~`scripts/verify-migrate.sh` 的 `asm_olc` 用例加速~~（2026-08-15 评估：
  输入已是 `R1.2k`（2k reads 缩小版），约 6 s 是双 k OLC 合理开销；脚本
  已历史化（仅旧 pgr 二进制可用），加速价值低，不做）。
- warning 清理（2026-08-15 部分完成）：`tadpole` 的
  `error_extension_pincer/tail` 未读字段已删；`Overlap` 的 dead_code 是
  bin/lib 双份 `libs`（cmd 混用 `crate::libs` 与 `anchr::`）的结构问题，
  非字段级清理，保留 `#[allow(dead_code)]`，待统一引用重构时一并解决。

## 技术债

- golden 数据体积：`tests/bbtools/Lambda/` 约 13 MB，可考虑外部数据源或
  生成脚本；
- 外部工具链依赖：流程命令强依赖 dazzler/hnsm（PATH），测试缺失时跳过
  可能掩盖回归（`project-understanding.md` §8.3）。
