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
| G37     | 7 | 13 | 187498 | 581339 | 121382 | 0 | 98.674 | 1.001 |
| MG1655  | 5 | 91 | 268281 | 4617679 | 112557 | 0 | 98.197 | 1.013 |
| DH5alpha| 13 | 105 | 258601 | 4496026 | 99473 | 0 | 97.800 | 1.003 |