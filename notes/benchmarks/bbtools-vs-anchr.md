# CLI benchmark: anchr vs BBTools 39.38

> 2026-08-13 随迁自 pgr（`benchmarks/bbtools-vs-pgr.md`）：fq 命令已迁至
> anchr，基准中的 fq 行改为 `anchr fq ...`；`kmer hist` 仍属 pgr，仅作参考。

Hyperfine 1.19.0 comparison of the anchr-native replacements against the
locally installed BBTools 39.38 (cbp package) on the Lambda test data
(40,000 reads, SRR5042715).

## Methodology

* `hyperfine --warmup 2 --runs 6` per pair; mean ± σ reported.
* anchr: parallel `target/release/anchr`; `anchr fq trim-adapter` runs with
  `--parallel 8` (matching BBTools); `anchr fq clump`, `anchr fq norm` (and pgr `kmer hist`)
  parallelize internally via rayon (logical CPU count, 32 on this host;
  `anchr fq norm` also accepts `--parallel`).
  BBTools: `threads=8`.
* Both sides write plain-text FASTQ (anchr does not gzip output); inputs are
  the same gzipped files (raw R1/R2, or the golden clumpify/trim/filter
  intermediates for the downstream steps).
* Parameters match the anchr trim pipeline (`k=31` clumpify, `k=23` trim,
  `k=27` filter, `k=31` khist, `min=3` bbnorm; `ordered=t`, fixed seeds).
* Updated 2026-08-10 after the `anchr fq norm` parallel + memory-bounded external
  path work and the `--threads` → `--parallel` rename.

## Results

| Step | BBTools 39.38 | anchr | Speedup |
|---|---:|---:|---:|
| clumpify / `anchr fq clump` | 211.5 ± 8.3 ms | 96.7 ± 4.7 ms | 2.19x |
| bbduk trim / `anchr fq trim-adapter` | 372.6 ± 20.6 ms | 84.5 ± 5.7 ms | 4.41x |
| bbduk filter / `anchr fq filter` | 192.2 ± 6.8 ms | 78.2 ± 10.0 ms | 2.46x |
| kmercountexact / pgr `kmer hist` (stays in pgr) | 2394.9 ± 345.7 ms | 686.9 ± 4.7 ms | 3.49x |
| repair / `anchr fq split` | 176.5 ± 2.9 ms | 21.7 ± 0.4 ms | 8.13x |
| reformat sample / `anchr fq sample` | 161.3 ± 2.2 ms | 24.2 ± 0.4 ms | 6.67x |
| bbnorm / `anchr fq norm` | 889.1 ± 46.7 ms | 228.7 ± 7.9 ms | 3.89x |

Sum of step means: BBTools ≈ 4.40 s, anchr ≈ 1.22 s (~3.6x end to end, without
the JVM/script startup and intermediate gzip files of the real pipeline).

## Output equivalence

`clump`, `trim`, `filter`, `split`, `sample`, and the khist/peaks text are
byte-identical between the two sides on these inputs. `anchr fq norm` keeps 39,846
reads vs bbnorm's 39,888 (exact canonical counts vs `bits=16` approximate
hash counts; see tests/bbtools/Lambda/README.md). The `--stats` files of
`anchr fq trim-adapter` match bbduk's 3-column `stats=` format byte for byte
(except the path-dependent `#File` line).

## Caveats

* `anchr fq trim-adapter` is the main worker-pool parallelized command (`--parallel`,
  default logical CPU count); on 50万-pair synthetic data it scaled from
  9.2 s (1 thread) to 1.4 s (8 threads), 6.6x, with ~15 MB peak memory.
* `anchr fq clump`, `anchr fq norm` (and pgr `kmer hist`) parallelize through rayon internally;
  the Lambda data is small, so their parallel gain is limited (clump even
  measured slightly slower than the earlier single-threaded run).
* `anchr fq norm` estimates its peak footprint from the input (48 B per base, gz
  expanded 8x) and switches to a memory-bounded external hash-bucket path
  (`--mem`, default 2g) when it exceeds the cap; on 50万-pair synthetic
  random data the external path ran ~8.5 s at a 256 MiB cap with ~210 MiB
  peak RSS, byte-identical to the in-memory path.
* `anchr fq split` / `anchr fq sample` are I/O-bound and remain single-threaded.
