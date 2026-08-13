#!/usr/bin/env bash

# unitig-bench.sh
# Benchmark unitig generation: `anchr asm unitig`, `anchr asm contig` and
# `anchr asm contig --no-bubbles` (anchr's own k-mer graph compaction /
# seeded traversal, with and without bubble popping) vs the raw bcalm and
# Bifrost tools, plus cuttlefish 2 (KMC3 + MPHF + DFA-state route), on the
# G37 simulated E. coli reads.
#
# All tools produce unitigs/contigs from solid k-mers, with these
# parameter/default differences (kept as each tool is normally used):
#   anchr asm unitig  --kmer 31 --min-count-seed 3 (like bcalm
#                      -abundance-min 3), min-contig-len default 0 (keep all
#                      unitigs, bcalm lossless compaction), single-threaded
#   anchr asm contig  --kmer 31 --no-bubbles (seeded traversal, parallel
#   --no-bubbles      paths kept separate), same defaults, single-threaded
#   anchr asm contig  --kmer 31 (seeded traversal, bubbles popped)
#   bcalm             -kmer-size 31 -abundance-min 3 -nb-cores N
#   Bifrost build     --kmer-length 31 --clip-tips --del-isolated -t N
#   cuttlefish build  -k 31 -t N --ref -c 3 (FASTA reference mode; note it
#                      filters (k+1)-mers while bcalm/anchr filter k-mers)
#
# Each configuration runs in a fresh mktemp directory. Wall time is
# measured with hyperfine; peak RSS with /usr/bin/time -v; output
# statistics with `pgr fa n50`.
#
# Usage:
#   bash scripts/unitig-bench.sh [small|medium|full] [runs]
#   G37=/path/to/g37 bash scripts/unitig-bench.sh full 5
#
# Requires in PATH: anchr (release preferred), bcalm, Bifrost, pgr,
# hyperfine, /usr/bin/time; cuttlefish via PATH or CUTTLEFISH=... .

set -euo pipefail

cd "$(dirname "$0")/.."

ANCHR="${ANCHR:-}"
if [ -z "$ANCHR" ] && [ -x "$PWD/target/release/anchr" ]; then
    ANCHR="$PWD/target/release/anchr"
else
    ANCHR="${ANCHR:-$PWD/target/debug/anchr}"
fi
CUTTLEFISH="${CUTTLEFISH:-$PWD/cuttlefish-2.2.0/bin/cuttlefish}"

G37="${G37:-$HOME/data/anchr/g37}"
SCALE="${1:-small}"
RUNS="${2:-3}"
K="${K:-31}"

THREADS=8

case "$SCALE" in
    small)  IN="$G37/4_down_sampling/Q0L0X40P000/pe.cor.fa" ;;
    medium) IN="$G37/4_down_sampling/Q0L0X80P000/pe.cor.fa" ;;
    full)   IN="$G37/2_illumina/merge/pe.cor.fa.gz" ;;
    *) echo "unknown scale: $SCALE (small|medium|full)" >&2; exit 1 ;;
esac

[ -e "$IN" ] || { echo "missing input: $IN" >&2; exit 1; }

# Prepare a plain (non-gz) input so all three tools read the same file.
WORK="$(mktemp -d /tmp/unitig-bench.XXXXXX)"
PLAIN="$WORK/pe.cor.fa"
case "$IN" in
    *.gz) gzip -dc "$IN" > "$PLAIN" ;;
    *)    cp "$IN" "$PLAIN" ;;
esac

# Runner: runs one configuration in a fresh dir and echoes that dir path.
RUNNER="$WORK/runner.sh"
cat > "$RUNNER" <<EOF
#!/usr/bin/env bash
set -euo pipefail
# \$1 = unitig|contig|contig-nb|bcalm|bifrost
cfg="\$1"
d="\$(mktemp -d /tmp/ub.XXXXXX)"
cd "\$d"
ln -s "$PLAIN" pe.cor.fa
case "\$cfg" in
    unitig)
        "$ANCHR" asm unitig pe.cor.fa -o unitigs.fa --kmer $K \\
            --min-count-seed 3 >/dev/null 2>&1
        ;;
    contig-nb)
        "$ANCHR" asm contig pe.cor.fa -o contigs.fa --kmer $K \\
            --no-bubbles >/dev/null 2>&1
        ;;
    contig)
        "$ANCHR" asm contig pe.cor.fa -o contigs.fa --kmer $K >/dev/null 2>&1
        ;;
    bcalm)
        bcalm -in pe.cor.fa -kmer-size $K -abundance-min 3 -verbose 0 \\
            -nb-cores $THREADS -out K$K >/dev/null 2>&1
        ;;
    bifrost)
        Bifrost build --input-seq-file pe.cor.fa --kmer-length $K --clip-tips \\
            --del-isolated --threads $THREADS --fasta --no-compress-out \\
            --output-file bf >/dev/null 2>&1
        ;;
    cuttlefish)
        ulimit -n 4096
        "$CUTTLEFISH" build -s pe.cor.fa -k $K -t $THREADS -o K$K -w . \\
            --ref -c 3 >/dev/null 2>&1
        ;;
    *)
        echo "unknown config: \$cfg" >&2
        exit 1
        ;;
esac
echo "\$d"
EOF
chmod +x "$RUNNER"

result_path() { # $1 = config
    case "$1" in
        bcalm)   echo "K$K.unitigs.fa" ;;
        bifrost) echo "bf.fasta" ;;
        cuttlefish) echo "K$K.fa" ;;
        unitig)  echo "unitigs.fa" ;;
        *)       echo "contigs.fa" ;;
    esac
}

stats() { # $1 = fasta path; prints "n50 sum avg esize count"
    if [ ! -s "$1" ]; then
        echo "MISSING"
        return
    fi
    pgr fa n50 -H -C -S -A -E "$1" | tr '\n' ' '
    echo
}

bench() { # $1 = config label, $2 = config id
    local label="$1" cfg="$2"
    local wall rss run_dir
    wall=$(hyperfine --warmup 1 --runs "$RUNS" -N --ignore-failure \
        "$RUNNER $cfg" 2>&1 | grep -E "Time \(mean" \
        | sed -E 's/\[User:.*$//; s/.*mean ± σ\):\s+//; s/[[:space:]]+$//')
    rss=$(/usr/bin/time -v "$RUNNER" "$cfg" 2>&1 | grep "Maximum resident" | awk '{print $6}')
    run_dir=$("$RUNNER" "$cfg")
    printf "%-8s %-16s %10s  %-24s  out: %s\n" \
        "$SCALE" "$label" "${rss:-?} KB" "$wall" \
        "$(stats "$run_dir/$(result_path "$cfg")")"
}

echo "scale=$SCALE  k=$K  threads=$THREADS  runs=$RUNS  input=$IN"
echo "note: -p binds counting (+ --dfa classification); default half(<=8), auto = all cores; walk is single-threaded/deterministic"
echo
printf "%-8s %-16s %10s  %-24s  %s\n" "scale" "config" "peak RSS" "wall (mean±σ)" "output (n50 sum avg esize count)"
printf "%-8s %-16s %10s  %-24s  %s\n" "-----" "------" "--------" "----------------" "-----------------------------"

bench "anchr unitig"        "unitig"
bench "anchr contig"        "contig"
bench "anchr contig no-bub" "contig-nb"
bench "bcalm"               "bcalm"
bench "Bifrost"             "bifrost"
bench "cuttlefish"          "cuttlefish"

echo
echo "workdir=$WORK"
