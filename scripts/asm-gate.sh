#!/usr/bin/env bash
# asm-gate.sh — layered quality gates for `anchr asm multik` regressions.
#
# Motivation (see notes/todo.md):
#   The full-chain quast gate (G37/MG1655/DH5alpha) is hours-long and
#   memory-heavy, so it only runs before merge — meanwhile N50/memory/
#   unitig-count regressions (run>=2 over-pruning, auto-k formula, K-table
#   caching OOM, olc --unitigs miss, cross-contig guardrail cost, ...) went
#   undetected by unit tests. These numeric regressions need a cheap inner
#   layer that reports early on every change, with the expensive full chain
#   as the final gate.
#
# Design (4 tiers):
#   L0 unit tests   cargo test (correctness / zero-panic) — not here.
#   L1 smoke        G37 MRX40P000 full multik (seconds, 2 GB): byte-diff
#                   unitig output vs golden + report N50/count. Catches
#                   "output must stay byte-identical" optimizations/refactors.
#   L2 single       G37 + MG1655 MRX40P000 multik->olc->extend (minutes):
#                   report count/N50/Total + peak RSS per dataset, diff vs
#                   baseline. Soft thresholds warn on N50 drop / unitig blowup /
#                   RSS over-limit. mis is NOT checked here (single-group mis
#                   is expected; multi-group anchor voting resolves it).
#   L3 full         per-dataset full chain multik->olc->extend->anchor ->
#                   cross-group merge -> quast (hours): 0 mis is a HARD red
#                   line; N50/GF drift is reported.
#
# Gate philosophy:
#   - mis = hard fail (authoritative per results/model_org.md); checked at L3.
#   - N50/unitig-count/peak-RSS = regression report (soft): large drift
#     warns, small drift is informational — do NOT auto-block legitimate
#     changes that pay a small cost to fix a chimera.
#
# Usage:
#   bash scripts/asm-gate.sh [smoke|single|full|all]
#   bash scripts/asm-gate.sh smoke --write   # (re)capture golden baseline
#   bash scripts/asm-gate.sh full <g37|mg1655|dh5alpha|all>
#   bash scripts/asm-gate.sh all             # smoke + single + full all
#
# Requires in PATH: anchr release binary, faops, quast.py. Baselines live in
#   scripts/asm_gate.md (single source of truth, sibling to this script).

set -euo pipefail

cd "$(dirname "$0")/.."
BIN="${ANCHR:-$PWD/target/release/anchr}"
QUAST="${QUAST:-/home/wangq/.cbp/bin/quast.py}"
REF="${REF:-$HOME/data/anchr/ref}"
GATE_MD="${GATE_MD:-$PWD/scripts/asm_gate.md}"

#--- datasets (override via env) ---
G37_DATA="${G37:-$HOME/data/anchr/g37}"
MG1655_DATA="${MG1655:-$HOME/data/anchr/mg1655}"
DH5_DATA="${DH5ALPHA:-$HOME/data/anchr/dh5alpha}"

#--- per-dataset MR group sets (assembly + gate merge closure) ---
# g37:    7 groups  MRX40P000-004, MRX80P000-001
# mg1655: gate uses 5 groups MRX40P000/P001/P002 + MRX80P000/P001
# dh5alpha: 13 groups MRX40P000-008, MRX80P000-003
G37_GROUPS="MRX40P000 MRX40P001 MRX40P002 MRX40P003 MRX40P004 MRX80P000 MRX80P001"
MG_GROUPS="MRX40P000 MRX40P001 MRX40P002 MRX80P000 MRX80P001"
MGDH_ALL_GROUPS="MRX40P000 MRX40P001 MRX40P002 MRX40P003 MRX40P004 MRX40P005 MRX40P006 MRX40P007 MRX40P008 MRX80P000 MRX80P001 MRX80P002 MRX80P003"

#--- soft threshold multipliers (L2/L3 soft report) ---
# N50 below baseline*N50_DROP -> warn; unitig count above baseline*COUNT_BLOW -> warn;
# peak RSS above RSS_GB -> warn.
N50_DROP="${N50_DROP:-0.6}"
COUNT_BLOW="${COUNT_BLOW:-3}"
RSS_GB="${RSS_GB:-20}"

# Persistent work for L3 (so full-chain results are resumable/reusable across
# runs, like the /tmp/*_full dirs). L1/L2 use a throwaway temp dir.
FULL_WORK="${FULL_WORK:-/tmp/asm-gate-full}"
WORK="$(mktemp -d /tmp/asm-gate.XXXXXX)"
declare -a FAILS=()
declare -a WARNS=()
declare -a NOTES=()

warn() { WARNS+=("$*"); echo >&2 -e "  \033[33mWARN\033[0m $*"; }
note() { NOTES+=("$*"); echo >&2 -e "  \033[32mnote\033[0m $*"; }
fail() { FAILS+=("$*"); echo >&2 -e "  \033[31mFAIL\033[0m $*"; }
finish() { rm -rf "$WORK"; }

# Prints "count n50 total" (N50 & total from `faops n50 -H -S`; count = seqs).
count_n50_total() { # fasta
    [ -s "$1" ] || { echo "0 0 0"; return; }
    local c n s
    c=$(grep -c '^>' "$1")
    read -r n < <(faops n50 -H -N 50 "$1" 2>/dev/null)
    read -r s < <(faops n50 -H -N 50 -S "$1" 2>/dev/null | tail -1)
    echo "$c ${n:-0} ${s:-0}"
}

# Runs multik, capturing peak RSS (KB) in $PEAK_RSS_KB. Use /usr/bin/time -v.
run_multik_rss() { # reads outdir parallel
    local reads="$1" dir="$2" par="${3:-8}"
    mkdir -p "$dir"
    ( cd "$dir" && /usr/bin/time -v "$BIN" asm multik "$reads" \
        --all-masters -p "$par" -o unitigs_all.fasta >/dev/null 2>time.log )
    PEAK_RSS_KB=$(awk -F' *Maximum resident set size \\(kbytes\\): *' '/Maximum resident/{print $2}' "$dir/time.log" 2>/dev/null || echo 0)
    local rss_gb; rss_gb=$(awk -v k="${PEAK_RSS_KB:-0}" 'BEGIN{printf "%.1f", k/1048576}')
    local st; st=$(count_n50_total "$dir/unitigs_all.fasta")
    note "    multik $(basename "$dir"): count/N50/Total = $st, peak RSS = ${rss_gb} GB"
    if awk -v g="${rss_gb}" -v cap="${RSS_GB}" 'BEGIN{exit !(g>cap)}'; then
        warn "multik peak RSS ${rss_gb} GB exceeds soft cap ${RSS_GB} GB"
    fi
}

# Wraps run_multik_rss with the smoke-golden byte check.
run_multik_here() { # reads outdir parallel
    run_multik_rss "$1" "$2" "$3"
}

#============================================================================#
# L1 smoke — G37 MRX40P000 multik, byte-diff vs golden (seconds, ~2 GB)
#============================================================================#
run_smoke() {
    note "L1 smoke: G37 MRX40P000 multik (auto) — byte-diff vs golden"
    local reads="$G37_DATA/6_down_sampling/MRX40P000/pe.cor.fa"
    [ -e "$reads" ] || { fail "missing reads: $reads"; return; }
    run_multik_here "$reads" "$WORK/g37" 8

    local golden
    golden=$(grep -E 'golden-md5[[:space:]`]*[0-9a-f]{32}' "$GATE_MD" 2>/dev/null | grep -oE '[0-9a-f]{32}' | head -1 || true)

    local cur="$WORK/g37/unitigs_all.fasta"
    if [ "${1:-}" = "--write" ] || [ -z "$golden" ]; then
        local m; m=$(md5sum "$cur" | awk '{print $1}')
        fail "  golden-md5 not recorded in $GATE_MD (current md5 = $m). Re-run with smoke --write after updating the doc."
        return
    fi

    local curmd5; curmd5=$(md5sum "$cur" | awk '{print $1}')
    local stats; stats=$(count_n50_total "$cur")
    note "  current: count/N50/Total = $stats"
    if [ "$curmd5" == "$golden" ]; then
        note "  L1 byte-identical to golden ✓"
    else
        warn "L1 golden changed ($golden → $curmd5). If intended, re-capture with smoke --write."
    fi
}

#============================================================================#
# L2 single — G37 + MG1655 MRX40P000 multik->olc->extend + peak RSS (minutes)
#============================================================================#
run_single() {
    note "L2 single: G37 + MG1655 MRX40P000 multik->olc->extend"
    local gd="$WORK/g37_single" md="$WORK/mg_single"
    local greads="$G37_DATA/6_down_sampling/MRX40P000/pe.cor.fa"
    local mreads="$MG1655_DATA/6_down_sampling/MRX40P000/pe.cor.fa"

    run_multik_here "$greads" "$gd" 8
    run_multik_here "$mreads" "$md" 8

    for ent in "g37|$gd|$greads|37655|2745" "mg1655|$md|$mreads|40990|8562"; do
        local tag dir reads bn bcount
        IFS='|' read -r tag dir reads bn bcount <<<"$ent"
        ( cd "$dir" && "$BIN" asm olc --unitigs unitigs_all.fasta \
            --min-overlap 1000 --min-contig-len 200 -o unitigs.fasta >/dev/null 2>&1 )
        ( cd "$dir" && "$BIN" asm extend unitigs.fasta "$reads" \
            --min-len 1000 -o unitigs.ext.fasta >/dev/null 2>&1 )
        mv "$dir/unitigs.ext.fasta" "$dir/unitigs.fasta"
        local st; st=$(count_n50_total "$dir/unitigs.fasta")
        note "  $tag .final: count/N50/Total = $st"
        # soft check vs multik-stage baseline (unitigs_all N50 / count)
        local n50 bN50 c bc
        read -r c n50 s <<<"$(count_n50_total "$dir/unitigs_all.fasta")"
        bN50="$bn"; bc="$bcount"
        if awk -v cur="$n50" -v base="$bN50" -v m="$N50_DROP" 'BEGIN{cur+=0;base+=0;exit !(base>0 && cur < base*m)}'; then
            warn "$tag unitigs_all N50 $n50 < ${bN50}*${N50_DROP} (soft)"
        fi
        if awk -v cur="$c" -v base="$bc" -v m="$COUNT_BLOW" 'BEGIN{exit !(base>0 && cur > base*m)}'; then
            warn "$tag unitigs_all count $c > ${bc}*${COUNT_BLOW} (soft)"
        fi
    done
}

#============================================================================#
# L3 full — per-dataset full chain + cross-group merge + quast (hours)
#============================================================================#
chain_group() { # org data groups workdir maxjobs
    local org="$1" data="$2" out="$3" maxjobs="${4:-2}"
    mkdir -p "$out/6_unitigs_multik"
    local pids=()
    for g in ${5}; do
        local d="$out/6_unitigs_multik/$g"
        [ -e "$d/anchor/anchor.fasta" ] && { note "  $org/$g reused"; continue; }
        mkdir -p "$d"
        local reads="$data/6_down_sampling/$g/pe.cor.fa"
        [ -e "$reads" ] || { fail "$org/$g: missing reads"; continue; }
        (
            cd "$d"
            "$BIN" asm multik "$reads" --all-masters -p 16 -o unitigs_all.fasta
            "$BIN" asm olc --unitigs unitigs_all.fasta \
                --min-overlap 1000 --min-contig-len 200 -o unitigs.fasta
            "$BIN" asm extend unitigs.fasta "$reads" --min-len 1000 -o unitigs.ext.fasta
            mv unitigs.ext.fasta unitigs.fasta
            mkdir -p anchor
            "$BIN" asm anchor unitigs.fasta "$reads" --uscale 2 --lscale 3 \
                -p 4 -o anchor/anchor.fasta
        ) >/dev/null 2>&1 &
        pids+=($!)
        while [ "${#pids[@]}" -ge "$maxjobs" ]; do
            for i in "${!pids[@]}"; do kill -0 "${pids[$i]}" 2>/dev/null || unset 'pids[$i]'; done
            pids=("${pids[@]}"); [ "${#pids[@]}" -lt "$maxjobs" ] && break; sleep 3
        done
    done
    wait
}

merge_and_quast() { # org out ref groups
    local org="$1" out="$2" ref="$3" groups="${4:-}"
    local merge; merge="$out/7_merge_mr_unitigs_multik"
    mkdir -p "$merge"
    for g in $groups; do
        [ -s "$out/6_unitigs_multik/$g/anchor/anchor.fasta" ] && echo "$out/6_unitigs_multik/$g/anchor/anchor.fasta"
    done | sort -r > "$merge/anchors.list"
    local n; n=$(grep -c . "$merge/anchors.list")
    [ "$n" -ge 2 ] || { fail "$org: not enough anchors ($n) to merge"; return; }
    "$BIN" asm olc --unitigs --list-files "$merge/anchors.list" \
        --min-overlap 1000 --min-contig-len 1000 \
        -o "$merge/anchor.merge.fasta" >/dev/null 2>&1
    ( cd "$merge" && "$QUAST" --threads 8 --min-contig 10 -r "$ref" \
        -o quast_$org anchor.merge.fasta >/dev/null 2>&1 )
    local r="$merge/quast_$org/report.tsv"
    [ -s "$r" ] || { fail "$org: quast failed"; return; }
    # quast report.tsv is "<metric>\t<value>" (NF=2) in this version.
    local n50 gf mis total contigs
    n50=$(awk -F'\t' '$1=="N50"{print $2}' "$r")
    gf=$(awk -F'\t' '$1=="Genome fraction (%)"{print $2}' "$r")
    mis=$(awk -F'\t' '$1=="# misassemblies"{print $2}' "$r")
    total=$(awk -F'\t' '$1=="Total length (>= 0 bp)"{print $2}' "$r")
    contigs=$(awk -F'\t' '$1=="# contigs (>= 0 bp)"{print $2}' "$r")
    note "  $org L3: #contigs=$contigs Total=$total N50=$n50 mis=$mis GF=${gf}%"
    if [ "${mis:-x}" = "0" ]; then note "  $org mis = 0 ✓"; else fail "$org: misassemblies = $mis (RED LINE, expect 0)"; fi
}

run_full() {
    local which="${1:-all}"
    case "$which" in
        all) run_full g37; run_full mg1655; run_full dh5alpha ;;
        g37)
            note "L3 full: G37 (7 groups)"
            chain_group g37 "$G37_DATA" "$FULL_WORK/g37" 2 "$G37_GROUPS"
            merge_and_quast g37 "$FULL_WORK/g37" "$REF/g37/genome.fa" "$G37_GROUPS" ;;
        mg1655)
            note "L3 full: MG1655 (gate 5 groups of 13)"
            chain_group mg1655 "$MG1655_DATA" "$FULL_WORK/mg1655" 2 "$MGDH_ALL_GROUPS"
            merge_and_quast mg1655 "$FULL_WORK/mg1655" "$REF/mg1655/genome.fa" "$MG_GROUPS" ;;
        dh5alpha)
            note "L3 full: DH5alpha (13 groups, no cv)"
            chain_group dh5alpha "$DH5_DATA" "$FULL_WORK/dh5alpha" 2 "$MGDH_ALL_GROUPS"
            merge_and_quast dh5alpha "$FULL_WORK/dh5alpha" "$REF/dh5alpha/genome.fa" "$MGDH_ALL_GROUPS" ;;
        *) fail "unknown dataset '$which'"; return ;;
    esac
}

runner() {
    local lvl="$1"; shift
    case "$lvl" in
        smoke)  run_smoke "$@" ;;
        single) run_single "$@" ;;
        full)   run_full "$@" ;;
        all)    run_smoke; run_single; run_full all ;;
        *) echo "usage: $0 [smoke|single|full <ds>|all]" >&2; exit 2 ;;
    esac
}

trap 'finish' EXIT
LEVEL="${1:-smoke}"
shift || true
runner "$LEVEL" "$@"

echo
if [ "${#FAILS[@]}" -gt 0 ]; then
    echo "GATE FAILED:"
    printf '  - %s\n' "${FAILS[@]}"
    exit 1
else
    echo "GATE PASSED (${#WARNS[@]} warnings, ${#NOTES[@]} notes)"
    [ "${#WARNS[@]}" -gt 0 ] && { echo "warnings:"; printf '  - %s\n' "${WARNS[@]}"; }
    exit 0
fi