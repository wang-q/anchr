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
#   L0 unit tests      cargo test (correctness / zero-panic) — not here.
#   L1 smoke           G37 MRX40P000 full multik (seconds, <1 GB): byte-diff
#                      unitig output vs golden + report N50/count. Catches
#                      "output must stay byte-identical" optimizations/refactors.
#   L2 single          G37 + MG1655 40x group, full multik->olc->extend
#                      (minutes): report N50/count across the chain. mis is
#                      checked only at full chain (single-group mis is expected;
#                      multi-group anchor voting resolves it, per model_org.md).
#   L3 full            G37/MG1655/DH5alpha full-chain merge + quast (hours):
#                      final gate before merge.
#
# Gate philosophy:
#   - mis = hard fail (authoritative per results/model_org.md); checked at L3.
#   - N50/unitig-count/peak-RSS = regression report (soft): large drift
#     warns, small drift is informational — do NOT auto-block legitimate
#     changes that pay a small cost to fix a chimera.
#
# Usage:
#   bash scripts/asm-gate.sh [smoke|single|full]
#   bash scripts/asm-gate.sh smoke --write   # (re)capture golden baseline
#   G37=... MG1655=... bash scripts/asm-gate.sh full
#
# Requires in PATH: anchr release binary, faops, quast.py. Baselines live in
#   results/asm_gate.md (single source of truth).

set -euo pipefail

cd "$(dirname "$0")/.."
BIN="${ANCHR:-$PWD/target/release/anchr}"
QUAST="${QUAST:-/home/wangq/.cbp/bin/quast.py}"

#--- datasets (override via env) ---
G37_DATA="${G37:-$HOME/data/anchr/g37}"
MG1655_DATA="${MG1655:-$HOME/data/anchr/mg1655}"
DH5_DATA="${DH5ALPHA:-$HOME/data/anchr/dh5alpha}"
G37_REF="${G37:-$HOME/data/anchr/ref/g37/genome.fa}"
MG1655_REF="${MG1655:-$HOME/data/anchr/ref/mg1655/genome.fa}"
DH5_REF="${DH5ALPHA:-$HOME/data/anchr/ref/dh5alpha/genome.fa}"

#--- baseline tables ---
# Authoritative baselines (full chain, quast --min-contig 10), single source
# of truth = results/asm_gate.md. Keep in sync.
#   G37 7g : 0 mis / N50 121382 / GF 98.674 / Total 581339
#   MG 5g  : 0 mis / N50 112557 / GF 98.197 / Total 4617679
#   DH 13g : 0 mis / N50  99473 / GF 97.800 / Total 4496026

WORK="$(mktemp -d /tmp/asm-gate.XXXXXX)"
trap 'rm -rf "$WORK"' EXIT
declare -a FAILS=()
declare -a WARNS=()
declare -a NOTES=()

warn() { WARNS+=("$*"); echo >&2 -e "  \033[33mWARN\033[0m $*"; }
note() { NOTES+=("$*"); echo >&2 -e "  \033[32mnote\033[0m $*"; }
fail() { FAILS+=("$*"); echo >&2 -e "  \033[31mFAIL\033[0m $*"; }

# Prints "count n50 total" (N50 & total from `faops n50 -H -S`; count = seqs).
count_n50_total() { # fasta
    [ -s "$1" ] || { echo "0 0 0"; return; }
    local c n s
    c=$(grep -c '^>' "$1")
    read -r n < <(faops n50 -H -N 50 "$1" 2>/dev/null)
    read -r s < <(faops n50 -H -N 50 -S "$1" 2>/dev/null | tail -1)
    echo "$c ${n:-0} ${s:-0}"
}

run_multik() { # reads workdir parallel
    (cd "$2" && "$BIN" asm multik "$1" --all-masters -p "${3:-16}" -o unitigs_all.fasta >/dev/null 2>&1)
}

#============================================================================#
# L1 smoke — G37 MRX40P000 multik, byte-diff vs golden (seconds, <1 GB)
#============================================================================#
run_smoke() {
    note "L1 smoke: G37 MRX40P000 multik (auto) — byte-diff vs golden"
    local reads="$G37_DATA/6_down_sampling/MRX40P000/pe.cor.fa"
    [ -e "$reads" ] || { fail "missing reads: $reads"; return; }
    run_multik "$reads" "$WORK" 8
    local cur="$WORK/unitigs_all.fasta"

    # Read golden md5 from results/asm_gate.md (any line containing golden-md5).
    local golden
    golden=$(grep -E 'golden-md5[[:space:]`]*[0-9a-f]{32}' results/asm_gate.md 2>/dev/null | grep -oE '[0-9a-f]{32}' | head -1 || true)

    if [ "${1:-}" = "--write" ] || [ -z "$golden" ]; then
        local m; m=$(md5sum "$cur" | awk '{print $1}')
        local st; st=$(count_n50_total "$cur")
        note "  captured: count/n50/total = $st, md5 = $m"
        fail "  baseline not yet recorded; run with --write, or add golden-md5 to results/asm_gate.md"
        return
    fi

    local curmd5; curmd5=$(md5sum "$cur" | awk '{print $1}')
    local stats; stats=$(count_n50_total "$cur")
    note "  current: count/n50/total = $stats, md5 = $curmd5"
    if [ "$curmd5" == "$golden" ]; then
        note "  L1 byte-identical to golden ✓"
    else
        warn "L1 golden changed ($golden → $curmd5). If intended, re-capture with smoke --write."
    fi
}

#============================================================================#
# L2 single — G37 + MG1655 40x group, multik->olc->extend (minutes)
#============================================================================#
run_single() {
    note "L2 single: G37 + MG1655 MRX40P000 multik->olc->extend"
    local gdir="$WORK/g37" mdir="$WORK/mg"
    mkdir -p "$gdir" "$mdir"
    local greads="$G37_DATA/6_down_sampling/MRX40P000/pe.cor.fa"
    local mreads="$MG1655_DATA/6_down_sampling/MRX40P000/pe.cor.fa"

    run_multik "$greads" "$gdir" 8
    run_multik "$mreads" "$mdir" 8

    for d in "$gdir" "$mdir"; do
        local reads; reads=$([[ "$d" == "$gdir" ]] && echo "$greads" || echo "$mreads")
        (cd "$d" && "$BIN" asm olc --unitigs unitigs_all.fasta \
            --min-overlap 1000 --min-contig-len 200 -o unitigs.fasta >/dev/null 2>&1)
        (cd "$d" && "$BIN" asm extend unitigs.fasta "$reads" \
            --min-len 1000 -o unitigs.ext.fasta >/dev/null 2>&1)
        mv "$d/unitigs.ext.fasta" "$d/unitigs.fasta"
        note "  $(basename "$d") chain: unitigs_all $(count_n50_total "$d/unitigs_all.fasta")"
        note "                     final  $(count_n50_total "$d/unitigs.fasta")"
    done
}

#============================================================================#
# L3 full — full-chain merge + quast for each dataset (hours)
#============================================================================#
run_full() {
    note "L3 full: G37 / MG1655 / DH5alpha full-chain gates (hours)"
    warn "L3 requires the per-group multik->olc->extend->anchor runs; see /tmp/run_full.sh"
    fail "L3 not yet automated end-to-end in this script; run the /tmp/*_full chains per results/model_org.md"
}

LEVEL="${1:-smoke}"
shift || true
case "$LEVEL" in
    smoke)  run_smoke "$@" ;;
    single) run_single "$@" ;;
    full)   run_full "$@" ;;
    *) echo "usage: $0 [smoke|single|full]" >&2; exit 2 ;;
esac

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