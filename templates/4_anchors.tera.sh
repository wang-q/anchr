{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 4_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX]"

DIR_PREFIX=${1:-"4_unitigs_multik"}

# Plain bash loops replace the former `parallel -j 2 ... :::` argument grid
# (nested double-quoted templates need \$ / \" escaping that breaks
# silently). Jobs run one at a time in the foreground, each anchr asm
# anchor with the full opt.parallel thread budget.
for Q in 0 {{ opt.qual }}; do
    for L in 0 {{ opt.len }}; do
        for X in {{ opt.cov }}; do
            for P in $(printf "%03d " {0..{{ opt.splitp }}}); do
                (
    if [ ! -e 4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz ]; then
        exit;
    fi

    echo >&2 "==> ${DIR_PREFIX}/Q${Q}L${L}X${X}P${P}"
    if [ -e ${DIR_PREFIX}/Q${Q}L${L}X${X}P${P}/anchor.fasta ]; then
        echo >&2 '    anchor.fasta already presents'
        exit;
    fi

    if [ ! -s ${DIR_PREFIX}/Q${Q}L${L}X${X}P${P}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta does not exist or is empty'
        exit;
    fi

    mkdir -p ${DIR_PREFIX}/Q${Q}L${L}X${X}P${P}
    cd ${DIR_PREFIX}/Q${Q}L${L}X${X}P${P}

    anchr asm anchor \
        unitigs.fasta \
        ../../4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz \
        --mincov 5 --mscale 3 \
        --lscale {{ opt.lscale }} \
        --uscale {{ opt.uscale }} \
        -p {{ opt.parallel }} \
        --stats anchor.stats.tsv \
        -o anchor.fasta

    echo >&2
                )
            done
        done
    done
done
