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

{% set parallel2 = opt.parallel | int / 2 -%}
{% set parallel2 = parallel2 | round(method="floor") -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
parallel --no-run-if-empty --linebuffer -k -j 2 "
    if [ ! -e 4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz ]; then
        exit;
    fi

    echo >&2 '==> ${DIR_PREFIX}/Q{1}L{2}X{3}P{4}'
    if [ -e ${DIR_PREFIX}/Q{1}L{2}X{3}P{4}/anchor.fasta ]; then
        echo >&2 '    anchor.fasta already presents'
        exit;
    fi

    if [ ! -s ${DIR_PREFIX}/Q{1}L{2}X{3}P{4}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta does not exist or is empty'
        exit;
    fi

    mkdir -p ${DIR_PREFIX}/Q{1}L{2}X{3}P{4}
    cd ${DIR_PREFIX}/Q{1}L{2}X{3}P{4}

    anchr asm anchor \
        unitigs.fasta \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
        --mincov 5 --mscale 3 \
        --lscale {{ opt.lscale }} \
        --uscale {{ opt.uscale }} \
        -p {{ parallel2 }} \
        --stats anchor.stats.tsv \
        -o anchor.fasta

    echo >&2
    " ::: 0 {{ opt.qual }} ::: 0 {{ opt.len }} ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
