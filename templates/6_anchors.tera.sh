{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 6_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX]"

DIR_PREFIX=${1:-"6_unitigs_multik"}

{% set parallel2 = opt.parallel | int / 2 -%}
{% set parallel2 = parallel2 | round(method="floor") -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
parallel --no-run-if-empty --linebuffer -k -j 2 "
    if [ ! -e 6_down_sampling/MRX{1}P{2}/pe.cor.fa ]; then
        exit;
    fi

    echo >&2 '==> ${DIR_PREFIX}/MRX{1}P{2}'
    if [ -e ${DIR_PREFIX}/MRX{1}P{2}/anchor.fasta ]; then
        echo >&2 '    anchor.fasta already presents'
        exit;
    fi

    if [ ! -s ${DIR_PREFIX}/MRX{1}P{2}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta does not exist or is empty'
        exit;
    fi

    mkdir -p ${DIR_PREFIX}/MRX{1}P{2}
    cd ${DIR_PREFIX}/MRX{1}P{2}

    anchr asm anchor \
        unitigs.fasta \
        ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa \
        --mincov 5 --mscale 3 \
        --lscale {{ opt.lscale }} \
        --uscale {{ opt.uscale }} \
        -p {{ parallel2 }} \
        --stats anchor.stats.tsv \
        -o anchor.fasta

    echo >&2
    " ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
