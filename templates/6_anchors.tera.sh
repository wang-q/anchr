{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 6_anchors.sh

{% set parallel2 = opt.parallel | int / 2 -%}
{% set parallel2 = parallel2 | round(method="floor") -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
parallel --no-run-if-empty --linebuffer -k -j 2 "
    if [ ! -e 6_down_sampling/MRX{1}P{2}/pe.cor.fa ]; then
        exit;
    fi

    echo >&2 '==> 6_unitigs/MRX{1}P{2}'
    if [ -e 6_unitigs/MRX{1}P{2}/anchor/anchor.fasta ]; then
        echo >&2 '    anchor.fasta already presents'
        exit;
    fi

    if [ ! -s 6_unitigs/MRX{1}P{2}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta does not exist or is empty'
        exit;
    fi

    if [ -d 6_unitigs/MRX{1}P{2}/anchor ]; then
        rm -fr 6_unitigs/MRX{1}P{2}/anchor
    fi
    mkdir -p 6_unitigs/MRX{1}P{2}/anchor
    cd 6_unitigs/MRX{1}P{2}/anchor

    anchr anchors \
        ../unitigs.fasta \
        ../pe.cor.fa \
        -p {{ parallel2 }} \
        -o anchors.sh
    bash anchors.sh

    echo >&2
    " ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
