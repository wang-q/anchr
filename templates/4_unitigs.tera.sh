{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn {{ outname }}

parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ ! -e 4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa ]; then
        exit;
    fi

    echo >&2 '==> 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}'
    if [ -e 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta already presents'
        exit;
    fi

    mkdir -p 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}
    cd 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}

    anchr unitigs \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/env.json \
        -u {{ unitigger }} \
        -p {{ opt.parallel }} \
        --kmer '31 41 51 61 71 81' \
        -o unitigs.sh
    bash unitigs.sh

    echo >&2
    " ::: 0 {{ opt.qual }} ::: 0 {{ opt.len }} ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
