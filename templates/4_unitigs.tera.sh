{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn {{ outname }}

parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ ! -e 4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz ]; then
        exit;
    fi

    echo >&2 '==> 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}'
    if [ -e 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta already presents'
        exit;
    fi

    mkdir -p 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}
    cd 4_unitigs_{{ unitigger }}/Q{1}L{2}X{3}P{4}

{% if unitigger == "bcalm" %}    # external bcalm unitigs per k, merged across k with the modern OLC step
    for K in 31 41 51 61 71 81; do
        bcalm \
            -in ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
            -kmer-size ${K} -abundance-min 3 -verbose 0 \
            -nb-cores {{ opt.parallel }} \
            -out K${K}
        mv K${K}.unitigs.fa unitigs_K${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-contig-len 1000 \
        -o unitigs.fasta
{% else %}    anchr asm multik \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
        -k 31,41,51,61,71,81 \
        -p {{ opt.parallel }} \
        -o unitigs.fasta
{% endif %}

    echo >&2
    " ::: 0 {{ opt.qual }} ::: 0 {{ opt.len }} ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
