{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn {{ outname }}

parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ ! -e 6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz ]; then
        exit;
    fi

    echo >&2 '==> 6_unitigs_{{ unitigger }}/MRX{1}P{2}'
    if [ -e 6_unitigs_{{ unitigger }}/MRX{1}P{2}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta already presents'
        exit;
    fi

    mkdir -p 6_unitigs_{{ unitigger }}/MRX{1}P{2}
    cd 6_unitigs_{{ unitigger }}/MRX{1}P{2}

{% if unitigger == "bcalm" %}    # external bcalm unitigs per k, merged across k with the modern OLC step
    for K in 31 41 51 61 71 81; do
        bcalm \
            -in ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
            -kmer-size ${K} -abundance-min 3 -verbose 0 \
            -nb-cores {{ opt.parallel }} \
            -out K${K}
        mv K${K}.unitigs.fa unitigs_K${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-contig-len 1000 \
        -o unitigs.fasta
{% elif unitigger == "unitig" %}    # in-house BCALM-semantics unitigs per k (asm unitig), merged across k
    for K in 31 41 51 61 71 81; do
        anchr asm unitig \
            ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
            -k ${K} \
            -p {{ opt.parallel }} \
            -o unitigs_K${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-contig-len 1000 \
        -o unitigs.fasta
{% else %}    anchr asm multik \
        ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
        -k 31,41,51,61,71,81 \
        -p {{ opt.parallel }} \
        -o unitigs.fasta
{% endif %}

    echo >&2
    " ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
