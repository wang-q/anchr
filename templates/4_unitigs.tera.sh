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

{% set parallel2 = opt.parallel | int / 2 -%}
{% set parallel2 = parallel2 | round(method="floor") -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
{# Per-master k range: 150 bp reads cap at 91 (k near read length fragments) #}
KS="31 41 51 61 71 81 91"
{% if unitigger == "bcalm" %}    # external bcalm unitigs per k, merged across k with the modern OLC step
    for K in ${KS}; do
        bcalm \
            -in ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
            -kmer-size \${K} -abundance-min 3 -verbose 0 \
            -nb-cores {{ opt.parallel }} \
            -out K\${K}
        mv K\${K}.unitigs.fa unitigs_K\${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% elif unitigger == "unitig" %}    # in-house BCALM-semantics unitigs per k (asm unitig), merged across k
    for K in ${KS}; do
        anchr asm unitig \
            ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
            -k \${K} \
            -p {{ opt.parallel }} \
            -o unitigs_K\${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% else %}    # per-master multik: every k builds its own skeleton (larger ks
    # validate it), masters run in parallel, then merged across masters
    for K in ${KS}; do
        K_LIST=\"\"
        for J in ${KS}; do
            if [ \${J} -ge \${K} ]; then
                K_LIST=\"\${K_LIST}\${K_LIST:+,}\${J}\"
            fi
        done
        (
            anchr asm multik \
                ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
                -k \${K_LIST} \
                -p {{ parallel2 }} \
                -o unitigs_K\${K}.fasta
        ) &
    done
    wait

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa.gz \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% endif %}

    echo >&2
    " ::: 0 {{ opt.qual }} ::: 0 {{ opt.len }} ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
