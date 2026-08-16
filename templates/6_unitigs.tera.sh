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

{% set parallel2 = opt.parallel | int / 2 -%}
{% set parallel2 = parallel2 | round(method="floor") -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
{# Per-master k range: merged reads (~450 bp) support masters up to 160.
   K128 crosses low-complexity gaps that 121 leaves broken; K160 adds a
   long-unitig layer that lifts the final N50 (G37 MR: 83.8K->121.5K,
   Dup 1.000, 0 mis after the indel-tolerant near-duplicate merge;
   see notes/benchmarks/g37-megahit-spades.md P2). k>=192 fragments on
   450 bp reads, so the list stops at 160. bcalm rejects k>127, so its
   branch keeps the old cap. #}
KS="31 41 51 61 71 81 101 121 128 160"
KS_BCALM="31 41 51 61 71 81 101 121"
{% if unitigger == "bcalm" %}    # external bcalm unitigs per k, merged across k with the modern OLC step
    for K in ${KS_BCALM}; do
        bcalm \
            -in ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
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
        ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% elif unitigger == "unitig" %}    # in-house BCALM-semantics unitigs per k (asm unitig), merged across k
    for K in ${KS}; do
        anchr asm unitig \
            ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
            -k \${K} \
            -p {{ opt.parallel }} \
            -o unitigs_K\${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% else %}    # per-master multik: every k builds its own skeleton (larger ks
    # validate it), masters run in parallel, then merged across masters
    # Keep short (200..1000 bp) reads-supported fragments: on merged reads
    # they cover low-complexity gap regions without chimeras (G37: GF
    # 98.869->99.083%, 0 mis; bcalm/unitig branches keep 1000 because their
    # raw unitigs do produce chimeric short fragments).
    for K in ${KS}; do
        K_LIST=\"\"
        for J in ${KS}; do
            if [ \${J} -ge \${K} ]; then
                K_LIST=\"\${K_LIST}\${K_LIST:+,}\${J}\"
            fi
        done
        (
            anchr asm multik \
                ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
                -k \${K_LIST} \
                -p {{ parallel2 }} \
                -o unitigs_K\${K}.fasta
        ) &
    done
    wait

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 200 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% endif %}

    echo >&2
    " ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
