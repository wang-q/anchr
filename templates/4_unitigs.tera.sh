{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn {{ outname }}

# Plain bash loops over the qual x len x cov x part grid replace the former
# `parallel ... :::` argument grid: nested double-quoted templates need
# \$ / \" escaping that breaks silently (bash -n cannot see inside the
# template string), while for-loops stay checkable. Jobs run one at a
# time, each with the full opt.parallel thread budget.
for Q in 0 {{ opt.qual }}; do
    for L in 0 {{ opt.len }}; do
        for X in {{ opt.cov }}; do
            for P in $(printf "%03d " {0..{{ opt.splitp }}}); do
                (
    if [ ! -e 4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz ]; then
        exit;
    fi

    echo >&2 "==> 4_unitigs_{{ unitigger }}/Q${Q}L${L}X${X}P${P}"
    if [ -e 4_unitigs_{{ unitigger }}/Q${Q}L${L}X${X}P${P}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta already presents'
        exit;
    fi

    mkdir -p 4_unitigs_{{ unitigger }}/Q${Q}L${L}X${X}P${P}
    cd 4_unitigs_{{ unitigger }}/Q${Q}L${L}X${X}P${P}

{# Per-master k range: 150 bp reads cap at 91 (k near read length fragments) #}
KS="31 41 51 61 71 81 91"
{% if unitigger == "bcalm" %}    # external bcalm unitigs per k, merged across k with the modern OLC step
    for K in ${KS}; do
        bcalm \
            -in ../../4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz \
            -kmer-size ${K} -abundance-min 3 -verbose 0 \
            -nb-cores {{ opt.parallel }} \
            -out K${K}
        mv K${K}.unitigs.fa unitigs_K${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% elif unitigger == "unitig" %}    # in-house BCALM-semantics unitigs per k (asm unitig), merged across k
    for K in ${KS}; do
        anchr asm unitig \
            ../../4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz \
            -k ${K} \
            -p {{ opt.parallel }} \
            -o unitigs_K${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% else %}    # single-invocation multi-master multik: every k in KS builds its own
    # skeleton validated by the larger ks (k-major order, the reads count
    # at each k is built once and shared by every master). Replaces the
    # per-master loop that re-counted the reads once per (master, round).
    # No guide: 5-group anchor voting makes it quality-neutral (MG1655
    # 0 mis with and without) while running ~2x faster.
    anchr asm multik \
        ../../4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz \
        -k $(echo ${KS} | tr ' ' ',') \
        --all-masters \
        -p {{ opt.parallel }} \
        -o unitigs_all.fasta

    anchr asm olc --unitigs unitigs_all.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% endif %}

    echo >&2
                )
            done
        done
    done
done
