{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn {{ outname }}

# Plain bash loops over the cov x part grid replace the former
# `parallel ... ::: cov ::: parts` argument grid: nested double-quoted
# templates need \$ / \" escaping that breaks silently (bash -n cannot see
# inside the template string), while for-loops stay checkable. Jobs run
# one at a time, each with the full opt.parallel thread budget.
for X in {{ opt.cov }}; do
    for P in $(printf "%03d " {0..{{ opt.splitp }}}); do
        (
    if [ ! -e 6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz ]; then
        exit;
    fi

    echo >&2 "==> 6_unitigs_{{ unitigger }}/MRX${X}P${P}"
    if [ -e 6_unitigs_{{ unitigger }}/MRX${X}P${P}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta already presents'
        exit;
    fi

    mkdir -p 6_unitigs_{{ unitigger }}/MRX${X}P${P}
    cd 6_unitigs_{{ unitigger }}/MRX${X}P${P}

{# The multik branch derives its master-k list from the read-length N50
   (anchr asm multik --print-ks) instead of hard-coding values, so the
   pipeline adapts to any read length. The static KS below is only the
   fallback for the unitig branch (per-k independent unitigs) and the
   bcalm cap (bcalm rejects k>127). #}
KS="31 41 51 61 71 81 101 121 128 160 192"
KS_BCALM="31 41 51 61 71 81 101 121"
{% if unitigger == "bcalm" %}    # external bcalm unitigs per k, merged across k with the modern OLC step
    for K in ${KS_BCALM}; do
        bcalm \
            -in ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz \
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
        ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% elif unitigger == "unitig" %}    # in-house BCALM-semantics unitigs per k (asm unitig), merged across k
    for K in ${KS}; do
        anchr asm unitig \
            ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz \
            -k ${K} \
            -p {{ opt.parallel }} \
            -o unitigs_K${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 1000 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% else %}    # per-master multik: every k builds its own skeleton (larger ks
    # validate it). The master-k list comes from the read-length N50 (no
    # hard-coded k values); the first master runs first and its unitigs
    # guide the rest (megahit seq2sdbg --contig guidance), so high-k
    # masters do not fragment on reads with thin high-k support. Every
    # master validates against every third later k (fixed spacing rule:
    # intermediate rounds re-count with little added signal).
    # Keep short (200..1000 bp) reads-supported fragments: on merged reads
    # they cover low-complexity gap regions without chimeras (G37: GF
    # 98.869->99.083%, 0 mis; bcalm/unitig branches keep 1000 because their
    # raw unitigs do produce chimeric short fragments).
    KS=$(anchr asm multik ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz --print-ks)
    K0=$(echo ${KS} | awk '{print $1}')
    # The first master is the serial prefix of this stage and only needs to
    # produce the guide skeleton: at most two validation rounds (median +
    # largest k, the strictest check) instead of the full every-third-k
    # chain. Later masters keep the dense rule — their output feeds the
    # final OLC merge.
    V0=$(echo ${KS} | awk 'NF>=5{printf "%s,%s", $(int(NF/2)+1), $NF} NF>=2&&NF<5{printf "%s", $NF}')
    if [ -n "${V0}" ]; then K0_LIST="${K0},${V0}"; else K0_LIST="${K0}"; fi
    anchr asm multik \
        ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz \
        -k ${K0_LIST} \
        -p {{ opt.parallel }} \
        -o unitigs_K${K0}.fasta

    # Masters run one at a time, each with the full opt.parallel budget.
    for K in $(echo ${KS} | awk '{for(i=2;i<=NF;i++) print $i}'); do
        VK=$(echo ${KS} | awk -v k=${K} '{for(i=1;i<=NF;i++) if($i>k && (i-1)%3==0) printf "%s,", $i}' | sed 's/,$//')
        if [ -n "${VK}" ]; then K_LIST="${K},${VK}"; else K_LIST="${K}"; fi
        anchr asm multik \
            ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz \
            --guide-contigs unitigs_K${K0}.fasta \
            -k ${K_LIST} \
            -p {{ opt.parallel }} \
            -o unitigs_K${K}.fasta
    done

    anchr asm olc --unitigs unitigs_K*.fasta \
        --min-overlap 1000 \
        --min-contig-len 200 \
        -o unitigs.fasta

    anchr asm extend unitigs.fasta \
        ../../6_down_sampling/MRX${X}P${P}/pe.cor.fa.gz \
        --min-len 1000 \
        -o unitigs.ext.fasta
    mv unitigs.ext.fasta unitigs.fasta
{% endif %}

    echo >&2
        )
    done
done
