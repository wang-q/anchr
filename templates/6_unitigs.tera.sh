{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn {{ outname }}

parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ ! -e 6_down_sampling/MRX{1}P{2}/pe.cor.fa ]; then
        exit;
    fi

    echo >&2 '==> 6_unitigs_multik/MRX{1}P{2}'
    if [ -e 6_unitigs_multik/MRX{1}P{2}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta already presents'
        exit;
    fi

    mkdir -p 6_unitigs_multik/MRX{1}P{2}
    cd 6_unitigs_multik/MRX{1}P{2}

    anchr asm multik \
        ../../6_down_sampling/MRX{1}P{2}/pe.cor.fa \
        -k 31,41,51,61,71,81 \
        -p {{ opt.parallel }} \
        -o unitigs.fasta

    echo >&2
    " ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
