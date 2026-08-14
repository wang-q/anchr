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

    echo >&2 '==> 4_unitigs_multik/Q{1}L{2}X{3}P{4}'
    if [ -e 4_unitigs_multik/Q{1}L{2}X{3}P{4}/unitigs.fasta ]; then
        echo >&2 '    unitigs.fasta already presents'
        exit;
    fi

    mkdir -p 4_unitigs_multik/Q{1}L{2}X{3}P{4}
    cd 4_unitigs_multik/Q{1}L{2}X{3}P{4}

    anchr asm multik \
        ../../4_down_sampling/Q{1}L{2}X{3}P{4}/pe.cor.fa \
        -k 31,41,51,61,71,81 \
        -p {{ opt.parallel }} \
        -o unitigs.fasta

    echo >&2
    " ::: 0 {{ opt.qual }} ::: 0 {{ opt.len }} ::: {{ opt.cov }} ::: $(printf "%03d " {0..{{ opt.splitp }}})
