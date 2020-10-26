{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 4_down_sampling.sh

parallel --no-run-if-empty --linebuffer -k -j 2 "
    if [ ! -e 2_illumina/Q{1}L{2}/pe.cor.fa.gz ]; then
        exit;
    fi
    echo '==> Q{1}L{2}X{3}'

    if [ -d 4_down_sampling/Q{1}L{2}X{3} ]; then
        echo '    Skip'
        exit
    fi

    # shortcut if cov == all
    if [[ {3} == 'all' ]]; then
        mkdir -p 4_down_sampling/Q{1}L{2}XallP000
        cd 4_down_sampling/Q{1}L{2}XallP000
        gzip -dcf ../../2_illumina/Q{1}L{2}/pe.cor.fa.gz > pe.cor.fa
        cp ../../2_illumina/Q{1}L{2}/env.json .
        exit;
    fi

    # actual sampling
    mkdir -p 4_down_sampling/Q{1}L{2}X{3}
    faops split-about -l 0 -e \
        2_illumina/Q{1}L{2}/pe.cor.fa.gz \
        \$(( {{ opt.genome }} * {3} )) \
        4_down_sampling/Q{1}L{2}X{3}

    MAX_SERIAL=\$(
        cat 2_illumina/Q{1}L{2}/env.json |
            jq '.SUM_OUT | tonumber | . / {{ opt.genome }} / {3} | floor | . - 1'
    )
    MAX_SERIAL=\$(( \${MAX_SERIAL} < {{ opt.splitp }} ? \${MAX_SERIAL} : {{ opt.splitp }} ))

    for i in \$( seq 0 1 \${MAX_SERIAL} ); do
        P=\$( printf '%03d' \${i})
        printf \"  * Part: %s\n\" \${P}

        mkdir -p \"4_down_sampling/Q{1}L{2}X{3}P\${P}\"

        mv  \"4_down_sampling/Q{1}L{2}X{3}/\${P}.fa\" \
            \"4_down_sampling/Q{1}L{2}X{3}P\${P}/pe.cor.fa\"
        cp 2_illumina/Q{1}L{2}/env.json \"4_down_sampling/Q{1}L{2}X{3}P\${P}\"
    done

    " ::: 0 {{ opt.qual }} ::: 0 {{ opt.len }} ::: {{ opt.cov }}
