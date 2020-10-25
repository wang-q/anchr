{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 6_down_sampling.sh

if [ ! -e 2_illumina/merge/pe.cor.fa.gz ]; then
    echo >&2 "2_illumina/merge/pe.cor.fa.gz not presents"
    exit;
fi

parallel --no-run-if-empty --linebuffer -k -j 2 "
    echo '==> MRX{}'

    if [ -d 6_downSampling/MRX{} ]; then
        echo '    Skip'
        exit
    fi

    # shortcut if cov2 == all
    if [[ {} == 'all' ]]; then
        mkdir -p 6_downSampling/MRXallP000
        cd 6_downSampling/MRXallP000
        gzip -d -c ../../2_illumina/merge/pe.cor.fa.gz > pe.cor.fa
        cp ../../2_illumina/merge/env.json .
        exit;
    fi

    # actual sampling
    mkdir -p 6_downSampling/MRX{}
    faops split-about -l 0 -e \
        2_illumina/merge/pe.cor.fa.gz \
        \$(( {{ opt.genome }} * {} )) \
        6_downSampling/MRX{}

    MAX_SERIAL=\$(
        cat 2_illumina/merge/env.json \
            | jq '.SUM_OUT | tonumber | . / {{ opt.genome }} / {} | floor | . - 1'
    )
    MAX_SERIAL=\$(( \${MAX_SERIAL} < {{ opt.splitp }} ? \${MAX_SERIAL} : {{ opt.splitp }} ))

    for i in \$( seq 0 1 \${MAX_SERIAL} ); do
        P=\$( printf '%03d' \${i})
        printf \"  * Part: %s\n\" \${P}

        mkdir -p \"6_downSampling/MRX{}P\${P}\"

        mv  \"6_downSampling/MRX{}/\${P}.fa\" \
            \"6_downSampling/MRX{}P\${P}/pe.cor.fa\"
        cp 2_illumina/merge/env.json \"6_downSampling/MRX{}P\${P}\"
    done

    " ::: {{ opt.cov }}
