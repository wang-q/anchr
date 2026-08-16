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

# Plain bash loops replace the former `parallel ... :::` argument grid
# (nested double-quoted templates need \$ / \" escaping that breaks
# silently). Jobs run one at a time in the foreground, each pigz with the
# full opt.parallel thread budget.
for X in {{ opt.cov }}; do
    (
    echo "==> MRX${X}"

    if [ -d 6_down_sampling/MRX${X} ]; then
        echo '    Skip'
        exit
    fi

    # shortcut if cov2 == all
    if [[ ${X} == 'all' ]]; then
        mkdir -p 6_down_sampling/MRXallP000
        cd 6_down_sampling/MRXallP000
        gzip -dcf ../../2_illumina/merge/pe.cor.fa.gz |
            pigz -p {{ opt.parallel }} > pe.cor.fa.gz
        cp ../../2_illumina/merge/env.json .
        exit;
    fi

    # actual sampling
    mkdir -p 6_down_sampling/MRX${X}
    pgr fa split about -e -c $(( {{ opt.genome }} * ${X} )) \
        2_illumina/merge/pe.cor.fa.gz \
        -o 6_down_sampling/MRX${X}

    MAX_SERIAL=$(
        cat 2_illumina/merge/env.json |
            jq ".SUM_OUT | tonumber | . / {{ opt.genome }} / ${X} | floor | . - 1"
    )
    MAX_SERIAL=$(( ${MAX_SERIAL} < {{ opt.splitp }} ? ${MAX_SERIAL} : {{ opt.splitp }} ))

    for i in $( seq 0 1 ${MAX_SERIAL} ); do
        P=$( printf '%03d' ${i})
        printf "  * Part: %s\n" ${P}

        mkdir -p "6_down_sampling/MRX${X}P${P}"

        mv  "6_down_sampling/MRX${X}/${P}.fa" \
            "6_down_sampling/MRX${X}P${P}/pe.cor.fa"
        pigz -p {{ opt.parallel }} "6_down_sampling/MRX${X}P${P}/pe.cor.fa"
        cp 2_illumina/merge/env.json "6_down_sampling/MRX${X}P${P}"
    done
    )
done
