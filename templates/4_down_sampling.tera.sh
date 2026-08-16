{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 4_down_sampling.sh

# Plain bash loops replace the former `parallel -j 2 ... :::` argument grid
# (nested double-quoted templates need \$ / \" escaping that breaks
# silently). Jobs run one at a time in the foreground, each pigz with the
# full opt.parallel thread budget.
for Q in 0 {{ opt.qual }}; do
    for L in 0 {{ opt.len }}; do
        for X in {{ opt.cov }}; do
            (
    if [ ! -e 2_illumina/Q${Q}L${L}/pe.cor.fa.gz ]; then
        exit;
    fi
    echo "==> Q${Q}L${L}X${X}"

    if [ -d 4_down_sampling/Q${Q}L${L}X${X} ]; then
        echo '    Skip'
        exit
    fi

    # shortcut if cov == all
    if [[ ${X} == 'all' ]]; then
        mkdir -p 4_down_sampling/Q${Q}L${L}XallP000
        cd 4_down_sampling/Q${Q}L${L}XallP000
        gzip -dcf ../../2_illumina/Q${Q}L${L}/pe.cor.fa.gz |
            pigz -p {{ opt.parallel }} > pe.cor.fa.gz
        cp ../../2_illumina/Q${Q}L${L}/env.json .
        exit;
    fi

    # actual sampling
    mkdir -p 4_down_sampling/Q${Q}L${L}X${X}
    pgr fa split about -e -c $(( {{ opt.genome }} * ${X} )) \
        2_illumina/Q${Q}L${L}/pe.cor.fa.gz \
        -o 4_down_sampling/Q${Q}L${L}X${X}

    MAX_SERIAL=$(
        cat 2_illumina/Q${Q}L${L}/env.json |
            jq ".SUM_OUT | tonumber | . / {{ opt.genome }} / ${X} | floor | . - 1"
    )
    MAX_SERIAL=$(( ${MAX_SERIAL} < {{ opt.splitp }} ? ${MAX_SERIAL} : {{ opt.splitp }} ))

    for i in $( seq 0 1 ${MAX_SERIAL} ); do
        P=$( printf '%03d' ${i})
        printf "  * Part: %s\n" ${P}

        mkdir -p "4_down_sampling/Q${Q}L${L}X${X}P${P}"

        mv  "4_down_sampling/Q${Q}L${L}X${X}/${P}.fa" \
            "4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa"
        pigz -p {{ opt.parallel }} "4_down_sampling/Q${Q}L${L}X${X}P${P}/pe.cor.fa"
        cp 2_illumina/Q${Q}L${L}/env.json "4_down_sampling/Q${Q}L${L}X${X}P${P}"
    done
        )
    done
    done
done
