{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 2_sga_preqc.sh

mkdir -p 2_illumina/sga_preqc
cd 2_illumina/sga_preqc

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -e ${PREFIX}.preqc.pdf ]; then
        continue;
    fi

    sga preprocess \
{% if opt.se == "0" -%}
        ../${PREFIX}1.fq.gz \
{% else -%}
        ../${PREFIX}1.fq.gz ../${PREFIX}2.fq.gz \
        --pe-mode 1 \
{% endif -%}
        -o ${PREFIX}.pp.fq

    sga index -a ropebwt -t {{ opt.parallel }} ${PREFIX}.pp.fq

    sga preqc -t {{ opt.parallel }} ${PREFIX}.pp.fq > ${PREFIX}.preqc_output

    sga-preqc-report.py ${PREFIX}.preqc_output -o ${PREFIX}.preqc

{% if opt.sgastats == "1" -%}
    sga stats -t {{ opt.parallel }} -n {{ opt.reads }} ${PREFIX}.pp.fq > ${PREFIX}.stats.txt
{% endif -%}

    find . -type f -name "${PREFIX}.pp.*" |
        parallel --no-run-if-empty -j 1 rm

done

{% if opt.sgastats == "1" -%}
echo -e "Table: statSgaStats\n" > statSgaStats.md
printf "| %s | %s | %s | %s |\n" \
    "Library" "incorrectBases" "perfectReads" "overlapDepth" \
    >> statSgaStats.md
printf "|:--|--:|--:|--:|\n" >> statSgaStats.md

# sga stats
#*** Stats:
#380308 out of 149120670 bases are potentially incorrect (0.002550)
#797208 reads out of 1000000 are perfect (0.797208)
#Mean overlap depth: 356.41
for PREFIX in R S T; do
    if [ ! -e ${PREFIX}.stats.txt ]; then
        continue;
    fi

    printf "| %s " "${PREFIX}" >> statSgaStats.md
    cat ${PREFIX}.stats.txt |
        perl -nl -e '
            BEGIN { our $stat = { }; };

            m{potentially incorrect \(([\d\.]+)\)} and $stat->{incorrectBases} = $1;
            m{perfect \(([\d\.]+)\)} and $stat->{perfectReads} = $1;
            m{overlap depth: ([\d\.]+)} and $stat->{overlapDepth} = $1;

            END {
                printf qq{| %.2f%% | %.2f%% | %s |\n},
                    $stat->{incorrectBases} * 100,
                    $stat->{perfectReads} * 100,
                    $stat->{overlapDepth};
            }
            ' \
        >> statSgaStats.md
done
{% endif -%}

cat statSgaStats.md

mv statSgaStats.md ../../
