{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_stat_reads.sh

if [ -e statReads.md ]; then
    log_debug "statReads.md presents";
    exit;
fi

echo -e "Table: statReads\n" > statReads.md
printf "| %s | %s | %s | %s |\n" \
    "Name" "N50" "Sum" "#" \
    >> statReads.md
printf "|:--|--:|--:|--:|\n" >> statReads.md

if [ -e 1_genome/genome.fa ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "Genome";   faops n50 -H -S -C 1_genome/genome.fa;) >> statReads.md
fi
if [ -e 1_genome/paralogs.fa ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "Paralogs"; faops n50 -H -S -C 1_genome/paralogs.fa;) >> statReads.md
fi

for PREFIX in R S T; do
    if [ -e 2_illumina/${PREFIX}1.fq.gz ]; then
        printf "| %s | %s | %s | %s |\n" \
            $(echo "Illumina.${PREFIX}"; stat_format 2_illumina/${PREFIX}1.fq.gz {% if opt.se == "0" %}2_illumina/${PREFIX}2.fq.gz{% endif %};) >> statReads.md
    fi
    if [ -e 2_illumina/trim/${PREFIX}1.fq.gz ]; then
        printf "| %s | %s | %s | %s |\n" \
            $(echo "trim.${PREFIX}"; stat_format 2_illumina/trim/${PREFIX}1.fq.gz {% if opt.se == "0" %}2_illumina/trim/${PREFIX}2.fq.gz 2_illumina/trim/${PREFIX}s.fq.gz{% endif %};) >> statReads.md
    fi

    parallel --no-run-if-empty -k -j 2 "
        stat_format () {
            echo \$(faops n50 -H -N 50 -S -C \$@) \
                | perl -nla -MNumber::Format -e '
                    printf qq(%d\t%s\t%d\n), \$F[0], Number::Format::format_bytes(\$F[1], base => 1000,), \$F[2];
                '
        }

        if [ ! -e 2_illumina/Q{1}L{2}/${PREFIX}1.fq.gz ]; then
            exit;
        fi

        printf \"| %s | %s | %s | %s |\n\" \
            \$(
                echo Q{1}L{2};
{% if opt.se == "0" %}
                stat_format \
                    2_illumina/Q{1}L{2}/${PREFIX}1.fq.gz \
                    2_illumina/Q{1}L{2}/${PREFIX}2.fq.gz \
                    2_illumina/Q{1}L{2}/${PREFIX}s.fq.gz;
{% else %}
                stat_format \
                    2_illumina/Q{1}L{2}/${PREFIX}1.fq.gz;
{% endif %}
            )
        " ::: {{ opt.qual }} ::: {{ opt.len }} \
        >> statReads.md
done

cat statReads.md
