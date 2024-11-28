{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_stat_reads.sh

cd 2_illumina

if [ -e statReads.tsv ]; then
    log_debug "statReads.tsv presents"
else

printf "%s\t%s\t%s\t%s\n" \
    "Name" "N50" "Sum" "#" \
    > statReads.tsv

for NAME in genome paralogs; do
    if [ -e ../1_genome/${NAME}.fa ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(echo "${NAME}"; stat_format ../1_genome/${NAME}.fa;)
    fi
done \
    >> statReads.tsv

if [ -e ../1_genome/repetitive.fa ]; then
    printf "%s\t%s\t%s\t%s\n" \
        $(echo "${NAME}"; stat_format ../1_genome/repetitive/repetitive.fa;)
fi \
    >> statReads.tsv

for PREFIX in R S T; do
    if [ -e ${PREFIX}1.fq.gz ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(echo "Illumina.${PREFIX}"; stat_format ${PREFIX}1.fq.gz {% if opt.se == "0" %}${PREFIX}2.fq.gz{% endif %};)
    fi
    if [ -e trim/${PREFIX}1.fq.gz ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(echo "trim.${PREFIX}"; stat_format trim/${PREFIX}1.fq.gz {% if opt.se == "0" %}trim/${PREFIX}2.fq.gz trim/${PREFIX}s.fq.gz{% endif %};)
    fi
done \
    >> statReads.tsv

for PREFIX in R S T; do
    for Q in 0 {{ opt.qual }}; do
        for L in 0 {{ opt.len }}; do
            if [ ! -e Q${Q}L${L}/${PREFIX}1.fq.gz ]; then
                continue
            fi

            printf "%s\t%s\t%s\t%s\n" \
                $(
                    echo Q${Q}L${L};
{% if opt.se == "0" %}
                    stat_format \
                        Q${Q}L${L}/${PREFIX}1.fq.gz \
                        Q${Q}L${L}/${PREFIX}2.fq.gz \
                        Q${Q}L${L}/${PREFIX}s.fq.gz;
{% else %}
                    stat_format \
                        Q${Q}L${L}/${PREFIX}1.fq.gz;
{% endif %}
                )
        done
    done
done \
    >> statReads.tsv

fi # end of statReads

cat statReads.tsv |
    rgr md stdin --right 2-4 \
    > statReads.md

echo -e "\nTable: statReads\n" >> statReads.md

cat statReads.md
mv statReads.md ${BASH_DIR}/../9_markdown
