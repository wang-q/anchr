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

if [ -e ../1_genome/repetitive/repetitive.fa ]; then
    printf "%s\t%s\t%s\t%s\n" \
        $(echo "repetitive"; stat_format ../1_genome/repetitive/repetitive.fa;)
fi \
    >> statReads.tsv

for PREFIX in R S T; do
    if [ -e ${PREFIX}1.fq.gz ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(echo "Illumina.${PREFIX}"; stat_format_fq ${PREFIX}1.fq.gz {% if opt.se == "0" %}${PREFIX}2.fq.gz{% endif %};)
    fi
    if [ -e trim/${PREFIX}1.fq ] || [ -e trim/${PREFIX}1.fq.gz ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(
                echo "trim.${PREFIX}";
{% if opt.se == "0" %}
                FQ1=trim/${PREFIX}1.fq; [ -e ${FQ1} ] || FQ1=trim/${PREFIX}1.fq.gz
                FQ2=trim/${PREFIX}2.fq; [ -e ${FQ2} ] || FQ2=trim/${PREFIX}2.fq.gz
                FQS=trim/${PREFIX}s.fq; [ -e ${FQS} ] || FQS=trim/${PREFIX}s.fq.gz
                stat_format_fq ${FQ1} ${FQ2} ${FQS};
{% else %}
                FQ1=trim/${PREFIX}1.fq; [ -e ${FQ1} ] || FQ1=trim/${PREFIX}1.fq.gz
                stat_format_fq ${FQ1};
{% endif %}
            )
    fi
done \
    >> statReads.tsv

for PREFIX in R S T; do
    for Q in 0 {{ opt.qual }}; do
        for L in 0 {{ opt.len }}; do
            FQ1=Q${Q}L${L}/${PREFIX}1.fq
            [ -e ${FQ1} ] || FQ1=Q${Q}L${L}/${PREFIX}1.fq.gz
            if [ ! -e ${FQ1} ]; then
                continue
            fi

            printf "%s\t%s\t%s\t%s\n" \
                $(
                    echo Q${Q}L${L};
{% if opt.se == "0" %}
                    FQ2=Q${Q}L${L}/${PREFIX}2.fq; [ -e ${FQ2} ] || FQ2=Q${Q}L${L}/${PREFIX}2.fq.gz
                    FQS=Q${Q}L${L}/${PREFIX}s.fq; [ -e ${FQS} ] || FQS=Q${Q}L${L}/${PREFIX}s.fq.gz
                    stat_format_fq \
                        ${FQ1} \
                        ${FQ2} \
                        ${FQS};
{% else %}
                    stat_format_fq ${FQ1};
{% endif %}
                )
        done
    done
done \
    >> statReads.tsv

fi # end of statReads

tva to md statReads.tsv --right 2-4 -o statReads.md
echo -e "\nTable: statReads\n" >> statReads.md

cat statReads.md
mkdir -p ${BASH_DIR}/../9_markdown
mv statReads.md ${BASH_DIR}/../9_markdown
