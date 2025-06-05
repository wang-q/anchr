{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# run
#----------------------------#
log_warn 2_trim.sh

mkdir -p 2_illumina/trim
cd 2_illumina/trim

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -e ${PREFIX}1.fq.gz ]; then
        log_debug "2_illumina/trim/${PREFIX}1.fq.gz presents"
        continue;
    fi

    anchr trim \
        {{ opt.trim }} \
        --qual "{{ opt.qual }}" \
        --len "{{ opt.len }}" \
    {% if opt.filter -%}
        --filter "{{ opt.filter }}" \
    {% endif -%}
    {% if opt.sample != "0" -%}
    {% if opt.genome != "0" -%}
        --sample $(( {{ opt.genome }} * {{ opt.sample }} )) \
    {% endif -%}
    {% endif -%}
        --parallel {{ opt.parallel }}{% if opt.xmx != "0" %} --xmx {{ opt.xmx }}{% endif %} \
        ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
        --prefix ${PREFIX} \
        -o trim.sh
    bash trim.sh

    log_info "stats of all .fq.gz files"

    if [ ! -e statTrimReads.tsv ]; then
        printf "%s\t%s\t%s\t%s\n" \
            "Name" "N50" "Sum" "#" \
            > statTrimReads.tsv
    fi

    for NAME in clumpify filteredbytile highpass sample trim filter ${PREFIX}1 ${PREFIX}2 ${PREFIX}s; do
        if [ ! -e ${NAME}.fq.gz ]; then
            continue;
        fi

        printf "%s\t%s\t%s\t%s\n" \
            $(echo ${NAME}; stat_format_fq ${NAME}.fq.gz;) >> statTrimReads.tsv
    done

    log_info "clear unneeded .fq.gz files"
    for NAME in temp clumpify filteredbytile highpass sample trim filter; do
        if [ -e ${NAME}.fq.gz ]; then
            rm ${NAME}.fq.gz
        fi
    done
done

cat statTrimReads.tsv |
    rgr md stdin --right 2-4 \
    > statTrimReads.md

echo -e "\nTable: statTrimReads\n" >> statTrimReads.md

for PREFIX in R S T; do
    if [ ! -s statTrimReads.md ]; then
        continue;
    fi

    if [ -e ${PREFIX}.trim.stats.txt ]; then
        echo >> statTrimReads.md
        echo '```text' >> statTrimReads.md
        echo "#${PREFIX}.trim" >> statTrimReads.md
        cat ${PREFIX}.trim.stats.txt |
            perl -nla -F"\t" -e '
                /^#(Matched|Name)/ and print and next;
                /^#/ and next;
                $F[2] =~ m{([\d.]+)} and $1 > 0.1 and print;
            ' \
            >> statTrimReads.md
        echo '```' >> statTrimReads.md
    fi

    if [ -e ${PREFIX}.filter.stats.txt ]; then
        echo >> statTrimReads.md
        echo '```text' >> statTrimReads.md
        echo "#${PREFIX}.filter" >> statTrimReads.md
        cat ${PREFIX}.filter.stats.txt |
            perl -nla -F"\t" -e '
                /^#(Matched|Name)/ and print and next;
                /^#/ and next;
                $F[2] =~ m{([\d.]+)} and $1 > 0.01 and print;
            ' \
            >> statTrimReads.md
        echo '```' >> statTrimReads.md
    fi

#    if [ -e ${PREFIX}.peaks.txt ]; then
#        echo >> statTrimReads.md
#        echo '```text' >> statTrimReads.md
#        echo "#${PREFIX}.peaks" >> statTrimReads.md
#        cat ${PREFIX}.peaks.txt |
#            grep "^#" \
#            >> statTrimReads.md
#        echo '```' >> statTrimReads.md
#    fi
done

if [ -s statTrimReads.md ]; then
    cat statTrimReads.md
    mv statTrimReads.md ${BASH_DIR}/../9_markdown
fi

cd ${BASH_DIR}/..
cd 2_illumina

parallel --no-run-if-empty --linebuffer -k -j 2 "
    ln -fs ./trim/Q{1}L{2}/ ./Q{1}L{2}
    " ::: {{ opt.qual }} ::: {{ opt.len }}
ln -fs ./trim ./Q0L0

log_info Done.

exit 0
