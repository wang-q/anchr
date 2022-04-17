{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 2_quorum.sh

for Q in 0 {{ opt.qual }}; do
    for L in 0 {{ opt.len }}; do
        cd ${BASH_DIR}

        if [ ! -d 2_illumina/Q${Q}L${L} ]; then
            continue;
        fi

        if [ -e 2_illumina/Q${Q}L${L}/pe.cor.fa.gz ]; then
            log_debug "2_illumina/Q${Q}L${L}/pe.cor.fa.gz presents"
            continue;
        fi

        START_TIME=$(date +%s)

        cd 2_illumina/Q${Q}L${L}

        for PREFIX in R S T; do
            if [ ! -e ${PREFIX}1.fq.gz ]; then
                continue;
            fi

            if [ -e ${PREFIX}.cor.fa.gz ]; then
                echo >&2 "    ${PREFIX}.cor.fa.gz exists"
                continue;
            fi

            log_info "Qual-Len: Q${Q}L${L}.${PREFIX}"

            anchr quorum \
                ${PREFIX}1.fq.gz \
{% if opt.se == "0" -%}
                ${PREFIX}2.fq.gz \
                $(
                    if [ -s ${PREFIX}s.fq.gz ]; then
                        echo ${PREFIX}s.fq.gz;
                    fi
                ) \
{% endif -%}
                -p {{ opt.parallel }} \
                --prefix ${PREFIX} \
                -o quorum.sh
            bash quorum.sh

            log_info "statQuorum.${PREFIX}"

            SUM_IN=$( cat env.json | jq '.SUM_IN | tonumber' )
            SUM_OUT=$( cat env.json | jq '.SUM_OUT | tonumber' )
            EST_G=$( cat env.json | jq '.ESTIMATED_GENOME_SIZE | tonumber' )
            SECS=$( cat env.json | jq '.RUNTIME | tonumber' )

            printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s |\n" \
                "Q${Q}L${L}.${PREFIX}" \
                $( perl -e "printf qq(%.1f), ${SUM_IN} / {{ opt.genome }};" ) \
                $( perl -e "printf qq(%.1f), ${SUM_OUT} / {{ opt.genome }};" ) \
                $( perl -e "printf qq(%.2f%%), (1 - ${SUM_OUT} / ${SUM_IN}) * 100;" ) \
                $( cat env.json | jq '.KMER' ) \
                $( perl -MNumber::Format -e "print Number::Format::format_bytes({{ opt.genome }}, base => 1000,);" ) \
                $( perl -MNumber::Format -e "print Number::Format::format_bytes(${EST_G}, base => 1000,);" ) \
                $( perl -e "printf qq(%.2f), ${EST_G} / {{ opt.genome }}" ) \
                $( printf "%d:%02d'%02d''\n" $((${SECS}/3600)) $((${SECS}%3600/60)) $((${SECS}%60)) ) \
                > statQuorum.${PREFIX}

        done

        log_info "Combine Q${Q}L${L} .cor.fa.gz files"
        if [ -e S1.fq.gz ]; then
            gzip -d -c [RST].cor.fa.gz |
                awk '{
                    OFS="\t"; \
                    getline seq; \
                    getline name2; \
                    getline seq2; \
                    print $0,seq,name2,seq2}' |
                tsv-sample |
                awk '{OFS="\n"; print $1,$2,$3,$4}' \
                > pe.cor.fa
            pigz -p {{ opt.parallel }} pe.cor.fa
            rm [RST].cor.fa.gz
        else
            mv R.cor.fa.gz pe.cor.fa.gz
        fi

        rm env.json
        log_debug "Reads stats with faops"
        SUM_OUT=$( faops n50 -H -N 0 -S pe.cor.fa.gz )
        save SUM_OUT

        save START_TIME
        END_TIME=$(date +%s)
        save END_TIME
        RUNTIME=$((END_TIME-START_TIME))
        save RUNTIME

        # save genome size
        ESTIMATED_GENOME_SIZE={{ opt.genome }}
        save ESTIMATED_GENOME_SIZE

    done
done

cd ${BASH_DIR}/2_illumina

if [ -e Q0L0/statQuorum.R ]; then
    echo -e "Table: statQuorum\n" > statQuorum.md
    printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s |\n" \
        "Name" "CovIn" "CovOut" "Discard%" \
        "Kmer" "RealG" "EstG" "Est/Real" \
        "RunTime" \
        >> statQuorum.md
    printf "|:--|--:|--:|--:|--:|--:|--:|--:|--:|\n" \
        >> statQuorum.md

    for PREFIX in R S T; do
        for Q in 0 {{ opt.qual }}; do
            for L in 0 {{ opt.len }}; do
                if [ -e Q${Q}L${L}/statQuorum.${PREFIX} ]; then
                    cat Q${Q}L${L}/statQuorum.${PREFIX} >> statQuorum.md;
                fi
            done
        done
    done

    cat statQuorum.md
    mv statQuorum.md ../
fi

log_info Done.

exit 0
