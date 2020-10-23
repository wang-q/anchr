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
            log_info "    faops interleave"

            # Create .cor.fa.gz
            faops interleave \
                -p pe \
                ${PREFIX}1.fq.gz \
{% if opt.se == "0" -%}
                ${PREFIX}2.fq.gz \
{% endif -%}
                > ${PREFIX}.interleave.fa

            if [ -e ${PREFIX}s.fq.gz ]; then
                faops interleave \
                    -p se \
                    ${PREFIX}s.fq.gz \
                    >> ${PREFIX}.interleave.fa
            fi

            # Shuffle interleaved reads.
            log_info Shuffle interleaved reads.
            cat ${PREFIX}.interleave.fa |
                awk '{
                    OFS="\t"; \
                    getline seq; \
                    getline name2; \
                    getline seq2; \
                    print $0,seq,name2,seq2}' |
                shuf |
                awk '{OFS="\n"; print $1,$2,$3,$4}' \
                > ${PREFIX}.cor.fa
            rm ${PREFIX}.interleave.fa
            pigz -p {{ opt.parallel }} ${PREFIX}.cor.fa

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
                shuf |
                awk '{OFS="\n"; print $1,$2,$3,$4}' \
                > pe.cor.fa
            pigz -p {{ opt.parallel }} pe.cor.fa
            rm [RST].cor.fa.gz
        else
            mv R.cor.fa.gz pe.cor.fa.gz
        fi
    done
done

log_info Done.

exit 0
