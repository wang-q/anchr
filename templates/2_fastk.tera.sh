{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 2_fastk.sh

mkdir -p 2_illumina/fastk
cd 2_illumina/fastk

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -d ${PREFIX}-GeneScope-21 ]; then
        continue;
    fi

    for KMER in 21 51 81; do
        log_info "PREFIX: ${PREFIX}; KMER: ${KMER}"

        log_info "FastK"
        FastK -v -T{{ opt.parallel }} -t1 -k${KMER} \
            ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
            -NTable-${KMER}

        log_info "GeneScope"
        Histex -G Table-${KMER} |
            Rscript ../../0_script/genescopefk.R -k ${KMER} -p 1 -o ${PREFIX}-GeneScope-${KMER}

        KatGC -T{{ opt.parallel }} -x1.9 -s Table-${KMER} ${PREFIX}-Merqury-KatGC-${KMER}

        Fastrm Table-${KMER}
    done
done

for PREFIX in R S T; do
    for KMER in 21 51 81; do
        if [ ! -e ${PREFIX}-GeneScope-${KMER}/summary.txt ]; then
            continue
        fi

        COV=$(
            cat ${PREFIX}-GeneScope-${KMER}/model.txt |
                grep '^kmercov' |
                tr -s ' ' '\t' |
                cut -f 2 |
                perl -nl -e 'printf qq(%.1f\n), $_'
        )

        cat ${PREFIX}-GeneScope-${KMER}/summary.txt |
            sed '1,6 d' |
            sed "1 s/^/K\t/" |
            sed "2 s/^/${PREFIX}.${KMER}\t/" |
            sed "3,7 s/^/\t/" |
            perl -nlp -e 's/\s{2,}/\t/g; s/\s+$//g;' |
            perl -nla -F'\t' -e '
                @fields = map {/\bNA\b/ ? q() : $_ } @F;        # Remove NA fields
                $fields[2] = q() if $fields[2] eq $fields[3];   # Remove identical fields
                print join qq(\t), @fields
            '

        printf "\tKmer Cov\t\t${COV}\n"
    done
done |
    keep-header -- grep -v 'property' \
    > statFastK.tsv

cat statFastK.tsv |
    mlr --itsv --omd cat |
    perl -nlp -e '$. == 2 and $_ = q(|:---|:---|---:|---:|)' \
    > statFastK.md

echo -e "\nTable: statFastK\n" >> statFastK.md

cat statFastK.md
mv statFastK.md ${BASH_DIR}/../9_markdown
