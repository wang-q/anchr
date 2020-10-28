{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 2_kat.sh

mkdir -p 2_illumina/kat
cd 2_illumina/kat

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -e ${PREFIX}-gcp-21.dist_analysis.json ]; then
        continue;
    fi

    for KMER in 21 31 41 51 61 71; do
        log_info "PREFIX: ${PREFIX}; KMER: ${KMER}"

        kat gcp \
            -t {{ opt.parallel }} -m ${KMER} \
            ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
            -o ${PREFIX}-gcp-${KMER}
    done
done

find . -type f -name "*.mx" | parallel --no-run-if-empty -j 1 rm

echo -e "Table: statKAT\n" > statKAT.md

for PREFIX in R S T; do
    find . -type f -name "${PREFIX}-gcp*.dist_analysis.json" |
        sort |
        xargs cat |
        sed 's/%//g' |
        jq "{
            k: (\"${PREFIX}.\" + (.coverage.k | tostring)),
            mean_freq: .coverage.mean_freq,
            est_genome_size: .coverage.est_genome_size,
            est_het_rate: .coverage.est_het_rate,
            mean_gc: .gc.mean_gc,
        }"
done |
    mlr --ijson --otsv cat |
    perl -nla -F"\t" -e '
        /mean_gc/ and print and next;
        $F[3] = sprintf q(%.4f), $F[3];
        $F[4] = sprintf q(%.2f), $F[4];
        print join qq(\t), @F;
    ' |
    mlr --itsv --omd cat \
    >> statKAT.md

cat statKAT.md
mv statKAT.md ../../

exit 0
