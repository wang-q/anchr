{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 2_insert_size.sh


mkdir -p 2_illumina/insert_size
cd 2_illumina/insert_size

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -e ${PREFIX}.ihist.tadpole.txt ]; then
        continue;
    fi

    tadpole.sh \
        in=../${PREFIX}1.fq.gz \
        in2=../${PREFIX}2.fq.gz \
        out=${PREFIX}.tadpole.contig.fasta \
        threads={{ opt.parallel }} \
        overwrite {% if opt.prefilter == "1" %}prefilter={{ opt.prefilter }}{% endif %}

    cat ${PREFIX}.tadpole.contig.fasta |
        faops dazz -l 0 -p T stdin stdout \
        > ${PREFIX}.tadpole.contig.fa

    bbmap.sh \
        in=../${PREFIX}1.fq.gz \
        in2=../${PREFIX}2.fq.gz \
        out=${PREFIX}.tadpole.sam.gz \
        ref=${PREFIX}.tadpole.contig.fa \
        threads={{ opt.parallel }} \
        pairedonly \
        reads=1000000 \
        nodisk overwrite

    reformat.sh \
        in=${PREFIX}.tadpole.sam.gz \
        ihist=${PREFIX}.ihist.tadpole.txt \
        overwrite

    picard SortSam \
        -I ${PREFIX}.tadpole.sam.gz \
        -O ${PREFIX}.tadpole.sort.bam \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY LENIENT

    picard CollectInsertSizeMetrics \
        -I ${PREFIX}.tadpole.sort.bam \
        -O ${PREFIX}.insert_size.tadpole.txt \
        --Histogram_FILE ${PREFIX}.insert_size.tadpole.pdf

    if [ -e ../../1_genome/genome.fa ]; then
        bbmap.sh \
            in=../${PREFIX}1.fq.gz \
            in2=../${PREFIX}2.fq.gz \
            out=${PREFIX}.genome.sam.gz \
            ref=../../1_genome/genome.fa \
            threads={{ opt.parallel }} \
            maxindel=0 strictmaxindel \
            reads=1000000 \
            nodisk overwrite

        reformat.sh \
            in=${PREFIX}.genome.sam.gz \
            ihist=${PREFIX}.ihist.genome.txt \
            overwrite

        picard SortSam \
            -I ${PREFIX}.genome.sam.gz \
            -O ${PREFIX}.genome.sort.bam \
            --SORT_ORDER coordinate

        picard CollectInsertSizeMetrics \
            -I ${PREFIX}.genome.sort.bam \
            -O ${PREFIX}.insert_size.genome.txt \
            --Histogram_FILE ${PREFIX}.insert_size.genome.pdf
    fi

    find . -name "${PREFIX}.*.sam.gz" -or -name "${PREFIX}.*.sort.bam" |
        parallel --no-run-if-empty -j 1 rm
done

echo -e "Table: statInsertSize\n" > statInsertSize.md
printf "| %s | %s | %s | %s | %s |\n" \
    "Group" "Mean" "Median" "STDev" "PercentOfPairs/PairOrientation" \
    >> statInsertSize.md
printf "|:--|--:|--:|--:|--:|\n" >> statInsertSize.md

# bbtools reformat.sh
#Mean	339.868
#Median	312
#Mode	251
#STDev	134.676
#PercentOfPairs	36.247
for PREFIX in R S T; do
    for G in genome tadpole; do
        if [ ! -e ${PREFIX}.ihist.${G}.txt ]; then
            continue;
        fi

        printf "| %s " "${PREFIX}.${G}.bbtools" >> statInsertSize.md
        cat ${PREFIX}.ihist.${G}.txt \
            | perl -nla -e '
                BEGIN { our $stat = { }; };

                m{\#(Mean|Median|STDev|PercentOfPairs)} or next;
                $stat->{$1} = $F[1];

                END {
                    printf qq{| %.1f | %s | %.1f | %.2f%% |\n},
                        $stat->{Mean},
                        $stat->{Median},
                        $stat->{STDev},
                        $stat->{PercentOfPairs};
                }
                ' \
            >> statInsertSize.md
    done
done

# picard CollectInsertSizeMetrics
#MEDIAN_INSERT_SIZE	MODE_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	PAIR_ORIENTATION	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_95_PERCENT	WIDTH_OF_99_PERCENT	SAMPLE	LIBRARY	READ_GROUP
#296	287	14	92	501	294.892521	21.587526	1611331	FR	7	11	17	23	29	35	41	49	63	81	145
for PREFIX in R S T; do
    for G in genome tadpole; do
        if [ ! -e ${PREFIX}.insert_size.${G}.txt ]; then
            continue;
        fi

        cat ${PREFIX}.insert_size.${G}.txt \
            | GROUP="${PREFIX}.${G}" perl -nla -F"\t" -e '
                next if @F < 9;
                next unless /^\d/;
                printf qq{| %s | %.1f | %s | %.1f | %s |\n},
                    qq{$ENV{GROUP}.picard},
                    $F[5],
                    $F[0],
                    $F[6],
                    $F[8];
                ' \
            >> statInsertSize.md
    done
done

cat statInsertSize.md

mv statInsertSize.md ../../
