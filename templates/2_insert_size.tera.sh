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

    anchr asm contig \
        ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
        -o ${PREFIX}.tadpole.contig.fasta

    anchr asm map \
        ${PREFIX}.tadpole.contig.fasta \
        ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
        --paired \
        --max-reads {{ opt.reads }} \
        --outm ${PREFIX}.tadpole.sam

    anchr sam ihist \
        ${PREFIX}.tadpole.sam \
        -o ${PREFIX}.ihist.tadpole.txt

    picard SortSam \
        -I ${PREFIX}.tadpole.sam \
        -O ${PREFIX}.tadpole.sort.bam \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY LENIENT

    picard CollectInsertSizeMetrics \
        -I ${PREFIX}.tadpole.sort.bam \
        -O ${PREFIX}.insert_size.tadpole.txt \
        --Histogram_FILE ${PREFIX}.insert_size.tadpole.pdf

    if [ -e ../../1_genome/genome.fa ]; then
        anchr asm map \
            ../../1_genome/genome.fa \
            ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
            --paired \
            --max-reads {{ opt.reads }} \
            --outm ${PREFIX}.genome.sam

        anchr sam ihist \
            ${PREFIX}.genome.sam \
            -o ${PREFIX}.ihist.genome.txt

        picard SortSam \
            -I ${PREFIX}.genome.sam \
            -O ${PREFIX}.genome.sort.bam \
            --SORT_ORDER coordinate

        picard CollectInsertSizeMetrics \
            -I ${PREFIX}.genome.sort.bam \
            -O ${PREFIX}.insert_size.genome.txt \
            --Histogram_FILE ${PREFIX}.insert_size.genome.pdf
    fi

    find . -name "${PREFIX}.*.sam" -or -name "${PREFIX}.*.sort.bam" |
        parallel --no-run-if-empty -j 1 rm
done

printf "%s\t%s\t%s\t%s\t%s\n" \
    "Group" "Mean" "Median" "STDev" "Pairs%/Orientation" \
    > statInsertSize.tsv

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

        cat ${PREFIX}.ihist.${G}.txt |
            GROUP="${PREFIX}.${G}" perl -nla -e '
                BEGIN { our $stat = { }; };

                m{\#(Mean|Median|STDev|PercentOfPairs)} or next;
                $stat->{$1} = $F[1];

                END {
                    printf qq(%s\t%.1f\t%s\t%.1f\t%.2f%%\n),
                        qq($ENV{GROUP}.bbtools),
                        $stat->{Mean},
                        $stat->{Median},
                        $stat->{STDev},
                        $stat->{PercentOfPairs};
                }
                '
    done
done \
    >> statInsertSize.tsv

# picard CollectInsertSizeMetrics
#MEDIAN_INSERT_SIZE	MODE_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	PAIR_ORIENTATION	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_95_PERCENT	WIDTH_OF_99_PERCENT	SAMPLE	LIBRARY	READ_GROUP
#296	287	14	92	501	294.892521	21.587526	1611331	FR	7	11	17	23	29	35	41	49	63	81	145
for PREFIX in R S T; do
    for G in genome tadpole; do
        if [ ! -e ${PREFIX}.insert_size.${G}.txt ]; then
            continue;
        fi

        cat ${PREFIX}.insert_size.${G}.txt |
            GROUP="${PREFIX}.${G}" perl -nla -F"\t" -e '
                next if @F < 9;
                next unless /^\d/;
                printf qq(%s\t%.1f\t%s\t%.1f\t%s\n),
                    qq($ENV{GROUP}.picard),
                    $F[5],
                    $F[0],
                    $F[6],
                    $F[8];
                '
    done
done \
    >> statInsertSize.tsv

cat statInsertSize.tsv |
    tva to md stdin --right 2-5 \
    > statInsertSize.md

echo -e "\nTable: statInsertSize\n" >> statInsertSize.md

cat statInsertSize.md
mkdir -p ${BASH_DIR}/../9_markdown
mv statInsertSize.md ${BASH_DIR}/../9_markdown
