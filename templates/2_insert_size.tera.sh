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

    if [ -e ${PREFIX}.ihist.contig.txt ]; then
        continue;
    fi

    anchr asm contig \
        ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
        -o ${PREFIX}.contig.fasta

    anchr asm map \
        ${PREFIX}.contig.fasta \
        ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
        --paired \
        --max-reads {{ opt.reads }} \
        --outm ${PREFIX}.contig.sam

    anchr sam ihist \
        ${PREFIX}.contig.sam \
        -o ${PREFIX}.ihist.contig.txt

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
    fi

    find . -name "${PREFIX}.*.sam" |
        parallel --no-run-if-empty -j 1 rm
done

printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "Group" "Mean" "Median" "STDev" "Pairs%" "Orientation" \
    > statInsertSize.tsv

# anchr sam ihist
#Mean	339.868
#Median	312
#Mode	251
#STDev	134.676
#PercentOfPairs	36.247
for PREFIX in R S T; do
    for G in genome contig; do
        if [ ! -e ${PREFIX}.ihist.${G}.txt ]; then
            continue;
        fi

        cat ${PREFIX}.ihist.${G}.txt |
            GROUP="${PREFIX}.${G}" perl -nla -e '
                BEGIN { our $stat = { }; };

                m{\#(Mean|Median|STDev|PercentOfPairs|Orientation)} or next;
                $stat->{$1} = $F[1];

                END {
                    printf qq(%s\t%.1f\t%s\t%.1f\t%.2f%%\t%s\n),
                        qq($ENV{GROUP}),
                        $stat->{Mean},
                        $stat->{Median},
                        $stat->{STDev},
                        $stat->{PercentOfPairs},
                        $stat->{Orientation};
                }
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
