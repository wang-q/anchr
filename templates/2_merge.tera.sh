{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# run
#----------------------------#
log_warn 2_merge.sh

if [ -e 2_illumina/merge/pe.cor.fa.gz ]; then
    log_debug "2_illumina/merge/pe.cor.fa.gz presents"
    exit;
fi

mkdir -p 2_illumina/merge
cd 2_illumina/merge

START_TIME=$(date +%s)
save START_TIME

NUM_THREADS={{ opt.parallel }}
save NUM_THREADS

# save genome size
ESTIMATED_GENOME_SIZE={{ opt.genome }}
save ESTIMATED_GENOME_SIZE

for PREFIX in R S T; do
    if [ ! -e ../trim/${PREFIX}1.fq.gz ]; then
        continue;
    fi

    # rotate ascii characters https://en.wikipedia.org/wiki/ROT13
    PREFIXM=$(echo ${PREFIX} | tr 'A-Z' 'V-ZA-U')   # M N O
    PREFIXU=$(echo ${PREFIX} | tr 'A-Z' 'D-ZA-C')   # U V W

    if [ -e ${PREFIXM}1.fq ] || [ -e ${PREFIXM}1.fq.gz ]; then
        log_debug "2_illumina/merge/${PREFIXM}1.fq(.gz) presents"
        continue;
    fi

    anchr mergeread \
        ../trim/${PREFIX}1.fq.gz ../trim/${PREFIX}2.fq.gz ../trim/${PREFIX}s.fq.gz \
        --parallel {{ opt.parallel }} \
        --prefixm ${PREFIXM} \
        --prefixu ${PREFIXU} \
        -o mergeread.sh
    bash mergeread.sh

    # Create .cor.fa.gz
    anchr fq interleave \
        --name-prefix unmerged \
        ${PREFIXU}1.fq \
        ${PREFIXU}2.fq \
        > ${PREFIXM}.interleave.fa

    anchr fq interleave \
        --name-prefix single \
        ${PREFIXU}s.fq \
        >> ${PREFIXM}.interleave.fa

    anchr fq interleave \
        --name-prefix merged \
        ${PREFIXM}1.fq \
        >> ${PREFIXM}.interleave.fa

    # Shuffle interleaved read pairs.
    # .interleave.fa is FASTA (4 lines per pair: name1/seq1 then name2/seq2),
    # so group 4 lines per row to shuffle whole pairs.
    log_info Shuffle interleaved read pairs.
    cat ${PREFIXM}.interleave.fa |
        awk '{
            OFS="\t"; \
            getline seq; \
            getline name2; \
            getline seq2; \
            print $0,seq,name2,seq2}' |
        tsv-sample |
        awk '{OFS="\n"; print $1,$2,$3,$4}' \
        > ${PREFIXM}.cor.fa
    rm ${PREFIXM}.interleave.fa
    pigz -p {{ opt.parallel }} ${PREFIXM}.cor.fa

    log_info "stats of all .fq files"
    if [ ! -e statMergeReads.tsv ]; then
        printf "%s\t%s\t%s\t%s\n" \
            "Name" "N50" "Sum" "#" \
            > statMergeReads.tsv
    fi

    for NAME in clumped ecco eccc ecct extended merged.raw unmerged.raw unmerged.trim ${PREFIXM}1 ${PREFIXU}1 ${PREFIXU}2 ${PREFIXU}s; do
        if [ ! -e ${NAME}.fq ]; then
            continue
        fi
        printf "%s\t%s\t%s\t%s\n" \
            $(echo ${NAME}; stat_format_fq ${NAME}.fq;) >> statMergeReads.tsv
    done

    printf "%s\t%s\t%s\t%s\n" \
        $(echo ${PREFIXM}.cor; stat_format ${PREFIXM}.cor.fa.gz;) >> statMergeReads.tsv

    log_info "stats of insert sizes"
    if [ ! -e statMergeInsert.tsv ]; then
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "Group" "Mean" "Median" "STDev" "Pairs%" \
            > statMergeInsert.tsv
    fi

    #Mean	339.868
    #Median	312
    #Mode	251
    #STDev	134.676
    #PercentOfPairs	36.247
    for NAME in ${PREFIXM}.ihist.merge1.txt ${PREFIXM}.ihist.merge.txt; do
        printf "| %s " ${NAME} >> statMergeReads.md
        cat ${NAME} |
            GROUP="${NAME}" perl -nla -e '
                BEGIN { our $stat = { }; };

                m{\#(Mean|Median|STDev|PercentOfPairs)} or next;
                $stat->{$1} = $F[1];

                END {
                    printf qq(%s\t%.1f\t%s\t%.1f\t%.2f%%\n),
                        qq($ENV{GROUP}),
                        $stat->{Mean},
                        $stat->{Median},
                        $stat->{STDev},
                        $stat->{PercentOfPairs};
                }
                ' \
            >> statMergeInsert.tsv
    done

    log_info "clear unneeded .fq files"
    for NAME in temp clumped ecco eccc ecct extended merged.raw unmerged.raw unmerged.trim; do
        if [ -e ${NAME}.fq ]; then
            rm ${NAME}.fq
        fi
    done

    log_info "compress kept merge outputs"
    for NAME in M1 U1 U2 Us; do
        if [ -e ${NAME}.fq ]; then
            pigz -p {{ opt.parallel }} ${NAME}.fq
        fi
    done

done

log_debug "Combine .cor.fa.gz files"
if [ -e ../S1.fq.gz ]; then
    gzip -d -c [MNO].cor.fa.gz |
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
    rm [MNO].cor.fa.gz
else
    mv M.cor.fa.gz pe.cor.fa.gz
fi

log_debug "Reads stats with pgr fa"
SUM_OUT=$( pgr fa n50 -H -N 0 -S pe.cor.fa.gz )
save SUM_OUT

cat statMergeReads.tsv |
    tva to md stdin --right 2-4 \
    > statMergeReads.md

echo -e "\nTable: statMergeReads\n" >> statMergeReads.md

cat statMergeReads.md
mkdir -p ${BASH_DIR}/../9_markdown
mv statMergeReads.md ${BASH_DIR}/../9_markdown

cat statMergeInsert.tsv |
    tva to md stdin --right 2-5 \
    > statMergeInsert.md

echo -e "\nTable: statMergeInsert\n" >> statMergeInsert.md

cat statMergeInsert.md
mkdir -p ${BASH_DIR}/../9_markdown
mv statMergeInsert.md ${BASH_DIR}/../9_markdown

END_TIME=$(date +%s)
save END_TIME
RUNTIME=$((END_TIME-START_TIME))
save RUNTIME

log_info Done.

exit 0
