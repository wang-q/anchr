{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 3_bowtie.sh

if [ ! -e 1_genome/genome.fa ]; then
    log_info "1_genome/genome.fa does not exist"
    exit;
fi

mkdir -p 3_bowtie
cd 3_bowtie

ln -fs ../1_genome/genome.fa genome.fa

# bowtie2 index
if [ ! -e genome.fa.rev.1.bt2 ]; then
    bowtie2-build --threads {{ opt.parallel }} genome.fa genome.fa
fi

# chr.sizes
if [ ! -e chr.sizes ]; then
    faops size genome.fa > chr.sizes
fi

#----------------------------#
# Mapping
#----------------------------#
{% set parallel2 = opt.parallel | int - 3 -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
if [ ! -e R.sort.bai ]; then
    if [ -f ../2_illumina/trim/{{ opt.bowtie }}/pe.cor.fa.gz ]; then
        gzip -dcf \
            ../2_illumina/trim/{{ opt.bowtie }}/pe.cor.fa.gz
    elif [ -f ../2_illumina/trim/{{ opt.bowtie }}/R1.fq.gz ]; then
        gzip -dcf \
            ../2_illumina/trim/{{ opt.bowtie }}/R1.fq.gz \
            ../2_illumina/trim/{{ opt.bowtie }}/R2.fq.gz \
            ../2_illumina/trim/{{ opt.bowtie }}/Rs.fq.gz |
            faops filter -l 0 stdin stdout # ignore QUAL
    fi |
        bowtie2 -p {{ parallel2 }} --very-fast -t \
            -x genome.fa \
            -f -U /dev/stdin \
            2> >(tee bowtie.R.log >&2) |
        picard CleanSam \
            --INPUT /dev/stdin \
            --OUTPUT /dev/stdout \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
        picard SortSam \
            --INPUT /dev/stdin \
            --OUTPUT R.sort.bam \
            --SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 1

    picard BuildBamIndex \
        --INPUT R.sort.bam \
        --VALIDATION_STRINGENCY LENIENT
fi

#----------------------------#
# Depth
#----------------------------#
if [ ! -f R.mosdepth.summary.txt ]; then
    if hash mosdepth 2>/dev/null; then
        # depth
        mosdepth R R.sort.bam

        # covered
        gzip -dcf R.per-base.bed.gz |
            perl -nla -F"\t" -e '
                $F[3] == 0 and next;
                $start = $F[1] + 1;
                $end = $F[2];
                if ($start == $F[2]) {
                    print qq($F[0]:$start);
                }
                else {
                    print qq($F[0]:$start-$end);
                }
            ' |
            spanr cover stdin -o covered.yml

        spanr stat chr.sizes covered.yml -o stdout |
            grep -v "^all" |
            sed 's/^chr/chrom/' |
            sed 's/,size/,covLength/' |
            sed 's/,coverage/,covRate/' |
            sed 's/,/\t/g' \
            > coverage.tsv

        # join
        cat coverage.tsv |
            tsv-join -H --filter-file R.mosdepth.summary.txt \
                --key-fields chrom --append-fields 3-6 \
            > join.tsv

    fi
fi

log_info Done.

exit 0
