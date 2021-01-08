{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 3_bwa.sh

if [ ! -e 1_genome/genome.fa ]; then
    log_info "1_genome/genome.fa does not exist"
    exit;
fi

mkdir -p 3_bwa
cd 3_bwa

ln -fs ../1_genome/genome.fa genome.fa

# bwa index
if [ ! -e genome.fa.bwt ]; then
    bwa index genome.fa
fi

# faidx
if [ ! -e genome.fa.fai ]; then
    samtools faidx genome.fa
fi

# dict
if [ ! -e genome.dict ]; then
    picard CreateSequenceDictionary --REFERENCE genome.fa
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
SAMPLE=$(readlinkf .. | xargs dirname)
if [ ! -e R.sort.bai ]; then
    bwa mem -t {{ parallel2 }} \
        -M -K 100000000 -v 3 -Y \
        genome.fa \
        ../2_illumina/{{ opt.bwa }}/R1.fq.gz \
        ../2_illumina/{{ opt.bwa }}/R2.fq.gz \
        2> >(tee R.bwa.log >&2) |
        samtools view -F 4 -u - | # Remove unmapped reads, write uncompressed BAM output
        picard AddOrReplaceReadGroups \
            --INPUT /dev/stdin \
            --OUTPUT /dev/stdout \
            --RGLB ${SAMPLE} \
            --RGPL ILLUMINA \
            --RGPU ${SAMPLE} \
            --RGSM ${SAMPLE} \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
        picard CleanSam \
            --INPUT /dev/stdin \
            --OUTPUT /dev/stdout \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
        picard FixMateInformation \
            --INPUT /dev/stdin \
            --OUTPUT R.mate.bam \
            --SORT_ORDER queryname \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 1

    picard MarkDuplicates \
        --INPUT R.mate.bam \
        --OUTPUT /dev/stdout \
        --METRICS_FILE R.dedup.metrics \
        --ASSUME_SORT_ORDER queryname \
        --REMOVE_DUPLICATES true \
        --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
    picard SortSam \
        --INPUT /dev/stdin \
        --OUTPUT R.sort.bam \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 1

    picard BuildBamIndex \
        --INPUT R.sort.bam \
        --VALIDATION_STRINGENCY LENIENT

    picard CollectWgsMetrics \
        --INPUT R.sort.bam \
        --REFERENCE_SEQUENCE genome.fa \
        --OUTPUT R.wgs.metrics \
        --VALIDATION_STRINGENCY SILENT
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
