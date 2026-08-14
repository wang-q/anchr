{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
START_TIME=$(date +%s)

#----------------------------#
# Interleave reads（s-filter 单文件输入）
#----------------------------#
log_info 'Processing pe and/or se library reads'

{% if args | length == 2 -%}
anchr fq interleave \
    --fq --name-prefix pe \
    '{{ args.0 }}' \
    '{{ args.1 }}' \
    > 'pe.renamed.fastq'
{% elif args | length == 3 -%}
anchr fq interleave \
    --fq --name-prefix pe \
    '{{ args.0 }}' \
    '{{ args.1 }}' \
    > 'pe.renamed.fastq'

anchr fq interleave \
    --fq --name-prefix se \
    '{{ args.2 }}' \
    > 'se.renamed.fastq'
cat 'se.renamed.fastq' >> 'pe.renamed.fastq'
{% else -%}
anchr fq interleave \
    --fq --name-prefix pe \
    '{{ args.0 }}' \
    > 'pe.renamed.fastq'
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# Self-filter reads（quorum 替代，2026-08-15）
#----------------------------#
log_info "Self-filter reads with anchr fq s-filter"
anchr fq s-filter \
    pe.renamed.fastq \
    -k 24 --good 3 --anchor-count 4 --min-count 3 \
    -o {{ opt.prefix }}.cor.fastq \
    --discard-file {{ opt.prefix }}.discard.fastq

# Discarded read names（statQuorum 用；s-filter 输出没有 :sub:/trunc 标记）
awk 'NR%4==1' {{ opt.prefix }}.discard.fastq |
    sed 's/^@//' \
    > {{ opt.prefix }}.discard.lst

# FASTA 输出（下游依赖 pe.cor.fa.gz）
anchr fq to-fa {{ opt.prefix }}.cor.fastq \
    -o {{ opt.prefix }}.cor.fa

#----------------------------#
# Estimating genome size（jellyfish 替代）
#----------------------------#
log_info Estimating genome size.

{% if opt.estsize == 'auto' -%}
ESTIMATED_GENOME_SIZE=$(
    pgr kmer gsize {{ opt.prefix }}.cor.fastq -k 31 |
        awk '/^genome_size/ {print $2}'
)
save ESTIMATED_GENOME_SIZE
log_debug "Estimated genome size: $ESTIMATED_GENOME_SIZE"
{% else -%}
ESTIMATED_GENOME_SIZE={{ opt.estsize }}
save ESTIMATED_GENOME_SIZE
log_debug "You set ESTIMATED_GENOME_SIZE of $ESTIMATED_GENOME_SIZE"
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# Reads stats
#----------------------------#
log_debug "Reads stats with pgr fa"
SUM_IN=$( pgr fa n50 -H -N 0 -S pe.renamed.fastq )
save SUM_IN
SUM_OUT=$( pgr fa n50 -H -N 0 -S {{ opt.prefix }}.cor.fastq )
save SUM_OUT

# s-filter 的 k-mer 长度固定 24（quorum mer length 对齐）
KMER=24
save KMER

#----------------------------#
# Shuffle interleaved reads.
#----------------------------#
log_info Shuffle interleaved reads.
cat {{ opt.prefix }}.cor.fa |
    awk '{
        OFS="\t"; \
        getline seq; \
        print $0,seq}' |
    tsv-sample |
    awk '{OFS="\n"; print $1,$2}' \
    > {{ opt.prefix }}.cor.fa.tmp
mv {{ opt.prefix }}.cor.fa.tmp {{ opt.prefix }}.cor.fa
pigz -p {{ opt.parallel }} {{ opt.prefix }}.cor.fa

#----------------------------#
# Done.
#----------------------------#
find . -type f -name "pe.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "se.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "{{ opt.prefix }}.cor.fastq" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "{{ opt.prefix }}.discard.fastq" | parallel --no-run-if-empty -j 1 rm

save START_TIME

END_TIME=$(date +%s)
save END_TIME

RUNTIME=$((END_TIME-START_TIME))
save RUNTIME

log_info Done.

exit 0
