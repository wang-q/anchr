{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 8_megahit.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 DIR_READS"

DIR_READS=${1:-"2_illumina/trim"}

# Convert to abs path
DIR_READS="$(cd "$(dirname "$DIR_READS")"; pwd)/$(basename "$DIR_READS")"

if [ -e 8_megahit/anchor.fasta ]; then
    log_info "8_megahit/anchor.fasta presents"
    exit;
fi

#----------------------------#
# spades
#----------------------------#
if [ -e 8_megahit/final.contigs.fa ]; then
    log_info "8_megahit/final.contigs.fa presents"
else
    log_info "Run megahit"

    megahit \
        -t {{ opt.parallel }} \
        --k-list 31,41,51,61,71,81 \
        --12 ${DIR_READS}/pe.cor.fa.gz \
        --min-count 3 \
        -o 8_megahit

    log_info "Clear intermediate files"
    find . -type d -path "*8_megahit/*" -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

mkdir -p 8_megahit
cd 8_megahit

anchr asm anchor \
    final.contigs.fa \
    ${DIR_READS}/pe.cor.fa.gz \
    --mincov 5 --mscale 3 \
    --lscale {{ opt.lscale }} \
    --uscale {{ opt.uscale }} \
    -p {{ opt.parallel }} \
    --stats anchor.stats.tsv \
    -o anchor.fasta

exit 0;
