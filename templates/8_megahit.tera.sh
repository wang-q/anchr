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

if [ -e 8_megahit/anchor/anchor.fasta ]; then
    log_info "8_megahit/anchor/anchor.fasta presents"
    exit;
fi

#----------------------------#
# spades
#----------------------------#
if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    log_info "8_megahit/megahit.non-contained.fasta presents"
else
    log_info "Run megahit"

    megahit \
        -t {{ opt.parallel }} \
        --k-list 31,41,51,61,71,81 \
        --12 ${DIR_READS}/pe.cor.fa.gz \
        --min-count 3 \
        -o 8_megahit

    anchr contained \
        8_megahit/final.contigs.fa \
        --len 1000 --idt 0.98 --ratio 0.99999 --parallel 16 \
        -o stdout |
        hnsm filter -a 1000 stdin -o 8_megahit/megahit.non-contained.fasta

    log_info "Clear intermediate files"
    find . -type d -path "*8_megahit/*" -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

mkdir -p 8_megahit/anchor
cd 8_megahit/anchor

anchr anchors \
    ../megahit.non-contained.fasta \
    ${DIR_READS}/pe.cor.fa.gz \
    --readl {{ opt.readl }} \
    --uscale {{ opt.uscale }} \
    --lscale {{ opt.lscale }} \
    -p {{ opt.parallel }} \
    --ratio 0.98 \
    -o anchors.sh
bash anchors.sh

exit 0;
