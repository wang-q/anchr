{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 8_mr_megahit.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0"

if [ -e 8_mr_megahit/anchor.fasta ]; then
    log_info "8_mr_megahit/anchor.fasta presents"
    exit;
fi

#----------------------------#
# megahit
#----------------------------#
if [ -e 8_mr_megahit/final.contigs.fa ]; then
    log_info "8_mr_megahit/final.contigs.fa presents"
else
    log_info "Run megahit"

    megahit \
        -t {{ opt.parallel }} \
        --k-min 45 --k-max 225 --k-step 26 \
        --12 ${BASH_DIR}/../2_illumina/merge/pe.cor.fa.gz \
        --min-count 3 \
        -o 8_mr_megahit

    log_info "Clear intermediate files"
    find . -type d -path "*8_mr_megahit/*" -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

mkdir -p 8_mr_megahit
cd 8_mr_megahit

anchr asm anchor \
    final.contigs.fa \
    ${BASH_DIR}/../2_illumina/merge/pe.cor.fa.gz \
    --mincov 5 --mscale 3 \
    --lscale {{ opt.lscale }} \
    --uscale {{ opt.uscale }} \
    -p {{ opt.parallel }} \
    --stats anchor.stats.tsv \
    -o anchor.fasta

exit 0;
