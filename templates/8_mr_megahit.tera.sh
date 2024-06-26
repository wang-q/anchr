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

if [ -e 8_mr_megahit/anchor/anchor.fasta ]; then
    log_info "8_mr_megahit/anchor/anchor.fasta presents"
    exit;
fi

#----------------------------#
# megahit
#----------------------------#
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    log_info "8_mr_megahit/megahit.non-contained.fasta presents"
else
    log_info "Run megahit"

    megahit \
        -t {{ opt.parallel }} \
        --k-min 45 --k-max 225 --k-step 26 \
        --12 ${BASH_DIR}/2_illumina/merge/pe.cor.fa.gz \
        --min-count 3 \
        -o 8_mr_megahit

    anchr contained \
        8_mr_megahit/final.contigs.fa \
        --len 1000 --idt 0.98 --ratio 0.99999 --parallel {{ opt.parallel }}  \
        -o stdout |
        faops filter -a 1000 -l 0 stdin 8_mr_megahit/megahit.non-contained.fasta

    log_info "Clear intermediate files"
    find . -type d -path "*8_mr_megahit/*" -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

mkdir -p 8_mr_megahit/anchor
cd 8_mr_megahit/anchor

anchr anchors \
    ../megahit.non-contained.fasta \
    ${BASH_DIR}/2_illumina/merge/pe.cor.fa.gz \
    --readl {{ opt.readl }} \
    --uscale {{ opt.uscale }} \
    --lscale {{ opt.lscale }} \
    -p {{ opt.parallel }} \
    --ratio 0.98 \
    -o anchors.sh
bash anchors.sh

find . -type f -name "pe.anchor.fa" | xargs rm
find . -type f -name "anchor.*.fasta" | xargs rm

exit 0;
