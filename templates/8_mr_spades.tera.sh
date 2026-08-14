{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 8_mr_spades.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0"

if [ -s 8_mr_spades/anchor.fasta ]; then
    log_info "8_mr_spades/anchor.fasta presents"
    exit;
fi

#----------------------------#
# spades
#----------------------------#
if [ -e 8_mr_spades/contigs.fasta ]; then
    log_info "8_mr_spades/contigs.fasta presents"
else
    log_info "Run spades"

    mkdir -p 8_mr_spades
    cd 8_mr_spades

    mkdir -p re-pair
    # pe.cor.fa.gz is fully shuffled by 2_merge (tsv-sample), so use the
    # name-indexed "repair" (rp) mode; "fint" only fixes partially-broken
    # interleaving and cannot pair fully shuffled reads.
    pgr fa filter --min-len 60 ${BASH_DIR}/../2_illumina/merge/pe.cor.fa.gz |
        repair.sh \
            in=stdin.fa \
            out=re-pair/R1.fa \
            out2=re-pair/R2.fa \
            outs=re-pair/Rs.fa \
            threads={{ opt.parallel }} \
            repair overwrite

    # spades seems ignore non-properly paired reads
    spades.py \
        -t {{ opt.parallel }} \
        --only-assembler \
        -k 25,55,95,125 \
        -1 re-pair/R1.fa \
        -2 re-pair/R2.fa \
        -s re-pair/Rs.fa \
        -o .

    log_info "Clear intermediate files"
    find . -type d -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

cd ${BASH_DIR}/..
mkdir -p 8_mr_spades
cd 8_mr_spades

anchr asm anchor \
    contigs.fasta \
    ${BASH_DIR}/../2_illumina/merge/pe.cor.fa.gz \
    --mincov 5 --mscale 3 \
    --lscale {{ opt.lscale }} \
    --uscale {{ opt.uscale }} \
    -p {{ opt.parallel }} \
    --stats anchor.stats.tsv \
    -o anchor.fasta

exit 0;
