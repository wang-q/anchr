{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 8_spades.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 DIR_READS"

DIR_READS=${1:-"2_illumina/trim"}

# Convert to abs path
DIR_READS="$(cd "$(dirname "$DIR_READS")"; pwd)/$(basename "$DIR_READS")"

if [ -e 8_spades/anchor/anchor.fasta ]; then
    log_info "8_spades/anchor/anchor.fasta presents"
    exit;
fi

#----------------------------#
# spades
#----------------------------#
if [ -e 8_spades/spades.non-contained.fasta ]; then
    log_info "8_spades/spades.non-contained.fasta presents"
else
    log_info "Run spades"

    mkdir -p 8_spades
    cd 8_spades

    mkdir -p re-pair
    hnsm filter -a 60 ${DIR_READS}/pe.cor.fa.gz |
        repair.sh \
            in=stdin.fa \
            out=re-pair/R1.fa \
            out2=re-pair/R2.fa \
            outs=re-pair/Rs.fa \
            threads={{ opt.parallel }} \
            fint overwrite

    # spades seems ignore non-properly paired reads
    spades.py \
        -t {{ opt.parallel }} \
        --only-assembler \
        -k 21,33,55,77 \
        -1 re-pair/R1.fa \
        -2 re-pair/R2.fa \
        -o .

    anchr contained \
        contigs.fasta \
        --len 1000 --idt 0.98 --ratio 0.99999 --parallel {{ opt.parallel }} \
        -o stdout |
        hnsm filter -a 1000 stdin -o spades.non-contained.fasta

    log_info "Clear intermediate files"
    find . -type d -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

cd ${BASH_DIR}/..
mkdir -p 8_spades/anchor
cd 8_spades/anchor

anchr anchors \
    ../spades.non-contained.fasta \
    ${DIR_READS}/pe.cor.fa.gz \
    --readl {{ opt.readl }} \
    --uscale {{ opt.uscale }} \
    --lscale {{ opt.lscale }} \
    -p {{ opt.parallel }} \
    --ratio 0.98 \
    -o anchors.sh
bash anchors.sh

exit 0;
