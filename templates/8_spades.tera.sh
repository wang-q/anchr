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

if [ -e 8_spades/anchor.fasta ]; then
    log_info "8_spades/anchor.fasta presents"
    exit;
fi

#----------------------------#
# spades
#----------------------------#
if [ -e 8_spades/contigs.fasta ]; then
    log_info "8_spades/contigs.fasta presents"
else
    log_info "Run spades"

    mkdir -p 8_spades
    cd 8_spades

    mkdir -p re-pair
    pgr fa filter --min-len 60 ${DIR_READS}/pe.cor.fa.gz |
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

    log_info "Clear intermediate files"
    find . -type d -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

cd ${BASH_DIR}/..
mkdir -p 8_spades
cd 8_spades

anchr asm anchor \
    contigs.fasta \
    ${DIR_READS}/pe.cor.fa.gz \
    --mincov 5 --mscale 3 \
    --lscale {{ opt.lscale }} \
    --uscale {{ opt.uscale }} \
    -p {{ opt.parallel }} \
    --stats anchor.stats.tsv \
    -o anchor.fasta

exit 0;
