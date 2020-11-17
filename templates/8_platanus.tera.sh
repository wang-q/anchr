{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 8_platanus.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 DIR_READS"

DIR_READS=${1:-"2_illumina/trim"}

DIR_READS="$(cd "$(dirname "$DIR_READS")"; pwd)/$(basename "$DIR_READS")"

if [ -e 8_platanus/anchor/anchor.fasta ]; then
    log_info "8_platanus/anchor/anchor.fasta presents"
    exit;
fi

#----------------------------#
# platanus
#----------------------------#
mkdir -p 8_platanus
cd 8_platanus

if [ -e platanus.non-contained.fasta ]; then
    log_info "platanus.non-contained.fasta presents"
else
    log_info "Run platanus"

    faops filter -l 0 -a 60 ${DIR_READS}/pe.cor.fa.gz stdout |
        repair.sh \
            in=stdin.fa \
            out=pe.fa \
            outs=se.fa \
            threads={{ opt.parallel }} \
            fint overwrite

    if [ -s pe.fa ]; then
        platanus assemble -t {{ opt.parallel }} -m 100 \
            -f pe.fa \
            $(
                if [ -s se.fa ]; then
                    echo se.fa
                fi
            ) \
            2>&1 | tee ass_log.txt

        platanus scaffold -t {{ opt.parallel }} \
            -c out_contig.fa -b out_contigBubble.fa \
            -ip1 pe.fa \
            2>&1 | tee sca_log.txt

        platanus gap_close -t {{ opt.parallel }} \
            -c out_scaffold.fa \
            -ip1 pe.fa \
            2>&1 | tee gap_log.txt

        dazz contained \
            out_gapClosed.fa \
            --len 1000 --idt 0.98 --proportion 0.99999 --parallel {{ opt.parallel }} \
            -o stdout |
            faops filter -a 1000 -l 0 stdin platanus.non-contained.fasta
    else
        platanus assemble -t {{ opt.parallel }} -m 100 \
            -f se.fa \
            2>&1 | tee ass_log.txt

        dazz contained \
            out_contig.fa out_contigBubble.fa \
            --len 1000 --idt 0.98 --proportion 0.99999 --parallel {{ opt.parallel }} \
            -o stdout |
            faops filter -a 1000 -l 0 stdin platanus.non-contained.fasta
    fi

    log_info "Clear intermediate files"
    find . -type f -name "[ps]e.fa" | parallel --no-run-if-empty -j 1 rm
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

cd ${BASH_DIR}
mkdir -p 8_platanus/anchor
cd 8_platanus/anchor

anchr anchors \
    ../platanus.non-contained.fasta \
    ${DIR_READS}/pe.cor.fa.gz \
    --readl {{ opt.readl }} \
    --uscale {{ opt.uscale }} \
    --lscale {{ opt.lscale }} \
    -p {{ opt.parallel }} \
    --ratio 0.98 \
    -o anchors.sh
bash anchors.sh

exit 0;
