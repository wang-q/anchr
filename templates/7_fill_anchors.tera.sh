{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 7_fill_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 FILE_ANCHOR FILE_LONG GAP_COV"

if [ "$#" -lt 2 ]; then
    echo >&2 "$USAGE"
    exit 1
fi

FILE_ANCHOR=$1
FILE_LONG=$2
GAP_COV=${3:-2}

if [ -e 7_fill_anchors/contig.fasta ]; then
    echo >&2 "7_fill_anchors/contig.fasta presents"
    exit;
fi

#----------------------------#
# fill anchors
#----------------------------#
mkdir -p 7_fill_anchors

log_info "overlap: between anchor-long"

dazz overlap2 \
    --parallel {{ opt.parallel }} \
    ${FILE_ANCHOR} \
    ${FILE_LONG} \
    -d 7_fill_anchors \
    -b 50 --len 1000 --idt 0.995 --all

cd 7_fill_anchors

CONTIG_COUNT=$(faops n50 -H -N 0 -C anchor.fasta)
log_debug "contigs: ${CONTIG_COUNT}"

log_info "group: anchor-long"
rm -fr group
dazz group \
    anchorLong.db \
    anchorLong.ovlp.tsv \
    --parallel {{ opt.parallel }} \
    --keep \
    --range "1-${CONTIG_COUNT}" --len 1000 --idt 0.995 --max {{ opt.fillmax }} -c ${GAP_COV}

log_info "Processing each groups"
{% set parallel2 = opt.parallel | int / 2 -%}
{% set parallel2 = parallel2 | round(method="floor") -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
cat group/groups.txt |
    parallel --no-run-if-empty --linebuffer -k -j {{ parallel2 }} '
        echo {};
        dazz orient \
            --len 1000 --idt 0.995 \
            group/{}.anchor.fasta \
            group/{}.long.fasta \
            -r group/{}.restrict.tsv \
            -o group/{}.strand.fasta;

        dazz overlap --len 1000 --idt 0.995 --all \
            group/{}.strand.fasta \
            -o stdout |
            ovlpr restrict \
                stdin group/{}.restrict.tsv \
                -o group/{}.ovlp.tsv;

        dazz layout \
            group/{}.strand.fasta \
            group/{}.ovlp.tsv \
            group/{}.relation.tsv \
            -o group/{}.contig.fasta
    '

log_info "Build contigs"
cat \
   group/non_grouped.fasta \
   group/*.contig.fasta |
   faops filter -l 0 -a 1000 stdin contig.fasta

log_info Done.

exit 0
