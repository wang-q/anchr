{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_quast.sh
{% set unitiggers = opt.unitigger | split(pat=" ") -%}

QUAST_TARGET=
QUAST_LABEL=

if [ -e 1_genome/genome.fa ]; then
    QUAST_TARGET+=" -R 1_genome/genome.fa "
fi

{% for u in unitiggers -%}
if [ -e 7_merge_unitigs_{{ u }}/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_unitigs_{{ u }}/anchor.merge.fasta "
    QUAST_LABEL+="merge_{{ u }},"
fi
{# Keep a blank line #}
{% if opt.merge == "1" and opt.se == "0" -%}
if [ -e 7_merge_mr_unitigs_{{ u }}/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_mr_unitigs_{{ u }}/anchor.merge.fasta "
    QUAST_LABEL+="merge_mr_{{ u }},"
fi
{# Keep a blank line #}
{% endif -%}
{% endfor -%}
if [ -e 7_merge_anchors/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_anchors/anchor.merge.fasta "
    QUAST_LABEL+="merge_anchors,"
fi

if [ -e 8_spades/contigs.fasta ]; then
    QUAST_TARGET+=" 8_spades/contigs.fasta "
    QUAST_LABEL+="spades,"
fi
if [ -e 8_mr_spades/contigs.fasta ]; then
    QUAST_TARGET+=" 8_mr_spades/contigs.fasta "
    QUAST_LABEL+="mr_spades,"
fi

if [ -e 8_megahit/final.contigs.fa ]; then
    QUAST_TARGET+=" 8_megahit/final.contigs.fa "
    QUAST_LABEL+="megahit,"
fi
if [ -e 8_mr_megahit/final.contigs.fa ]; then
    QUAST_TARGET+=" 8_mr_megahit/final.contigs.fa "
    QUAST_LABEL+="mr_megahit,"
fi

if [ -e 1_genome/paralogs.fa ]; then
    QUAST_TARGET+=" 1_genome/paralogs.fa "
    QUAST_LABEL+="paralogs,"
fi

if [ -e 1_genome/repetitive/repetitive.fa ]; then
    QUAST_TARGET+=" 1_genome/repetitive/repetitive.fa "
    QUAST_LABEL+="repetitive,"
fi

QUAST_LABEL=$( echo "${QUAST_LABEL}" | sed 's/,$//' )

rm -fr 9_quast
quast.py --threads {{ opt.parallel }} --min-contig 100 \
    ${QUAST_TARGET} \
    --label ${QUAST_LABEL} \
    -o 9_quast
