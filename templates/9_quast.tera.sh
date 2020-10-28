{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_quast.sh

QUAST_TARGET=
QUAST_LABEL=

if [ -e 1_genome/genome.fa ]; then
    QUAST_TARGET+=" -R 1_genome/genome.fa "
fi

if [ -e 7_merge_unitigs_anchors/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_unitigs_anchors/anchor.merge.fasta "
    QUAST_LABEL+="merge_unitigs,"
fi

if [ -e 7_merge_mr_unitigs_anchors/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_mr_unitigs_anchors/anchor.merge.fasta "
    QUAST_LABEL+="merge_mr_unitigs,"
fi

if [ -e 7_merge_anchors/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_anchors/anchor.merge.fasta "
    QUAST_LABEL+="merge_anchors,"
fi

if [ -e 7_glue_anchors/contig.fasta ]; then
    QUAST_TARGET+=" 7_glue_anchors/contig.fasta "
    QUAST_LABEL+="glue_anchors,"
fi
if [ -e 7_fill_anchors/contig.fasta ]; then
    QUAST_TARGET+=" 7_fill_anchors/contig.fasta "
    QUAST_LABEL+="fill_anchors,"
fi

if [ -e 8_spades/spades.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_spades/spades.non-contained.fasta "
    QUAST_LABEL+="spades,"
fi
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_mr_spades/spades.non-contained.fasta "
    QUAST_LABEL+="mr_spades,"
fi

if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_megahit/megahit.non-contained.fasta "
    QUAST_LABEL+="megahit,"
fi
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_mr_megahit/megahit.non-contained.fasta "
    QUAST_LABEL+="mr_megahit,"
fi

if [ -e 8_platanus/platanus.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_platanus/platanus.non-contained.fasta "
    QUAST_LABEL+="platanus,"
fi

if [ -e 1_genome/paralogs.fa ]; then
    QUAST_TARGET+=" 1_genome/paralogs.fa "
    QUAST_LABEL+="paralogs,"
fi

QUAST_LABEL=$( echo "${QUAST_LABEL}" | sed 's/,$//' )

rm -fr 9_quast
quast --no-check --threads {{ opt.parallel }} \
    ${QUAST_TARGET} \
    --label ${QUAST_LABEL} \
    -o 9_quast
