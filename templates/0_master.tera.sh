{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 0_master.sh

#----------------------------#
# Illumina QC
#----------------------------#
if [ -e 2_fastqc.sh ]; then
    bash 2_fastqc.sh;
fi
if [ -e 2_kmergenie.sh ]; then
    bash 2_kmergenie.sh;
fi

if [ -e 2_insert_size.sh ]; then
    bash 2_insert_size.sh;
fi

if [ -e 2_sga_preqc.sh ]; then
    bash 2_sga_preqc.sh;
fi

#----------------------------#
# trim reads
#----------------------------#
if [ -e 2_trim.sh ]; then
    bash 2_trim.sh;
fi

if [ -e 9_stat_reads.sh ]; then
    bash 9_stat_reads.sh;
fi

#----------------------------#
# merge reads
#----------------------------#
if [ -e 2_merge.sh ]; then
    bash 2_merge.sh;
fi

#----------------------------#
# quorum
#----------------------------#
if [ -e 2_quorum.sh ]; then
    bash 2_quorum.sh;
fi

#----------------------------#
# down sampling trimmed reads; build unitigs and anchors
#----------------------------#
if [ -e 4_down_sampling.sh ]; then
    bash 4_down_sampling.sh;
fi

if [ -e 4_unitigs.sh ]; then
    bash 4_unitigs.sh;
fi
if [ -e 4_anchors.sh ]; then
    bash 4_anchors.sh;
fi
if [ -e 9_stat_anchors.sh ]; then
    bash 9_stat_anchors.sh 4_unitigs statUnitigsAnchors.md
fi

#----------------------------#
# down sampling merged reads
#----------------------------#
if [ -e 6_down_sampling.sh ]; then
    bash 6_down_sampling.sh
fi

if [ -e 6_unitigs.sh ]; then
    bash 6_unitigs.sh;
fi
if [ -e 6_anchors.sh ]; then
    bash 6_anchors.sh;
fi
if [ -e 9_stat_mr_anchors.sh ]; then
    bash 9_stat_mr_anchors.sh 6_unitigs statMRUnitigsAnchors.md
fi

#----------------------------#
# merge anchors
#----------------------------#
if [ -e 7_merge_anchors.sh ]; then
    bash 7_merge_anchors.sh 4_unitigs 7_merge_unitigs_anchors;
fi

{% if opt.merge == "1" -%}
if [ -e 7_merge_anchors.sh ]; then
    bash 7_merge_anchors.sh 6_unitigs 7_merge_mr_unitigs_anchors
fi
{% endif -%}

if [ -e 7_merge_anchors.sh ]; then
    bash 7_merge_anchors.sh 7_merge 7_merge_anchors;
fi

if [ -e 9_stat_merge_anchors.sh ]; then
    bash 9_stat_merge_anchors.sh
fi

#----------------------------#
# spades, megahit and platanus
#----------------------------#
if [ -e 8_spades.sh ]; then
    bash 8_spades.sh;
fi
if [ -e 8_spades_MR.sh ]; then
    bash 8_spades_MR.sh;
fi
if [ -e 8_megahit.sh ]; then
    bash 8_megahit.sh;
fi
if [ -e 8_megahit_MR.sh ]; then
    bash 8_megahit_MR.sh;
fi
if [ -e 8_platanus.sh ]; then
    bash 8_platanus.sh;
fi

if [ -e 9_statOtherAnchors.sh ]; then
    bash 9_statOtherAnchors.sh;
fi

#----------------------------#
# final stats
#----------------------------#
if [ -e 9_statFinal.sh ]; then
    bash 9_statFinal.sh;
fi
if [ -e 9_quast.sh ]; then
    bash 9_quast.sh;
fi
