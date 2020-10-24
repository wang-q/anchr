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

if [ -e 2_insertSize.sh ]; then
    bash 2_insertSize.sh;
fi

if [ -e 2_sgaPreQC.sh ]; then
    bash 2_sgaPreQC.sh;
fi

#----------------------------#
# trim reads
#----------------------------#
if [ -e 2_trim.sh ]; then
    bash 2_trim.sh;
fi

if [ -e 9_statReads.sh ]; then
    bash 9_statReads.sh;
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
# down sampling, k-unitigs and anchors
#----------------------------#
if [ -e 4_downSampling.sh ]; then
    bash 4_downSampling.sh;
fi

if [ -e 4_kunitigs.sh ]; then
    bash 4_kunitigs.sh;
fi
if [ -e 4_anchors.sh ]; then
    bash 4_anchors.sh;
fi
if [ -e 9_statAnchors.sh ]; then
    bash 9_statAnchors.sh 4_kunitigs statKunitigsAnchors.md
fi

#----------------------------#
# down sampling mergereads
#----------------------------#
if [ -e 6_downSampling.sh ]; then
    bash 6_downSampling.sh
fi

if [ -e 6_kunitigs.sh ]; then
    bash 6_kunitigs.sh;
fi
if [ -e 6_anchors.sh ]; then
    bash 6_anchors.sh;
fi
if [ -e 9_statMRAnchors.sh ]; then
    bash 9_statMRAnchors.sh 6_kunitigs statMRKunitigsAnchors.md
fi

#----------------------------#
# merge anchors
#----------------------------#
if [ -e 7_mergeAnchors.sh ]; then
    bash 7_mergeAnchors.sh 4_kunitigs 7_mergeKunitigsAnchors;
fi

{% if opt.merge -%}
if [ -e 7_mergeAnchors.sh ]; then
    bash 7_mergeAnchors.sh 6_kunitigs 7_mergeMRKunitigsAnchors
fi
{% endif -%}

if [ -e 7_mergeAnchors.sh ]; then
    bash 7_mergeAnchors.sh 7_merge 7_mergeAnchors;
fi

if [ -e 7_mergeAnchors.sh ]; then
    bash 9_statMergeAnchors.sh
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
