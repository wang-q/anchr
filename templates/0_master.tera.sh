{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 0_master.sh
{# Keep a blank line #}
#----------------------------#
# Illumina QC
#----------------------------#
if [ -e 0_script/2_fastqc.sh ]; then
    bash 0_script/2_fastqc.sh
fi

if [ -e 0_script/2_insert_size.sh ]; then
    bash 0_script/2_insert_size.sh
fi

if [ -e 0_script/2_fastk.sh ]; then
    bash 0_script/2_fastk.sh
fi

#----------------------------#
# trim reads
#----------------------------#
if [ -e 0_script/2_trim.sh ]; then
    bash 0_script/2_trim.sh
fi

if [ -e 0_script/9_stat_reads.sh ]; then
    bash 0_script/9_stat_reads.sh
fi

#----------------------------#
# merge reads
#----------------------------#
if [ -e 0_script/2_merge.sh ]; then
    bash 0_script/2_merge.sh
fi

#----------------------------#
# quorum
#----------------------------#
if [ -e 0_script/2_quorum.sh ]; then
    bash 0_script/2_quorum.sh
fi

#----------------------------#
# mapping
#----------------------------#
if [ -e 0_script/3_bwa.sh ]; then
    bash 0_script/3_bwa.sh
fi
if [ -e 0_script/3_gatk.sh ]; then
    bash 0_script/3_gatk.sh
fi

#----------------------------#
# down sampling trimmed reads; build unitigs and anchors
#----------------------------#
if [ -e 0_script/4_down_sampling.sh ]; then
    bash 0_script/4_down_sampling.sh
fi

if [ -e 0_script/4_unitigs_multik.sh ]; then
    bash 0_script/4_unitigs_multik.sh
fi
if [ -e 0_script/4_anchors.sh ]; then
    bash 0_script/4_anchors.sh 4_unitigs_multik
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_anchors.sh 4_unitigs_multik statUnitigsMultik.md
fi
{# Keep a blank line #}
{% if opt.merge == "1" and opt.se == "0" -%}
#----------------------------#
# down sampling merged reads
#----------------------------#
if [ -e 0_script/6_down_sampling.sh ]; then
    bash 0_script/6_down_sampling.sh
fi

if [ -e 0_script/6_unitigs_multik.sh ]; then
    bash 0_script/6_unitigs_multik.sh
fi
if [ -e 0_script/6_anchors.sh ]; then
    bash 0_script/6_anchors.sh 6_unitigs_multik
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_mr_anchors.sh 6_unitigs_multik statMRUnitigsMultik.md
fi
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# merge anchors
#----------------------------#
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 4_unitigs_multik 7_merge_unitigs_multik
fi
{# Keep a blank line #}
{% if opt.merge == "1" and opt.se == "0" -%}
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 6_unitigs_multik 7_merge_mr_unitigs_multik
fi
{% endif -%}
{# Keep a blank line #}
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 7_merge 7_merge_anchors
fi

if [ -e 0_script/9_stat_merge_anchors.sh ]; then
    bash 0_script/9_stat_merge_anchors.sh
fi

#----------------------------#
# spades, megahit
#----------------------------#
if [ -e 0_script/8_spades.sh ]; then
    bash 0_script/8_spades.sh
fi
if [ -e 0_script/8_mr_spades.sh ]; then
    bash 0_script/8_mr_spades.sh
fi
if [ -e 0_script/8_megahit.sh ]; then
    bash 0_script/8_megahit.sh
fi
if [ -e 0_script/8_mr_megahit.sh ]; then
    bash 0_script/8_mr_megahit.sh
fi

if [ -e 0_script/9_stat_other_anchors.sh ]; then
    bash 0_script/9_stat_other_anchors.sh
fi

{# Keep a blank line #}
#----------------------------#
# final stats
#----------------------------#
if [ -e 0_script/9_stat_final.sh ]; then
    bash 0_script/9_stat_final.sh
fi
if [ -e 0_script/9_quast.sh ]; then
    bash 0_script/9_quast.sh
fi
