{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 0_master.sh
{% set unitiggers = opt.unitigger | split(pat=" ") -%}
{# Keep a blank line #}
#----------------------------#
# Illumina QC
#----------------------------#
if [ -e 2_fastqc.sh ]; then
    bash 2_fastqc.sh;
fi

if [ -e 2_insert_size.sh ]; then
    bash 2_insert_size.sh;
fi

if [ -e 2_kat.sh ]; then
    bash 2_kat.sh;
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
# mapping
#----------------------------#
if [ -e 3_bwa.sh ]; then
    bash 3_bwa.sh;
fi
if [ -e 3_gatk.sh ]; then
    bash 3_gatk.sh;
fi

#----------------------------#
# down sampling trimmed reads; build unitigs and anchors
#----------------------------#
if [ -e 4_down_sampling.sh ]; then
    bash 4_down_sampling.sh;
fi

{% for u in unitiggers -%}
if [ -e 4_unitigs_{{ u }}.sh ]; then
    bash 4_unitigs_{{ u }}.sh;
fi
if [ -e 4_anchors.sh ]; then
    bash 4_anchors.sh 4_unitigs_{{ u }};
fi
if [ -e 9_stat_anchors.sh ]; then
    bash 9_stat_anchors.sh 4_unitigs_{{ u }} statUnitigs{{ u | title }}.md
fi
{% endfor -%}
{# Keep a blank line #}
{% if opt.merge == "1" and opt.se == "0" -%}
#----------------------------#
# down sampling merged reads
#----------------------------#
if [ -e 6_down_sampling.sh ]; then
    bash 6_down_sampling.sh
fi

{% for u in unitiggers -%}
if [ -e 6_unitigs_{{ u }}.sh ]; then
    bash 6_unitigs_{{ u }}.sh;
fi
if [ -e 6_anchors.sh ]; then
    bash 6_anchors.sh 6_unitigs_{{ u }};
fi
if [ -e 9_stat_anchors.sh ]; then
    bash 9_stat_mr_anchors.sh 6_unitigs_{{ u }} statMRUnitigs{{ u | title }}.md
fi
{% endfor -%}
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# merge anchors
#----------------------------#
{% for u in unitiggers -%}
if [ -e 7_merge_anchors.sh ]; then
    bash 7_merge_anchors.sh 4_unitigs_{{ u }} 7_merge_unitigs_{{ u }};
fi
{% endfor -%}
{# Keep a blank line #}
{% if opt.merge == "1" and opt.se == "0" -%}
{% for u in unitiggers -%}
if [ -e 7_merge_anchors.sh ]; then
    bash 7_merge_anchors.sh 6_unitigs_{{ u }} 7_merge_mr_unitigs_{{ u }}
fi
{% endfor -%}
{% endif -%}
{# Keep a blank line #}
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
if [ -e 8_mr_spades.sh ]; then
    bash 8_mr_spades.sh;
fi
if [ -e 8_megahit.sh ]; then
    bash 8_megahit.sh;
fi
if [ -e 8_mr_megahit.sh ]; then
    bash 8_mr_megahit.sh;
fi
if [ -e 8_platanus.sh ]; then
    bash 8_platanus.sh;
fi

if [ -e 9_stat_other_anchors.sh ]; then
    bash 9_stat_other_anchors.sh;
fi

#----------------------------#
# extend anchors
#----------------------------#
{% if opt.extend == "1" -%}
rm -fr 7_extend_anchors
mkdir -p 7_extend_anchors
cat \
    8_spades/spades.non-contained.fasta \
    8_megahit/megahit.non-contained.fasta \
    8_platanus/platanus.non-contained.fasta \
{% if opt.merge == "1" and opt.se == "0" -%}
    8_mr_spades/spades.non-contained.fasta \
    8_mr_megahit/megahit.non-contained.fasta \
{% endif -%}
    | faops dazz -a -l 0 stdin stdout \
    | faops filter -a 1000 -l 0 stdin 7_extend_anchors/contigs.2GS.fasta

if [ -e 7_glue_anchors.sh ]; then
    bash 7_glue_anchors.sh 7_merge_anchors/anchor.merge.fasta 7_extend_anchors/contigs.2GS.fasta 3;
fi
if [ -e 7_fill_anchors.sh ]; then
    bash 7_fill_anchors.sh 7_glue_anchors/contig.fasta 7_extend_anchors/contigs.2GS.fasta 3;
fi
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# final stats
#----------------------------#
if [ -e 9_stat_final.sh ]; then
    bash 9_stat_final.sh;
fi
if [ -e 9_quast.sh ]; then
    bash 9_quast.sh;
fi
