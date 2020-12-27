{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 0_bsub.sh
{% set unitiggers = opt.unitigger | split(pat=" ") -%}
{# Keep a blank line #}
#----------------------------#
# Illumina QC
#----------------------------#
if [ -e 2_fastqc.sh ]; then
    bsub -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-2_fastqc" \
        "bash 2_fastqc.sh"
fi

if [ -e 2_insert_size.sh ]; then
    bsub -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-2_insert_size" \
        "bash 2_insert_size.sh"
fi

if [ -e 2_kat.sh ]; then
    bsub -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-2_kat" \
        "bash 2_kat.sh"
fi

#----------------------------#
# trim reads
#----------------------------#
bsub -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-2_trim" \
    "bash 2_trim.sh"

bsub -w "ended(${BASE_NAME}-2_trim)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-9_stat_reads" \
    "bash 9_stat_reads.sh"

if [ -e 3_bowtie.sh ]; then
    bsub  -w "ended(${BASE_NAME}-2_trim)" \
        -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-3_bowtie" \
        "bash 3_bowtie.sh"
fi

{% if opt.merge == "1" and opt.se == "0" -%}
#----------------------------#
# merge reads
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_trim)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-2_merge" \
    "bash 2_merge.sh"
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# quorum
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_trim)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-2_quorum" \
    "bash 2_quorum.sh"

#----------------------------#
# down sampling trimmed reads; build unitigs and anchors
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_quorum)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-4_down_sampling" \
    "bash 4_down_sampling.sh"

{% for u in unitiggers -%}
bsub -w "ended(${BASE_NAME}-4_down_sampling)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-4_unitigs_{{ u }}" \
    "
    bash 4_unitigs_{{ u }}.sh
    bash 4_anchors.sh 4_unitigs_{{ u }}
    bash 9_stat_anchors.sh 4_unitigs_{{ u }} statUnitigs{{ u | title }}.md
    "

{% endfor -%}
{# Keep a blank line #}
{% if opt.merge == "1" and opt.se == "0" -%}
#----------------------------#
# down sampling merged reads
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_merge)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-6_down_sampling" \
    "bash 6_down_sampling.sh"

{% for u in unitiggers -%}
bsub -w "ended(${BASE_NAME}-6_down_sampling)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-6_unitigs_{{ u }}" \
    "
    bash 6_unitigs_{{ u }}.sh
    bash 6_anchors.sh 6_unitigs_{{ u }}
    bash 9_stat_mr_anchors.sh 6_unitigs_{{ u }} statMRUnitigs{{ u | title }}.md
    "

{% endfor -%}
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# merge anchors
#----------------------------#
{% for u in unitiggers -%}
bsub -w "ended(${BASE_NAME}-4_unitigs_{{ u }})" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-7_merge_anchors_4_unitigs_{{ u }}" \
    "bash 7_merge_anchors.sh 4_unitigs_{{ u }} 7_merge_unitigs_{{ u }}"
{% endfor -%}
{# Keep a blank line #}
{% if opt.merge == "1" and opt.se == "0" -%}
{% for u in unitiggers -%}
bsub -w "ended(${BASE_NAME}-6_unitigs_{{ u }})" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-7_merge_anchors_6_unitigs_{{ u }}" \
    "bash 7_merge_anchors.sh 6_unitigs_{{ u }} 7_merge_mr_unitigs_{{ u }}"
{% endfor -%}
{% endif -%}
{# Keep a blank line #}
bsub -w "ended(${BASE_NAME}-2_quorum) {% for u in unitiggers %}&& ended(${BASE_NAME}-7_merge_anchors_4_unitigs_{{ u }}){% endfor %} {% if opt.merge == "1" and opt.se == "0" %}{% for u in unitiggers %}&& ended(${BASE_NAME}-7_merge_anchors_4_unitigs_{{ u }}){% endfor %}{% endif %}" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-7_merge_anchors" \
    "bash 7_merge_anchors.sh 7_merge 7_merge_anchors"
bsub -w "ended(${BASE_NAME}-7_merge_anchors)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-9_stat_merge_anchors" \
    "bash 9_stat_merge_anchors.sh"

#----------------------------#
# spades, megahit and platanus
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_quorum)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-8_spades" \
    "bash 8_spades.sh"

bsub -w "ended(${BASE_NAME}-2_quorum)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-8_megahit" \
    "bash 8_megahit.sh"

bsub -w "ended(${BASE_NAME}-2_quorum)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-8_platanus" \
    "bash 8_platanus.sh"

{% if opt.merge == "1" and opt.se == "0" -%}
bsub -w "ended(${BASE_NAME}-2_merge)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-8_mr_spades" \
    "bash 8_mr_spades.sh"
bsub -w "ended(${BASE_NAME}-2_merge)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-8_mr_megahit" \
    "bash 8_mr_megahit.sh"
{% endif -%}
{# Keep a blank line #}
bsub -w "ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-8_megahit) && ended(${BASE_NAME}-8_platanus) {% if opt.merge == "1" and opt.se == "0" %}&& ended(${BASE_NAME}-8_mr_spades) && ended(${BASE_NAME}-8_mr_megahit){% endif %}" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-9_stat_other_anchors" \
    "bash 9_stat_other_anchors.sh"

#----------------------------#
# extend anchors
#----------------------------#
{% if opt.extend == "1" -%}
bsub -w "ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-8_megahit) && ended(${BASE_NAME}-8_platanus) {% if opt.merge == "1" and opt.se == "0" %}&& ended(${BASE_NAME}-8_mr_spades)&& ended(${BASE_NAME}-8_mr_megahit){% endif %}" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-contigs_2GS" \
    '
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
    '

bsub -w "ended(${BASE_NAME}-7_merge_anchors) && ended(${BASE_NAME}-contigs_2GS)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-7_glue_anchors" \
    "bash 7_glue_anchors.sh 7_merge_anchors/anchor.merge.fasta 7_extend_anchors/contigs.2GS.fasta 3"

bsub -w "ended(${BASE_NAME}-7_glue_anchors)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-7_fill_anchors" \
    "bash 7_fill_anchors.sh 7_glue_anchors/contig.fasta 7_extend_anchors/contigs.2GS.fasta 3"
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# final stats
#----------------------------#
bsub -w "ended(${BASE_NAME}-7_merge_anchors) && ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-8_platanus) {% if opt.extend == "1" %}&& ended(${BASE_NAME}-7_fill_anchors){% endif %}" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-9_stat_final" \
    "bash 9_stat_final.sh"

bsub -w "ended(${BASE_NAME}-7_merge_anchors) && ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-8_platanus) {% if opt.extend == "1" %}&& ended(${BASE_NAME}-7_fill_anchors){% endif %}" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-9_quast" \
    "bash 9_quast.sh"

bsub -w "ended(${BASE_NAME}-9_stat_final) && ended(${BASE_NAME}-9_quast)" \
    -q {{ opt.queue }} -n {{ opt.parallel }} -J "${BASE_NAME}-0_cleanup" \
    "bash 0_cleanup.sh"
