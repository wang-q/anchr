{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 3_gatk.sh

if [ ! -e 1_genome/genome.fa ]; then
    log_info "1_genome/genome.fa does not exist"
    exit;
fi

if [ ! -e 3_bwa/R.sort.bai ]; then
    log_info "3_bwa/R.sort.bai does not exist"
    exit;
fi

mkdir -p 3_gatk
cd 3_gatk

{% set parallel2 = opt.parallel | int - 3 -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}

#----------------------------#
# Mutect2
#----------------------------#
#https://github.com/gatk-workflows/gatk4-mitochondria-pipeline/blob/master/tasks/align-and-call.wdl
gatk --java-options "-Xmx{{ opt.xmx }}" \
    Mutect2 \
    --native-pair-hmm-threads {{ parallel2 }} \
    -R ../3_bwa/genome.fa \
    -I ../3_bwa/R.sort.bam \
    --max-reads-per-alignment-start 100 \
    --max-mnp-distance 0 \
    --annotation StrandBiasBySample \
    --mitochondria-mode \
    -O R.raw.vcf

#----------------------------#
# Filter
#----------------------------#
gatk --java-options "-Xmx{{ opt.xmx }}" \
    FilterMutectCalls \
    -R ../3_bwa/genome.fa \
    -V R.raw.vcf \
    --max-alt-allele-count 4 \
    --mitochondria-mode \
    -O R.filtered.vcf

log_info Done.

exit 0
