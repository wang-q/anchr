{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 3_freebayes.sh

if [ ! -e 1_genome/genome.fa ]; then
    log_info "1_genome/genome.fa does not exist"
    exit;
fi

if [ ! -e 3_bwa/R.sort.bai ]; then
    log_info "3_bwa/R.sort.bai does not exist"
    exit;
fi

mkdir -p 3_freebayes
cd 3_freebayes

#----------------------------#
# freebayes
#----------------------------#
# https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/non-dip/tutorial.html
if [ ! -e R.raw.vcf ]; then
    freebayes \
        -f ../3_bwa/genome.fa \
        ../3_bwa/R.sort.bam \
        --min-alternate-count 5 \
        --min-alternate-fraction 0.01 \
        > R.raw.vcf
fi

#----------------------------#
# Filter
#----------------------------#

# QUAL > 1: removes really bad sites
# QUAL / AO > 10: additional contribution of each obs should be 10 log units (~ Q10 per read)
# SAF > 0 & SAR > 0: reads on both strands
# RPR > 1 & RPL > 1: at least two reads “balanced” to each side of the site
if [ ! -e R.filtered.vcf ]; then
    cat R.raw.vcf |
        vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
        > R.filtered.vcf
fi

log_info Done.

exit 0
