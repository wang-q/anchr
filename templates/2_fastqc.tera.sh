{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 2_fastqc.sh

mkdir -p 2_illumina/fastqc
cd 2_illumina/fastqc

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ ! -e ${PREFIX}1_fastqc.html ]; then
        fastqc -t {{ opt.parallel }} \
            ../${PREFIX}1.fq.gz{% if opt.se == "0" %} ../${PREFIX}2.fq.gz{% endif %} \
            -o .
    fi
done

log_info Done.

exit 0
