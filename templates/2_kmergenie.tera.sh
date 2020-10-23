{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
mkdir -p 2_illumina/kmergenie
cd 2_illumina/kmergenie

{% set parallel2 = opt.parallel | int / 2 -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ ! -e ${PREFIX}1.dat.pdf ]; then
        parallel --no-run-if-empty --linebuffer -k -j 2 "
            kmergenie -l 21 -k 121 -s 10 -t {{ parallel2 }} --one-pass ../{}.fq.gz -o {}
            " ::: ${PREFIX}1 {% if opt.se == "0" %}${PREFIX}2{% endif %}
    fi
done

log_info Done.

exit 0
