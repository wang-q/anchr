{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
rm -f temp.fq.gz;

#----------------------------#
# Pipeline
#----------------------------#
# from bbmap/bbmap/pipelines/assemblyPipeline.sh

# Reorder reads for speed of subsequent phases
log_info "clumpify"
clumpify.sh \
    in={{ args.0 }} \
{% if args.1 -%}
    in2={{ args.1 }} \
{% endif -%}
    out=clumped.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
    dedupe dupesubs=0 \
    overwrite
{% if args.2 -%}
log_info "clumpify SE"
clumpify.sh \
    in={{ args.2 }} \
    out=clumpeds.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
    dedupe dupesubs=0 \
    overwrite
cat clumpeds.fq.gz >> clumped.fq.gz
rm clumpeds.fq.gz
{% endif -%}
rm -f temp.fq.gz; ln -s clumped.fq.gz temp.fq.gz

{% set ecphases = opt.ecphase | split(pat=" ") -%}
{% for ecphase in ecphases -%}
log_info Error-correct phase {{ ecphase }}

{% if ecphase == "1" -%}
# Error-correct phase 1
# error-correct via overlap
bbmerge.sh \
    in=temp.fq.gz out=ecco.fq.gz \
    ihist={{ opt.prefixm }}.ihist.merge1.txt \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
{% if opt.prefilter -%}
    prefilter={{ opt.prefilter }} \
{% endif -%}
    ecco mix vstrict overwrite
rm temp.fq.gz; ln -s ecco.fq.gz temp.fq.gz
{% endif -%}

{% if ecphase == "2" -%}
# Error-correct phase 2
clumpify.sh \
    in=temp.fq.gz out=eccc.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
    passes=4 ecc unpair repair overwrite
rm temp.fq.gz; ln -s eccc.fq.gz temp.fq.gz
{% endif -%}

{% if ecphase == "3" -%}
# Error-correct phase 3
# Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
# For large genomes, tadpole and bbmerge (during the "Merge" phase) may need the flag
# "prefilter=1" or "prefilter=2" to avoid running out of memory.
# "prefilter" makes these take twice as long though so don't use it if you have enough memory.
tadpole.sh \
    in=temp.fq.gz out=ecct.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
{% if opt.prefilter -%}
    prefilter={{ opt.prefilter }} \
{% endif -%}
    ecc tossjunk tossdepth=2 tossuncorrectable overwrite
rm temp.fq.gz; ln -s ecct.fq.gz temp.fq.gz
{% endif -%}

{% endfor -%}
{# Keep a blank line #}
log_info "Read extension"
tadpole.sh \
    in=temp.fq.gz out=extended.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
{% if opt.prefilter -%}
    prefilter={{ opt.prefilter }} \
{% endif -%}
    mode=extend el=20 er=20 k=62 overwrite
rm temp.fq.gz; ln -s extended.fq.gz temp.fq.gz

log_info "Read merging"
bbmerge-auto.sh \
    in=temp.fq.gz out=merged.raw.fq.gz outu=unmerged.raw.fq.gz \
    ihist={{ opt.prefixm }}.ihist.merge.txt \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
{% if opt.prefilter -%}
    prefilter={{ opt.prefilter }} \
{% endif -%}
    strict k=81 extend2=80 rem overwrite

log_info "Dedupe merged reads"
clumpify.sh \
    in=merged.raw.fq.gz \
    out={{ opt.prefixm }}1.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
    dedupe dupesubs=0 \
    overwrite

log_info "Quality-trim the unmerged reads"
bbduk.sh \
    in=unmerged.raw.fq.gz out=unmerged.trim.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
    qtrim=r trimq={{ opt.qual }} minlen={{ opt.len }} overwrite

# Separates unmerged reads
repair.sh \
    in=unmerged.trim.fq.gz \
    out={{ opt.prefixu }}1.fq.gz \
    out2={{ opt.prefixu }}2.fq.gz \
    outs={{ opt.prefixu }}s.fq.gz \
    threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %} \
    repair overwrite

#----------------------------#
# Done.
#----------------------------#
log_info Done.

exit 0
