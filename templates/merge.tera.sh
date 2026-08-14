{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
rm -f temp.fq;

#----------------------------#
# Pipeline
#----------------------------#
# from bbmap/bbmap/pipelines/assemblyPipeline.sh

# Reorder reads for speed of subsequent phases
log_info "clump with anchr fq clump"
anchr fq clump \
    {{ args.0 }} \
{% if args.1 -%}
    {{ args.1 }} \
{% endif -%}
    -o clumped.fq \
    --dedupe \
    --parallel {{ opt.parallel }}
{% if args.2 -%}
log_info "clump SE with anchr fq clump"
anchr fq clump \
    {{ args.2 }} \
    -o clumpeds.fq \
    --dedupe \
    --parallel {{ opt.parallel }}
cat clumpeds.fq >> clumped.fq
rm clumpeds.fq
{% endif -%}
rm -f temp.fq; ln -s clumped.fq temp.fq

{% set ecphases = opt.ecphase | split(pat=" ") -%}
{% for ecphase in ecphases -%}
log_info Error-correct phase {{ ecphase }}

{% if ecphase == "1" -%}
# Error-correct phase 1
# error-correct via overlap
anchr fq ec-overlap \
    temp.fq \
    --no-make-vector --vstrict \
    --ihist {{ opt.prefixm }}.ihist.merge1.txt \
    -o ecco.fq
rm temp.fq; ln -s ecco.fq temp.fq
{% endif -%}

{% if ecphase == "2" -%}
# Error-correct phase 2
# 跳过：clumpify ecc（phase 2）与 phase 3（tadpole ecc）冗余且常卡住，
# anchr 侧不实现（见 design/fq-merge-replace.md §3.3）；temp 链不变。
{% endif -%}

{% if ecphase == "3" -%}
# Error-correct phase 3
# Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
anchr fq ec-kmer \
    temp.fq \
    --toss-junk --toss-depth 2 --toss-uncorrectable \
    -o ecct.fq
rm temp.fq; ln -s ecct.fq temp.fq
{% endif -%}

{% endfor -%}
{# Keep a blank line #}
log_info "Read extension"
anchr fq extend \
    temp.fq \
    -k 62 --el 20 --er 20 \
    -o extended.fq
rm temp.fq; ln -s extended.fq temp.fq

log_info "Read merging"
anchr fq merge \
    temp.fq \
    --no-make-vector --strict --extend2 80 --rem \
    --ihist {{ opt.prefixm }}.ihist.merge.txt \
    -o merged.raw.fq \
    --outu unmerged.raw.fq

log_info "Dedupe merged reads"
anchr fq clump \
    merged.raw.fq \
    -o {{ opt.prefixm }}1.fq \
    --dedupe \
    --parallel {{ opt.parallel }}

log_info "Quality-trim the unmerged reads"
anchr fq clean \
    unmerged.raw.fq \
    --trim-quality {{ opt.qual }} \
    --minlen {{ opt.len }} \
    -o unmerged.trim.fq

# Separates unmerged reads
anchr fq split \
    unmerged.trim.fq \
    -o {{ opt.prefixu }}1.fq \
    --outfile-2 {{ opt.prefixu }}2.fq \
    --outfile-single {{ opt.prefixu }}s.fq

#----------------------------#
# Done.
#----------------------------#
log_info Done.

exit 0
