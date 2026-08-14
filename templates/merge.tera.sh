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

# Error-correct: overlap
log_info "Error-correct: overlap"
anchr fq ec-overlap \
    temp.fq \
    --no-make-vector --vstrict \
    --ihist {{ opt.prefixm }}.ihist.merge1.txt \
    -o ecco.fq
rm temp.fq; ln -s ecco.fq temp.fq
rm -f clumped.fq

# Error-correct: kmer graph
log_info "Error-correct: kmer"
anchr fq ec-kmer \
    temp.fq \
    --toss-junk --toss-depth 2 --toss-uncorrectable \
    -o ecct.fq
rm temp.fq; ln -s ecct.fq temp.fq
rm -f clumped.fq ecco.fq
{# Keep a blank line #}
log_info "Read extension"
anchr fq extend \
    temp.fq \
    -k 62 --el 20 --er 20 \
-o extended.fq
rm temp.fq; ln -s extended.fq temp.fq
rm -f clumped.fq ecco.fq ecct.fq

log_info "Read merging"
anchr fq merge \
    temp.fq \
    --no-make-vector --strict --extend2 80 --rem \
    --ihist {{ opt.prefixm }}.ihist.merge.txt \
    -o merged.raw.fq \
    --outu unmerged.raw.fq
rm -f extended.fq

log_info "Dedupe merged reads"
anchr fq clump \
    merged.raw.fq \
    -o {{ opt.prefixm }}1.fq \
    --dedupe \
    --parallel {{ opt.parallel }}
rm -f merged.raw.fq

log_info "Quality-trim the unmerged reads"
anchr fq clean \
    unmerged.raw.fq \
    --trim-quality {{ opt.qual }} \
    --minlen {{ opt.len }} \
    -o unmerged.trim.fq
rm -f unmerged.raw.fq

# Separates unmerged reads
anchr fq split \
    unmerged.trim.fq \
    -o {{ opt.prefixu }}1.fq \
    --outfile-2 {{ opt.prefixu }}2.fq \
    --outfile-single {{ opt.prefixu }}s.fq
rm -f unmerged.trim.fq

#----------------------------#
# Done.
#----------------------------#
log_info Done.

exit 0
