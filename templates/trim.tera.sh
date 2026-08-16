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
# As we're going to precess reads from different sources, don't dedupe here.
# 1. dedupe, Remove duplicate reads.
# 2. optical, mark or remove optical duplicates only. Normal Illumina names needed.
log_info "clump with anchr fq clump"
if [ ! -e clumpify.fq ]; then
    anchr fq clump \
        {{ args.0 }} \
{% if args.1 -%}
        {{ args.1 }} \
{% endif -%}
        -o clumpify.fq \
{% if opt.dedupe == "1" -%}
        --dedupe \
{% endif -%}
        --parallel {{ opt.parallel }}
fi
rm -f temp.fq; ln -s clumpify.fq temp.fq
printf "%s\t%s\t%s\t%s\n" \
    $(echo clumpify; stat_format_fq clumpify.fq;) >> statTrimReads.tsv

{% if opt.cutoff != "0" -%}
# Remove reads without high depth kmer
log_info "kmer cutoff with anchr fq norm"
if [ ! -e highpass.fq ]; then
    anchr fq norm \
        temp.fq \
        -o highpass.fq \
        --min {{ opt.cutoff }} \
        --parallel {{ opt.parallel }}
fi
rm temp.fq; ln -s highpass.fq temp.fq
printf "%s\t%s\t%s\t%s\n" \
    $(echo highpass; stat_format_fq highpass.fq;) >> statTrimReads.tsv
rm -f clumpify.fq
{% endif -%}
{# Keep a blank line #}
{% if opt.sample != "0" -%}
# Down sampling reads. 300x is fine
log_info "sample with anchr fq sample"
if [ ! -e sample.fq ]; then
    anchr fq sample \
        temp.fq \
        -o sample.fq \
        --bases {{ opt.sample }}
fi
rm temp.fq; ln -s sample.fq temp.fq
printf "%s\t%s\t%s\t%s\n" \
    $(echo sample; stat_format_fq sample.fq;) >> statTrimReads.tsv
rm -f clumpify.fq highpass.fq
{% endif -%}
{# Keep a blank line #}
# Trim 5' adapters and discard reads with Ns
# Use anchr fq clean to quality and length trim the Illumina reads and remove adapter sequences
# 1. ftm = 5, right trim read length to a multiple of 5
# 2. k = 23, Kmer length used for finding contaminants
# 3. ktrim=r, Trim reads to remove bases matching reference kmers to the right
# 4. mink=7, look for shorter kmers at read tips down to 7 bps
# 5. hdist=1, hamming distance for query kmers
# 6. tbo, trim adapters based on where paired reads overlap
# 7. tpe, when kmer right-trimming, trim both reads to the minimum length of either
# 8. qtrim=r, trim read right ends to remove bases with low quality
# 9. trimq=15, regions with average quality below 15 will be trimmed.
# 10. minlen=60, reads shorter than 60 bps after trimming will be discarded.
log_info "trim with anchr fq clean"
if [ ! -e trim.fq ]; then
    anchr fq clean \
        temp.fq \
        --ref {{ opt.adapter }} \
        --k {{ opt.trimk }} --min-k 11 --hamming-distance 1 \
        --trim-quality {{ opt.trimq }} \
        --minlen {% set lens = opt.len | split(pat=" ") %}{{ lens.0}} \
        --max-ns 0 --force-trim-mod 5 \
        --stats {{ opt.prefix }}.trim.stats.txt \
        -o trim.fq
fi
rm temp.fq; ln -s trim.fq temp.fq
printf "%s\t%s\t%s\t%s\n" \
    $(echo trim; stat_format_fq trim.fq;) >> statTrimReads.tsv
rm -f clumpify.fq highpass.fq sample.fq

# Remove synthetic artifacts, spike-ins and 3' adapters by kmer-matching.
log_info "filter with anchr fq filter"
if [ ! -e filter.fq ]; then
    cat {% set fs = opt.filter | split(pat=" ") %}{% for filter in fs %}{% if filter == "adapter" %}{{ opt.adapter }} {% endif %}{% if filter == "artifact" %}{{ opt.artifact }} {% endif %}{% endfor %}> filter.ref.fa
    anchr fq filter \
        temp.fq \
        --ref filter.ref.fa \
        --k {{ opt.matchk }} \
        --stats {{ opt.prefix }}.filter.stats.txt \
        -o filter.fq
fi
rm temp.fq; ln -s filter.fq temp.fq
printf "%s\t%s\t%s\t%s\n" \
    $(echo filter; stat_format_fq filter.fq;) >> statTrimReads.tsv
rm -f trim.fq

log_info "kmer histogram and peaks with pgr kmer hist"
if [ ! -e peaks.final.txt ]; then
    pgr kmer hist \
        temp.fq \
        -k {{ opt.cutk }} \
        --khist-text {{ opt.prefix }}.khist.txt \
        --peaks {{ opt.prefix }}.peaks.txt \
        -o {{ opt.prefix }}.hist
fi

# Revert to normal pair-end fastq files
log_info "re-pair with anchr fq split"
if [ ! -e {{ opt.prefix }}1.trim.fq.gz ]; then
{% if args.1 -%}
    anchr fq split \
        temp.fq \
        -o {{ opt.prefix }}1.fq \
        --outfile-2 {{ opt.prefix }}2.fq \
        --outfile-single {{ opt.prefix }}s.fq
{% else -%}
    cp -L temp.fq {{ opt.prefix }}1.fq
{% endif -%}
fi
printf "%s\t%s\t%s\t%s\n" \
    $(echo {{ opt.prefix }}1; stat_format_fq {{ opt.prefix }}1.fq;) >> statTrimReads.tsv
printf "%s\t%s\t%s\t%s\n" \
    $(echo {{ opt.prefix }}2; stat_format_fq {{ opt.prefix }}2.fq;) >> statTrimReads.tsv
printf "%s\t%s\t%s\t%s\n" \
    $(echo {{ opt.prefix }}s; stat_format_fq {{ opt.prefix }}s.fq;) >> statTrimReads.tsv
rm -f filter.fq

#----------------------------#
# Sickle
#----------------------------#
log_info "sickle ::: Qual {{ opt.qual }} ::: Len {{ opt.len }}"
# Plain bash loops replace the former `parallel -j 2 ... :::` argument grid
# (nested double-quoted templates need \$ / \" escaping that breaks
# silently). Jobs run one at a time in the foreground, each pigz with the
# full opt.parallel thread budget.
for Q in {{ opt.qual }}; do
    for L in {{ opt.len }}; do
        (
    mkdir -p Q${Q}L${L}
    cd Q${Q}L${L}

    printf '==> Qual-Len: %s\n' "Q${Q}L${L}"
    if [ -e {{ opt.prefix }}1.fq ]; then
        echo '    {{ opt.prefix }}1.fq already presents'
        exit;
    fi

{% if args.1 -%}
    anchr fq trim-qual \
        -q ${Q} \
        -l ${L} \
        ../{{ opt.prefix }}1.fq \
        ../{{ opt.prefix }}2.fq \
        -o {{ opt.prefix }}1.fq \
        --outfile-2 {{ opt.prefix }}2.fq \
        --outfile-single {{ opt.prefix }}s.fq
    anchr fq trim-qual \
        -q ${Q} \
        -l ${L} \
        ../{{ opt.prefix }}s.fq \
        -o {{ opt.prefix }}s.temp.fq
    cat {{ opt.prefix }}s.temp.fq >> {{ opt.prefix }}s.fq
    rm {{ opt.prefix }}s.temp.fq
{% else -%}
    anchr fq trim-qual \
        -q ${Q} \
        -l ${L} \
        ../{{ opt.prefix }}1.fq \
        -o {{ opt.prefix }}1.fq
{% endif -%}

    pigz -p {{ opt.parallel }} *.fq
        )
    done
done

exit 0
