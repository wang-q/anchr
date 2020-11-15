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
# As we're going to precess reads from different sources, don't dedupe here.
# 1. dedupe, Remove duplicate reads.
# 2. optical, mark or remove optical duplicates only. Normal Illumina names needed.
log_info "clumpify"
if [ ! -e clumpify.fq.gz ]; then
    clumpify.sh \
        in={{ args.0 }} \
{% if args.1 -%}
        in2={{ args.1 }} \
{% endif -%}
        out=clumpify.fq.gz \
{% if opt.dedupe == "1" -%}
        dedupe dupesubs=0 \
{% endif -%}
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
fi
rm -f temp.fq.gz; ln -s clumpify.fq.gz temp.fq.gz

{% if opt.tile == "1" -%}
# Remove low-quality reads by positions in flowcell
log_info "filteredbytile"
if [ ! -e filteredbytile.fq.gz ]; then
    filterbytile.sh \
        in=temp.fq.gz \
        out=filteredbytile.fq.gz \
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
fi
rm temp.fq.gz; ln -s filteredbytile.fq.gz temp.fq.gz
{% endif -%}
{# Keep a blank line #}
{% if opt.cutoff != "0" -%}
# Remove reads without high depth kmer
log_info "kmer cutoff with bbnorm.sh"
if [ ! -e highpass.fq.gz ]; then
    bbnorm.sh \
        in=temp.fq.gz \
        out=highpass.fq.gz \
        passes=1 bits=16 min={{ opt.cutoff }} target=9999999 \
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
fi
rm temp.fq.gz; ln -s highpass.fq.gz temp.fq.gz
{% endif -%}
{# Keep a blank line #}
{% if opt.sample != "0" -%}
# Down sampling reads. 300x is fine
log_info "sample with reformat.sh"
if [ ! -e sample.fq.gz ]; then
    reformat.sh \
        in=temp.fq.gz \
        out=sample.fq.gz \
        samplebasestarget={{ opt.sample }} \
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
fi
rm temp.fq.gz; ln -s sample.fq.gz temp.fq.gz
{% endif -%}
{# Keep a blank line #}
# Trim 5' adapters and discard reads with Ns
# Use bbduk.sh to quality and length trim the Illumina reads and remove adapter sequences
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
log_info "trim with bbduk.sh"
if [ ! -e trim.fq.gz ]; then
    bbduk.sh \
        in=temp.fq.gz \
        out=trim.fq.gz \
        ref={{ opt.adapter }} \
        maxns=0 ktrim=r k={{ opt.trimk }} mink=11 hdist=1 tbo tpe \
        minlen={% set lens = opt.len | split(pat=" ") %}{{ lens.0}} qtrim=r trimq={{ opt.trimq }} ftm=5 \
        stats={{ opt.prefix }}.trim.stats.txt overwrite \
        tossbrokenreads=t \
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
fi
rm temp.fq.gz; ln -s trim.fq.gz temp.fq.gz

# Remove synthetic artifacts, spike-ins and 3' adapters by kmer-matching.
log_info "filter with bbduk.sh"
if [ ! -e filter.fq.gz ]; then
    bbduk.sh \
        in=temp.fq.gz \
        out=filter.fq.gz \
        ref={% set fs = opt.filter | split(pat=" ") %}{% for filter in fs %}{% if filter == "adapter" %}{{ opt.adapter }},{% endif %}{% if filter == "artifact" %}{{ opt.artifact }},{% endif %}{% endfor %} \
        k={{ opt.matchk }} cardinality \
        stats={{ opt.prefix }}.filter.stats.txt overwrite \
        tossbrokenreads=t \
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
fi
rm temp.fq.gz; ln -s filter.fq.gz temp.fq.gz

log_info "kmer histogram and peaks"
if [ ! -e peaks.final.txt ]; then
    kmercountexact.sh \
        in=temp.fq.gz \
        khist={{ opt.prefix }}.khist.txt peaks={{ opt.prefix }}.peaks.txt k={{ opt.cutk }} \
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
fi

# Revert to normal pair-end fastq files
log_info "re-pair with repair.sh"
if [ ! -e {{ opt.prefix }}1.trim.fq.gz ]; then
{% if args.1 -%}
    repair.sh \
        in=temp.fq.gz \
        out={{ opt.prefix }}1.fq.gz \
        out2={{ opt.prefix }}2.fq.gz \
        outs={{ opt.prefix }}s.fq.gz \
        repair \
        threads={{ opt.parallel }}{% if opt.xmx != "0" %} -Xmx{{ opt.xmx }}{% endif %}
{% else -%}
    cp -L temp.fq.gz {{ opt.prefix }}1.fq.gz
{% endif -%}
fi

#----------------------------#
# Sickle
#----------------------------#
log_info "sickle ::: Qual {{ opt.qual }} ::: Len {{ opt.len }}"
parallel --no-run-if-empty --linebuffer -k -j 2 "
    mkdir -p Q{1}L{2}
    cd Q{1}L{2}

    printf '==> Qual-Len: %s\n'  Q{1}L{2}
    if [ -e {{ opt.prefix }}1.fq.gz ]; then
        echo '    {{ opt.prefix }}1.fq.gz already presents'
        exit;
    fi

{% if args.1 -%}
    sickle pe \
        -t sanger \
        -q {1} \
        -l {2} \
        -f ../{{ opt.prefix }}1.fq.gz \
        -r ../{{ opt.prefix }}2.fq.gz \
        -o {{ opt.prefix }}1.fq \
        -p {{ opt.prefix }}2.fq \
        -s {{ opt.prefix }}s.fq
    sickle se \
        -t sanger \
        -q {1} \
        -l {2} \
        -f ../{{ opt.prefix }}s.fq.gz \
        -o {{ opt.prefix }}s.temp.fq
    cat {{ opt.prefix }}s.temp.fq >> {{ opt.prefix }}s.fq
    rm {{ opt.prefix }}s.temp.fq
{% else -%}
    sickle se \
        -t sanger \
        -q {1} \
        -l {2} \
        -f ../{{ opt.prefix }}1.fq.gz \
        -o {{ opt.prefix }}1.fq
{% endif -%}

    pigz *.fq
    " ::: {{ opt.qual }} ::: {{ opt.len }}

exit 0
