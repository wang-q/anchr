{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 0_cleanup.sh
{% set unitiggers = opt.unitigger | split(pat=" ") -%}
{# Keep a blank line #}
# Illumina
parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ -e 2_illumina/{}.fq.gz ]; then
        rm 2_illumina/{}.fq.gz;
        touch 2_illumina/{}.fq.gz;
    fi
    " ::: clumpify filteredbytile sample trim filter

# insertSize
rm -f 2_illumina/insert_size/*tadpole.contig.fa*

# bwa
find 3_bwa -type f -name "genome.fa*"        | parallel --no-run-if-empty -j 1 rm
find 3_bwa -type f -name "*mate.ba[mi]"      | parallel --no-run-if-empty -j 1 rm
find 3_bwa -type f -name "*.per-base.bed.gz" | parallel --no-run-if-empty -j 1 rm

# quorum
find 2_illumina -type f -name "quorum_mer_db.jf" | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "k_u_hash_0"       | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "*.tmp"            | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "pe.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "se.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "pe.cor.sub.fa"    | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "pe.cor.log"       | parallel --no-run-if-empty -j 1 rm

# down sampling
{% for u in unitiggers -%}
find . -type f -path "*4_unitigs_{{ u }}/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_{{ u }}/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_{{ u }}/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm

find . -type f -path "*6_unitigs_{{ u }}/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_{{ u }}/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_{{ u }}/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm
{% endfor -%}
{# Keep a blank line #}
# tempdir
find . -type d -name "\?" | xargs rm -fr

# anchorLong and anchorFill
find . -type d -name "group"         -path "*7_anchor*" | parallel --no-run-if-empty -j 1 rm -fr
find . -type f -name "long.fasta"    -path "*7_anchor*" | parallel --no-run-if-empty -j 1 rm
find . -type f -name ".anchorLong.*" -path "*7_anchor*" | parallel --no-run-if-empty -j 1 rm

# spades
find . -type d -path "*8_spades/*" -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr

# platanus
find . -type f -path "*8_platanus/*" -name "[ps]e.fa" | parallel --no-run-if-empty -j 1 rm

# quast
find . -type d -name "nucmer_output" | parallel --no-run-if-empty -j 1 rm -fr
find . -type f -path "*contigs_reports/*" -name "*.stdout*" -or -name "*.stderr*" | parallel --no-run-if-empty -j 1 rm

# LSF outputs and dumps
find . -type f -name "output.*" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "core.*"   | parallel --no-run-if-empty -j 1 rm

# cat all .md
if [ -e statInsertSize.md ]; then
    echo;
    cat statInsertSize.md;
    echo;
fi
if [ -e statKAT.md ]; then
    echo;
    cat statKAT.md;
    echo;
fi
if [ -e statSgaStats.md ]; then
    echo;
    cat statSgaStats.md;
    echo;
fi
if [ -e statReads.md ]; then
    echo;
    cat statReads.md;
    echo;
fi
if [ -e statTrimReads.md ]; then
    echo;
    cat statTrimReads.md;
    echo;
fi
if [ -e statMergeReads.md ]; then
    echo;
    cat statMergeReads.md;
    echo;
fi
if [ -e statQuorum.md ]; then
    echo;
    cat statQuorum.md;
    echo;
fi
if [ -e statAnchors.md ]; then
    echo;
    cat statAnchors.md;
    echo;
fi
{% for u in unitiggers -%}
if [ -e statUnitigs{{ u | title }}.md ]; then
    echo;
    cat statUnitigs{{ u | title }}.md
    echo;
fi
if [ -e statMRUnitigs{{ u | title }}.md ]; then
    echo;
    cat statMRUnitigs{{ u | title }}.md;
    echo;
fi
{% endfor -%}
if [ -e statMergeAnchors.md ]; then
    echo;
    cat statMergeAnchors.md;
    echo;
fi
if [ -e statOtherAnchors.md ]; then
    echo;
    cat statOtherAnchors.md;
    echo;
fi
if [ -e statFinal.md ]; then
    echo;
    cat statFinal.md;
    echo;
fi
