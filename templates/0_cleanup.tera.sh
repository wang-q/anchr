{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 0_cleanup.sh

# Illumina
parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ -e 2_illumina/{}.fq.gz ]; then
        rm 2_illumina/{}.fq.gz;
        touch 2_illumina/{}.fq.gz;
    fi
    " ::: clumpify filteredbytile sample trim filter

# insertSize
rm -f 2_illumina/insertSize/*tadpole.contig.fasta

# quorum
find 2_illumina -type f -name "quorum_mer_db.jf" | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "k_u_hash_0"       | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "*.tmp"            | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "pe.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "se.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "pe.cor.sub.fa"    | parallel --no-run-if-empty -j 1 rm
find 2_illumina -type f -name "pe.cor.log"       | parallel --no-run-if-empty -j 1 rm

# down sampling
rm -fr 4_down_sampling/
find . -type f -path "*4_unitigs/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_tadpole/*" -name "unitigs_K*.fasta"   | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_tadpole/*/anchor*" -name "basecov.txt"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_tadpole/*/anchor*" -name "*.sam"        | parallel --no-run-if-empty -j 1 rm

rm -fr 6_down_sampling
find . -type f -path "*6_unitigs/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_tadpole/*" -name "unitigs_K*.fasta"   | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_tadpole/*/anchor*" -name "basecov.txt"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_tadpole/*/anchor*" -name "*.sam"        | parallel --no-run-if-empty -j 1 rm

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
if [ -e statUnitigsAnchors.md ]; then
    echo;
    cat statUnitigsAnchors.md;
    echo;
fi
if [ -e statTadpoleAnchors.md ]; then
    echo;
    cat statTadpoleAnchors.md;
    echo;
fi
if [ -e statMRUnitigsAnchors.md ]; then
    echo;
    cat statMRUnitigsAnchors.md;
    echo;
fi
if [ -e statMRTadpoleAnchors.md ]; then
    echo;
    cat statMRTadpoleAnchors.md;
    echo;
fi
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
