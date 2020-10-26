{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 0_real_clean.sh

# Illumina
rm -f 2_illumina/Q*

parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ -e 2_illumina/{1}.{2}.fq.gz ]; then
        rm 2_illumina/{1}.{2}.fq.gz;
    fi
    " ::: R1 R2 Rs ::: uniq shuffle sample bbduk clean

rm -fr 2_illumina/trim/
rm -fr 2_illumina/merge/

# down sampling
rm -fr 4_down_sampling
rm -fr 4_unitigs*
rm -fr 4_tadpole*

rm -fr 6_down_sampling
rm -fr 6_unitigs*
rm -fr 6_tadpole*

# mergeAnchors, anchorLong and anchorFill
rm -fr 7_merge*
rm -fr 7_anchor*
rm -fr 7_fillAnchor

# spades, platanus, and megahit
rm -fr 8_spades*
rm -fr 8_platanus*
rm -fr 8_megahit*

# quast
rm -fr 9_quast*

# tempdir
find . -type d -name "\?" | parallel --no-run-if-empty -j 1 rm -fr

# LSF outputs and dumps
find . -type f -name "output.*" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "core.*"   | parallel --no-run-if-empty -j 1 rm

# .md
rm *.md

# bash
rm *.sh
