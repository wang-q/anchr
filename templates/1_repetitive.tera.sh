{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 1_repetitive.sh

mkdir -p 1_genome/repetitive
cd 1_genome/repetitive

if [ ! -s repetitive.fa ]; then
    pgr fa size ../genome.fa > chr.sizes

    pgr rept s-kmer ../genome.fa -k 21 --fill-kmer 2 --min-len 100 \
        --fill-fragment 10 \
        -o repetitive.json

    pgr runlist convert repetitive.json \
        > region.txt

    pgr fa range ../genome.fa -r region.txt |
        pgr fa filter -N -d --min-len 100 stdin \
        > repetitive.fa

    pgr runlist stat chr.sizes repetitive.json \
        > statRepetitive.tsv
fi

cat statRepetitive.tsv |
    rgr md stdin --num \
    > statRepetitive.md

echo -e "\nTable: statRepetitive\n" >> statRepetitive.md

cat statRepetitive.md
mkdir -p ${BASH_DIR}/../9_markdown
mv statRepetitive.md ${BASH_DIR}/../9_markdown
