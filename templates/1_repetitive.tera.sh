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

    FastK -v -p -k21 ../genome

    cat chr.sizes |
        number-lines |
        parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
            Profex ../genome {1} |
                sed "1,2 d" |
                perl -nl -e '\''/(\d+).+(\d+)/ and printf qq(%s\t%s\n), $1 + 1, $2'\'' |
                tsv-filter --ge 2:2 |
                cut -f 1 |
                sed "s/^/{2}:/" |
                spanr cover stdin |
                spanr span --op fill -n 2 stdin |
                spanr span --op excise -n 100 stdin |
                spanr span --op fill -n 10 stdin
        ' |
        spanr combine stdin \
        > repetitive.json

    Fastrm ../genome

    spanr convert repetitive.json > region.txt
    pgr fa range ../genome.fa -r region.txt |
        pgr fa filter -N -d --min-len 100 stdin \
        > repetitive.fa

    spanr stat chr.sizes repetitive.json |
        tr ',' '\t' \
        > statRepetitive.tsv
fi

cat statRepetitive.tsv |
    rgr md stdin --num \
    > statRepetitive.md

echo -e "\nTable: statRepetitive\n" >> statRepetitive.md

cat statRepetitive.md
mkdir -p ${BASH_DIR}/../9_markdown
mv statRepetitive.md ${BASH_DIR}/../9_markdown
