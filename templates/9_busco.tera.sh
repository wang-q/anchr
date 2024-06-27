{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_busco.sh
{% set unitiggers = opt.unitigger | split(pat=" ") -%}

ARRAY=()

if [ -e 1_genome/genome.fa ]; then
    ARRAY+=('Genome::1_genome/genome.fa')
fi
if [ -e 1_genome/paralogs.fa ]; then
    ARRAY+=('Paralogs::1_genome/paralogs.fa')
fi
if [ -e 1_genome/repetitive.fa ]; then
    ARRAY+=('repetitive::1_genome/repetitive.fa')
fi

{% for u in unitiggers -%}
if [ -e 7_merge_unitigs_{{ u }}/anchor.merge.fasta ]; then
    ARRAY+=('merge_{{ u }}::7_merge_unitigs_{{ u }}/anchor.merge.fasta')
fi
{% endfor -%}
{# Keep a blank line #}
{% for u in unitiggers -%}
if [ -e 7_merge_mr_unitigs_{{ u }}/anchor.merge.fasta ]; then
    ARRAY+=('merge_mr_{{ u }}::7_merge_mr_unitigs_{{ u }}/anchor.merge.fasta')
fi
{% endfor -%}
{# Keep a blank line #}
if [ -e 7_merge_anchors/anchor.merge.fasta ]; then
    ARRAY+=('merge_anchors::7_merge_anchors/anchor.merge.fasta')
fi

if [ -e 7_glue_anchors/contig.fasta ]; then
    ARRAY+=('glue_anchors::7_glue_anchors/contig.fasta')
fi
if [ -e 7_fill_anchors/contig.fasta ]; then
    ARRAY+=('fill_anchors::7_fill_anchors/contig.fasta')
fi

if [ -e 8_spades/spades.non-contained.fasta ]; then
    ARRAY+=('spades::8_spades/spades.non-contained.fasta')
fi
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    ARRAY+=('mr_spades::8_mr_spades/spades.non-contained.fasta')
fi

if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    ARRAY+=('megahit::8_megahit/megahit.non-contained.fasta')
fi
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    ARRAY+=('mr_megahit::8_mr_megahit/megahit.non-contained.fasta')
fi

mkdir -p 9_busco

for item in "${ARRAY[@]}" ; do
    NAME="${item%%::*}"
    FILE="${item##*::}"

    echo "==> ${NAME}"

    cat ${FILE} |
        sed "s/\//_/g" > tmp.fasta

    busco \
        -i tmp.fasta \
        -o ${NAME} \
        --out_path 9_busco \
        --auto-lineage-prok \
        -m genome --cpu {{ opt.parallel }}

    rm tmp.fasta

    # Clean unused directories
    find 9_busco/${NAME} -maxdepth 1 -mindepth 1 -type d |
        grep -v "run_" |
        xargs rm -fr

    find 9_busco/${NAME} -maxdepth 2 -mindepth 2 -type d |
        xargs rm -fr

done

# Save lineages
find 9_busco/ -maxdepth 2 -mindepth 2 -type d |
    parallel -j 1 'echo {/}' |
    grep "run_" |
    sort -u \
    > 9_busco/lineages.txt

for L in $(cat 9_busco/lineages.txt); do
    echo "Table: statBusco ${L}"
    echo
    echo '| NAME | C | S | D | F | M | Total |'
    echo '|:--|--:|--:|--:|--:|--:|--:|'

    for item in "${ARRAY[@]}" ; do
        NAME="${item%%::*}"

        if [ -f 9_busco/${NAME}/${L}/short_summary.txt ]; then
            cat 9_busco/${NAME}/${L}/short_summary.txt |
                NAME=${NAME} perl -n -MJSON -e '
                    if (m/(\d+)\s*Complete BUSCOs/){
                        $C = $1;
                    }
                    if (m/(\d+)\s*Complete and single-copy BUSCOs/){
                        $S = $1;
                    }
                    if (m/(\d+)\s*Complete and duplicated BUSCOs/){
                        $D = $1;
                    }
                    if (m/(\d+)\s*Fragmented BUSCOs/){
                        $F = $1;
                    }
                    if (m/(\d+)\s*Missing BUSCOs/){
                        $M = $1;
                    }
                    if (m/(\d+)\s*Total BUSCO groups searched/){
                        $T = $1;
                    }
                    END{
                        printf "| %s | %d | %d | %d | %d | %d | %d |\n",
                            $ENV{NAME},
                            $C || 0,
                            $S || 0,
                            $D || 0,
                            $F || 0,
                            $M || 0,
                            $T || 0;
                    }
                '
        fi
    done

    echo
    echo

done \
    > statBusco.md

find . -type f -name "busco*log" | xargs rm
