{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_statFinal.sh

tempfile=$(mktemp /tmp/stat_final_XXXXXXXX)
trap 'rm -f "$tempfile"' EXIT

printf "%s\t" \
    "Name" "N50" "Sum" "#" |
    sed 's/\t$/\n/' \
    > ${tempfile}

# genome
if [ -e 1_genome/genome.fa ]; then
    printf "%s\t" "Genome" $( stat_format 1_genome/genome.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 1_genome/paralogs.fa ]; then
    printf "%s\t" "Paralogs" $( stat_format 1_genome/paralogs.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 1_genome/repetitive/repetitive.fa ]; then
    printf "%s\t" "repetitive" $( stat_format 1_genome/repetitive/repetitive.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# anchors
for D in 7_merge_anchors; do
    if [ -e ${D}/anchor.merge.fasta ]; then
        printf "%s\t" "${D}.anchors" $( stat_format ${D}/anchor.merge.fasta ) |
            sed 's/\t$/\n/' \
            >> ${tempfile}
    fi
    if [ -e ${D}/others.non-contained.fasta ]; then
        printf "%s\t" "${D}.others" $( stat_format ${D}/others.non-contained.fasta ) |
            sed 's/\t$/\n/' \
            >> ${tempfile}
    fi
done

# extended anchors
if [ -e 7_glue_anchors/contig.fasta ]; then
    printf "%s\t" "glue_anchors" $( stat_format 7_glue_anchors/contig.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 7_fill_anchors/contig.fasta ]; then
    printf "%s\t" "fill_anchors" $( stat_format 7_fill_anchors/contig.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# spades
if [ -e 8_spades/contigs.fasta ]; then
    printf "%s\t" "spades.contig" $( stat_format 8_spades/contigs.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_spades/scaffolds.fasta ]; then
    printf "%s\t" "spades.scaffold" $( stat_format 8_spades/scaffolds.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_spades/spades.non-contained.fasta ]; then
    printf "%s\t" "spades.non-contained" $( stat_format 8_spades/spades.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# mr_spades
if [ -e 8_mr_spades/contigs.fasta ]; then
    printf "%s\t" "mr_spades.contig" $( stat_format 8_mr_spades/contigs.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_mr_spades/scaffolds.fasta ]; then
    printf "%s\t" "mr_spades.scaffold" $( stat_format 8_mr_spades/scaffolds.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    printf "%s\t" "mr_spades.non-contained" $( stat_format 8_mr_spades/spades.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# megahit
if [ -e 8_megahit/final.contigs.fa ]; then
    printf "%s\t" "megahit.contig" $( stat_format 8_megahit/final.contigs.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    printf "%s\t" "megahit.non-contained" $( stat_format 8_megahit/megahit.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# mr_megahit
if [ -e 8_mr_megahit/final.contigs.fa ]; then
    printf "%s\t" "mr_megahit.contig" $( stat_format 8_mr_megahit/final.contigs.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    printf "%s\t" "mr_megahit.non-contained" $( stat_format 8_mr_megahit/megahit.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

rgr md ${tempfile} --right 2-4 -o statFinal.md
echo -e "Table: statFinal\n" >> statFinal.md

cat statFinal.md
mv statFinal.md ${BASH_DIR}/../9_markdown
