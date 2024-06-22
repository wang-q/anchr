{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_statFinal.sh

echo -e "Table: statFinal\n" > statFinal.md
printf "| %s | %s | %s | %s |\n" \
    "Name" "N50" "Sum" "#" \
    >> statFinal.md
printf "|:--|--:|--:|--:|\n" >> statFinal.md

# genome
if [ -e 1_genome/genome.fa ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "Genome";   faops n50 -H -S -C 1_genome/genome.fa;) >> statFinal.md
fi
if [ -e 1_genome/paralogs.fa ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "Paralogs"; faops n50 -H -S -C 1_genome/paralogs.fa;) >> statFinal.md
fi
if [ -e 1_genome/repetitives.fa ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "Repetitives"; faops n50 -H -S -C 1_genome/repetitives.fa;) >> statFinal.md
fi

# anchors
for D in 7_merge_anchors; do
    if [ -e ${D}/anchor.merge.fasta ]; then
        printf "| %s | %s | %s | %s |\n" \
            $(echo "${D}.anchors";   faops n50 -H -S -C ${D}/anchor.merge.fasta;) >> statFinal.md
    fi
    if [ -e ${D}/others.non-contained.fasta ]; then
        printf "| %s | %s | %s | %s |\n" \
            $(echo "${D}.others";   faops n50 -H -S -C ${D}/others.non-contained.fasta;) >> statFinal.md
    fi
done

# extended anchors
if [ -e 7_glue_anchors/contig.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "glue_anchors"; faops n50 -H -S -C 7_glue_anchors/contig.fasta;) >> statFinal.md
fi
if [ -e 7_fill_anchors/contig.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "fill_anchors"; faops n50 -H -S -C 7_fill_anchors/contig.fasta;) >> statFinal.md
fi

# spades
if [ -e 8_spades/contigs.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "spades.contig"; faops n50 -H -S -C 8_spades/contigs.fasta;) >> statFinal.md
fi
if [ -e 8_spades/scaffolds.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "spades.scaffold"; faops n50 -H -S -C 8_spades/scaffolds.fasta;) >> statFinal.md
fi
if [ -e 8_spades/spades.non-contained.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "spades.non-contained"; faops n50 -H -S -C 8_spades/spades.non-contained.fasta;) >> statFinal.md
fi

# mr_spades
if [ -e 8_mr_spades/contigs.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "mr_spades.contig"; faops n50 -H -S -C 8_mr_spades/contigs.fasta;) >> statFinal.md
fi
if [ -e 8_mr_spades/scaffolds.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "mr_spades.scaffold"; faops n50 -H -S -C 8_mr_spades/scaffolds.fasta;) >> statFinal.md
fi
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "mr_spades.non-contained"; faops n50 -H -S -C 8_mr_spades/spades.non-contained.fasta;) >> statFinal.md
fi

# megahit
if [ -e 8_megahit/final.contigs.fa ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "megahit.contig"; faops n50 -H -S -C 8_megahit/final.contigs.fa;) >> statFinal.md
fi
if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "megahit.non-contained"; faops n50 -H -S -C 8_megahit/megahit.non-contained.fasta;) >> statFinal.md
fi

# mr_megahit
if [ -e 8_mr_megahit/final.contigs.fa ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "mr_megahit.contig"; faops n50 -H -S -C 8_mr_megahit/final.contigs.fa;) >> statFinal.md
fi
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    printf "| %s | %s | %s | %s |\n" \
        $(echo "mr_megahit.non-contained"; faops n50 -H -S -C 8_mr_megahit/megahit.non-contained.fasta;) >> statFinal.md
fi

cat statFinal.md
