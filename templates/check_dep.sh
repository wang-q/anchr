#!/usr/bin/env bash

# Check external dependencies

#----------------------------#
# common
#----------------------------#
hash anchr 2>/dev/null || {
    echo >&2 "anchr is required but it's not installed (cargo install --path .)."
}

hash pgr 2>/dev/null || {
    echo >&2 "pgr is required but it's not installed (cargo install --path ../pgr)."
}

hash parallel 2>/dev/null || {
    echo >&2 "parallel is required but it's not installed."
}

hash jq 2>/dev/null || {
    echo >&2 "jq is required but it's not installed."
}

hash pigz 2>/dev/null || {
    echo >&2 "pigz is required but it's not installed."
}

hash hnsm 2>/dev/null || {
    echo >&2 "hnsm is required but it's not installed."
}

hash faops 2>/dev/null || {
    echo >&2 "faops is required but it's not installed."
}

hash quast.py 2>/dev/null || {
    echo >&2 "quast is required but it's not installed."
}

#----------------------------#
# QC
#----------------------------#
hash fastqc 2>/dev/null || {
    echo >&2 "fastqc is optional (reference for fq qc) but it's not installed."
}

hash picard 2>/dev/null || {
    echo >&2 "picard is legacy (replaced by asm map + SAM TLEN) but it's not installed."
}

#----------------------------#
# trim, merge, and quorum
#----------------------------#
hash bbduk.sh 2>/dev/null || {
    echo >&2 "bbtools is legacy (replaced by anchr fq) but it's not installed."
}

hash sickle 2>/dev/null || {
    echo >&2 "sickle is legacy (replaced by anchr fq trim-qual) but it's not installed."
}

hash tsv-sample 2>/dev/null || {
    echo >&2 "tsv-sample is required but it's not installed."
}

hash jellyfish 2>/dev/null || {
    echo >&2 "jellyfish is legacy (replaced by pgr kmer hist) but it's not installed."
}

hash quorum 2>/dev/null || {
    echo >&2 "quorum is legacy (replaced by anchr fq s-filter) but it's not installed."
}

hash masurca 2>/dev/null || {
    echo >&2 "masurca is optional but it's not installed."
}

perl -MNumber::Format -e "1" 2>/dev/null || {
    echo >&2 "Number::Format is required but it's not installed."
    echo >&2 "Install with cpanm: cpanm Number::Format"
}

#----------------------------#
# mapping
#----------------------------#
hash bwa 2>/dev/null || {
    echo >&2 "bwa is required but it's not installed."
}

hash samtools 2>/dev/null || {
    echo >&2 "samtools is required but it's not installed."
}

hash gatk 2>/dev/null || {
    echo >&2 "gatk is required but it's not installed."
    echo >&2 "Install with homebrew: brew install brewsci/bio/gatk"
}

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    hash mosdepth 2>/dev/null || {
        echo >&2 "mosdepth is required but it's not installed."
    }
fi

#----------------------------#
# unitigs
#----------------------------#
hash bcalm 2>/dev/null || {
    echo >&2 "bcalm is legacy (replaced by asm unitig) but it's not installed."
}

hash Bifrost 2>/dev/null || {
    echo >&2 "bifrost is legacy (replaced by asm unitig) but it's not installed."
}

hash fasta2DB 2>/dev/null || {
    echo >&2 "DAZZ_DB is required but it's not installed."
}

hash daligner 2>/dev/null || {
    echo >&2 "daligner is required but it's not installed."
}

hash dazz 2>/dev/null || {
    echo >&2 "dazz is required but it's not installed."
    echo >&2 "Install with cpanm: cpanm App::Dazz"
}

#----------------------------#
# anchors
#----------------------------#
hash spanr 2>/dev/null || {
    echo >&2 "spanr is required but it's not installed."
    echo >&2 "Install with cbp: cbp install intspan"
}

hash fasops 2>/dev/null || {
    echo >&2 "fasops is required but it's not installed."
    echo >&2 "Install with cpanm: cpanm App::Fasops"
}

#----------------------------#
# group anchors
#----------------------------#
#hash dot 2>/dev/null || {
#    echo >&2 "GraphViz is required but it's not installed.";
#    echo >&2 "Install with homebrew: brew install graphviz";
#}

hash poa 2>/dev/null || {
    echo >&2 "poa is required but it's not installed."
}

#perl -MGraphViz -e "1" 2>/dev/null || {
#    echo >&2 "GraphViz is required but it's not installed.";
#    echo >&2 "Install with cpanm: cpanm GraphViz";
#}

perl -MAlignDB::IntSpan -e "1" 2>/dev/null || {
    echo >&2 "AlignDB::IntSpan is required but it's not installed."
    echo >&2 "Install with cpanm: cpanm AlignDB::IntSpan"
}

perl -MGraph -e "1" 2>/dev/null || {
    echo >&2 "Graph is required but it's not installed."
    echo >&2 "Install with cpanm: cpanm Graph"
}

#----------------------------#
# sort_on_ref.sh
#----------------------------#
hash sparsemem 2>/dev/null || {
    echo >&2 "sparsemem is required but it's not installed."
    echo >&2 "Install with homebrew: brew install wang-q/tap/sparsemem"
}

echo >&2 OK
