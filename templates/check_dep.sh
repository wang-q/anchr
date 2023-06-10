#!/usr/bin/env bash

# Check external dependencies

#----------------------------#
# common
#----------------------------#
hash parallel 2>/dev/null || {
    echo >&2 "parallel is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install parallel"
}

hash mlr 2>/dev/null || {
    echo >&2 "miller is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install miller"
}

hash jq 2>/dev/null || {
    echo >&2 "jq is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install jq"
}

hash pigz 2>/dev/null || {
    echo >&2 "pigz is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install pigz"
}

hash faops 2>/dev/null || {
    echo >&2 "faops is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install wang-q/tap/faops"
}

#----------------------------#
# QC
#----------------------------#
hash fastqc 2>/dev/null || {
    echo >&2 "fastqc is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install fastqc"
}

hash picard 2>/dev/null || {
    echo >&2 "picard is required but it's not installed."
    echo >&2 "    Install with homebrew: brew install picard-tools"
}

hash kat 2>/dev/null || {
    echo >&2 "KAT is optional but it's not installed."
    echo >&2 "    Install with homebrew: brew install brewsci/bio/kat"
    echo >&2 "                           pip3 install tabulate"
}

#----------------------------#
# trim, merge, and quorum
#----------------------------#
hash bbduk.sh 2>/dev/null || {
    echo >&2 "bbtools is required but it's not installed."
    echo >&2 "  Install with homebrew: brew install wang-q/tap/bbtools@37.77"
}

hash sickle 2>/dev/null || {
    echo >&2 "sickle is required but it's not installed."
    echo >&2 "  Install with homebrew: brew install sickle"
}

hash tsv-sample 2>/dev/null || {
    echo >&2 "tsv-sample is required but it's not installed."
    echo >&2 "Install with homebrew: brew install wang-q/tap/tsv-utils"
}

hash jellyfish 2>/dev/null || {
    echo >&2 "jellyfish is required but it's not installed."
    echo >&2 "Install with homebrew: brew install jellyfish"
}

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    hash quorum 2>/dev/null || {
        echo >&2 "quorum is required but it's not installed."
        echo >&2 "Install with homebrew: brew install wang-q/tap/quorum@1.1.2"
    }
    hash masurca 2>/dev/null || {
        echo >&2 "masurca is optional but it's not installed."
        echo >&2 "Install with homebrew: brew install brewsci/bio/masurca"
    }
fi

perl -MNumber::Format -e "1" 2>/dev/null || {
    echo >&2 "Number::Format is required but it's not installed."
    echo >&2 "Install with cpanm: cpanm Number::Format"
}

#----------------------------#
# mapping
#----------------------------#
hash bwa 2>/dev/null || {
    echo >&2 "bwa is required but it's not installed."
    echo >&2 "Install with homebrew: brew install bwa"
}

hash samtools 2>/dev/null || {
    echo >&2 "samtools is required but it's not installed."
    echo >&2 "Install with homebrew: brew install samtools"
}

hash gatk 2>/dev/null || {
    echo >&2 "gatk is required but it's not installed."
    echo >&2 "Install with homebrew: brew install brewsci/bio/gatk"
}

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    hash mosdepth 2>/dev/null || {
        echo >&2 "mosdepth is required but it's not installed."
        echo >&2 "Install with homebrew: brew install brewsci/bio/mosdepth"
    }
fi

#----------------------------#
# unitigs
#----------------------------#
hash bcalm 2>/dev/null || {
    echo >&2 "bcalm is required but it's not installed."
    echo >&2 "Install with homebrew: brew install brewsci/bio/bcalm"
}

hash fasta2DB 2>/dev/null || {
    echo >&2 "DAZZ_DB is required but it's not installed."
    echo >&2 "Install with homebrew: brew install --HEAD wang-q/tap/dazz_db"
}

hash daligner 2>/dev/null || {
    echo >&2 "daligner is required but it's not installed."
    echo >&2 "Install with homebrew: brew install --HEAD wang-q/tap/daligner"
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
    echo >&2 "Install with homebrew: brew install wang-q/tap/intspan"
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
    echo >&2 "Install with homebrew: brew install homebrew/science/poa"
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
