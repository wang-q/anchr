#!/usr/bin/env bash

# Check external dependencies

#----------------------------#
# common
#----------------------------#
hash parallel 2>/dev/null || {
    echo >&2 "parallel is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install parallel";
    exit 1;
}

hash mlr 2>/dev/null || {
    echo >&2 "miller is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install miller";
    exit 1;
}

hash jq 2>/dev/null || {
    echo >&2 "jq is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install jq";
    exit 1;
}

hash pigz 2>/dev/null || {
    echo >&2 "pigz is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install pigz";
    exit 1;
}

hash faops 2>/dev/null || {
    echo >&2 "faops is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install wang-q/tap/faops";
    exit 1;
}

#----------------------------#
# QC
#----------------------------#
hash fastqc 2>/dev/null || {
    echo >&2 "fastqc is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install fastqc";
    exit 1;
}

hash picard 2>/dev/null || {
    echo >&2 "picard is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install picard-tools";
    exit 1;
}

hash kat 2>/dev/null || {
    echo >&2 "picard is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install brewsci/bio/kat";
    echo >&2 "                       pip3 install tabulate";
    exit 1;
}

#----------------------------#
# trim, merge, and quorum
#----------------------------#
hash bbduk.sh 2>/dev/null || {
    echo >&2 "bbtools is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install brewsci/bio/bbtools";
    exit 1;
}

hash sickle 2>/dev/null || {
    echo >&2 "sickle is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install sickle";
    exit 1;
}

hash tsv-sample 2>/dev/null || {
    echo >&2 "tsv-sample is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install wang-q/tap/tsv-utils";
    exit 1;
}

hash jellyfish 2>/dev/null || {
    echo >&2 "jellyfish is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install brewsci/bio/jellyfish";
    exit 1;
}

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    hash masurca 2>/dev/null || {
        echo >&2 "masurca is required but it's not installed.";
        echo >&2 "Install with homebrew: brew install brewsci/bio/masurca";
        exit 1;
    }
fi

perl -MNumber::Format -e "1" 2>/dev/null || {
    echo >&2 "Number::Format is required but it's not installed.";
    echo >&2 "Install with cpanm: cpanm Number::Format";
    exit 1;
}

#----------------------------#
# unitigs
#----------------------------#
hash bcalm 2>/dev/null || {
    echo >&2 "bcalm is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install brewsci/bio/bcalm";
    exit 1;
}

hash fasta2DB 2>/dev/null || {
    echo >&2 "DAZZ_DB is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install wang-q/tap/dazz_db@20201008";
    exit 1;
}

hash daligner 2>/dev/null || {
    echo >&2 "daligner is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install wang-q/tap/daligner@20201008";
    exit 1;
}

hash dazz 2>/dev/null || {
    echo >&2 "dazz is required but it's not installed.";
    echo >&2 "Install with cpanm: cpanm App::dazz";
    exit 1;
}

#----------------------------#
# anchors
#----------------------------#
hash spanr 2>/dev/null || {
    echo >&2 "spanr is required but it's not installed.";
    echo >&2 "Install with cargo: cargo install intspan";
    exit 1;
}

hash fasops 2>/dev/null || {
    echo >&2 "fasops is required but it's not installed.";
    echo >&2 "Install with cpanm: cpanm App::Fasops";
    exit 1;
}

#----------------------------#
# group anchors
#----------------------------#
#hash dot 2>/dev/null || {
#    echo >&2 "GraphViz is required but it's not installed.";
#    echo >&2 "Install with homebrew: brew install graphviz";
#    exit 1;
#}

hash poa 2>/dev/null || {
    echo >&2 "poa is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install homebrew/science/poa";
    exit 1;
}

#perl -MGraphViz -e "1" 2>/dev/null || {
#    echo >&2 "GraphViz is required but it's not installed.";
#    echo >&2 "Install with cpanm: cpanm GraphViz";
#    exit 1;
#}

perl -MAlignDB::IntSpan -e "1" 2>/dev/null || {
    echo >&2 "AlignDB::IntSpan is required but it's not installed.";
    echo >&2 "Install with cpanm: cpanm AlignDB::IntSpan";
    exit 1;
}

perl -MGraph -e "1" 2>/dev/null || {
    echo >&2 "Graph is required but it's not installed.";
    echo >&2 "Install with cpanm: cpanm Graph";
    exit 1;
}

#----------------------------#
# sort_on_ref.sh
#----------------------------#
hash sparsemem 2>/dev/null || {
    echo >&2 "sparsemem is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install wang-q/tap/sparsemem";
    exit 1;
}

echo OK
