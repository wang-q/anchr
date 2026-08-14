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

hash quast.py 2>/dev/null || {
    echo >&2 "quast is required but it's not installed."
}

hash tsv-sample 2>/dev/null || {
    echo >&2 "tsv-sample is required but it's not installed."
}

perl -MNumber::Format -e "1" 2>/dev/null || {
    echo >&2 "Number::Format is required but it's not installed."
    echo >&2 "Install with cpanm: cpanm Number::Format"
}

#----------------------------#
# QC / insert-size / merge (optional references and retained picard)
#----------------------------#
hash fastqc 2>/dev/null || {
    echo >&2 "fastqc is optional (reference for fq qc) but it's not installed."
}

hash picard 2>/dev/null || {
    echo >&2 "picard is required (2_insert_size / 3_bwa) but it's not installed."
}

#----------------------------#
# variant detection (3_*)
#----------------------------#
hash bwa 2>/dev/null || {
    echo >&2 "bwa is required but it's not installed."
}

hash samtools 2>/dev/null || {
    echo >&2 "samtools is required but it's not installed."
}

hash freebayes 2>/dev/null || {
    echo >&2 "freebayes is required but it's not installed."
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
# optional references
#----------------------------#
hash spades.py 2>/dev/null || {
    echo >&2 "spades is optional (reference assembler) but it's not installed."
}

hash megahit 2>/dev/null || {
    echo >&2 "megahit is optional (reference assembler) but it's not installed."
}

echo >&2 OK
