#!/usr/bin/env bash

brew tap brewsci/bio

check_install() {
    if brew list --versions "$1" >/dev/null; then
        echo "$1 already installed"
    else
        brew install "$1"
    fi
}

for package in gatk; do
    check_install brewsci/bio/${package}
done

cbp install openjdk
cbp install parallel jq pigz tsv-utils tva
cbp install fastqc picard
cbp install bwa samtools freebayes mosdepth
cbp install quast spades megahit

cpanm Number::Format

exit 0
