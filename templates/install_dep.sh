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
cbp install parallel jq pigz hnsm
cbp install fastqc picard
cbp install bbtools sickle tsv-utils
cbp install jellyfish quorum superreads
cbp install bwa samtools freebayes mosdepth # gatk
cbp install bcalm bifrost
cbp install faops intspan
cbp install spoa
cbp install dazzdb daligner

exit 0
