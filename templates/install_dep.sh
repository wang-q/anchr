#!/usr/bin/env bash

brew tap brewsci/bio
brew tap wang-q/tap

check_install() {
    if brew list --versions "$1" >/dev/null; then
        echo "$1 already installed"
    else
        brew install "$1"
    fi
}

for package in bbtools jellyfish; do
    check_install ${package}
done

for package in gatk; do
    check_install brewsci/bio/${package}
done

cbp install openjdk jq parallel pigz
cbp install fastqc sickle bwa samtools picard
cbp install spoa
cbp install bcalm bifrost
cbp install tsv-utils faops intspan
cbp install dazzdb daligner

# linux
cbp install mosdepth

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    check_install brewsci/bio/masurca
    check_install wang-q/tap/quorum@1.1.2
fi

exit 0
