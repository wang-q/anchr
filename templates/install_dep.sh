#!/usr/bin/env bash

brew tap brewsci/bio
brew tap brewsci/science
brew tap wang-q/tap

check_install () {
    if brew list --versions "$1" > /dev/null; then
        echo "$1 already installed"
    else
        brew install "$1";
    fi
}

for package in openjdk jq parallel pigz; do
    check_install ${package}
done

for package in fastqc sickle bwa samtools picard-tools; do
    check_install ${package}
done

for package in jellyfish bcalm; do
    check_install brewsci/bio/${package};
done

# shellcheck disable=SC2043
for package in poa; do
    check_install brewsci/science/${package};
done

for package in tsv-utils faops sparsemem dazz_db@20201008 daligner@20201008 bbtools@37.77 intspan; do
    check_install wang-q/tap/${package};
done

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    check_install brewsci/bio/masurca
    check_install brewsci/bio/mosdepth
fi

exit 0
