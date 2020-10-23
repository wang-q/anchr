# Anchr

![Publish](https://github.com/wang-q/Anchr/workflows/Publish/badge.svg)
![Build](https://github.com/wang-q/Anchr/workflows/Build/badge.svg)

Anchr - **A**ssembler of **N**-free **CHR**omosomes

## INSTALL

Current release: 0.1.0

```bash
cargo install --force --path .

```

## SYNOPSIS

```
$ Anchr help

```

## RUNTIME DEPENDENCIES

* Command line tools managed by `Linuxbrew`

```bash
brew install r perl
brew install parallel wget pigz
brew install datamash mlr

brew tap wang-q/tap
brew install wang-q/tap/tsv-utils wang-q/tap/intspan

```

## EXAMPLES

```bash

```

## PACKAGE

```bash
tar cvfz anchr.$(date +"%Y%m%d").tar.gz \
    $(git ls-files | grep -v "results/")

RESULT=LUAD
pandoc results/${RESULT}.md \
    --standalone \
    -t latex \
    --pdf-engine xelatex \
    -N \
    -V fontsize=10pt \
    -V mainfont="Charter" \
    -V monofont="Fira Mono" \
    -V geometry:"top=2cm, bottom=2cm, left=1.5cm, right=1.5cm" \
    -V geometry:a4paper \
    --highlight-style pygments \
    -o ${RESULT}.$(date +"%Y%m%d").pdf

cloc --md $(git ls-files | grep -v "results/" | grep -v ".tsv")

```

## AUTHOR

Qiang Wang <wang-q@outlook.com>

## LICENSE

Copyright by Qiang Wang.

MIT

Written by Qiang Wang <wang-q@outlook.com>, 2020.
