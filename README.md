# Anchr

![Publish](https://github.com/wang-q/Anchr/workflows/Publish/badge.svg)
![Build](https://github.com/wang-q/Anchr/workflows/Build/badge.svg)

Anchr - **A**ssembler of **N**-free **CHR**omosomes

## INSTALL

Current release: 0.1.5

```bash
# Via cargo
cargo install --force --path .

# Compiled binary
curl -fsSL https://github.com/wang-q/Anchr/releases/download/v0.1.5/Anchr-x86_64-unknown-linux-musl.tar.gz |
  tar xvz
cp target/x86_64-unknown-linux-musl/release/Anchr ~/bin

```

## SYNOPSIS

```
$ Anchr help
anchr 0.1.1-alpha.0
wang-q <wang-q@outlook.com>
Anchr - Assembler of N-free CHRomosomes

USAGE:
    Anchr [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    dep         Dependencies
    ena         ENA scripts
    help        Prints this message or the help of the given subcommand(s)
    merge       Merge Illumina PE reads with bbtools
    quorum      Run quorum to discard bad reads
    template    Creates Bash scripts
    trim        Trim Illumina PE/SE fastq files

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
