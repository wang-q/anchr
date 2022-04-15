# Anchr

[![Publish](https://github.com/wang-q/anchr/actions/workflows/publish.yml/badge.svg)](https://github.com/wang-q/anchr/actions)
[![Build](https://github.com/wang-q/anchr/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/anchr/actions)
[![Codecov](https://img.shields.io/codecov/c/github/wang-q/anchr/master.svg)](https://codecov.io/github/wang-q/anchr?branch=main)
[![Lines of code](https://tokei.rs/b1/github/wang-q/anchr?category=code)](https://github.com//wang-q/anchr)

Anchr - the **A**ssembler of **N**-free **CHR**omosomes

## INSTALL

Current release: 0.3.17

```shell
# Via cargo
cargo install --force --path .

cargo install --git https://github.com/wang-q/anchr --branch main

# Static binary for Linux
mkdir -p ${HOME}/bin
curl -fsSL $(
    curl -fsSL https://api.github.com/repos/wang-q/anchr/releases/latest |
        jq -r '.assets[] | select(.name == "anchr-x86_64-unknown-linux-musl.tar.gz").browser_download_url'
    ) |
    tar xvz
cp target/x86_64-unknown-linux-musl/release/anchr ${HOME}/bin
rm -fr target

```

## SYNOPSIS

```text
$ anchr help
anchr 0.3.16-alpha.1
wang-q <wang-q@outlook.com>
Anchr - the Assembler of N-free CHRomosomes

USAGE:
    anchr [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    anchors     Select anchors (proper covered regions) from contigs
    dep         Dependencies
    ena         ENA scripts
    help        Prints this message or the help of the given subcommand(s)
    merge       Merge Illumina PE reads with bbtools
    quorum      Run quorum to discard bad reads
    template    Creates Bash scripts
    trim        Trim Illumina PE/SE fastq files
    unitigs     Create unitigs from trimmed/merged reads

```

## RUNTIME DEPENDENCIES

* Command line tools managed by `Linuxbrew`

```bash
brew install perl cpanminus
brew install parallel wget pigz
brew install datamash miller

# cite parallel
# parallel --will-cite

brew tap wang-q/tap
brew install wang-q/tap/tsv-utils wang-q/tap/intspan

# Myer's dazzler wrapper
cpanm --installdeps App::Dazz
cpanm -nq App::Dazz

anchr dep install | bash
anchr dep check | bash

# Optional
# assembly quality assessment. https://github.com/ablab/quast/issues/140
brew install brewsci/bio/quast --HEAD

# Optional
# leading assemblers
brew install spades
spades.py --test
brew install brewsci/bio/megahit
brew install wang-q/tap/platanus

```

## EXAMPLES

### Fetching data

Data soruce: *E. coli* virus Lambda

```shell
# ena
mkdir -p ~/data/ena
cd ~/data/ena

cat << EOF > source.csv
SRX2365802,Lambda,HiSeq 2500
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -j 4 -x 4 -s 2 -c --file-allocation=none -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

# sampling reads as test materials
seqtk sample -s 23 SRR5042715_1.fastq.gz 20000 | pigz > R1.fq.gz
seqtk sample -s 23 SRR5042715_2.fastq.gz 20000 | pigz > R2.fq.gz

```

| name   | srx        | platform | layout | ilength | srr        | spot     | base  |
|:-------|:-----------|:---------|:-------|:--------|:-----------|:---------|:------|
| Lambda | SRX2365802 | ILLUMINA | PAIRED |         | SRR5042715 | 16540237 | 3.33G |

### Individual subcommands

```shell
mkdir -p ~/data/anchr_test
cd ~/data/anchr_test

# Lambda
for F in R1.fq.gz R2.fq.gz; do
    1>&2 echo ${F}
    curl -fsSLO "https://raw.githubusercontent.com/wang-q/anchr/main/tests/Lambda/${F}"
done

# trim
mkdir -p trim
pushd trim

anchr trim \
    ../R1.fq.gz ../R2.fq.gz \
    -q 25 -l 60 \
    -o stdout |
    bash
popd

# merge
mkdir -p merge
pushd merge

anchr merge \
    ../trim/R1.fq.gz ../trim/R2.fq.gz ../trim/Rs.fq.gz \
    --ecphase "1 2 3" \
    --parallel 4 \
    -o stdout |
    bash
popd

# quorum
pushd trim
anchr quorum \
    R1.fq.gz R2.fq.gz \
    -o stdout |
    bash
popd

pushd trim/Q25L60
anchr quorum \
    R1.fq.gz R2.fq.gz Rs.fq.gz \
    -o stdout |
    bash
popd

# unitigs
gzip -dcf trim/pe.cor.fa.gz > trim/pe.cor.fa

mkdir -p superreads
pushd superreads

anchr unitigs \
    ../trim/pe.cor.fa ../trim/env.json \
    --kmer "31 41 51 61 71 81" \
    --parallel 4 \
    -o unitigs.sh
bash unitigs.sh
popd

# unitigs - tadpole
mkdir -p tadpole
pushd tadpole

anchr unitigs \
    ../trim/pe.cor.fa ../trim/env.json \
    -u tadpole \
    --kmer "31 41 51 61 71 81" \
    --parallel 4 \
    -o unitigs.sh
bash unitigs.sh
popd

# unitigs - bcalm
mkdir -p bcalm
pushd bcalm

anchr unitigs \
    ../trim/pe.cor.fa ../trim/env.json \
    -u bcalm \
    --kmer "31 41 51 61 71 81" \
    --parallel 4 \
    -o unitigs.sh
bash unitigs.sh
popd

# anchors
mkdir -p bcalm/anchors
pushd bcalm/anchors

anchr anchors \
    ../unitigs.fasta \
    ../pe.cor.fa \
    --readl 150 \
    --keepedge \
    -p 4 \
    -o anchors.sh
bash anchors.sh
popd

```

### `anchr template`

With a conventional directory structure, `anchr template` creates all scripts from reads QC to
assembly evaluations.

* E. coli
    * [*Escherichia* virus Lambda](results/e_coli.md#escherichia-virus-lambda)
    * [*E.
      coli* str. K-12 substr. MG1655](results/e_coli.md#escherichia-coli-str-k-12-substr-mg1655)
    * [*E.
      coli* str. K-12 substr. DH5alpha](results/e_coli.md#escherichia-coli-str-k-12-substr-dh5alpha)

* FDA-ARGOS bacteria
    * [Ca_jej_jejuni_NCTC_11168_ATCC_700819](results/fda_argos.md#ca_jej_jejuni_nctc_11168_atcc_700819)
    * [Clostridio_dif_630](results/fda_argos.md#clostridio_dif_630)
    * [Co_dip_NCTC_13129](results/fda_argos.md#co_dip_nctc_13129)
    * [Fr_tul_tularensis_SCHU_S4](results/fda_argos.md#fr_tul_tularensis_schu_s4)

* Yeast
    * [*Saccharomyces cerevisiae* S288c](results/yeast.md#saccharomyces-cerevisiae-s288c)

## AUTHOR

Qiang Wang <wang-q@outlook.com>

## LICENSE

Copyright by Qiang Wang.

MIT

Written by Qiang Wang <wang-q@outlook.com>, 2020.
