# Anchr

[![Publish](https://github.com/wang-q/anchr/actions/workflows/publish.yml/badge.svg)](https://github.com/wang-q/anchr/actions)
[![Build](https://github.com/wang-q/anchr/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/anchr/actions)
[![Codecov](https://img.shields.io/codecov/c/github/wang-q/anchr/master.svg)](https://codecov.io/github/wang-q/anchr?branch=main)
[![Lines of code](https://tokei.rs/b1/github/wang-q/anchr?category=code)](https://github.com//wang-q/anchr)

Anchr - the **A**ssembler of **N**-free **CHR**omosomes

## INSTALL

Current release: 0.3.18

```shell
# Via cargo
cargo install --path . --force #--offline

cargo install --git https://github.com/wang-q/anchr --branch main

# test
cargo test -- --test-threads=1

## Static binary for Linux
#mkdir -p ${HOME}/bin
#curl -fsSL $(
#    curl -fsSL https://api.github.com/repos/wang-q/anchr/releases/latest |
#        jq -r '.assets[] | select(.name == "anchr-x86_64-unknown-linux-musl.tar.gz").browser_download_url'
#    ) |
#    tar xvz
#cp target/x86_64-unknown-linux-musl/release/anchr ${HOME}/bin
#rm -fr target

# build under WSL 2
mkdir -p /tmp/cargo
export CARGO_TARGET_DIR=/tmp/cargo
cargo build

# build for CentOS 7
# rustup target add x86_64-unknown-linux-gnu
# pip3 install cargo-zigbuild
cargo zigbuild --target x86_64-unknown-linux-gnu.2.17 --release
ll $CARGO_TARGET_DIR/x86_64-unknown-linux-gnu/release/

```

## SYNOPSIS

```text
$ anchr help
Anchr - the Assembler of N-free CHRomosomes

Usage: anchr [COMMAND]

Commands:
  anchors    Select anchors (proper covered regions) from contigs
  contained  Discard contained unitigs
  covered    Covered regions from .ovlp.tsv files
  dazzname   Rename FASTA records for dazz_db
  dep        Dependencies
  ena        ENA scripts
  merge      Merge overlapped unitigs
  mergeread  Merge Illumina PE reads with bbtools
  orient     Orient overlapped sequences to the same strand
  overlap    Detect overlaps by daligner
  paf2ovlp   Convert minimap .paf to overlaps
  quorum     Run quorum to discard bad reads
  restrict   Restrict overlaps to known pairs
  show2ovlp  Convert LAshow outputs to overlaps
  template   Creates Bash scripts
  trim       Trim Illumina PE/SE fastq files
  unitigs    Create unitigs from trimmed/merged reads
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```

## RUNTIME DEPENDENCIES

* Command line tools managed by `Homebrew`

```bash
brew install perl cpanminus
brew install r
brew install parallel wget pigz
brew install miller prettier

# cite parallel
# parallel --citation
# will cite

brew tap wang-q/tap
brew install wang-q/tap/tsv-utils wang-q/tap/intspan

brew install --HEAD wang-q/tap/dazz_db
brew install --HEAD wang-q/tap/daligner

anchr dep install | bash
anchr dep check | bash

# Optional: fastk
brew install --HEAD wang-q/tap/fastk
brew install --HEAD wang-q/tap/merquryfk
brew install --HEAD wang-q/tap/fastga

brew install binutils
brew link binutils --force
parallel -j 1 -k --line-buffer '
    Rscript -e '\'' if (!requireNamespace("{}", quietly = FALSE)) { install.packages("{}", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN") } '\''
    ' ::: \
        argparse minpack.lm \
        ggplot2 scales viridis

# Optional
# quast - assembly quality assessment
curl -LO https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
tar xvfz quast-5.2.0.tar.gz
cd quast-5.2.0
python3 ./setup.py install

quast.py --test

# Optional: leading assemblers
brew install brewsci/bio/megahit

brew install spades
spades.py --test

```

## EXAMPLES

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

# mergeread
mkdir -p merge
pushd merge

anchr mergeread \
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

### Fetching data

* Data source: *Mycoplasma genitalium* G37

```shell
mkdir -p ~/data/anchr/g37/ena
cd ~/data/anchr/g37/ena

cat << EOF > source.csv
ERX452667,G37,MiSeq
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

rgr md ena_info.tsv --fmt

aria2c -j 4 -x 4 -s 2 -c --file-allocation=none -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

## sampling reads as test materials
#seqtk sample -s 23 SRR5042715_1.fastq.gz 20000 | pigz > R1.fq.gz
#seqtk sample -s 23 SRR5042715_2.fastq.gz 20000 | pigz > R2.fq.gz

```

| name | srx       | platform | layout | ilength | srr       |   spots | bases  |
|------|-----------|----------|--------|--------:|-----------|--------:|--------|
| G37  | ERX452667 | ILLUMINA | PAIRED |     447 | ERR486835 | 680,644 | 97.37M |

### `anchr template`

With a conventional directory structure, `anchr template` creates all scripts from reads QC to
assembly evaluations.

* Model bacteria
    * [*Mycoplasma genitalium* G37](results/model.md#mycoplasma-genitalium-g37)
    * [*E. coli* str. K-12 substr. MG1655](results/model.md#e-coli-str-k-12-substr-mg1655)
    * [*E. coli* str. K-12 substr. DH5alpha](results/model.md#e-coli-str-k-12-substr-dh5alpha)

* FDA-ARGOS bacteria
    * [Ca_jej_jejuni_NCTC_11168_ATCC_700819](results/fda_argos.md#ca_jej_jejuni_nctc_11168_atcc_700819)
    * [Clostridio_dif_630](results/fda_argos.md#clostridio_dif_630)
    * [Co_dip_NCTC_13129](results/fda_argos.md#co_dip_nctc_13129)
    * [Fr_tul_tularensis_SCHU_S4](results/fda_argos.md#fr_tul_tularensis_schu_s4)

* Yeast
    * [*Saccharomyces cerevisiae* S288c](results/yeast.md#saccharomyces-cerevisiae-s288c)

### Overlaps - Standalone

```shell
anchr dazzname tests/ovlpr/1_4.anchor.fasta -o stdout

anchr show2ovlp tests/ovlpr/1_4.show.txt tests/ovlpr/1_4.replace.tsv --orig

anchr paf2ovlp tests/ovlpr/1_4.pac.paf

echo "tests/ovlpr/1_4.anchor.fasta;tests/ovlpr/1_4.pac.fasta" |
    parallel --colsep ";" -j 1 "
        minimap2 -cx asm20 {1} {2} |
            anchr paf2ovlp stdin |
            tsv-sort
        minimap2 -cx asm20 {2} {1} |
            anchr paf2ovlp stdin |
            tsv-sort
    " |
    anchr covered stdin --mean

anchr covered tests/ovlpr/1_4.pac.paf.ovlp.tsv
anchr covered tests/ovlpr/11_2.long.paf --paf
anchr covered tests/ovlpr/1_4.pac.paf.ovlp.tsv --base
anchr covered tests/ovlpr/1_4.pac.paf.ovlp.tsv --mean

anchr restrict tests/ovlpr/1_4.ovlp.tsv tests/ovlpr/1_4.restrict.tsv

```

### Overlaps - Daligner pipelines

```shell
anchr overlap tests/ovlpr/1_4.pac.fasta
anchr overlap tests/ovlpr/1_4.pac.fasta --idt 0.8 --len 2500 --serial

anchr orient tests/ovlpr/1_4.anchor.fasta tests/ovlpr/1_4.pac.fasta
anchr orient tests/ovlpr/1_4.anchor.fasta tests/ovlpr/1_4.pac.fasta -r tests/ovlpr/1_4.2.restrict.tsv

anchr contained tests/ovlpr/contained.fasta

cargo run --bin anchr merge tests/ovlpr/merge.fasta -o test.fasta

```

## AUTHOR

Qiang Wang <wang-q@outlook.com>

## LICENSE

Copyright by Qiang Wang.

MIT

Written by Qiang Wang <wang-q@outlook.com>, 2020.
