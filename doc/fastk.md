# FastK, GENESCOPE.FK, MERQURY.FK, and FastGA

<!-- TOC -->
* [FastK, GENESCOPE.FK, MERQURY.FK, and FastGA](#fastk-genescopefk-merquryfk-and-fastga)
  * [Installation](#installation)
  * [GeneScope](#genescope)
  * [Merqury](#merqury)
  * [Find repetitive regions](#find-repetitive-regions)
  * [FastGA](#fastga)
  * [Paralogs](#paralogs)
<!-- TOC -->

[FastK](https://github.com/thegenemyers/FASTK)

[GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)

[MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)

## Installation

```shell
brew install --HEAD wang-q/tap/fastk
brew install --HEAD wang-q/tap/merquryfk

# GeneScope.FK is bundled

parallel -j 1 -k --line-buffer '
    Rscript -e '\'' if (!requireNamespace("{}", quietly = FALSE)) { install.packages("{}", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN") } '\''
    ' ::: \
        argparse minpack.lm \
        ggplot2 scales viridis

```

## GeneScope

```shell
cd ~/data/anchr/g37/ena/

time FastK -v -t1 -k21 *.fastq.gz -NTable-21
#real    0m1.984s
#user    0m3.318s
#sys     0m0.341s

time FastK -v -t1 -k51 *.fastq.gz -NTable-51
#real    0m2.082s
#user    0m3.555s
#sys     0m0.507s

Histex -G Table-21 | Rscript <(gzip -dcf ~/Scripts/anchr/templates/genescopefk.R.gz) -k 21 -p 1 -o GeneScope-21

Histex -G Table-51 | Rscript <(gzip -dcf ~/Scripts/anchr/templates/genescopefk.R.gz) -k 51 -p 1 -o GeneScope-51

# disk usages
ll |
    grep Table |
    tr -s ' ' '\t' |
    cut -f 5,9 |
    sed 's/\*$//' |
    ( echo -e 'Size\tName' && cat ) |
    rgr md stdin --fmt

Fastrm Table-21 Table-51

# Reports
for K in 21 51; do
    COV=$(
        cat GeneScope-${K}/model.txt |
            grep '^kmercov' |
            tr -s ' ' '\t' |
            cut -f 2 |
            perl -nl -e 'printf qq{%.1f\n}, $_'
    )

    cat GeneScope-${K}/summary.txt |
        sed '1,6 d' |
        sed '1 s/^/K\t/' |
        sed "2 s/^/${K}\t/" |
        sed "3,7 s/^/\t/" |
        perl -nlp -e 's/\s{2,}/\t/g; s/\s+$//g;' |
        perl -nla -F'\t' -e '
            @fields = map {/\bNA\b/ ? q{} : $_ } @F;        # Remove NA fields
            $fields[2] = q{} if $fields[2] eq $fields[3];   # Remove identical fields
            print join qq{\t}, @fields
        '

    printf "\tKmer Cov\t\t${COV}\n"
done |
    keep-header -- grep -v '^K' |
    rgr md stdin --right 3-4

```

|       Size | Name             |
|-----------:|------------------|
|    262,164 | Table-21.hist    |
|    524,304 | Table-21.ktab    |
|  4,382,748 | .Table-21.ktab.1 |
|  4,457,262 | .Table-21.ktab.2 |
|  4,480,296 | .Table-21.ktab.3 |
|  4,378,950 | .Table-21.ktab.4 |
|    262,164 | Table-51.hist    |
|    524,304 | Table-51.ktab    |
| 12,613,366 | .Table-51.ktab.1 |
| 10,755,393 | .Table-51.ktab.2 |
| 11,851,397 | .Table-51.ktab.3 |
| 11,481,742 | .Table-51.ktab.4 |

| K  | property              |        min |        max |
|----|-----------------------|-----------:|-----------:|
| 21 | Homozygous (a)        |            |       100% |
|    | Genome Haploid Length |            | 577,872 bp |
|    | Genome Repeat Length  |   3,638 bp |   3,643 bp |
|    | Genome Unique Length  | 573,811 bp | 574,653 bp |
|    | Model Fit             |   92.6091% |   93.3366% |
|    | Read Error Rate       |            |  0.137192% |
|    | Kmer Cov              |            |      148.8 |
| 51 | Homozygous (a)        |            |       100% |
|    | Genome Haploid Length |            | 578,025 bp |
|    | Genome Repeat Length  |            |       0 bp |
|    | Genome Unique Length  |            | 578,025 bp |
|    | Model Fit             |   95.6385% |   95.7255% |
|    | Read Error Rate       |            | 0.0942288% |
|    | Kmer Cov              |            |      112.2 |

## Merqury

> these tables must be produced with the option -t1 set

```shell
cd ~/data/anchr/mg1655/2_illumina

mkdir "Merqury"

# KatGC
time FastK -v -t1 -k21 R*.fq.gz -NTable-21
#real    0m37.733s
#user    1m7.119s
#sys     0m7.661s

time FastK -v -t1 -k51 R*.fq.gz -NTable-51
#real    0m40.730s
#user    1m17.317s
#sys     0m12.324s

KatGC -T4 -x1.9 -s Table-21 Merqury/KatGC-21
KatGC -T4 -x1.9 -s Table-51 Merqury/KatGC-51

Fastrm Table-21 Table-51

# KatComp
FastK -v -t1 -k21 R1.fq.gz -NR1-21
FastK -v -t1 -k21 R2.fq.gz -NR2-21

KatComp -T4 -s R1-21 R2-21 Merqury/KatComp-21

Fastrm R1-21 R2-21

```

## Find repetitive regions

```shell
cd ~/data/anchr/mg1655/1_genome/

hnsm size genome.fa > chr.sizes

FastK -v -p -k21 genome

cat chr.sizes |
    number-lines |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        Profex genome {1} |
            sed "1,2 d" |
            perl -nl -e '\''/(\d+).+(\d+)/ and print qq{$1\t$2}'\'' |
            tsv-filter --ge 2:2 |
            cut -f 1 |
            sed "s/^/{2}:/" |
            spanr cover stdin |
            spanr span --op fill -n 2 stdin |
            spanr span --op excise -n 100 stdin |
            spanr span --op fill -n 10 stdin
    ' |
    spanr combine stdin \
    > repetitive.json

Fastrm genome

spanr convert repetitive.json > region.txt
hnsm range genome.fa -r region.txt |
    hnsm filter -N -d -a 100 stdin \
    > repetitive.fa

spanr stat chr.sizes repetitive.json
#chr,chrLength,size,coverage
#NC_000913,4641652,91989,0.0198
#all,4641652,91989,0.0198

hnsm size repetitive.fa | tsv-summarize --sum 2
#91989

hnsm size paralogs.fa | tsv-summarize --sum 2
#187300

```

## FastGA

```shell
mkdir -p ~/data/anchr/test
cd ~/data/anchr/test

curl -LO https://raw.githubusercontent.com/thegenemyers/FASTGA/refs/heads/main/EXAMPLE/HAP1.fasta.gz
curl -LO https://raw.githubusercontent.com/thegenemyers/FASTGA/refs/heads/main/EXAMPLE/HAP2.fasta.gz

FastGA -vk -1:H1vH2 HAP1 HAP2

ALNshow H1vH2 | more

ALNtoPAF H1vH2 | more

ALNtoPSL H1vH2 | more

ALNshow -a H1vH2 | more

ALNplot H1vH2 -p -W300 -l5000 -i.8 -n0

```

## Paralogs

```shell
cd ~/data/anchr/mg1655/1_genome/

time FastGA -v -l100 -i0.7 -1:self genome.fa

ALNchain -v -s1000 -z100 -g20 self

ALNtoPAF self |
    anchr paf2ovlp stdin |
    tsv-summarize --sum 3

ALNtoPAF self.chain |
    anchr paf2ovlp stdin |
    tsv-summarize --sum 3

ALNplot self.chain -p -W500 -n0

ALNtoPAF self.chain |
    anchr paf2ovlp stdin |
    tsv-filter --ge 3:1000 |
    perl -nla -F"\t" -e '
        # f
        if ($F[4] == 0) {
            printf qq(%s(+):%d-%d\t), $F[0], $F[5], $F[6]
        } else {
            printf qq(%s(-):%d-%d\t), $F[0], $F[6], $F[5]
        }

        # g
        if ($F[8] == 0) {
            printf qq(%s(+):%d-%d\n), $F[1], $F[9], $F[10]
        } else {
            printf qq(%s(-):%d-%d\n), $F[1], $F[10], $F[9]
        }
    ' |
    linkr sort stdin |
    linkr clean stdin \
    > self.link.tsv

rgr merge self.link.tsv -c 0.95 -o self.merge.tsv

linkr clean self.link.tsv -r self.merge.tsv --bundle 500 |
    linkr connect stdin |
    linkr filter stdin -r 0.8 -o self.filter.tsv

cat self.filter.tsv |
    tr '\t' '\n' |
    spanr cover stdin -o paralogs.json

rm self.*.tsv

spanr convert paralogs.json > region.txt
hnsm range genome.fa -r region.txt |
    hnsm filter -N -d -a 100 stdin \
    > paralogs.fa

spanr stat chr.sizes paralogs.json
#chr,chrLength,size,coverage
#NC_000913,4641652,260353,0.0561
#all,4641652,260353,0.0561

```
