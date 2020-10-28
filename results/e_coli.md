# Assemble genomes of model organisms by Anchr

[TOC levels=1-3]: # ""

- [Assemble genomes of model organisms by Anchr](#assemble-genomes-of-model-organisms-by-anchr)
- [More tools on downloading and preprocessing data](#more-tools-on-downloading-and-preprocessing-data)
  - [Extra external executables](#extra-external-executables)
  - [Other leading assemblers](#other-leading-assemblers)
- [*Escherichia* virus Lambda](#escherichia-virus-lambda)
  - [lambda: reference](#lambda-reference)
  - [lambda: download](#lambda-download)
  - [lambda: template](#lambda-template)
  - [lambda: run](#lambda-run)
- [*Escherichia coli* str. K-12 substr. MG1655](#escherichia-coli-str-k-12-substr-mg1655)
  - [mg1655: reference](#mg1655-reference)
  - [mg1655: download](#mg1655-download)
  - [mg1655: template](#mg1655-template)
  - [mg1655: run](#mg1655-run)
- [*Escherichia coli* str. K-12 substr. DH5α](#escherichia-coli-str-k-12-substr-dh5α)
  - [dh5alpha: reference](#dh5alpha-reference)
  - [dh5alpha: download](#dh5alpha-download)
  - [dh5alpha: template](#dh5alpha-template)
  - [dh5alpha: run](#dh5alpha-run)


# More tools on downloading and preprocessing data

## Extra external executables

```shell script
brew install aria2 curl                     # downloading tools
brew install miller

brew tap brewsci/bio
brew tap brewsci/science

brew install openblas                       # numpy
brew install python
brew install --HEAD quast         # assembly quality assessment. https://github.com/ablab/quast/issues/140
quast --test                                # may recompile the bundled nucmer

##brew install r
#brew install ntcard
#brew install wang-q/tap/kmergenie@1.7051

brew install --ignore-dependencies picard-tools

brew install kat
pip3 install tabulate

#kat comp -t 4 -n -o R 2_illumina/R1.fq.gz 2_illumina/R2.fq.gz -m 51
#kat plot spectra-mx -i -o R-spectra R-main.mx

```

## Other leading assemblers

```shell script
brew install spades
brew install megahit
brew install wang-q/tap/platanus

```

# *Escherichia* virus Lambda

* Taxonomy ID: [10710](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10710)
* Ref. Assembly: [GCF_000840245.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000840245.1/)

## lambda: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/ref
cd ~/data/anchr/ref

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/ \
    lambda/

```

```shell script
mkdir -p ~/data/anchr/lambda/1_genome
cd ~/data/anchr/lambda/1_genome

find ~/data/anchr/ref/lambda/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    xargs gzip -dcf |
    faops filter -N -s stdin genome.fa

touch paralogs.fa

```

## lambda: download

```shell script
mkdir -p ~/data/anchr/lambda/ena
cd ~/data/anchr/lambda/ena

cat << EOF > source.csv
SRX2365802,Lambda,HiSeq 2500
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

# sampling reads as test materials
seqtk sample -s 23 SRR5042715_1.fastq.gz 100000 | pigz > R1.fq.gz
seqtk sample -s 23 SRR5042715_2.fastq.gz 100000 | pigz > R2.fq.gz

```

| name   | srx        | platform | layout | ilength | srr        | spot     | base  |
|:-------|:-----------|:---------|:-------|:--------|:-----------|:---------|:------|
| Lambda | SRX2365802 | ILLUMINA | PAIRED |         | SRR5042715 | 16540237 | 3.33G |

* Illumina

```shell script
cd ~/data/anchr/lambda

mkdir -p 2_illumina
cd 2_illumina

seqtk sample -s 23 ../ena/SRR5042715_1.fastq.gz 100000 | pigz > R1.fq.gz
seqtk sample -s 23 ../ena/SRR5042715_2.fastq.gz 100000 | pigz > R2.fq.gz

```

## lambda: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=lambda

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 48502 \
    --parallel 4 \
    --xmx 6g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 5 --cutk 31" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    --cov "40 80" \
    --statp 2 

```

## lambda: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=lambda

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

bash 0_master.sh
#bash 0_cleanup.sh

```

Table: statInsertSize

| Group             |  Mean | Median | STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|------:|-------------------------------:|
| R.genome.bbtools  | 430.6 |    413 | 367.3 |                         96.29% |
| R.tadpole.bbtools | 418.5 |    407 | 101.3 |                         82.25% |
| R.genome.picard   | 424.6 |    413 | 103.1 |                             FR |
| R.tadpole.picard  | 418.6 |    407 | 101.3 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 334       | 46251           | 0            | 51.29   |
| R.31 | 297       | 46451           | 0            | 50.64   |
| R.41 | 260       | 46286           | 0            | 50.53   |
| R.51 | 217       | 46045           | 0            | 50.65   |
| R.61 | 180       | 46218           | 0            | 50.77   |
| R.71 | 143       | 46402           | 0            | 50.77   |


Table: statReads

| Name       |   N50 |    Sum |      # |
|:-----------|------:|-------:|-------:|
| Genome     | 48502 |  48502 |      1 |
| Paralogs   |     0 |      0 |      0 |
| Illumina.R |   108 |  21.6M | 200000 |
| trim.R     |   105 | 18.82M | 181790 |
| Q20L60     |   105 |  18.6M | 179679 |
| Q25L60     |   105 | 17.94M | 175041 |
| Q30L60     |   105 | 15.99M | 161943 |


Table: statTrimReads

| Name     | N50 |    Sum |      # |
|:---------|----:|-------:|-------:|
| clumpify | 108 | 21.56M | 199606 |
| highpass | 108 |  21.5M | 199078 |
| trim     | 105 | 18.82M | 181790 |
| filter   | 105 | 18.82M | 181790 |
| R1       | 105 |  9.48M |  90895 |
| R2       | 105 |  9.34M |  90895 |
| Rs       |   0 |      0 |      0 |


```text
#R.trim
#Matched	3807	1.91232%
#Name	Reads	ReadsPct
Reverse_adapter	1879	0.94385%
TruSeq_Universal_Adapter	1793	0.90065%
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	235058
#error_kmers	185792
#genomic_kmers	49266
#main_peak	269
#genome_size_in_peaks	2577144
#genome_size	2634418
#haploid_genome_size	48785
#fold_coverage	5
#haploid_fold_coverage	269
#ploidy	54
#het_rate	0.00080
#percent_repeat_in_peaks	0.000
#percent_repeat	98.128
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |    Sum |      # |
|:--------------|----:|-------:|-------:|
| clumped       | 105 | 18.82M | 181786 |
| ecco          | 105 | 18.82M | 181786 |
| eccc          | 105 | 18.82M | 181786 |
| ecct          | 105 | 18.77M | 181266 |
| extended      | 145 | 25.97M | 181266 |
| merged.raw    | 389 | 15.02M |  40052 |
| unmerged.raw  | 145 | 14.42M | 101162 |
| unmerged.trim | 145 | 14.42M | 101154 |
| M1            | 389 | 15.01M |  40032 |
| U1            | 145 |  7.27M |  50577 |
| U2            | 145 |  7.14M |  50577 |
| Us            |   0 |      0 |      0 |
| M.cor         | 267 | 29.46M | 181218 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 163.4 |    175 |  31.9 |          0.44% |
| M.ihist.merge.txt  | 374.9 |    383 |  47.4 |         44.19% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |   EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|-------:|---------:|----------:|
| Q0L0.R   | 388.0 |  370.5 |    4.50% | "73" | 48.5K | 49.26K |     1.02 | 0:00'13'' |
| Q20L60.R | 383.5 |  368.4 |    3.94% | "73" | 48.5K | 48.97K |     1.01 | 0:00'15'' |
| Q25L60.R | 370.0 |  358.0 |    3.24% | "73" | 48.5K | 48.62K |     1.00 | 0:00'13'' |
| Q30L60.R | 329.9 |  321.7 |    2.50% | "73" | 48.5K | 48.52K |     1.00 | 0:00'12'' |


Table: statUnitigsAnchors.md

| Name          | CovCor | Mapped% | N50Anchor |    Sum |  # | N50Others |   Sum |  # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|-------:|---:|----------:|------:|---:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.55% |     35403 | 46.88K |  3 |      1088 | 1.94K | 19 |   40.0 |  3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q0L0X40P001   |   40.0 |  99.14% |     30673 | 47.51K |  4 |        88 |   383 | 12 |   40.0 |  4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q0L0X40P002   |   40.0 |  99.60% |     30484 | 48.13K |  4 |        91 |   620 | 12 |   40.0 |  4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q0L0X80P000   |   80.0 |  97.35% |     30637 | 46.46K |  5 |       137 |   672 | 14 |   80.0 |  7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:03 |   0:00:11 |
| Q0L0X80P001   |   80.0 |  98.09% |     30513 | 47.23K |  6 |        71 |   496 | 15 |   79.0 |  9.0 |  17.3 | 158.0 | "31,41,51,61,71,81" |   0:00:04 |   0:00:10 |
| Q0L0X80P002   |   80.0 |  98.56% |     10585 | 47.43K |  6 |       191 |   697 | 14 |   82.0 | 10.0 |  17.3 | 164.0 | "31,41,51,61,71,81" |   0:00:04 |   0:00:10 |
| Q20L60X40P000 |   40.0 |  99.90% |     35649 | 48.11K |  3 |       190 |   779 | 14 |   40.0 |  4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q20L60X40P001 |   40.0 |  99.83% |     40312 | 49.37K |  3 |        51 |   538 | 12 |   39.0 |  2.5 |  10.5 |  69.8 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q20L60X40P002 |   40.0 |  99.83% |     12925 |  48.6K |  4 |        69 |   911 | 18 |   40.0 |  3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q20L60X80P000 |   80.0 |  99.83% |     35702 | 48.14K |  3 |       178 |   559 |  8 |   80.0 |  7.5 |  19.2 | 153.8 | "31,41,51,61,71,81" |   0:00:04 |   0:00:09 |
| Q20L60X80P001 |   80.0 |  99.52% |     26903 | 48.01K |  5 |       109 |   582 | 12 |   81.0 |  8.0 |  19.0 | 157.5 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q20L60X80P002 |   80.0 |  99.61% |     11448 | 48.31K |  6 |        64 |   493 | 14 |   80.0 |  7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:04 |   0:00:10 |
| Q25L60X40P000 |   40.0 |  99.48% |     35412 | 47.77K |  3 |        97 |   590 |  9 |   41.0 |  2.0 |  11.7 |  70.5 | "31,41,51,61,71,81" |   0:00:02 |   0:00:10 |
| Q25L60X40P001 |   40.0 |  99.84% |     30597 |    48K |  4 |        79 |   569 | 11 |   41.0 |  1.5 |  12.2 |  68.2 | "31,41,51,61,71,81" |   0:00:02 |   0:00:09 |
| Q25L60X40P002 |   40.0 |  99.87% |     35732 | 48.21K |  3 |        62 |   550 | 14 |   41.0 |  5.0 |   8.7 |  82.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q25L60X80P000 |   80.0 |  99.55% |     14140 | 48.79K |  4 |       120 |   746 | 12 |   81.0 |  7.0 |  20.0 | 153.0 | "31,41,51,61,71,81" |   0:00:04 |   0:00:09 |
| Q25L60X80P001 |   80.0 |  99.91% |     35725 | 48.22K |  3 |       167 |   581 | 10 |   81.0 | 10.0 |  17.0 | 162.0 | "31,41,51,61,71,81" |   0:00:04 |   0:00:10 |
| Q25L60X80P002 |   80.0 |  99.84% |     21453 | 48.31K |  4 |       121 |   462 |  9 |   79.0 |  9.0 |  17.3 | 158.0 | "31,41,51,61,71,81" |   0:00:04 |   0:00:10 |
| Q30L60X40P000 |   40.0 |  99.95% |      2784 |  12.5K |  4 |        77 |   346 |  7 |   41.0 |  1.0 |  12.7 |  66.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q30L60X40P001 |   40.0 |  99.91% |     35374 |  47.8K |  2 |       316 | 1.03K |  8 |   40.0 |  4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q30L60X40P002 |   40.0 |  99.91% |     31880 | 48.01K |  4 |       292 |   517 |  2 |   41.0 |  1.5 |  12.2 |  68.2 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| Q30L60X80P000 |   80.0 |  99.91% |     35641 | 47.96K |  3 |       174 |   588 |  7 |   80.0 |  4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:00:04 |   0:00:09 |
| Q30L60X80P001 |   80.0 |  98.91% |     35694 | 47.31K |  3 |       188 |   521 |  8 |   83.0 |  9.0 |  18.7 | 165.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:10 |
| Q30L60X80P002 |   80.0 |  99.92% |     35724 | 48.19K |  3 |       164 |   559 |  8 |   82.0 |  9.0 |  18.3 | 163.5 | "31,41,51,61,71,81" |   0:00:04 |   0:00:10 |


Table: statMRUnitigsAnchors.md

| Name      | CovCor | Mapped% | N50Anchor |    Sum |  # | N50Others | Sum |  # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|-------:|---:|----------:|----:|---:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  99.55% |     35413 | 48.01K |  2 |       205 | 582 |  4 |   40.0 | 2.5 |  10.8 |  71.2 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| MRX40P001 |   40.0 |  99.23% |     35536 | 47.97K |  3 |       198 | 703 |  6 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:03 |   0:00:09 |
| MRX40P002 |   40.0 |  99.64% |     35802 | 48.32K |  2 |       143 | 270 |  3 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:04 |   0:00:09 |
| MRX80P000 |   80.0 |  99.04% |     35412 | 47.88K |  3 |       202 | 583 |  6 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:04 |   0:00:10 |
| MRX80P001 |   80.0 |  99.15% |     35647 | 48.15K |  3 |       158 | 306 |  5 |   80.0 | 8.0 |  18.7 | 156.0 | "31,41,51,61,71,81" |   0:00:04 |   0:00:09 |
| MRX80P002 |   80.0 |  99.32% |     30819 | 48.31K |  4 |        86 | 441 |  7 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:04 |   0:00:09 |


Table: statMergeAnchors.md

| Name                       | Mapped% | N50Anchor |    Sum |  # | N50Others |   Sum |  # | median | MAD | lower | upper | RunTimeAN |
|:---------------------------|--------:|----------:|-------:|---:|----------:|------:|---:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors            |   0.00% |     48514 | 48.51K |  1 |         0 |     0 |  0 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_mr_unitigs_anchors |   0.00% |     48514 | 48.51K |  1 |         0 |     0 |  0 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_unitigs_anchors    |   0.00% |     48348 | 48.35K |  1 |      1236 | 1.24K |  1 |    0.0 | 0.0 |   0.0 |   0.0 |           |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |    Sum |  # | N50Others | Sum |  # | median | MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|-------:|---:|----------:|----:|---:|-------:|----:|------:|------:|----------:|
| 8_spades     |  99.83% |     48024 | 48.02K |  1 |       267 | 498 |  2 |  370.0 | 1.0 | 122.3 | 559.5 |   0:00:10 |
| 8_mr_spades  |  99.57% |     26356 | 48.08K |  2 |       429 | 543 |  3 |  603.0 | 1.0 | 200.0 | 909.0 |   0:00:10 |
| 8_megahit    |  99.83% |     48024 | 48.02K |  1 |       267 | 495 |  2 |  370.0 | 1.0 | 122.3 | 559.5 |   0:00:10 |
| 8_mr_megahit |  99.74% |     48072 | 48.07K |  1 |       242 | 442 |  2 |  605.0 | 0.0 | 201.7 | 907.5 |   0:00:10 |


Table: statFinal

| Name                     |   N50 |   Sum |  # |
|:-------------------------|------:|------:|---:|
| Genome                   | 48502 | 48502 |  1 |
| Paralogs                 |     0 |     0 |  0 |
| 7_merge_anchors.anchors  | 48514 | 48514 |  1 |
| 7_merge_anchors.others   |     0 |     0 |  0 |
| spades.contig            | 48522 | 48600 |  2 |
| spades.scaffold          | 48522 | 48600 |  2 |
| spades.non-contained     | 48522 | 48522 |  1 |
| spades_MR.contig         | 48627 | 48753 |  2 |
| spades_MR.scaffold       | 48627 | 48753 |  2 |
| spades_MR.non-contained  | 48627 | 48627 |  1 |
| megahit.contig           | 48519 | 48519 |  1 |
| megahit.non-contained    | 48519 | 48519 |  1 |
| megahit_MR.contig        | 48514 | 48514 |  1 |
| megahit_MR.non-contained | 48514 | 48514 |  1 |
| platanus.non-contained   |     0 |     0 |  0 |


# *Escherichia coli* str. K-12 substr. MG1655

* Taxonomy ID: [511145](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=511145)
* Genome: INSDC [U00096.3](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3)
* Assembly: [GCF_000005845.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2)
* Proportion of paralogs (> 1000 bp): 0.0325

## mg1655: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/ref
cd ~/data/anchr/ref

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/ \
    mg1655/

```

```shell script
mkdir -p ~/data/anchr/mg1655/1_genome
cd ~/data/anchr/mg1655/1_genome

find ~/data/anchr/ref/mg1655/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    xargs gzip -dcf |
    faops filter -N -s stdin genome.fa

cat ~/data/anchr/paralogs/model/Results/mg1655/mg1655.multi.fas |
    faops filter -N -d stdin stdout \
    > paralogs.fa

```

## mg1655: download

```shell script
cd ~/data/anchr/mg1655

mkdir -p ena
cd ena

cat << EOF > source.csv
ERX008638,mg1655,PE100, GA IIx
EOF

ena_info.pl -v source.csv > ena_info.yml
ena_prep.pl ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name   | srx       | platform | layout | ilength | srr       | spot     | base  |
|:-------|:----------|:---------|:-------|:--------|:----------|:---------|:------|
| mg1655 | ERX008638 | ILLUMINA | PAIRED | 311     | ERR022075 | 22720100 | 4.27G |


* Illumina

```shell script
cd ~/data/anchr/mg1655

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/ERR022075_1.fastq.gz R1.fq.gz
ln -s ../ena/ERR022075_2.fastq.gz R2.fq.gz

#aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz
#aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz
#
#ln -s MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz R1.fq.gz
#ln -s MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz R2.fq.gz

```

## mg1655: template

* Rsync to hpcc

```shell script
rsync -avP \
    --exclude="p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz" \
    ~/data/anchr/mg1655/ \
    wangq@202.119.37.251:data/anchr/mg1655

# rsync -avP wangq@202.119.37.251:data/anchr/mg1655/ ~/data/anchr/mg1655

```

* template

  Filling tiles: Trouble parsing header ERR022075.3334391 EAS600_70:5:18:8880:11520/1

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4641652 \
    --parallel 24 \
    --xmx 10g \
    --queue mpi \
    \
    --fastqc \
    --kmergenie \
    \
    --trim "--dedupe --cutoff 5 --cutk 31" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3"

```

## mg1655: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/mergereads statReads.md 

bash 0_bsub.sh
#bsub -q largemem -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

```


# *Escherichia coli* str. K-12 substr. DH5α

* Taxonomy ID: [83333](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=83333)
* Genome: [CP017100](https://www.ncbi.nlm.nih.gov/nuccore/CP017100)
* Assembly: [GCF_001723515.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_001723515.1)
* Proportion of paralogs (> 1000 bp): 0.0342

## dh5alpha: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/ref
cd ~/data/anchr/ref

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/001/723/505/GCF_001723505.1_ASM172350v1/ \
    dh5alpha/

```

```shell script
mkdir -p ~/data/anchr/dh5alpha/1_genome
cd ~/data/anchr/dh5alpha/1_genome

find ~/data/anchr/ref/dh5alpha/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    xargs gzip -dcf |
    faops filter -N -s stdin genome.fa

cat ~/data/anchr/paralogs/e_coli/Results/dh5alpha/dh5alpha.multi.fas |
    faops filter -N -d stdin stdout \
    > paralogs.fa

```

## dh5alpha: download

```shell script
cd ~/data/anchr/dh5alpha

mkdir -p ena
cd ena

cat << EOF > source.csv
SRP251726,dh5alpha,HiSeq 2500 PE125
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name     | srx        | platform        | layout | ilength | srr         | spot    | base  |
|:---------|:-----------|:----------------|:-------|:--------|:------------|:--------|:------|
| dh5alpha | SRX7856678 | ILLUMINA        | PAIRED |         | SRR11245239 | 5881654 | 1.37G |
| dh5alpha | SRX7856679 | OXFORD_NANOPORE | SINGLE |         | SRR11245238 | 346489  | 3.35G |


* Illumina

```shell script
cd ~/data/anchr/dh5alpha

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR11245239_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR11245239_2.fastq.gz R2.fq.gz

```

## dh5alpha: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4583637 \
    --parallel 4 \
    --xmx 6g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 5 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    --cov "40 80" \
    --statp 2 \
    --redoanchors

```

## dh5alpha: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

bash 0_master.sh
#bash 0_cleanup.sh

```

Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 432.7 |    346 | 1241.2 |                         99.98% |
| R.tadpole.bbtools | 389.1 |    341 |  206.1 |                         95.12% |
| R.genome.picard   | 394.8 |    346 |  208.3 |                             FR |
| R.tadpole.picard  | 389.1 |    341 |  205.7 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 260       | 4346013         | 0.0814       | 51.62   |
| R.31 | 233       | 4360998         | 0.0563       | 51.54   |
| R.41 | 209       | 4433088         | 0.0607       | 51.50   |
| R.51 | 182       | 4413726         | 0.0150       | 51.47   |
| R.61 | 157       | 4456634         | 0.0332       | 51.45   |
| R.71 | 132       | 4495949         | 0.0194       | 51.44   |
| R.81 | 107       | 4475448         | 0.0160       | 51.42   |


Table: statReads

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 4583637 | 4583637 |        1 |
| Paralogs   |    1737 |  188158 |      111 |
| Illumina.R |     125 |   1.47G | 11763308 |
| trim.R     |     125 |   1.37G | 10962374 |
| Q25L60     |     125 |   1.25G | 10281036 |
| Q30L60     |     125 |   1.13G |  9405631 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 125 |   1.37G | 10970448 |
| highpass | 125 |   1.37G | 10966250 |
| trim     | 125 |   1.37G | 10962374 |
| filter   | 125 |   1.37G | 10962374 |
| R1       | 125 |    683M |  5481187 |
| R2       | 125 | 683.64M |  5481187 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	5620	0.05125%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	35781051
#error_kmers	31303388
#genomic_kmers	4477663
#main_peak	223
#genome_size_in_peaks	4561456
#genome_size	4571195
#haploid_genome_size	4571195
#fold_coverage	223
#haploid_fold_coverage	223
#ploidy	1
#percent_repeat_in_peaks	1.841
#percent_repeat	1.987
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |        # |
|:--------------|----:|--------:|---------:|
| clumped       | 125 |   1.37G | 10959556 |
| ecco          | 125 |   1.37G | 10959556 |
| eccc          | 125 |   1.37G | 10959556 |
| ecct          | 125 |   1.37G | 10952902 |
| extended      | 165 |    1.8G | 10952902 |
| merged.raw    | 343 |    1.1G |  3511427 |
| unmerged.raw  | 165 | 645.95M |  3930048 |
| unmerged.trim | 165 | 645.95M |  3930048 |
| M1            | 343 |   1.06G |  3404031 |
| U1            | 165 | 322.88M |  1965024 |
| U2            | 165 | 323.07M |  1965024 |
| Us            |   0 |       0 |        0 |
| M.cor         | 250 |   1.71G | 10738110 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 180.7 |    181 |  32.1 |         23.40% |
| M.ihist.merge.txt  | 312.5 |    310 |  87.7 |         64.12% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 298.2 |  262.2 |   12.07% | "63" | 4.58M |  4.6M |     1.00 | 0:08'26'' |
| Q25L60.R | 273.7 |  253.6 |    7.35% | "63" | 4.58M | 4.53M |     0.99 | 0:06'56'' |
| Q30L60.R | 246.5 |  232.4 |    5.71% | "63" | 4.58M | 4.52M |     0.99 | 0:06'23'' |


Table: statUnitigsAnchors.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.16% |     22066 | 4.36M | 355 |       561 | 63.96K |  762 |   39.0 | 2.0 |  11.0 |  67.5 | "31,41,51,61,71,81" |   0:02:34 |   0:00:33 |
| Q0L0X40P001   |   40.0 |  97.34% |     22460 | 4.42M | 308 |       863 | 69.54K |  702 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:02:34 |   0:00:34 |
| Q0L0X40P002   |   40.0 |  97.15% |     25244 | 4.37M | 300 |       761 | 57.69K |  653 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:02:34 |   0:00:28 |
| Q0L0X80P000   |   80.0 |  96.07% |     12127 | 4.41M | 551 |        42 | 57.69K | 1170 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:04:24 |   0:00:33 |
| Q0L0X80P001   |   80.0 |  95.80% |     12425 | 4.38M | 529 |        42 | 50.83K | 1111 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:04:25 |   0:00:32 |
| Q0L0X80P002   |   80.0 |  95.90% |     12642 | 4.38M | 547 |        45 | 62.69K | 1154 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:04:25 |   0:00:35 |
| Q25L60X40P000 |   40.0 |  97.70% |     33418 | 4.27M | 237 |       609 |  47.7K |  560 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:02:33 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  97.69% |     36000 | 4.43M | 238 |       721 | 51.07K |  550 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:02:35 |   0:00:30 |
| Q25L60X40P002 |   40.0 |  97.71% |     33538 | 4.37M | 219 |       390 | 39.38K |  505 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:02:35 |   0:00:29 |
| Q25L60X80P000 |   80.0 |  96.96% |     24644 | 4.36M | 320 |        47 | 41.02K |  724 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:04:20 |   0:00:33 |
| Q25L60X80P001 |   80.0 |  97.08% |     24735 | 4.36M | 309 |        47 | 36.44K |  679 |   79.0 | 3.0 |  23.3 | 132.0 | "31,41,51,61,71,81" |   0:04:20 |   0:00:32 |
| Q25L60X80P002 |   80.0 |  97.02% |     25603 | 4.44M | 309 |        44 | 36.25K |  686 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:04:21 |   0:00:33 |
| Q30L60X40P000 |   40.0 |  97.73% |     32737 |  4.2M | 234 |       793 | 46.28K |  502 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:02:33 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  97.84% |     31807 | 4.23M | 237 |       574 | 43.16K |  517 |   39.0 | 2.0 |  11.0 |  67.5 | "31,41,51,61,71,81" |   0:02:31 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  97.76% |     35160 | 4.24M | 219 |       908 | 53.22K |  535 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:02:35 |   0:00:29 |
| Q30L60X80P000 |   80.0 |  97.14% |     25893 | 4.44M | 294 |        42 | 30.48K |  643 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:04:17 |   0:00:33 |
| Q30L60X80P001 |   80.0 |  97.19% |     25578 | 4.35M | 290 |        48 | 36.98K |  648 |   79.0 | 3.0 |  23.3 | 132.0 | "31,41,51,61,71,81" |   0:04:18 |   0:00:37 |


Table: statMRUnitigsAnchors.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.60% |     59590 | 4.16M | 130 |        87 | 29.95K | 329 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:03:26 |   0:00:36 |
| MRX40P001 |   40.0 |  97.55% |     55454 | 4.09M | 129 |       129 | 40.57K | 336 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:03:28 |   0:00:36 |
| MRX40P002 |   40.0 |  97.62% |     57820 | 4.44M | 138 |       114 | 40.35K | 347 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:03:30 |   0:00:31 |
| MRX80P000 |   80.0 |  97.55% |     64083 | 4.44M | 121 |        81 |  27.6K | 323 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:06:08 |   0:00:34 |
| MRX80P001 |   80.0 |  97.51% |     63276 | 4.45M | 124 |        84 | 30.97K | 334 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:06:03 |   0:00:31 |
| MRX80P002 |   80.0 |  97.55% |     63612 | 4.27M | 119 |        85 | 31.49K | 326 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:06:03 |   0:00:35 |


Table: statMergeAnchors.md

| Name                       | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:---------------------------|--------:|----------:|------:|----:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors            |  97.02% |     73542 | 4.44M | 115 |     10851 | 246.22K | 70 |  263.0 | 11.0 |  76.7 | 444.0 |   0:00:45 |
| 7_merge_mr_unitigs_anchors |  98.86% |     73576 | 4.44M | 114 |     41158 |  77.41K | 26 |  262.0 | 12.5 |  74.8 | 449.2 |   0:01:24 |
| 7_merge_unitigs_anchors    |  98.95% |     73618 | 4.45M | 115 |      8683 | 194.92K | 63 |  263.0 | 13.0 |  74.7 | 453.0 |   0:01:26 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.45% |    112304 | 4.29M |  84 |      1036 | 24.41K | 157 |  263.0 | 12.0 |  75.7 | 448.5 |   0:00:47 |
| 8_mr_spades  |  98.94% |    125905 | 4.06M |  68 |       758 | 23.46K | 122 |  376.0 | 18.0 | 107.3 | 645.0 |   0:00:57 |
| 8_megahit    |  98.44% |     67313 | 4.46M | 119 |      1039 | 29.93K | 230 |  263.0 | 14.0 |  73.7 | 457.5 |   0:00:49 |
| 8_mr_megahit |  99.23% |    132514 | 4.49M |  79 |       705 | 34.04K | 144 |  376.0 | 18.0 | 107.3 | 645.0 |   0:00:53 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 4583637 | 4583637 |   1 |
| Paralogs                 |    1737 |  188158 | 111 |
| 7_merge_anchors.anchors  |   73542 | 4438968 | 115 |
| 7_merge_anchors.others   |   10851 |  246219 |  70 |
| spades.contig            |  125538 | 4514902 | 165 |
| spades.scaffold          |  143522 | 4515652 | 155 |
| spades.non-contained     |  132337 | 4493720 |  79 |
| spades_MR.contig         |  155808 | 4523183 |  89 |
| spades_MR.scaffold       |  203812 | 4523493 |  85 |
| spades_MR.non-contained  |  155808 | 4513386 |  61 |
| megahit.contig           |   82955 | 4508133 | 168 |
| megahit.non-contained    |   82955 | 4487915 | 116 |
| megahit_MR.contig        |  133730 | 4557848 | 158 |
| megahit_MR.non-contained |  133730 | 4522134 |  70 |
| platanus.non-contained   |       0 |       0 |   0 |

