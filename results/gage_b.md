# Assemble genomes from GAGE-B data sets

[TOC levels=1-3]: # ""

- [Assemble genomes from GAGE-B data sets](#assemble-genomes-from-gage-b-data-sets)
- [*Bacillus cereus* ATCC 10987](#bacillus-cereus-atcc-10987)
  - [Bcer_100x: reference](#bcer_100x-reference)
  - [Bcer_100x: download](#bcer_100x-download)
  - [Bcer_100x: template](#bcer_100x-template)
  - [Bcer_100x: run](#bcer_100x-run)
- [*Mycobacterium abscessus* 6G-0125-R](#mycobacterium-abscessus-6g-0125-r)
  - [Mabs_100x: reference](#mabs_100x-reference)
  - [Mabs_100x: download](#mabs_100x-download)
  - [Mabs_100x: template](#mabs_100x-template)
  - [Mabs_100x: run](#mabs_100x-run)
- [*Rhodobacter sphaeroides* 2.4.1](#rhodobacter-sphaeroides-241)
  - [Rsph_100x: reference](#rsph_100x-reference)
  - [Rsph_100x: download](#rsph_100x-download)
  - [Rsph_100x: template](#rsph_100x-template)
  - [Rsph_100x: run](#rsph_100x-run)
- [*Vibrio cholerae* CP1032(5)](#vibrio-cholerae-cp10325)
  - [Vcho_100x: reference](#vcho_100x-reference)
  - [Vcho_100x: download](#vcho_100x-download)
  - [Vcho_100x: template](#vcho_100x-template)
  - [Vcho_100x: run](#vcho_100x-run)
- [*Mycobacterium abscessus* 6G-0125-R Full](#mycobacterium-abscessus-6g-0125-r-full)
  - [Mabs_full: reference](#mabs_full-reference)
  - [Mabs_full: download](#mabs_full-download)
  - [Mabs_full: template](#mabs_full-template)
  - [Mabs_full: run](#mabs_full-run)
- [*Rhodobacter sphaeroides* 2.4.1 Full](#rhodobacter-sphaeroides-241-full)
  - [Rsph_full: reference](#rsph_full-reference)
  - [Rsph_full: download](#rsph_full-download)
  - [Rsph_full: template](#rsph_full-template)
  - [Rsph_full: run](#rsph_full-run)
- [*Vibrio cholerae* CP1032(5) Full](#vibrio-cholerae-cp10325-full)
  - [Vcho_full: reference](#vcho_full-reference)
  - [Vcho_full: download](#vcho_full-download)
  - [Vcho_full: template](#vcho_full-template)
  - [Vcho_full: run](#vcho_full-run)

* Rsync to hpcc

```shell script
for D in \
    Bcer_100x Mabs_100x Rsph_100x Vcho_100x \
    Mabs_full Rsph_full Vcho_full \
    ; do
    rsync -avP \
        ~/data/anchr/${D}/ \
        wangq@202.119.37.251:data/anchr/${D}
done

```


# *Bacillus cereus* ATCC 10987

## Bcer_100x: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Bcer_100x/1_genome
cd ~/data/anchr/Bcer_100x/1_genome

cp ~/data/anchr/ref/Bcer/genome.fa .
cp ~/data/anchr/ref/Bcer/paralogs.fa .
cp ~/data/anchr/ref/Bcer/repetitives.fa .

```

## Bcer_100x: download

* Illumina

```shell script
cd ~/data/anchr/Bcer_100x

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/datasets/B_cereus_MiSeq.tar.gz

# NOT gzipped tar
tar xvf B_cereus_MiSeq.tar.gz raw/frag_1__cov100x.fastq
tar xvf B_cereus_MiSeq.tar.gz raw/frag_2__cov100x.fastq

cat raw/frag_1__cov100x.fastq |
    pigz -p 8 -c \
    > R1.fq.gz
cat raw/frag_2__cov100x.fastq |
    pigz -p 8 -c \
    > R2.fq.gz

rm -fr raw

```

* GAGE-B assemblies

```shell script
cd ~/data/anchr/Bcer_100x

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/B_cereus_MiSeq.tar.gz

tar xvfz B_cereus_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz mira_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz sga_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz soap_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz spades_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz velvet_ctg.fasta

```

## Bcer_100x: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 5432652 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --tile --cutoff 10 --cutk 31" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 50 60 all" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

```

## Bcer_100x: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_100x

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Bcer_100x bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

rm -fr 9_quast_competitor
quast --threads 16 \
    -R 1_genome/genome.fa \
    8_competitor/abyss_ctg.fasta \
    8_competitor/cabog_ctg.fasta \
    8_competitor/mira_ctg.fasta \
    8_competitor/msrca_ctg.fasta \
    8_competitor/sga_ctg.fasta \
    8_competitor/soap_ctg.fasta \
    8_competitor/spades_ctg.fasta \
    8_competitor/velvet_ctg.fasta \
    7_merge_anchors/anchor.merge.fasta \
    7_merge_anchors/others.non-contained.fasta \
    1_genome/paralogs.fa \
    --label "abyss,cabog,mira,msrca,sga,soap,spades,velvet,merge,others,paralogs" \
    -o 9_quast_competitor

```


Table: statInsertSize

| Group             |  Mean | Median | STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|------:|-------------------------------:|
| R.genome.bbtools  | 589.6 |    583 | 726.8 |                         96.28% |
| R.tadpole.bbtools | 567.8 |    575 | 154.7 |                         86.55% |
| R.genome.picard   | 582.1 |    585 | 146.5 |                             FR |
| R.tadpole.picard  | 573.8 |    577 | 147.2 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 75        | 5368162         | 0.1207       | 38.81   |
| R.31 | 71        | 5319712         | 0.1500       | 38.55   |
| R.41 | 66        | 5307740         | 0.1333       | 38.38   |
| R.51 | 62        | 5224995         | 0.1404       | 38.25   |
| R.61 | 58        | 5195289         | 0.1313       | 38.14   |
| R.71 | 54        | 5251178         | 0.0978       | 38.05   |
| R.81 | 51        | 5188933         | 0.1113       | 37.97   |


Table: statReads

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 5224283 | 5432652 |       2 |
| Paralogs    |    2295 |  220468 |     101 |
| Repetitives |    2461 |  113050 |     173 |
| Illumina.R  |     251 | 481.02M | 2080000 |
| trim.R      |     250 | 404.63M | 1808658 |
| Q20L60      |     250 | 397.05M | 1759736 |
| Q25L60      |     250 | 379.81M | 1707367 |
| Q30L60      |     250 | 344.43M | 1612284 |


Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 480.99M | 2079856 |
| filteredbytile | 251 | 463.83M | 2007096 |
| highpass       | 251 | 460.12M | 1991060 |
| trim           | 250 | 404.68M | 1808850 |
| filter         | 250 | 404.63M | 1808658 |
| R1             | 250 | 209.41M |  904329 |
| R2             | 247 | 195.22M |  904329 |
| Rs             |   0 |       0 |       0 |


```text
#R.trim
#Matched	1242	0.06238%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	99	0.00547%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	13976982
#error_kmers	8762014
#genomic_kmers	5214968
#main_peak	68
#genome_size_in_peaks	5298354
#genome_size	5421175
#haploid_genome_size	5421175
#fold_coverage	68
#haploid_fold_coverage	68
#ploidy	1
#percent_repeat_in_peaks	1.576
#percent_repeat	2.021
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 250 | 404.63M | 1808656 |
| ecco          | 250 | 404.63M | 1808656 |
| eccc          | 250 | 404.63M | 1808656 |
| ecct          | 250 | 400.68M | 1787862 |
| extended      | 290 | 471.31M | 1787862 |
| merged.raw    | 586 | 316.63M |  584218 |
| unmerged.raw  | 285 | 150.01M |  619426 |
| unmerged.trim | 285 | 150.01M |  619402 |
| M1            | 586 |  316.6M |  584164 |
| U1            | 290 |   80.1M |  309701 |
| U2            | 270 |   69.9M |  309701 |
| Us            |   0 |       0 |       0 |
| M.cor         | 518 | 467.19M | 1787730 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 362.6 |    389 |  97.8 |         19.12% |
| M.ihist.merge.txt  | 542.0 |    564 | 120.0 |         65.35% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   |  74.5 |   64.8 |   12.98% | "127" | 5.43M | 5.35M |     0.98 | 0:00'54'' |
| Q20L60.R |  73.1 |   64.7 |   11.48% | "127" | 5.43M | 5.35M |     0.98 | 0:00'49'' |
| Q25L60.R |  69.9 |   63.8 |    8.71% | "127" | 5.43M | 5.34M |     0.98 | 0:00'48'' |
| Q30L60.R |  63.4 |   59.8 |    5.75% | "127" | 5.43M | 5.34M |     0.98 | 0:00'44'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  97.38% |     15742 | 5.18M | 511 |       312 | 80.85K | 1173 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:44 |
| Q0L0X50P000    |   50.0 |  97.38% |     16430 | 5.21M | 505 |       413 |  81.2K | 1174 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:45 |
| Q0L0X60P000    |   60.0 |  97.34% |     17079 | 5.24M | 491 |       179 | 64.71K | 1148 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:43 |
| Q0L0XallP000   |   64.8 |  97.34% |     17266 | 5.23M | 503 |       265 | 73.95K | 1171 |   64.0 | 7.0 |  16.7 | 156.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:44 |
| Q20L60X40P000  |   40.0 |  97.69% |     16913 | 5.11M | 491 |       347 | 86.29K | 1119 |   39.0 | 4.0 |  10.3 |  94.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:44 |
| Q20L60X50P000  |   50.0 |  97.61% |     19176 |  5.2M | 462 |       274 | 73.31K | 1086 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:45 |
| Q20L60X60P000  |   60.0 |  97.56% |     19570 | 5.21M | 454 |       180 | 62.27K | 1100 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:45 |
| Q20L60XallP000 |   64.7 |  97.52% |     19100 | 5.22M | 465 |       295 | 72.05K | 1105 |   64.0 | 7.0 |  16.7 | 156.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:46 |
| Q25L60X40P000  |   40.0 |  97.89% |     21378 |  4.9M | 387 |       176 | 63.55K |  999 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:45 |
| Q25L60X50P000  |   50.0 |  97.96% |     20784 | 5.07M | 412 |       210 | 62.88K | 1000 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:45 |
| Q25L60X60P000  |   60.0 |  97.97% |     21379 | 5.11M | 396 |       167 | 56.47K |  986 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:47 |
| Q25L60XallP000 |   63.8 |  97.98% |     20654 | 5.11M | 415 |       180 | 61.31K | 1015 |   63.0 | 7.0 |  16.3 | 154.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:49 |
| Q30L60X40P000  |   40.0 |  98.24% |     19888 | 4.86M | 415 |       538 |  80.3K |  959 |   39.0 | 4.0 |  10.3 |  94.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:43 |
| Q30L60X50P000  |   50.0 |  98.29% |     22000 | 5.02M | 386 |       448 | 65.02K |  957 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:46 |
| Q30L60XallP000 |   59.8 |  98.37% |     21751 | 4.92M | 369 |       337 | 58.02K |  960 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:44 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.50% |     32698 | 5.22M | 267 |       158 | 36.21K | 603 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:39 |
| MRX40P001  |   40.0 |  97.52% |     31702 | 5.23M | 270 |       132 | 30.74K | 609 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:40 |
| MRX50P000  |   50.0 |  97.47% |     31938 | 5.22M | 269 |       141 | 32.45K | 605 |   49.0 |  6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:39 |
| MRX60P000  |   60.0 |  97.45% |     32708 | 5.23M | 264 |       140 | 28.68K | 603 |   59.0 |  8.0 |  14.3 | 150.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:41 |
| MRXallP000 |   86.0 |  97.37% |     31806 |  5.3M | 274 |       135 | 26.39K | 625 |   84.0 | 11.0 |  20.7 | 212.0 | "31,41,51,61,71,81" |   0:01:59 |   0:00:41 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  97.52% |     17419 | 5.16M | 458 |       429 | 76.18K | 1113 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:41 |
| Q0L0X50P000    |   50.0 |  97.46% |     17821 | 5.22M | 469 |       353 | 74.02K | 1137 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:43 |
| Q0L0X60P000    |   60.0 |  97.38% |     17763 | 5.26M | 469 |       146 | 59.53K | 1160 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:43 |
| Q0L0XallP000   |   64.8 |  97.38% |     17680 | 5.26M | 485 |       354 | 72.59K | 1172 |   64.0 | 7.0 |  16.7 | 156.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:44 |
| Q20L60X40P000  |   40.0 |  97.56% |     20647 | 5.25M | 438 |       260 |  73.2K | 1126 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:45 |
| Q20L60X50P000  |   50.0 |  97.54% |     19360 | 5.21M | 458 |       416 |  77.7K | 1103 |   49.0 | 5.0 |  13.0 | 118.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:42 |
| Q20L60X60P000  |   60.0 |  97.58% |     19751 | 5.25M | 447 |       167 | 61.15K | 1111 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:43 |
| Q20L60XallP000 |   64.7 |  97.56% |     19551 | 5.22M | 452 |       180 | 66.43K | 1122 |   64.0 | 7.0 |  16.7 | 156.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:45 |
| Q25L60X40P000  |   40.0 |  97.70% |     21573 | 5.23M | 400 |       317 | 63.32K | 1041 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:43 |
| Q25L60X50P000  |   50.0 |  97.76% |     21425 | 5.26M | 419 |       414 | 64.54K | 1049 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:42 |
| Q25L60X60P000  |   60.0 |  97.80% |     21425 | 5.14M | 407 |       410 | 63.06K | 1014 |   60.0 | 7.0 |  15.3 | 148.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:46 |
| Q25L60XallP000 |   63.8 |  97.86% |     21287 | 5.18M | 412 |       166 | 58.13K | 1036 |   63.0 | 7.0 |  16.3 | 154.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:46 |
| Q30L60X40P000  |   40.0 |  98.08% |     19406 | 5.11M | 438 |       445 | 73.06K | 1077 |   40.0 | 4.0 |  10.7 |  96.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:42 |
| Q30L60X50P000  |   50.0 |  98.16% |     20654 | 5.08M | 422 |       491 | 75.51K | 1048 |   50.0 | 5.0 |  13.3 | 120.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:45 |
| Q30L60XallP000 |   59.8 |  98.18% |     21845 |  5.2M | 387 |       478 |  56.2K |  994 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:42 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.33% |     33935 | 5.07M | 253 |       148 | 28.72K | 522 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:36 |
| MRX40P001  |   40.0 |  97.38% |     32819 | 5.23M | 254 |       144 | 27.89K | 528 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:38 |
| MRX50P000  |   50.0 |  97.35% |     33200 | 5.07M | 257 |       141 | 27.71K | 542 |   49.0 |  6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:38 |
| MRX60P000  |   60.0 |  97.28% |     33944 | 5.23M | 258 |       139 | 25.62K | 537 |   59.0 |  7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:37 |
| MRXallP000 |   86.0 |  97.26% |     32712 | 5.31M | 265 |       152 | 25.86K | 554 |   84.0 | 10.0 |  21.3 | 208.0 | "31,41,51,61,71,81" |   0:01:56 |   0:00:39 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  98.26% |     19705 | 5.14M | 418 |       429 |  74.6K | 1000 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:41 |
| Q0L0X50P000    |   50.0 |  98.18% |     20498 | 5.22M | 432 |       252 | 64.34K | 1008 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:42 |
| Q0L0X60P000    |   60.0 |  98.22% |     20760 | 5.23M | 424 |       156 |  55.3K | 1032 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:42 |
| Q0L0XallP000   |   64.8 |  98.19% |     20422 | 5.22M | 442 |       168 | 60.03K | 1028 |   64.0 | 7.0 |  16.7 | 156.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:43 |
| Q20L60X40P000  |   40.0 |  97.97% |     21733 | 5.15M | 399 |       290 | 64.76K |  963 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:43 |
| Q20L60X50P000  |   50.0 |  98.00% |     20905 | 5.19M | 423 |       429 | 71.78K | 1018 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:44 |
| Q20L60X60P000  |   60.0 |  98.28% |     21397 | 5.23M | 413 |       199 | 58.22K |  987 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:42 |
| Q20L60XallP000 |   64.7 |  98.24% |     21397 | 5.24M | 411 |       141 | 50.58K |  983 |   64.0 | 8.0 |  16.0 | 160.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:43 |
| Q25L60X40P000  |   40.0 |  98.44% |     21980 | 5.12M | 380 |       405 | 63.03K |  966 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:43 |
| Q25L60X50P000  |   50.0 |  98.46% |     21767 |  5.1M | 391 |       405 | 64.19K |  975 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:44 |
| Q25L60X60P000  |   60.0 |  98.48% |     22437 | 5.16M | 382 |       160 | 53.96K |  950 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:45 |
| Q25L60XallP000 |   63.8 |  98.48% |     22390 | 5.16M | 384 |       151 |  50.2K |  952 |   63.0 | 8.0 |  15.7 | 158.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:46 |
| Q30L60X40P000  |   40.0 |  98.54% |     19986 | 5.07M | 426 |       478 | 75.57K | 1031 |   40.0 | 4.0 |  10.7 |  96.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:42 |
| Q30L60X50P000  |   50.0 |  98.58% |     22532 | 5.12M | 387 |       448 | 62.91K |  970 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:43 |
| Q30L60XallP000 |   59.8 |  98.58% |     22530 | 5.11M | 372 |       224 | 50.86K |  936 |   59.0 | 7.0 |  15.0 | 146.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:45 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.70% |     39818 | 5.23M | 245 |       175 | 30.26K | 503 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:36 |
| MRX40P001  |   40.0 |  97.42% |     34498 | 5.23M | 240 |       140 | 25.65K | 489 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:37 |
| MRX50P000  |   50.0 |  97.40% |     38973 | 5.23M | 249 |       152 | 27.39K | 517 |   49.0 |  6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:38 |
| MRX60P000  |   60.0 |  97.37% |     38978 | 5.23M | 245 |       140 | 24.37K | 520 |   59.0 |  8.0 |  14.3 | 150.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:40 |
| MRXallP000 |   86.0 |  97.34% |     34382 | 5.31M | 254 |       135 | 21.96K | 533 |   85.0 | 11.0 |  21.0 | 214.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:38 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  97.58% |     25734 | 5.29M | 363 |      1033 | 108.72K | 111 |   65.0 | 8.0 |  16.3 | 162.0 |   0:00:49 |
| 7_merge_mr_unitigs_bcalm      |  98.20% |     25394 | 5.15M | 359 |      1007 |  27.36K |  31 |   65.0 | 8.0 |  16.3 | 162.0 |   0:00:59 |
| 7_merge_mr_unitigs_superreads |  98.20% |     23695 | 5.13M | 373 |      1007 |  29.12K |  34 |   65.0 | 8.0 |  16.3 | 162.0 |   0:00:58 |
| 7_merge_mr_unitigs_tadpole    |  98.21% |     24808 |  5.1M | 351 |      1010 |  27.56K |  31 |   65.0 | 8.0 |  16.3 | 162.0 |   0:00:58 |
| 7_merge_unitigs_bcalm         |  98.12% |     23637 | 5.29M | 388 |      1031 |  87.33K |  89 |   65.0 | 8.0 |  16.3 | 162.0 |   0:00:59 |
| 7_merge_unitigs_superreads    |  97.84% |     23734 | 5.26M | 371 |      1042 |  83.87K |  87 |   64.0 | 8.0 |  16.0 | 160.0 |   0:00:57 |
| 7_merge_unitigs_tadpole       |  97.88% |     24004 | 5.25M | 361 |      1019 |   72.3K |  77 |   64.0 | 8.0 |  16.0 | 160.0 |   0:00:54 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |     Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|--------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.97% |     32723 |   1.32M |  71 |      1634 |  8.29K |  84 |   64.0 |  9.0 |  15.3 | 164.0 |   0:00:40 |
| 8_mr_spades  |  98.79% |     64637 |   4.66M | 134 |       446 | 13.74K | 229 |   85.0 | 12.0 |  20.3 | 218.0 |   0:00:40 |
| 8_megahit    |  98.70% |     28403 |    3.8M | 239 |       504 | 31.61K | 353 |   64.0 |  8.0 |  16.0 | 160.0 |   0:00:38 |
| 8_mr_megahit |  98.82% |     51534 |   4.87M | 159 |       416 | 19.42K | 296 |   85.0 | 11.0 |  21.0 | 214.0 |   0:00:39 |
| 8_platanus   |  97.39% |     27954 | 666.32K |  45 |        50 |    956 |  42 |   64.0 |  9.0 |  15.3 | 164.0 |   0:00:39 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 5224283 | 5432652 |   2 |
| Paralogs                 |    2295 |  220468 | 101 |
| Repetitives              |    2461 |  113050 | 173 |
| 7_merge_anchors.anchors  |   25734 | 5293709 | 363 |
| 7_merge_anchors.others   |    1033 |  108715 | 111 |
| glue_anchors             |   28761 | 5292797 | 331 |
| fill_anchors             |   69316 | 5306605 | 136 |
| spades.contig            |  207470 | 5366665 | 154 |
| spades.scaffold          |  285416 | 5367072 | 139 |
| spades.non-contained     |  207470 | 5349445 |  58 |
| mr_spades.contig         |  101571 | 5368070 | 127 |
| mr_spades.scaffold       |  284423 | 5374661 |  67 |
| mr_spades.non-contained  |  101571 | 5360578 | 102 |
| megahit.contig           |   59600 | 5360031 | 204 |
| megahit.non-contained    |   59600 | 5341675 | 159 |
| mr_megahit.contig        |   75019 | 5384488 | 170 |
| mr_megahit.non-contained |   75019 | 5369154 | 138 |
| platanus.contig          |   18988 | 5419998 | 658 |
| platanus.scaffold        |  284634 | 5396847 | 256 |
| platanus.non-contained   |  284634 | 5346766 |  39 |


# *Mycobacterium abscessus* 6G-0125-R

## Mabs_100x: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Mabs_100x/1_genome
cd ~/data/anchr/Mabs_100x/1_genome

cp ~/data/anchr/ref/Mabs/genome.fa .
cp ~/data/anchr/ref/Mabs/paralogs.fa .
cp ~/data/anchr/ref/Mabs/repetitives.fa .

```

## Mabs_100x: download

* Illumina

```shell script
cd ~/data/anchr/Mabs_100x

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/datasets/M_abscessus_MiSeq.tar.gz

# NOT gzipped tar
tar xvf M_abscessus_MiSeq.tar.gz raw/reads_1.fastq
tar xvf M_abscessus_MiSeq.tar.gz raw/reads_2.fastq

cat raw/reads_1.fastq |
    pigz -p 8 -c \
    > R1.fq.gz
cat raw/reads_2.fastq |
    pigz -p 8 -c \
    > R2.fq.gz

rm -fr raw

```

* GAGE-B assemblies

```shell script
cd ~/data/anchr/Mabs_100x

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/M_abscessus_MiSeq.tar.gz

tar xvfz M_abscessus_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz mira_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz sga_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz soap_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz spades_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz velvet_ctg.fasta

```

## Mabs_100x: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Mabs_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 5090491 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --tile --cutoff 5 --cutk 31" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 all" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    \
    --extend

```

## Mabs_100x: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Mabs_100x

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Mabs_100x bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

rm -fr 9_quast_competitor
quast --threads 16 \
    -R 1_genome/genome.fa \
    8_competitor/abyss_ctg.fasta \
    8_competitor/cabog_ctg.fasta \
    8_competitor/mira_ctg.fasta \
    8_competitor/msrca_ctg.fasta \
    8_competitor/sga_ctg.fasta \
    8_competitor/soap_ctg.fasta \
    8_competitor/spades_ctg.fasta \
    8_competitor/velvet_ctg.fasta \
    7_merge_anchors/anchor.merge.fasta \
    7_merge_anchors/others.non-contained.fasta \
    1_genome/paralogs.fa \
    --label "abyss,cabog,mira,msrca,sga,soap,spades,velvet,merge,others,paralogs" \
    -o 9_quast_competitor

```

Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 460.0 |    278 | 2532.2 |                         14.46% |
| R.tadpole.bbtools | 260.8 |    262 |   51.3 |                         69.95% |
| R.genome.picard   | 295.7 |    279 |   47.4 |                             FR |
| R.genome.picard   | 287.1 |    271 |   33.8 |                             RF |
| R.tadpole.picard  | 268.0 |    267 |   49.1 |                             FR |
| R.tadpole.picard  | 251.6 |    255 |   47.9 |                             RF |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 59        | 5595943         | 0.0000       | 59.42   |
| R.31 | 53        | 5711858         | 0.0000       | 59.51   |
| R.41 | 49        | 5611136         | 0.0000       | 59.84   |
| R.51 | 44        | 5640306         | 0.0000       | 60.06   |
| R.61 | 40        | 5660378         | 0.0000       | 60.29   |
| R.71 | 36        | 5682672         | 0.0000       | 60.46   |
| R.81 | 32        | 5667474         | 0.0117       | 60.50   |


Table: statReads

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 5067172 | 5090491 |       2 |
| Paralogs    |    1693 |   83291 |      53 |
| Repetitives |     192 |   15322 |      82 |
| Illumina.R  |     251 |    512M | 2039840 |
| trim.R      |     177 | 280.57M | 1704906 |
| Q20L60      |     178 | 271.98M | 1640727 |
| Q25L60      |     174 | 248.16M | 1545205 |
| Q30L60      |     165 | 204.44M | 1369930 |


Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 511.87M | 2039328 |
| filteredbytile | 251 | 487.06M | 1940498 |
| highpass       | 251 | 473.11M | 1884884 |
| trim           | 177 | 281.53M | 1709470 |
| filter         | 177 | 280.57M | 1704906 |
| R1             | 187 | 149.89M |  852453 |
| R2             | 167 | 130.68M |  852453 |
| Rs             |   0 |       0 |       0 |


```text
#R.trim
#Matched	1392511	73.87781%
#Name	Reads	ReadsPct
Reverse_adapter	721366	38.27111%
pcr_dimer	389747	20.67751%
TruSeq_Universal_Adapter	112545	5.97092%
PCR_Primers	98405	5.22075%
TruSeq_Adapter_Index_1_6	46355	2.45930%
Nextera_LMP_Read2_External_Adapter	14125	0.74938%
TruSeq_Adapter_Index_11	4993	0.26490%
```

```text
#R.filter
#Matched	4563	0.26692%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	4563	0.26692%
```

```text
#R.peaks
#k	31
#unique_kmers	12857688
#error_kmers	7911738
#genomic_kmers	4945950
#main_peak	45
#genome_size_in_peaks	4945242
#genome_size	5229974
#haploid_genome_size	5229974
#fold_coverage	45
#haploid_fold_coverage	45
#ploidy	1
#percent_repeat_in_peaks	0.028
#percent_repeat	1.915
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 177 | 278.96M | 1694298 |
| ecco          | 177 | 278.91M | 1694298 |
| eccc          | 177 | 278.91M | 1694298 |
| ecct          | 176 | 270.05M | 1644118 |
| extended      | 214 | 335.46M | 1644118 |
| merged.raw    | 235 | 189.27M |  813263 |
| unmerged.raw  | 210 |   3.27M |   17592 |
| unmerged.trim | 210 |   3.27M |   17586 |
| M1            | 235 | 181.08M |  778271 |
| U1            | 230 |   1.86M |    8793 |
| U2            | 188 |    1.4M |    8793 |
| Us            |   0 |       0 |       0 |
| M.cor         | 234 | 185.13M | 1574128 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 191.1 |    186 |  46.5 |         90.51% |
| M.ihist.merge.txt  | 232.7 |    226 |  51.5 |         98.93% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   |  55.1 |   44.5 |   19.27% | "45" | 5.09M | 5.15M |     1.01 | 0:00'38'' |
| Q20L60.R |  53.4 |   43.9 |   17.90% | "45" | 5.09M | 5.15M |     1.01 | 0:00'38'' |
| Q25L60.R |  48.8 |   41.5 |   14.89% | "43" | 5.09M | 5.14M |     1.01 | 0:00'36'' |
| Q30L60.R |  40.2 |   35.5 |   11.71% | "39" | 5.09M | 5.13M |     1.01 | 0:00'33'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  96.94% |      5398 | 4.72M | 1174 |       818 | 357.64K | 2521 |   28.0 | 2.0 |   8.0 |  64.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:36 |
| Q0L0X40P000    |   40.0 |  96.21% |      4930 | 4.72M | 1263 |       862 | 397.84K | 2570 |   38.0 | 2.0 |  11.3 |  84.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:39 |
| Q0L0XallP000   |   44.5 |  95.79% |      5488 | 4.84M | 1187 |       853 | 270.01K | 2487 |   42.0 | 3.0 |  12.0 |  96.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:37 |
| Q20L60X30P000  |   30.0 |  97.33% |      5264 |  4.7M | 1160 |       871 | 387.77K | 2591 |   28.0 | 2.0 |   8.0 |  64.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:37 |
| Q20L60X40P000  |   40.0 |  96.44% |      4778 | 4.64M | 1264 |       891 | 388.06K | 2518 |   38.0 | 2.0 |  11.3 |  84.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:38 |
| Q20L60XallP000 |   43.9 |  96.13% |      5687 | 4.82M | 1157 |       834 | 273.77K | 2400 |   42.0 | 3.0 |  12.0 |  96.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:36 |
| Q25L60X30P000  |   30.0 |  98.03% |      5474 | 4.38M | 1106 |       871 | 424.57K | 2289 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:38 |
| Q25L60X40P000  |   40.0 |  97.43% |      5256 | 4.59M | 1175 |       872 | 385.67K | 2383 |   38.0 | 2.0 |  11.3 |  84.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:39 |
| Q25L60XallP000 |   41.5 |  97.36% |      5743 | 4.73M | 1104 |       842 | 300.27K | 2313 |   40.0 | 3.0 |  11.3 |  92.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:39 |
| Q30L60X30P000  |   30.0 |  98.73% |      5582 | 4.24M | 1021 |       854 | 350.84K | 2230 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:40 |
| Q30L60XallP000 |   35.5 |  98.61% |      6431 | 4.52M |  961 |       784 | 304.58K | 2194 |   34.0 | 2.0 |  10.0 |  76.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:42 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  99.01% |     28337 |  4.9M | 305 |       231 | 41.17K | 580 |   29.0 | 1.0 |   9.0 |  62.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:35 |
| MRXallP000 |   36.4 |  98.91% |     27610 | 5.07M | 308 |       220 | 36.09K | 616 |   35.0 | 2.0 |  10.3 |  78.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:38 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  98.30% |      6079 |  4.6M | 1057 |       840 | 337.62K | 2855 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:42 |
| Q0L0X40P000    |   40.0 |  97.88% |      5791 | 4.77M | 1109 |       833 | 370.43K | 2905 |   38.0 | 2.0 |  11.3 |  84.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:43 |
| Q0L0XallP000   |   44.5 |  97.39% |      5968 | 4.87M | 1117 |       775 | 294.36K | 2888 |   43.0 | 3.0 |  12.3 |  98.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:42 |
| Q20L60X30P000  |   30.0 |  98.28% |      5749 | 4.69M | 1110 |       765 | 343.84K | 2925 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:40 |
| Q20L60X40P000  |   40.0 |  97.81% |      6741 | 4.86M | 1024 |       728 | 248.92K | 2867 |   38.0 | 3.0 |  10.7 |  88.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:44 |
| Q20L60XallP000 |   43.9 |  97.54% |      6515 | 4.88M | 1034 |       766 | 251.62K | 2824 |   42.0 | 3.0 |  12.0 |  96.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:43 |
| Q25L60X30P000  |   30.0 |  98.39% |      5961 | 4.57M | 1067 |       831 |  369.8K | 2929 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:41 |
| Q25L60X40P000  |   40.0 |  98.22% |      5832 | 4.65M | 1083 |       873 | 395.27K | 2903 |   38.0 | 2.0 |  11.3 |  84.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:44 |
| Q25L60XallP000 |   41.5 |  98.13% |      6283 |  4.8M | 1027 |       783 |  297.2K | 2871 |   40.0 | 3.0 |  11.3 |  92.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:43 |
| Q30L60X30P000  |   30.0 |  98.62% |      5811 | 4.65M | 1084 |       757 | 322.24K | 2919 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:41 |
| Q30L60XallP000 |   35.5 |  98.70% |      6713 | 4.71M |  983 |       771 | 292.82K | 2919 |   34.0 | 2.0 |  10.0 |  76.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:44 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  99.34% |     37727 | 4.81M | 242 |       195 | 35.09K | 501 |   29.0 | 1.0 |   9.0 |  62.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:39 |
| MRXallP000 |   36.4 |  99.24% |     38610 | 4.71M | 227 |       338 | 28.64K | 422 |   35.0 | 1.0 |  11.0 |  74.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:37 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  98.81% |      6887 | 4.08M | 849 |       883 | 287.82K | 1916 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:40 |
| Q0L0X40P000    |   40.0 |  98.55% |      7015 |  4.5M | 898 |       820 |  272.4K | 1988 |   39.0 | 2.0 |  11.7 |  86.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:43 |
| Q0L0XallP000   |   44.5 |  98.25% |      7236 | 4.69M | 929 |       781 | 261.78K | 2089 |   43.0 | 2.0 |  13.0 |  94.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:41 |
| Q20L60X30P000  |   30.0 |  98.88% |      6737 | 4.08M | 849 |       833 | 303.02K | 2045 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:42 |
| Q20L60X40P000  |   40.0 |  98.52% |      6541 | 4.27M | 903 |       789 | 290.48K | 2039 |   38.0 | 2.0 |  11.3 |  84.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:42 |
| Q20L60XallP000 |   43.9 |  98.29% |      6943 | 4.66M | 961 |       789 | 268.23K | 2072 |   42.0 | 2.0 |  12.7 |  92.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:41 |
| Q25L60X30P000  |   30.0 |  98.92% |      6800 |    4M | 838 |       852 | 319.53K | 2026 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:40 |
| Q25L60X40P000  |   40.0 |  98.74% |      6596 | 4.32M | 901 |       854 | 303.62K | 2091 |   38.0 | 2.0 |  11.3 |  84.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:43 |
| Q25L60XallP000 |   41.5 |  98.70% |      7261 | 4.58M | 890 |       791 | 259.81K | 2066 |   40.0 | 2.0 |  12.0 |  88.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:40 |
| Q30L60X30P000  |   30.0 |  99.13% |      6192 | 4.19M | 907 |       793 | 297.89K | 2219 |   29.0 | 2.0 |   8.3 |  66.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:42 |
| Q30L60XallP000 |   35.5 |  99.13% |      7536 | 4.54M | 857 |       771 | 270.78K | 2287 |   34.0 | 2.0 |  10.0 |  76.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:43 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  99.28% |     36670 | 4.51M | 224 |       402 | 27.74K | 353 |   29.0 | 1.0 |   9.0 |  62.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:36 |
| MRXallP000 |   36.4 |  99.28% |     43472 | 4.57M | 210 |       467 | 25.07K | 355 |   35.0 | 1.0 |  11.0 |  74.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |   0.00% |     56936 | 5.27M | 186 |      1410 |  247.2K | 193 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_mr_unitigs_bcalm      |   0.00% |     44166 | 4.92M | 207 |       934 |  17.16K |  18 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_mr_unitigs_superreads |   0.00% |     31629 |  5.1M | 266 |       911 |  21.14K |  21 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_mr_unitigs_tadpole    |   0.00% |     50753 | 4.72M | 194 |      1017 |  13.41K |  12 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_unitigs_bcalm         |   0.00% |     16327 | 5.15M | 524 |      1194 | 406.18K | 360 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_unitigs_superreads    |   0.00% |     11838 | 5.13M | 670 |      1183 | 515.21K | 461 |    0.0 | 0.0 |   0.0 |   0.0 |           |
| 7_merge_unitigs_tadpole       |   0.00% |     17015 | 5.11M | 510 |      1251 | 412.26K | 350 |    0.0 | 0.0 |   0.0 |   0.0 |           |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |     Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|--------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|----------:|
| 8_spades     |  99.32% |      4992 | 177.41K |  46 |      1161 | 21.48K |  64 |   43.0 | 1.0 |  13.7 |  90.0 |   0:00:36 |
| 8_mr_spades  |  99.35% |     45608 |   3.35M | 138 |       752 | 17.47K | 186 |   35.0 | 1.0 |  11.0 |  74.0 |   0:00:35 |
| 8_megahit    |  99.20% |     11326 | 748.29K | 105 |      1097 | 21.74K | 118 |   43.0 | 2.0 |  13.0 |  94.0 |   0:00:34 |
| 8_mr_megahit |  99.61% |     41860 |   3.44M | 152 |       752 | 16.86K | 184 |   35.0 | 1.0 |  11.0 |  74.0 |   0:00:36 |
| 8_platanus   |  98.30% |      9061 |   3.29M | 525 |       835 | 94.77K | 681 |   43.0 | 2.0 |  13.0 |  94.0 |   0:00:34 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 5067172 | 5090491 |   2 |
| Paralogs                 |    1693 |   83291 |  53 |
| Repetitives              |     192 |   15322 |  82 |
| 7_merge_anchors.anchors  |   56936 | 5270926 | 186 |
| 7_merge_anchors.others   |    1410 |  247200 | 193 |
| glue_anchors             |   69408 | 5265822 | 146 |
| fill_anchors             |  104416 | 5292049 | 105 |
| spades.contig            |  147271 | 5135031 |  91 |
| spades.scaffold          |  147271 | 5135221 |  87 |
| spades.non-contained     |  147271 | 5123159 |  64 |
| mr_spades.contig         |  115558 | 5134632 | 112 |
| mr_spades.scaffold       |  115585 | 5134841 | 109 |
| mr_spades.non-contained  |  115558 | 5121206 |  85 |
| megahit.contig           |  107304 | 5138484 | 143 |
| megahit.non-contained    |  107304 | 5118115 |  91 |
| mr_megahit.contig        |   95831 | 5140708 | 106 |
| mr_megahit.non-contained |   95831 | 5133066 |  91 |
| platanus.contig          |   16538 | 5191588 | 892 |
| platanus.scaffold        |   25843 | 5130043 | 407 |
| platanus.non-contained   |   25843 | 5098283 | 306 |


# *Rhodobacter sphaeroides* 2.4.1

## Rsph_100x: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Rsph_100x/1_genome
cd ~/data/anchr/Rsph_100x/1_genome

cp ~/data/anchr/ref/Rsph/genome.fa .
cp ~/data/anchr/ref/Rsph/paralogs.fa .
cp ~/data/anchr/ref/Rsph/repetitives.fa .

```

## Rsph_100x: download

* Illumina

```shell script
cd ~/data/anchr/Rsph_100x

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/datasets/R_sphaeroides_MiSeq.tar.gz

# NOT gzipped tar
tar xvf R_sphaeroides_MiSeq.tar.gz raw/insert_540_1__cov100x.fastq
tar xvf R_sphaeroides_MiSeq.tar.gz raw/insert_540_2__cov100x.fastq

cat raw/insert_540_1__cov100x.fastq |
    pigz -p 8 -c \
    > R1.fq.gz
cat raw/insert_540_2__cov100x.fastq |
    pigz -p 8 -c \
    > R2.fq.gz

rm -fr raw

```

* GAGE-B assemblies

```shell script
cd ~/data/anchr/Rsph_100x

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/R_sphaeroides_MiSeq.tar.gz

tar xvfz R_sphaeroides_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz mira_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz sga_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz soap_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz spades_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz velvet_ctg.fasta

```


## Rsph_100x: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4602977 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 3" \
    \
    --cov "30 all" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

```

## Rsph_100x: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_100x

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Rsph_100x bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

rm -fr 9_quast_competitor
quast --no-check --threads 16 \
    -R 1_genome/genome.fa \
    8_competitor/abyss_ctg.fasta \
    8_competitor/cabog_ctg.fasta \
    8_competitor/mira_ctg.fasta \
    8_competitor/msrca_ctg.fasta \
    8_competitor/sga_ctg.fasta \
    8_competitor/soap_ctg.fasta \
    8_competitor/spades_ctg.fasta \
    8_competitor/velvet_ctg.fasta \
    7_merge_anchors/anchor.merge.fasta \
    7_merge_anchors/others.non-contained.fasta \
    1_genome/paralogs.fa \
    --label "abyss,cabog,mira,msrca,sga,soap,spades,velvet,merge,others,paralogs" \
    -o 9_quast_competitor

```

Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 443.0 |    422 | 1014.7 |                         31.69% |
| R.tadpole.bbtools | 408.5 |    420 |  107.0 |                         66.01% |
| R.genome.picard   | 412.9 |    422 |   39.3 |                             FR |
| R.tadpole.picard  | 408.5 |    421 |   46.7 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 46        | 5066911         | 0.2691       | 67.75   |
| R.31 | 41        | 5000728         | 0.1951       | 68.11   |
| R.41 | 37        | 4992140         | 0.1372       | 68.22   |
| R.51 | 33        | 4997930         | 0.1016       | 68.27   |
| R.61 | 30        | 5016426         | 0.0838       | 68.29   |
| R.71 | 26        | 5019864         | 0.0737       | 68.30   |
| R.81 | 23        | 4975776         | 0.0812       | 68.30   |


Table: statReads

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 3188524 | 4602977 |       7 |
| Paralogs    |    2337 |  146789 |      66 |
| Repetitives |     572 |   57281 |     165 |
| Illumina.R  |     251 |  451.8M | 1800000 |
| trim.R      |       0 |       0 |       0 |
| Q20L60      |       0 |       0 |       0 |
| Q25L60      |       0 |       0 |       0 |
| Q30L60      |       0 |       0 |       0 |


Table: statTrimReads

| Name     | N50 |     Sum |       # |
|:---------|----:|--------:|--------:|
| clumpify | 251 | 447.53M | 1782988 |
| trim     | 148 |  200.1M | 1452702 |
| filter   | 148 |  200.1M | 1452702 |
| R1       | 164 | 100.14M |  655190 |
| R2       | 133 |  81.61M |  655190 |
| Rs       | 141 |  18.34M |  142322 |


```text
#R.trim
#Matched	113970	6.39208%
#Name	Reads	ReadsPct
Reverse_adapter	81598	4.57647%
pcr_dimer	14481	0.81218%
PCR_Primers	8081	0.45323%
TruSeq_Universal_Adapter	5665	0.31773%
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	8019347
#error_kmers	3640676
#genomic_kmers	4378671
#main_peak	32
#genome_size_in_peaks	4616686
#genome_size	5068031
#haploid_genome_size	5068031
#fold_coverage	32
#haploid_fold_coverage	32
#ploidy	1
#percent_repeat_in_peaks	5.610
#percent_repeat	11.140
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 148 | 200.09M | 1452575 |
| ecco          | 148 | 199.85M | 1452574 |
| ecct          | 148 | 198.72M | 1443939 |
| extended      | 186 | 255.79M | 1443939 |
| merged.raw    | 455 | 197.36M |  475339 |
| unmerged.raw  | 171 |   80.1M |  493260 |
| unmerged.trim | 171 |  80.08M |  492922 |
| M1            | 455 | 197.19M |  474957 |
| U1            | 172 |  19.95M |  123250 |
| U2            | 151 |  17.79M |  123250 |
| Us            | 182 |  42.34M |  246422 |
| M.cor         | 443 | 277.99M | 1689258 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 184.4 |    179 |  66.0 |         10.44% |
| M.ihist.merge.txt  | 415.2 |    452 |  89.1 |         65.84% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   |  43.5 |   38.7 |   11.10% | "37" |  4.6M | 4.55M |     0.99 | 0:00'34'' |
| Q20L60.R |  42.1 |   37.9 |    9.98% | "37" |  4.6M | 4.55M |     0.99 | 0:00'31'' |
| Q25L60.R |  36.8 |   34.9 |    5.03% | "35" |  4.6M | 4.54M |     0.99 | 0:00'30'' |
| Q30L60.R |  27.2 |   26.6 |    2.20% | "31" |  4.6M | 4.52M |     0.98 | 0:00'26'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  97.74% |     10118 | 3.67M | 582 |      2020 | 495.69K | 1503 |   27.0 | 2.0 |   7.7 |  62.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:39 |
| Q0L0XallP000   |   38.7 |  97.71% |     11561 | 3.87M | 546 |      2020 | 419.45K | 1413 |   35.0 | 3.0 |   9.7 |  82.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:39 |
| Q20L60X30P000  |   30.0 |  97.71% |      9712 | 3.83M | 634 |      1985 | 415.68K | 1441 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:39 |
| Q20L60XallP000 |   37.9 |  97.78% |     10648 | 3.78M | 571 |      2464 | 409.31K | 1321 |   35.0 | 3.0 |   9.7 |  82.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:37 |
| Q25L60X30P000  |   30.0 |  98.30% |      8044 | 3.45M | 658 |      2795 | 491.89K | 1504 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:37 |
| Q25L60XallP000 |   34.9 |  98.36% |      8774 |  3.7M | 638 |      3317 | 492.11K | 1433 |   32.0 | 3.0 |   8.7 |  76.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:37 |
| Q30L60XallP000 |   26.6 |  97.60% |      5721 | 3.68M | 918 |      2263 | 714.18K | 1981 |   25.0 | 3.0 |   6.3 |  62.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:34 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.39% |     16212 | 4.03M | 420 |      5224 | 250.42K | 862 |   27.0 | 2.0 |   7.7 |  62.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:32 |
| MRX30P001  |   30.0 |  97.46% |     16306 | 4.05M | 429 |      7573 | 269.31K | 846 |   27.0 | 2.0 |   7.7 |  62.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:33 |
| MRXallP000 |   60.4 |  97.34% |     18753 | 4.17M | 377 |      7573 | 206.69K | 801 |   55.0 | 5.0 |  15.0 | 130.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:34 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  97.28% |      8627 | 4.08M |  719 |      2421 | 457.07K | 2114 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:37 |
| Q0L0XallP000   |   38.7 |  97.53% |     10601 | 4.14M |  625 |      1885 | 394.71K | 1959 |   36.0 | 4.0 |   9.3 |  88.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:40 |
| Q20L60X30P000  |   30.0 |  97.07% |      8171 | 4.06M |  753 |      2815 | 500.95K | 2039 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:37 |
| Q20L60XallP000 |   37.9 |  97.55% |      9748 | 4.03M |  655 |      2815 | 520.14K | 1892 |   35.0 | 3.0 |   9.7 |  82.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:37 |
| Q25L60X30P000  |   30.0 |  97.08% |      6598 | 3.95M |  880 |      3212 | 563.22K | 2262 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:35 |
| Q25L60XallP000 |   34.9 |  97.60% |      7409 | 4.02M |  791 |      3598 | 551.72K | 2102 |   33.0 | 4.0 |   8.3 |  82.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:37 |
| Q30L60XallP000 |   26.6 |  93.56% |      3959 | 3.69M | 1149 |      2347 | 624.71K | 2814 |   25.0 | 3.0 |   6.3 |  62.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:33 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.36% |     15973 | 4.06M | 439 |      5247 | 247.15K |  988 |   27.0 | 2.0 |   7.7 |  62.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:32 |
| MRX30P001  |   30.0 |  97.37% |     16187 | 4.06M | 449 |      6101 | 277.33K | 1027 |   27.0 | 2.0 |   7.7 |  62.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:33 |
| MRXallP000 |   60.4 |  97.45% |     16781 | 4.18M | 412 |     11860 | 241.05K |  867 |   55.0 | 4.0 |  15.7 | 126.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:34 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  98.27% |     10049 | 4.02M |  631 |      3089 | 477.52K | 1673 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:38 |
| Q0L0XallP000   |   38.7 |  98.43% |     12008 | 4.08M |  560 |      2966 | 410.03K | 1521 |   36.0 | 4.0 |   9.3 |  88.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:40 |
| Q20L60X30P000  |   30.0 |  98.19% |      9524 | 3.98M |  671 |      2837 | 478.77K | 1713 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:37 |
| Q20L60XallP000 |   37.9 |  98.38% |     10305 | 3.87M |  583 |      3394 | 504.59K | 1545 |   35.0 | 3.0 |   9.7 |  82.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:40 |
| Q25L60X30P000  |   30.0 |  97.92% |      7316 | 3.87M |  794 |      3417 | 556.13K | 1926 |   28.0 | 3.0 |   7.3 |  68.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:36 |
| Q25L60XallP000 |   34.9 |  98.14% |      8009 |    4M |  736 |      4316 | 512.31K | 1780 |   33.0 | 4.0 |   8.3 |  82.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q30L60XallP000 |   26.6 |  96.22% |      4883 | 3.77M | 1018 |      2300 | 669.59K | 2439 |   25.0 | 3.0 |   6.3 |  62.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:33 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.50% |     16212 | 4.03M | 419 |      5247 | 253.17K | 902 |   27.0 | 2.0 |   7.7 |  62.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:35 |
| MRX30P001  |   30.0 |  97.50% |     16306 | 4.02M | 423 |      6101 | 273.78K | 916 |   27.0 | 2.0 |   7.7 |  62.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:33 |
| MRXallP000 |   60.4 |  97.59% |     17622 | 4.15M | 397 |     11860 | 238.64K | 825 |   55.0 | 4.0 |  15.7 | 126.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  87.16% |     12877 | 4.18M | 554 |      2439 | 973.23K | 490 |   35.0 | 3.0 |   9.7 |  82.0 |   0:00:41 |
| 7_merge_mr_unitigs_bcalm      |  88.11% |     12488 | 4.02M | 528 |      4249 | 371.26K | 154 |   35.0 | 3.0 |   9.7 |  82.0 |   0:00:41 |
| 7_merge_mr_unitigs_superreads |  87.03% |     10972 | 3.78M | 567 |      4204 | 405.32K | 163 |   35.0 | 2.0 |  10.3 |  78.0 |   0:00:38 |
| 7_merge_mr_unitigs_tadpole    |  87.19% |     10883 | 3.76M | 561 |      4201 | 411.28K | 174 |   35.0 | 2.0 |  10.3 |  78.0 |   0:00:40 |
| 7_merge_unitigs_bcalm         |  90.13% |     11210 | 4.16M | 620 |      2795 | 779.15K | 387 |   35.0 | 3.0 |   9.7 |  82.0 |   0:00:42 |
| 7_merge_unitigs_superreads    |  89.03% |     12466 | 4.11M | 560 |      2572 | 810.27K | 404 |   35.0 | 3.0 |   9.7 |  82.0 |   0:00:43 |
| 7_merge_unitigs_tadpole       |  89.63% |     12036 | 4.16M | 577 |      2664 | 719.09K | 371 |   35.0 | 3.0 |   9.7 |  82.0 |   0:00:42 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |     Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|--------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 8_spades     |  99.20% |     11574 | 368.63K |  49 |     37282 | 174.51K |  77 |   35.0 | 2.0 |  10.3 |  78.0 |   0:00:33 |
| 8_mr_spades  |  99.16% |     21420 |   2.57M | 206 |     19232 | 207.88K | 286 |   55.0 | 3.5 |  16.0 | 124.0 |   0:00:33 |
| 8_megahit    |  98.52% |     12400 |   1.64M | 220 |     11870 | 247.84K | 300 |   35.0 | 3.0 |   9.7 |  82.0 |   0:00:32 |
| 8_mr_megahit |  99.29% |     20129 |   4.08M | 342 |     16123 | 230.05K | 571 |   55.0 | 4.0 |  15.7 | 126.0 |   0:00:34 |
| 8_platanus   |  96.34% |     12457 | 735.28K | 102 |      6042 | 194.06K | 174 |   35.0 | 3.0 |   9.7 |  82.0 |   0:00:32 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 3188524 | 4602977 |    7 |
| Paralogs                 |    2337 |  146789 |   66 |
| Repetitives              |     572 |   57281 |  165 |
| 7_merge_anchors.anchors  |   12877 | 4179148 |  554 |
| 7_merge_anchors.others   |    2439 |  973227 |  490 |
| glue_anchors             |   14774 | 4175948 |  517 |
| fill_anchors             |   37984 | 4197543 |  210 |
| spades.contig            |  150729 | 4576779 |  136 |
| spades.scaffold          |  172916 | 4577123 |  131 |
| spades.non-contained     |  150729 | 4562257 |   71 |
| mr_spades.contig         |   55603 | 4565108 |  165 |
| mr_spades.scaffold       |   93957 | 4566245 |  117 |
| mr_spades.non-contained  |   55603 | 4555295 |  148 |
| megahit.contig           |   52830 | 4573210 |  248 |
| megahit.non-contained    |   52830 | 4540669 |  182 |
| mr_megahit.contig        |   32000 | 4575220 |  277 |
| mr_megahit.non-contained |   32000 | 4562392 |  252 |
| platanus.contig          |   15555 | 4617410 | 1657 |
| platanus.scaffold        |   85196 | 4561389 |  574 |
| platanus.non-contained   |   89576 | 4473094 |  135 |


# *Vibrio cholerae* CP1032(5)

* *Vibrio cholerae* O1 biovar El Tor str. N16961
  * Taxid: [243277](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=243277)
  * Assembly: [GCF_000006745.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000006745.1)
  * Proportion of paralogs (> 1000 bp): 0.0216
* *Vibrio cholerae* CP1032(5)
  * Assembly: [GCF_000279305.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000279305.1)

## Vcho_100x: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Vcho_100x/1_genome
cd ~/data/anchr/Vcho_100x/1_genome

cp ~/data/anchr/ref/Vcho/genome.fa .
cp ~/data/anchr/ref/Vcho/paralogs.fa .
cp ~/data/anchr/ref/Vcho/repetitives.fa .

```

## Vcho_100x: download

* Illumina

```shell script
cd ~/data/anchr/Vcho_100x

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/datasets/V_cholerae_MiSeq.tar.gz

# NOT gzipped tar
tar xvf V_cholerae_MiSeq.tar.gz raw/reads_1.fastq
tar xvf V_cholerae_MiSeq.tar.gz raw/reads_2.fastq

cat raw/reads_1.fastq |
    pigz -p 8 -c \
    > R1.fq.gz
cat raw/reads_2.fastq |
    pigz -p 8 -c \
    > R2.fq.gz

rm -fr raw

```

* GAGE-B assemblies

```shell script
cd ~/data/anchr/Vcho_100x

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/V_cholerae_MiSeq.tar.gz

tar xvfz V_cholerae_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz mira_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz sga_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz soap_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz spades_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz velvet_ctg.fasta

```


## Vcho_100x: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Vcho_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4033464 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --tile --cutoff 5 --cutk 31" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 50 all" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

```

## Vcho_100x: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Vcho_100x

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

#bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

rm -fr 9_quast_competitor
quast --threads 16 \
    -R 1_genome/genome.fa \
    8_competitor/abyss_ctg.fasta \
    8_competitor/cabog_ctg.fasta \
    8_competitor/mira_ctg.fasta \
    8_competitor/msrca_ctg.fasta \
    8_competitor/sga_ctg.fasta \
    8_competitor/soap_ctg.fasta \
    8_competitor/spades_ctg.fasta \
    8_competitor/velvet_ctg.fasta \
    7_merge_anchors/anchor.merge.fasta \
    7_merge_anchors/others.non-contained.fasta \
    1_genome/paralogs.fa \
    --label "abyss,cabog,mira,msrca,sga,soap,spades,velvet,merge,others,paralogs" \
    -o 9_quast_competitor

```

Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 395.7 |    275 | 1929.6 |                         16.74% |
| R.tadpole.bbtools | 268.3 |    266 |   53.1 |                         82.27% |
| R.genome.picard   | 294.0 |    277 |   48.0 |                             FR |
| R.genome.picard   | 280.2 |    268 |   29.0 |                             RF |
| R.tadpole.picard  | 271.9 |    268 |   48.1 |                             FR |
| R.tadpole.picard  | 260.7 |    262 |   44.9 |                             RF |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 66        | 4104185         | 0.0297       | 46.75   |
| R.31 | 61        | 4036245         | 0.1260       | 47.56   |
| R.41 | 56        | 4036627         | 0.0665       | 48.36   |
| R.51 | 52        | 4033882         | 0.0306       | 48.69   |
| R.61 | 47        | 4035875         | 0.0405       | 48.94   |
| R.71 | 43        | 4034383         | 0.0347       | 49.07   |
| R.81 | 39        | 4032895         | 0.0130       | 48.98   |


Table: statReads

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 2961149 | 4033464 |       2 |
| Paralogs    |    3424 |  119270 |      49 |
| Repetitives |    1070 |  120471 |     244 |
| Illumina.R  |     251 |    400M | 1593624 |
| trim.R      |     189 | 260.67M | 1438624 |
| Q20L60      |     189 | 256.24M | 1410440 |
| Q25L60      |     187 | 242.89M | 1359601 |
| Q30L60      |     181 | 215.72M | 1260864 |


Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 397.98M | 1585572 |
| filteredbytile | 251 | 377.77M | 1505070 |
| highpass       | 251 |  373.5M | 1488058 |
| trim           | 189 | 262.28M | 1446258 |
| filter         | 189 | 260.67M | 1438624 |
| R1             | 193 | 133.98M |  719312 |
| R2             | 184 | 126.68M |  719312 |
| Rs             |   0 |       0 |       0 |


```text
#R.trim
#Matched	1214411	81.61046%
#Name	Reads	ReadsPct
Reverse_adapter	588085	39.52030%
pcr_dimer	339985	22.84756%
PCR_Primers	174574	11.73167%
TruSeq_Adapter_Index_1_6	44985	3.02307%
TruSeq_Universal_Adapter	44820	3.01198%
Nextera_LMP_Read2_External_Adapter	18442	1.23933%
```

```text
#R.filter
#Matched	7632	0.52771%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	7632	0.52771%
```

```text
#R.peaks
#k	31
#unique_kmers	9683593
#error_kmers	5963436
#genomic_kmers	3720157
#main_peak	57
#genome_size_in_peaks	3782096
#genome_size	3928244
#haploid_genome_size	3928244
#fold_coverage	57
#haploid_fold_coverage	57
#ploidy	1
#percent_repeat_in_peaks	1.732
#percent_repeat	2.918
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 189 | 254.95M | 1403642 |
| ecco          | 189 | 254.92M | 1403642 |
| eccc          | 189 | 254.92M | 1403642 |
| ecct          | 189 | 252.86M | 1392188 |
| extended      | 228 | 308.31M | 1392188 |
| merged.raw    | 239 | 165.52M |  691254 |
| unmerged.raw  | 227 |   2.02M |    9680 |
| unmerged.trim | 227 |   2.02M |    9678 |
| M1            | 239 | 153.77M |  642588 |
| U1            | 240 |    1.1M |    4839 |
| U2            | 215 | 922.12K |    4839 |
| Us            |   0 |       0 |       0 |
| M.cor         | 238 | 156.43M | 1294854 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 197.6 |    192 |  44.6 |         91.25% |
| M.ihist.merge.txt  | 239.4 |    232 |  51.5 |         99.31% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   |  64.6 |   54.5 |   15.65% | "109" | 4.03M | 3.94M |     0.98 | 0:00'37'' |
| Q20L60.R |  63.5 |   54.2 |   14.71% | "111" | 4.03M | 3.92M |     0.97 | 0:00'37'' |
| Q25L60.R |  60.2 |   52.8 |   12.25% | "109" | 4.03M | 3.92M |     0.97 | 0:00'37'' |
| Q30L60.R |  53.5 |   48.4 |    9.55% | "105" | 4.03M | 3.92M |     0.97 | 0:00'33'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  93.60% |      7053 | 3.57M | 732 |      1045 | 215.89K | 1458 |   39.0 | 5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:34 |
| Q0L0X50P000    |   50.0 |  92.89% |      6346 |  3.6M | 797 |      1011 | 205.12K | 1538 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:36 |
| Q0L0XallP000   |   54.5 |  92.65% |      6309 | 3.62M | 803 |      1007 | 187.99K | 1542 |   54.0 | 7.0 |  13.3 | 136.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:35 |
| Q20L60X40P000  |   40.0 |  96.05% |     14275 | 3.21M | 391 |      1049 | 143.48K |  777 |   40.0 | 5.0 |  10.0 | 100.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:34 |
| Q20L60X50P000  |   50.0 |  95.47% |     12647 | 3.41M | 428 |      1068 | 145.35K |  841 |   50.0 | 7.0 |  12.0 | 128.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:35 |
| Q20L60XallP000 |   54.2 |  95.25% |     10755 | 3.41M | 464 |      1068 | 154.31K |  856 |   54.0 | 7.0 |  13.3 | 136.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:35 |
| Q25L60X40P000  |   40.0 |  96.30% |     13867 | 3.07M | 359 |      1191 | 139.25K |  668 |   40.0 | 6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:35 |
| Q25L60X50P000  |   50.0 |  96.12% |     14191 | 3.54M | 400 |      1137 | 144.71K |  763 |   50.0 | 7.0 |  12.0 | 128.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:35 |
| Q25L60XallP000 |   52.8 |  95.90% |     13983 | 3.46M | 413 |      1116 | 142.59K |  786 |   53.0 | 7.0 |  13.0 | 134.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:34 |
| Q30L60X40P000  |   40.0 |  96.91% |     14632 | 3.12M | 365 |      1166 | 159.11K |  757 |   40.0 | 6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:37 |
| Q30L60XallP000 |   48.4 |  96.47% |     14960 | 3.41M | 380 |      1118 | 138.33K |  729 |   48.0 | 6.0 |  12.0 | 120.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:34 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   38.8 |  96.72% |     31324 | 3.62M | 217 |      1043 | 80.71K | 412 |   39.0 | 6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:30 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  95.72% |      9378 | 3.61M | 577 |      1022 | 208.46K | 1449 |   40.0 | 5.0 |  10.0 | 100.0 | "31,41,51,61,71,81" |   0:01:53 |   0:00:44 |
| Q0L0X50P000    |   50.0 |  94.58% |      7873 | 3.63M | 670 |       978 | 214.82K | 1491 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:01:53 |   0:00:43 |
| Q0L0XallP000   |   54.5 |  94.19% |      7604 | 3.71M | 704 |      1000 | 184.65K | 1557 |   53.0 | 7.0 |  13.0 | 134.0 | "31,41,51,61,71,81" |   0:01:42 |   0:00:33 |
| Q20L60X40P000  |   40.0 |  97.35% |     17535 | 3.54M | 357 |      1002 | 159.24K |  986 |   40.0 | 5.0 |  10.0 | 100.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:37 |
| Q20L60X50P000  |   50.0 |  97.05% |     18390 | 3.48M | 328 |       998 | 146.38K |  845 |   50.0 | 7.0 |  12.0 | 128.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:37 |
| Q20L60XallP000 |   54.2 |  96.91% |     19916 | 3.64M | 333 |       965 | 141.03K |  769 |   54.0 | 8.0 |  12.7 | 140.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:38 |
| Q25L60X40P000  |   40.0 |  97.53% |     17539 | 3.42M | 329 |       859 |  141.1K |  904 |   40.0 | 5.0 |  10.0 | 100.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:38 |
| Q25L60X50P000  |   50.0 |  97.36% |     21909 | 3.55M | 308 |      1010 | 136.38K |  762 |   50.0 | 7.0 |  12.0 | 128.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:38 |
| Q25L60XallP000 |   52.8 |  97.25% |     20216 | 3.49M | 311 |      1044 | 139.31K |  731 |   53.0 | 8.0 |  12.3 | 138.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:38 |
| Q30L60X40P000  |   40.0 |  97.67% |     17442 | 3.39M | 340 |      1014 | 142.68K |  929 |   40.0 | 6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:40 |
| Q30L60XallP000 |   48.4 |  97.66% |     20308 | 3.49M | 308 |      1045 | 147.54K |  807 |   49.0 | 7.0 |  11.7 | 126.0 | "31,41,51,61,71,81" |   0:02:04 |   0:00:38 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |  Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|-----:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   38.8 |  96.97% |     36539 | 3.6M | 190 |      1010 | 71.66K | 425 |   39.0 | 7.5 |   8.0 | 108.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:37 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  96.43% |     12052 |  3.5M | 477 |      1046 | 188.69K | 1031 |   40.0 | 6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:33 |
| Q0L0X50P000    |   50.0 |  95.76% |      9696 | 3.54M | 552 |      1041 | 190.86K | 1090 |   49.0 | 6.0 |  12.3 | 122.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:34 |
| Q0L0XallP000   |   54.5 |  95.49% |      9417 | 3.57M | 585 |      1032 | 182.38K | 1149 |   54.0 | 7.0 |  13.3 | 136.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:32 |
| Q20L60X40P000  |   40.0 |  97.65% |     19379 | 3.27M | 314 |      1122 | 160.69K |  700 |   40.0 | 6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:37 |
| Q20L60X50P000  |   50.0 |  97.39% |     19985 | 3.14M | 275 |      1046 | 124.51K |  588 |   50.0 | 7.0 |  12.0 | 128.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:34 |
| Q20L60XallP000 |   54.2 |  97.33% |     21292 | 3.41M | 296 |      1046 | 125.51K |  575 |   54.0 | 8.0 |  12.7 | 140.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:35 |
| Q25L60X40P000  |   40.0 |  97.54% |     17757 | 3.06M | 292 |      1120 | 143.81K |  673 |   40.0 | 6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:34 |
| Q25L60X50P000  |   50.0 |  97.57% |     21951 | 3.24M | 278 |      1046 | 135.05K |  601 |   50.0 | 8.0 |  11.3 | 132.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:35 |
| Q25L60XallP000 |   52.8 |  97.45% |     20004 | 3.06M | 273 |      1121 | 129.57K |  554 |   53.0 | 8.0 |  12.3 | 138.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:35 |
| Q30L60X40P000  |   40.0 |  97.73% |     17528 | 3.09M | 309 |      1027 | 161.91K |  714 |   40.0 | 6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:36 |
| Q30L60XallP000 |   48.4 |  97.63% |     23684 | 3.36M | 254 |      1183 | 128.71K |  609 |   49.0 | 8.0 |  11.0 | 130.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:34 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   38.8 |  96.74% |     35883 | 3.41M | 189 |      1043 | 75.51K | 346 |   39.0 | 7.0 |   8.3 | 106.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  96.42% |     21645 | 3.75M | 300 |      1503 | 335.04K | 245 |   55.0 | 7.5 |  13.3 | 140.0 |   0:00:39 |
| 7_merge_mr_unitigs_bcalm      |  90.52% |     22047 | 3.54M | 272 |      1205 |  87.86K |  75 |   55.0 | 9.0 |  12.3 | 146.0 |   0:00:34 |
| 7_merge_mr_unitigs_superreads |  91.07% |     21954 | 3.57M | 278 |      1177 |  98.27K |  88 |   55.0 | 8.0 |  13.0 | 142.0 |   0:00:31 |
| 7_merge_mr_unitigs_tadpole    |  86.08% |     22601 | 3.36M | 251 |      1153 |   86.2K |  78 |   55.0 | 9.0 |  12.3 | 146.0 |   0:00:32 |
| 7_merge_unitigs_bcalm         |  97.07% |     19954 | 3.75M | 322 |      1286 |  269.1K | 219 |   55.0 | 7.0 |  13.7 | 138.0 |   0:00:44 |
| 7_merge_unitigs_superreads    |  96.94% |     18363 | 3.74M | 333 |      1491 |    271K | 201 |   55.0 | 8.0 |  13.0 | 142.0 |   0:00:46 |
| 7_merge_unitigs_tadpole       |  97.00% |     22251 | 3.77M | 293 |      1429 |  253.7K | 194 |   55.0 | 8.0 |  13.0 | 142.0 |   0:00:42 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.40% |     53389 | 1.73M |  99 |      1739 | 68.54K | 162 |   55.0 | 11.0 |  11.0 | 154.0 |   0:00:46 |
| 8_mr_spades  |  98.38% |     58758 | 3.44M | 161 |      1118 | 56.54K | 266 |   39.0 |  8.0 |   7.7 | 110.0 |   0:00:40 |
| 8_megahit    |  97.88% |     27853 | 2.16M | 155 |      1215 | 76.46K | 201 |   55.0 |  9.0 |  12.3 | 146.0 |   0:00:32 |
| 8_mr_megahit |  99.12% |     52927 | 3.08M | 143 |      1123 | 58.34K | 232 |   39.0 |  8.0 |   7.7 | 110.0 |   0:00:34 |
| 8_platanus   |  96.59% |     28111 | 3.69M | 237 |      1051 | 76.16K | 356 |   55.0 | 10.0 |  11.7 | 150.0 |   0:00:34 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 2961149 | 4033464 |   2 |
| Paralogs                 |    3424 |  119270 |  49 |
| Repetitives              |    1070 |  120471 | 244 |
| 7_merge_anchors.anchors  |   21645 | 3747558 | 300 |
| 7_merge_anchors.others   |    1503 |  335039 | 245 |
| glue_anchors             |   22517 | 3746705 | 279 |
| fill_anchors             |   81095 | 3776774 | 122 |
| spades.contig            |  127702 | 3945352 | 179 |
| spades.scaffold          |  127702 | 3945652 | 176 |
| spades.non-contained     |  127702 | 3918525 |  87 |
| mr_spades.contig         |   92782 | 3948994 | 232 |
| mr_spades.scaffold       |   94852 | 3949540 | 225 |
| mr_spades.non-contained  |   92782 | 3910606 | 127 |
| megahit.contig           |   83091 | 3940868 | 187 |
| megahit.non-contained    |   83091 | 3905347 | 117 |
| mr_megahit.contig        |  100511 | 3960332 | 161 |
| mr_megahit.non-contained |  100511 | 3932475 | 107 |
| platanus.contig          |   47238 | 3985145 | 555 |
| platanus.scaffold        |   55299 | 3932316 | 329 |
| platanus.non-contained   |   55299 | 3875871 | 171 |


# *Mycobacterium abscessus* 6G-0125-R Full

## Mabs_full: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Mabs_full/1_genome
cd ~/data/anchr/Mabs_full/1_genome

cp ~/data/anchr/ref/Mabs/genome.fa .
cp ~/data/anchr/ref/Mabs/paralogs.fa .
cp ~/data/anchr/ref/Mabs/repetitives.fa .

```

## Mabs_full: download

```shell script
cd ~/data/anchr/Mabs_full

mkdir -p ena
cd ena

cat << EOF > source.csv
SRX246890,Mabs_full,MiSeq PE250
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name      | srx       | platform | layout | ilength | srr       | spot    | base  |
|:----------|:----------|:---------|:-------|:--------|:----------|:--------|:------|
| Mabs_full | SRX246890 | ILLUMINA | PAIRED |         | SRR768269 | 4370570 | 2.04G |

* Illumina

```shell script
cd ~/data/anchr/Mabs_full

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR768269_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR768269_2.fastq.gz R2.fq.gz

```


## Mabs_full: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Mabs_full

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 5090491 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    \
    --extend

```

## Mabs_full: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Mabs_full

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 
# rm -fr 4_down_sampling 6_down_sampling

# BASE_NAME=Mabs_full bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
# bkill -J "${BASE_NAME}-*"

# bash 0_master.sh
# bash 0_cleanup.sh

```


# *Rhodobacter sphaeroides* 2.4.1 Full

## Rsph_full: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Rsph_full/1_genome
cd ~/data/anchr/Rsph_full/1_genome

cp ~/data/anchr/ref/Rsph/genome.fa .
cp ~/data/anchr/ref/Rsph/paralogs.fa .
cp ~/data/anchr/ref/Rsph/repetitives.fa .

```

## Rsph_full: download

```shell script
cd ~/data/anchr/Rsph_full

mkdir -p ena
cd ena

cat << EOF > source.csv
SRX160386,Rsph_full,MiSeq PE250
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name      | srx       | platform | layout | ilength | srr       | spot    | base  |
|:----------|:----------|:---------|:-------|:--------|:----------|:--------|:------|
| Rsph_full | SRX160386 | ILLUMINA | PAIRED | 540     | SRR522246 | 8440668 | 3.95G |

* Illumina

```shell script
cd ~/data/anchr/Rsph_full

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR522246_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR522246_2.fastq.gz R2.fq.gz

```


## Rsph_full: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_full

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4602977 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 50 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    \
    --extend

```

## Rsph_full: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_full

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 
# rm -fr 4_down_sampling 6_down_sampling

# BASE_NAME=Rsph_full bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
# bkill -J "${BASE_NAME}-*"

# bash 0_master.sh
# bash 0_cleanup.sh

```


# *Vibrio cholerae* CP1032(5) Full

## Vcho_full: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Vcho_full/1_genome
cd ~/data/anchr/Vcho_full/1_genome

cp ~/data/anchr/ref/Vcho/genome.fa .
cp ~/data/anchr/ref/Vcho/paralogs.fa .
cp ~/data/anchr/ref/Vcho/repetitives.fa .

```

## Vcho_full: download

```shell script
cd ~/data/anchr/Vcho_full

mkdir -p ena
cd ena

cat << EOF > source.csv
SRX247310,Vcho_full,MiSeq PE250
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name      | srx       | platform | layout | ilength | srr       | spot    | base  |
|:----------|:----------|:---------|:-------|:--------|:----------|:--------|:------|
| Vcho_full | SRX247310 | ILLUMINA | PAIRED |         | SRR769320 | 3510275 | 1.64G |

* Illumina

```shell script
cd ~/data/anchr/Vcho_full

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR769320_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR769320_2.fastq.gz R2.fq.gz

```


## Vcho_full: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Vcho_full

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4033464 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

```

## Vcho_full: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Vcho_full

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 
# rm -fr 4_down_sampling 6_down_sampling

# BASE_NAME=Vcho_full bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
# bkill -J "${BASE_NAME}-*"

# bash 0_master.sh
# bash 0_cleanup.sh

```


Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 397.5 |    275 | 1954.3 |                         16.64% |
| R.tadpole.bbtools | 273.9 |    268 |   56.3 |                         85.78% |
| R.genome.picard   | 293.8 |    277 |   47.8 |                             FR |
| R.genome.picard   | 280.5 |    268 |   29.3 |                             RF |
| R.tadpole.picard  | 275.2 |    270 |   46.1 |                             FR |
| R.tadpole.picard  | 268.0 |    267 |   42.2 |                             RF |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 275       | 3752720         | 0.1393       | 46.74   |
| R.31 | 256       | 3717826         | 0.1774       | 47.29   |
| R.41 | 238       | 3723151         | 0.1368       | 48.20   |
| R.51 | 222       | 3821403         | 0.1185       | 48.60   |
| R.61 | 207       | 3888451         | 0.0407       | 48.90   |
| R.71 | 188       | 3822080         | 0.0926       | 49.07   |
| R.81 | 171       | 3819756         | 0.0961       | 48.99   |


Table: statReads

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 2961149 | 4033464 |       2 |
| Paralogs    |    3424 |  119270 |      49 |
| Repetitives |    1070 |  120471 |     244 |
| Illumina.R  |     251 |   1.76G | 7020550 |
| trim.R      |     188 |   1.19G | 6603424 |
| Q20L60      |     188 |   1.17G | 6449509 |
| Q25L60      |     186 |    1.1G | 6192918 |
| Q30L60      |     180 | 968.84M | 5703824 |


Table: statTrimReads

| Name     | N50 |     Sum |       # |
|:---------|----:|--------:|--------:|
| clumpify | 251 |   1.73G | 6883006 |
| trim     | 188 |    1.2G | 6640516 |
| filter   | 188 |   1.19G | 6603424 |
| R1       | 192 | 611.65M | 3301712 |
| R2       | 183 | 576.91M | 3301712 |
| Rs       |   0 |       0 |       0 |


```text
#R.trim
#Matched	5590979	81.22874%
#Name	Reads	ReadsPct
Reverse_adapter	2713329	39.42070%
pcr_dimer	1554664	22.58699%
PCR_Primers	797469	11.58606%
TruSeq_Universal_Adapter	219469	3.18856%
TruSeq_Adapter_Index_1_6	203266	2.95316%
Nextera_LMP_Read2_External_Adapter	83758	1.21688%
```

```text
#R.filter
#Matched	37067	0.55819%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	37067	0.55819%
```

```text
#R.peaks
#k	31
#unique_kmers	34114702
#error_kmers	30343795
#genomic_kmers	3770907
#main_peak	241
#genome_size_in_peaks	3855193
#genome_size	3980124
#haploid_genome_size	3980124
#fold_coverage	241
#haploid_fold_coverage	241
#ploidy	1
#percent_repeat_in_peaks	2.193
#percent_repeat	2.870
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 190 |   1.09G | 6025130 |
| ecco          | 190 |   1.09G | 6025130 |
| ecct          | 190 |   1.08G | 5948138 |
| extended      | 228 |   1.32G | 5948138 |
| merged.raw    | 240 | 711.76M | 2951541 |
| unmerged.raw  | 225 |   9.16M |   45056 |
| unmerged.trim | 225 |   9.16M |   45046 |
| M1            | 239 |  525.6M | 2185202 |
| U1            | 238 |   4.99M |   22523 |
| U2            | 212 |   4.17M |   22523 |
| Us            |   0 |       0 |       0 |
| M.cor         | 239 | 536.95M | 4415450 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 198.8 |    193 |  44.7 |         94.38% |
| M.ihist.merge.txt  | 241.1 |    233 |  52.0 |         99.24% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   | 294.7 |  245.8 |   16.60% | "109" | 4.03M | 4.55M |     1.13 | 0:02'12'' |
| Q20L60.R | 288.9 |  244.3 |   15.43% | "109" | 4.03M |  4.5M |     1.11 | 0:02'06'' |
| Q25L60.R | 272.7 |  237.6 |   12.87% | "109" | 4.03M | 4.37M |     1.08 | 0:02'03'' |
| Q30L60.R | 240.3 |  215.6 |   10.28% | "105" | 4.03M | 4.14M |     1.03 | 0:02'01'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  79.17% |      2641 | 3.19M | 1358 |      1025 | 298.72K | 2820 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:34 |
| Q0L0X40P001   |   40.0 |  80.21% |      2865 | 3.23M | 1313 |      1091 | 297.33K | 2737 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  79.87% |      2776 | 3.21M | 1296 |      1045 | 288.96K | 2726 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:28 |
| Q0L0X80P000   |   80.0 |  60.69% |      1820 | 2.54M | 1433 |      1087 | 273.35K | 2989 |   72.0 |  9.0 |  18.0 | 180.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:28 |
| Q0L0X80P001   |   80.0 |  60.20% |      1842 | 2.55M | 1429 |      1070 | 240.65K | 2963 |   72.0 |  9.0 |  18.0 | 180.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q0L0X80P002   |   80.0 |  58.78% |      1834 | 2.47M | 1397 |      1101 | 256.95K | 2921 |   72.0 |  9.0 |  18.0 | 180.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q20L60X40P000 |   40.0 |  80.71% |      2872 | 3.22M | 1298 |      1033 | 311.45K | 2724 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:29 |
| Q20L60X40P001 |   40.0 |  80.77% |      2778 | 3.23M | 1323 |      1034 | 312.66K | 2816 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:29 |
| Q20L60X40P002 |   40.0 |  80.12% |      2667 | 3.22M | 1358 |      1048 | 308.13K | 2849 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:28 |
| Q20L60X80P000 |   80.0 |  61.93% |      1935 | 2.59M | 1393 |      1118 | 260.46K | 2880 |   73.0 |  9.0 |  18.3 | 182.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q20L60X80P001 |   80.0 |  61.79% |      1869 | 2.57M | 1439 |      1196 | 278.95K | 2994 |   73.0 |  9.0 |  18.3 | 182.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q20L60X80P002 |   80.0 |  62.12% |      1922 |  2.6M | 1423 |      1144 | 257.82K | 2941 |   73.0 |  9.0 |  18.3 | 182.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q25L60X40P000 |   40.0 |  83.07% |      2894 | 3.29M | 1286 |      1025 | 303.97K | 2665 |   38.0 |  5.0 |   9.3 |  96.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:29 |
| Q25L60X40P001 |   40.0 |  83.60% |      3039 | 3.35M | 1266 |      1077 | 289.25K | 2676 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:29 |
| Q25L60X40P002 |   40.0 |  82.77% |      3059 | 3.31M | 1294 |      1092 |  301.4K | 2729 |   37.0 |  5.0 |   9.0 |  94.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:29 |
| Q25L60X80P000 |   80.0 |  68.09% |      2091 | 2.81M | 1436 |      1114 | 272.75K | 2974 |   74.0 |  9.0 |  18.7 | 184.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q25L60X80P001 |   80.0 |  67.60% |      2018 | 2.79M | 1459 |      1140 | 261.52K | 3038 |   74.0 |  9.0 |  18.7 | 184.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:28 |
| Q30L60X40P000 |   40.0 |  93.22% |      5946 | 3.58M |  810 |      1193 | 250.76K | 1675 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| Q30L60X40P001 |   40.0 |  93.08% |      6212 | 3.59M |  813 |      1036 | 252.23K | 1689 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  93.11% |      6008 | 3.64M |  830 |      1076 | 242.49K | 1693 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:29 |
| Q30L60X80P000 |   80.0 |  88.30% |      4115 | 3.56M | 1088 |      1246 | 216.71K | 2187 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:31 |
| Q30L60X80P001 |   80.0 |  88.55% |      3849 | 3.56M | 1127 |      1085 | 202.32K | 2256 |   77.0 |  9.0 |  19.7 | 190.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:32 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.40% |     25584 | 3.16M | 229 |      1095 |  86.96K | 414 |   40.0 |  6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:33 |
| MRX40P001 |   40.0 |  96.69% |     26161 | 3.54M | 250 |      1060 | 117.98K | 488 |   40.0 |  6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:33 |
| MRX40P002 |   40.0 |  96.49% |     27829 | 3.32M | 231 |      1183 |  79.59K | 399 |   40.0 |  6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:32 |
| MRX80P000 |   80.0 |  95.40% |     19221 | 3.55M | 310 |      1067 |  86.97K | 595 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:33 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  94.94% |      9017 | 3.68M |  614 |       981 | 206.72K | 1583 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:38 |
| Q0L0X40P001   |   40.0 |  95.01% |      8612 | 3.63M |  600 |      1023 | 218.79K | 1561 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:39 |
| Q0L0X40P002   |   40.0 |  95.39% |      8845 | 3.65M |  624 |       932 | 184.84K | 1601 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:35 |
| Q0L0X80P000   |   80.0 |  88.79% |      4401 | 3.57M | 1061 |      1137 | 226.51K | 2169 |   77.0 |  9.0 |  19.7 | 190.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:33 |
| Q0L0X80P001   |   80.0 |  88.37% |      4211 | 3.58M | 1073 |      1069 |  213.3K | 2196 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:31 |
| Q0L0X80P002   |   80.0 |  87.70% |      4024 | 3.57M | 1122 |      1060 | 218.86K | 2299 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:31 |
| Q20L60X40P000 |   40.0 |  95.34% |      9007 | 3.59M |  585 |      1037 | 225.03K | 1536 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:34 |
| Q20L60X40P001 |   40.0 |  95.06% |      8221 | 3.68M |  637 |      1030 | 208.79K | 1578 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:34 |
| Q20L60X40P002 |   40.0 |  95.09% |      8891 | 3.63M |  635 |      1001 | 204.34K | 1620 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:35 |
| Q20L60X80P000 |   80.0 |  88.67% |      4373 | 3.58M | 1065 |      1106 | 215.34K | 2160 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:31 |
| Q20L60X80P001 |   80.0 |  88.87% |      4059 | 3.57M | 1082 |      1192 | 246.53K | 2275 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:33 |
| Q20L60X80P002 |   80.0 |  88.28% |      4104 | 3.57M | 1075 |      1102 | 213.34K | 2201 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:33 |
| Q25L60X40P000 |   40.0 |  95.85% |     10047 | 3.66M |  548 |      1059 | 201.63K | 1515 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:36 |
| Q25L60X40P001 |   40.0 |  95.23% |      9079 | 3.64M |  595 |      1033 | 199.32K | 1574 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:35 |
| Q25L60X40P002 |   40.0 |  95.37% |      8815 | 3.69M |  611 |      1034 | 225.07K | 1586 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:35 |
| Q25L60X80P000 |   80.0 |  90.20% |      4961 | 3.61M |  967 |      1161 | 219.43K | 2006 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:31 |
| Q25L60X80P001 |   80.0 |  89.40% |      4465 | 3.59M | 1026 |      1189 | 232.41K | 2116 |   77.0 | 10.0 |  19.0 | 194.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:32 |
| Q30L60X40P000 |   40.0 |  97.00% |     13248 | 3.58M |  450 |       969 | 177.23K | 1331 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:36 |
| Q30L60X40P001 |   40.0 |  97.02% |     12106 | 3.65M |  460 |       959 | 179.07K | 1260 |   40.0 |  6.0 |   9.3 | 104.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:37 |
| Q30L60X40P002 |   40.0 |  96.66% |     13237 | 3.62M |  455 |      1019 | 156.32K | 1202 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:35 |
| Q30L60X80P000 |   80.0 |  94.27% |      7674 | 3.68M |  673 |      1110 |  204.2K | 1402 |   79.0 | 10.0 |  19.7 | 198.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:34 |
| Q30L60X80P001 |   80.0 |  94.09% |      7463 | 3.69M |  692 |      1106 | 187.03K | 1407 |   79.0 | 10.5 |  19.3 | 200.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:33 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.28% |     35874 | 3.47M | 201 |      1045 | 92.61K | 454 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:01:17 |   0:01:15 |
| MRX40P001 |   40.0 |  97.09% |     33413 | 3.58M | 215 |      1006 | 89.56K | 428 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:01:12 |   0:01:14 |
| MRX40P002 |   40.0 |  97.28% |     33370 | 3.61M | 223 |      1007 | 89.51K | 476 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:35 |
| MRX80P000 |   80.0 |  96.36% |     33407 | 3.65M | 215 |      1098 | 67.16K | 385 |   80.0 | 14.0 |  17.3 | 216.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:33 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.82% |     11030 | 3.59M | 504 |      1105 | 203.82K | 1134 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:34 |
| Q0L0X40P001   |   40.0 |  96.14% |     11113 | 3.49M | 484 |      1331 | 244.28K | 1111 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  96.22% |     13291 | 3.66M | 440 |      1016 | 172.96K | 1069 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  94.63% |      7549 | 3.68M | 698 |      1176 | 219.05K | 1694 |   79.0 | 10.0 |  19.7 | 198.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:35 |
| Q0L0X80P001   |   80.0 |  94.34% |      7626 | 3.68M | 724 |      1181 | 181.31K | 1620 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:34 |
| Q0L0X80P002   |   80.0 |  94.37% |      6703 | 3.67M | 764 |      1260 | 232.29K | 1859 |   80.0 | 10.0 |  20.0 | 200.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:35 |
| Q20L60X40P000 |   40.0 |  96.27% |     11105 | 3.42M | 455 |      1092 | 211.39K | 1060 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:34 |
| Q20L60X40P001 |   40.0 |  96.00% |     11901 | 3.72M | 483 |      1006 | 171.77K | 1116 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:34 |
| Q20L60X40P002 |   40.0 |  96.12% |     12671 | 3.68M | 469 |      1025 | 181.35K | 1120 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:34 |
| Q20L60X80P000 |   80.0 |  94.47% |      6913 | 3.67M | 730 |      1166 | 211.19K | 1745 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| Q20L60X80P001 |   80.0 |  94.71% |      7062 | 3.69M | 748 |      1241 |  248.1K | 1804 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:36 |
| Q20L60X80P002 |   80.0 |  94.45% |      7526 | 3.68M | 742 |      1195 | 228.12K | 1718 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:36 |
| Q25L60X40P000 |   40.0 |  96.61% |     11799 | 3.59M | 455 |      1180 | 199.81K | 1042 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:34 |
| Q25L60X40P001 |   40.0 |  96.27% |     12655 | 3.61M | 442 |      1086 | 176.46K | 1071 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:34 |
| Q25L60X40P002 |   40.0 |  96.32% |     12311 | 3.63M | 494 |      1111 | 220.64K | 1134 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:34 |
| Q25L60X80P000 |   80.0 |  95.07% |      7913 | 3.62M | 647 |      1257 | 226.79K | 1611 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:35 |
| Q25L60X80P001 |   80.0 |  95.04% |      7505 | 3.66M | 700 |      1132 | 220.15K | 1704 |   80.0 | 10.0 |  20.0 | 200.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| Q30L60X40P000 |   40.0 |  97.29% |     15599 | 3.43M | 379 |      1045 | 187.84K |  924 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:35 |
| Q30L60X40P001 |   40.0 |  97.35% |     15788 | 3.41M | 367 |      1025 | 192.45K |  877 |   39.0 |  5.0 |   9.7 |  98.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:36 |
| Q30L60X40P002 |   40.0 |  97.23% |     15453 | 3.43M | 383 |      1109 | 167.64K |  888 |   39.0 |  5.5 |   9.3 | 100.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  96.40% |     11012 | 3.64M | 491 |      1116 | 197.59K | 1117 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:35 |
| Q30L60X80P001 |   80.0 |  96.36% |     11296 | 3.69M | 509 |      1216 | 160.97K | 1073 |   80.0 | 11.0 |  19.3 | 204.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:34 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.13% |     37057 | 3.13M | 171 |      1046 | 66.79K | 312 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:36 |
| MRX40P001 |   40.0 |  96.89% |     34195 | 3.24M | 181 |       975 | 69.41K | 325 |   39.0 |  6.5 |   8.7 | 104.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:35 |
| MRX40P002 |   40.0 |  96.81% |     38879 | 3.22M | 188 |      1110 | 74.45K | 315 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:33 |
| MRX80P000 |   80.0 |  96.80% |     39421 | 3.62M | 174 |       979 | 50.44K | 322 |   81.0 | 15.0 |  17.0 | 222.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:35 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  96.52% |     38204 |  3.8M | 222 |      2256 |   1.15M | 545 |  247.0 | 41.0 |  55.0 | 658.0 |   0:00:50 |
| 7_merge_mr_unitigs_bcalm      |  97.13% |     37076 | 3.78M | 229 |      1291 | 157.29K | 127 |  248.0 | 42.0 |  54.7 | 664.0 |   0:01:02 |
| 7_merge_mr_unitigs_superreads |  96.60% |     31083 | 3.75M | 241 |      1708 | 157.43K | 107 |  249.0 | 39.0 |  57.0 | 654.0 |   0:01:03 |
| 7_merge_mr_unitigs_tadpole    |  95.86% |     35753 | 3.61M | 221 |      1192 | 128.59K | 101 |  247.0 | 39.5 |  56.0 | 652.0 |   0:00:56 |
| 7_merge_unitigs_bcalm         |  98.49% |     35647 | 3.79M | 232 |      1823 |  856.8K | 513 |  247.0 | 39.0 |  56.3 | 650.0 |   0:01:30 |
| 7_merge_unitigs_superreads    |  98.21% |     30218 | 3.79M | 256 |      1865 | 902.26K | 503 |  241.0 | 33.0 |  58.3 | 614.0 |   0:01:32 |
| 7_merge_unitigs_tadpole       |  98.40% |     37674 |  3.8M | 222 |      2095 | 816.89K | 438 |  246.0 | 41.0 |  54.7 | 656.0 |   0:01:27 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  97.47% |    102499 | 1.72M |  66 |      1352 | 35.34K | 102 |  247.0 | 61.5 |  41.3 | 740.0 |   0:00:42 |
| 8_mr_spades  |  98.44% |     82666 | 2.02M |  84 |      1443 |    29K | 132 |  134.0 | 28.0 |  26.0 | 380.0 |   0:00:42 |
| 8_megahit    |  97.09% |     47666 | 1.97M | 101 |      1366 | 46.51K | 158 |  247.0 | 43.5 |  53.3 | 668.0 |   0:00:39 |
| 8_mr_megahit |  99.02% |     96551 | 1.09M |  54 |      1362 | 19.46K |  88 |  134.0 | 27.5 |  26.3 | 378.0 |   0:00:35 |
| 8_platanus   |  96.50% |     45123 | 3.54M | 182 |      1003 | 62.59K | 272 |  247.0 | 46.0 |  51.7 | 678.0 |   0:00:44 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 2961149 | 4033464 |    2 |
| Paralogs                 |    3424 |  119270 |   49 |
| Repetitives              |    1070 |  120471 |  244 |
| 7_merge_anchors.anchors  |   38204 | 3801489 |  222 |
| 7_merge_anchors.others   |    2256 | 1154079 |  545 |
| glue_anchors             |   47560 | 3800341 |  195 |
| fill_anchors             |  123134 | 3829423 |   93 |
| spades.contig            |  199415 | 4361080 | 1043 |
| spades.scaffold          |  259449 | 4361326 | 1040 |
| spades.non-contained     |  246373 | 3937748 |   67 |
| mr_spades.contig         |  197574 | 3966769 |  174 |
| mr_spades.scaffold       |  246404 | 3967170 |  168 |
| mr_spades.non-contained  |  197574 | 3932977 |   71 |
| megahit.contig           |  130066 | 4184034 |  662 |
| megahit.non-contained    |  130226 | 3920864 |   91 |
| mr_megahit.contig        |  355251 | 3975128 |  100 |
| mr_megahit.non-contained |  355251 | 3952556 |   48 |
| platanus.contig          |   59932 | 4005160 |  383 |
| platanus.scaffold        |   71841 | 3947598 |  253 |
| platanus.non-contained   |   71841 | 3910474 |  134 |

