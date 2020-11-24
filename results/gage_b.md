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
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 60 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    --redo \
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

Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 446.5 |    278 | 2423.1 |                         14.23% |
| R.tadpole.bbtools | 254.9 |    258 |   50.6 |                         67.36% |
| R.genome.picard   | 295.6 |    279 |   47.2 |                             FR |
| R.genome.picard   | 287.3 |    271 |   33.9 |                             RF |
| R.tadpole.picard  | 264.0 |    264 |   49.2 |                             FR |
| R.tadpole.picard  | 243.6 |    249 |   47.4 |                             RF |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 243       | 5067726         | 0.1772       | 59.12   |
| R.31 | 223       | 5281018         | 0.0000       | 59.24   |
| R.41 | 208       | 5290346         | 0.0000       | 59.68   |
| R.51 | 189       | 5322717         | 0.0000       | 59.96   |
| R.61 | 173       | 5430771         | 0.0000       | 60.23   |
| R.71 | 157       | 5625554         | 0.0000       | 60.41   |
| R.81 | 140       | 5513611         | 0.0000       | 60.42   |


Table: statReads

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 5067172 | 5090491 |       2 |
| Paralogs    |    1693 |   83291 |      53 |
| Repetitives |     192 |   15322 |      82 |
| Illumina.R  |     251 |   2.19G | 8741140 |
| trim.R      |     176 |   1.25G | 7610016 |
| Q20L60      |     177 |   1.21G | 7305465 |
| Q25L60      |     174 |    1.1G | 6854243 |
| Q30L60      |     164 | 896.06M | 6036227 |


Table: statTrimReads

| Name     | N50 |     Sum |       # |
|:---------|----:|--------:|--------:|
| clumpify | 251 |   2.19G | 8732340 |
| highpass | 251 |   2.12G | 8447442 |
| trim     | 176 |   1.25G | 7631352 |
| filter   | 176 |   1.25G | 7610016 |
| R1       | 186 | 667.87M | 3805008 |
| R2       | 166 | 580.17M | 3805008 |
| Rs       |   0 |       0 |       0 |


```text
#R.trim
#Matched	6206957	73.47736%
#Name	Reads	ReadsPct
Reverse_adapter	3222541	38.14813%
pcr_dimer	1717723	20.33424%
TruSeq_Universal_Adapter	512651	6.06871%
PCR_Primers	436127	5.16283%
TruSeq_Adapter_Index_1_6	208483	2.46800%
Nextera_LMP_Read2_External_Adapter	62460	0.73940%
TruSeq_Adapter_Index_11	23556	0.27885%
```

```text
#R.filter
#Matched	21327	0.27947%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	21327	0.27947%
```

```text
#R.peaks
#k	31
#unique_kmers	39624605
#error_kmers	34502033
#genomic_kmers	5122572
#main_peak	190
#genome_size_in_peaks	40125019
#genome_size	42519439
#haploid_genome_size	5314929
#fold_coverage	23
#haploid_fold_coverage	190
#ploidy	8
#het_rate	0.00001
#percent_repeat_in_peaks	0.126
#percent_repeat	87.943
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 177 |   1.22G | 7435420 |
| ecco          | 177 |   1.22G | 7435420 |
| eccc          | 177 |   1.22G | 7435420 |
| ecct          | 176 |   1.19G | 7236082 |
| extended      | 214 |   1.47G | 7236082 |
| merged.raw    | 236 | 835.97M | 3577142 |
| unmerged.raw  | 210 |  15.17M |   81798 |
| unmerged.trim | 210 |  15.17M |   81776 |
| M1            | 236 | 705.18M | 3020325 |
| U1            | 230 |   8.66M |   40888 |
| U2            | 186 |   6.51M |   40888 |
| Us            |   0 |       0 |       0 |
| M.cor         | 235 | 723.36M | 6122426 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 191.2 |    187 |  46.4 |         92.44% |
| M.ihist.merge.txt  | 233.7 |    227 |  51.7 |         98.87% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 245.2 |  197.7 |   19.34% | "45" | 5.09M | 5.57M |     1.09 | 0:02'13'' |
| Q20L60.R | 237.2 |  194.6 |   18.00% | "45" | 5.09M |  5.5M |     1.08 | 0:02'10'' |
| Q25L60.R | 215.5 |  182.3 |   15.41% | "43" | 5.09M | 5.23M |     1.03 | 0:02'01'' |
| Q30L60.R | 176.2 |  154.8 |   12.14% | "39" | 5.09M | 5.18M |     1.02 | 0:01'40'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  75.83% |      2351 | 3.99M | 1872 |      1003 | 278.96K | 4134 |   36.0 |  9.0 |   6.0 | 108.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:30 |
| Q0L0X40P001   |   40.0 |  76.18% |      2287 | 4.02M | 1921 |       946 | 275.76K | 4248 |   36.0 |  9.0 |   6.0 | 108.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:30 |
| Q0L0X40P002   |   40.0 |  76.15% |      2371 | 4.03M | 1885 |       808 | 258.29K | 4157 |   36.0 |  9.0 |   6.0 | 108.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| Q0L0X60P000   |   60.0 |  61.44% |      1833 | 3.23M | 1826 |      1013 | 303.83K | 3981 |   53.0 | 13.0 |   9.0 | 158.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:26 |
| Q0L0X60P001   |   60.0 |  61.76% |      1801 | 3.27M | 1851 |      1010 | 297.06K | 4012 |   53.0 | 14.0 |   8.3 | 162.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:27 |
| Q0L0X60P002   |   60.0 |  60.61% |      1827 |  3.2M | 1814 |      1010 | 295.19K | 3938 |   53.0 | 14.0 |   8.3 | 162.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:27 |
| Q0L0X80P000   |   80.0 |  48.36% |      1581 | 2.54M | 1604 |      1012 | 294.99K | 3512 |   70.0 | 18.0 |  11.3 | 212.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:25 |
| Q0L0X80P001   |   80.0 |  48.13% |      1599 | 2.52M | 1579 |      1016 | 307.73K | 3478 |   70.0 | 18.0 |  11.3 | 212.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:25 |
| Q20L60X40P000 |   40.0 |  78.68% |      2378 | 4.14M | 1914 |       938 | 265.82K | 4252 |   36.0 |  9.0 |   6.0 | 108.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:28 |
| Q20L60X40P001 |   40.0 |  78.45% |      2372 | 4.13M | 1910 |       848 | 280.27K | 4252 |   36.0 |  9.0 |   6.0 | 108.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| Q20L60X40P002 |   40.0 |  77.59% |      2396 | 4.07M | 1886 |      1003 | 296.89K | 4156 |   36.0 |  9.0 |   6.0 | 108.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:28 |
| Q20L60X60P000 |   60.0 |  65.35% |      1887 | 3.46M | 1906 |      1004 | 289.43K | 4129 |   53.0 | 14.0 |   8.3 | 162.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:27 |
| Q20L60X60P001 |   60.0 |  63.69% |      1864 | 3.36M | 1880 |      1012 | 308.17K | 4080 |   53.0 | 14.0 |   8.3 | 162.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:27 |
| Q20L60X60P002 |   60.0 |  64.96% |      1936 | 3.44M | 1886 |      1007 | 282.98K | 4100 |   53.0 | 14.0 |   8.3 | 162.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:27 |
| Q20L60X80P000 |   80.0 |  52.91% |      1649 | 2.79M | 1709 |      1010 | 299.56K | 3713 |   70.0 | 18.0 |  11.3 | 212.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:26 |
| Q20L60X80P001 |   80.0 |  51.84% |      1643 | 2.73M | 1685 |      1016 | 298.24K | 3658 |   70.0 | 18.0 |  11.3 | 212.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:25 |
| Q25L60X40P000 |   40.0 |  95.56% |      5784 | 4.94M | 1136 |        97 | 161.13K | 2773 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:31 |
| Q25L60X40P001 |   40.0 |  95.53% |      5763 | 4.96M | 1142 |        91 | 156.24K | 2761 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:32 |
| Q25L60X40P002 |   40.0 |  95.40% |      5763 | 4.95M | 1141 |        93 | 155.05K | 2765 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:32 |
| Q25L60X60P000 |   60.0 |  92.58% |      4390 | 4.83M | 1380 |        73 | 164.92K | 3131 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:31 |
| Q25L60X60P001 |   60.0 |  92.14% |      4454 | 4.79M | 1406 |       111 | 181.94K | 3205 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:31 |
| Q25L60X60P002 |   60.0 |  92.73% |      4541 | 4.81M | 1360 |       373 | 205.47K | 3126 |   57.0 | 13.0 |  10.3 | 166.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:32 |
| Q25L60X80P000 |   80.0 |  89.15% |      3529 | 4.64M | 1594 |       457 | 224.34K | 3505 |   75.0 | 17.0 |  13.7 | 218.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:33 |
| Q25L60X80P001 |   80.0 |  89.49% |      3636 | 4.66M | 1580 |       461 | 228.65K | 3553 |   75.0 | 17.0 |  13.7 | 218.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:32 |
| Q30L60X40P000 |   40.0 |  97.84% |      8548 | 5.03M |  859 |        71 | 118.79K | 2423 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:37 |
| Q30L60X40P001 |   40.0 |  97.61% |      8182 | 5.03M |  895 |        90 | 128.65K | 2487 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:36 |
| Q30L60X40P002 |   40.0 |  97.86% |      8861 | 5.02M |  846 |        66 | 121.73K | 2518 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:36 |
| Q30L60X60P000 |   60.0 |  96.47% |      6339 | 4.96M | 1093 |       139 | 155.38K | 2746 |   57.0 | 13.0 |  10.3 | 166.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:34 |
| Q30L60X60P001 |   60.0 |  96.13% |      6137 | 4.96M | 1090 |       133 | 158.28K | 2759 |   57.0 | 13.0 |  10.3 | 166.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:35 |
| Q30L60X80P000 |   80.0 |  94.68% |      5060 | 4.89M | 1275 |       230 | 166.46K | 3040 |   76.0 | 17.0 |  14.0 | 220.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:34 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.27% |     23723 | 5.09M | 344 |        65 | 42.48K |  813 |   39.0 |  7.0 |   8.3 | 106.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:33 |
| MRX40P001 |   40.0 |  98.11% |     20152 | 5.09M | 388 |        73 | 48.58K |  905 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:31 |
| MRX40P002 |   40.0 |  98.15% |     23868 | 5.09M | 372 |        72 | 47.15K |  863 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:31 |
| MRX60P000 |   60.0 |  97.18% |     13184 |  5.1M | 584 |        71 | 55.64K | 1284 |   58.0 |  9.0 |  13.3 | 152.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:32 |
| MRX60P001 |   60.0 |  96.84% |     12246 | 5.07M | 625 |        70 | 58.35K | 1368 |   58.0 |  9.0 |  13.3 | 152.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:34 |
| MRX80P000 |   80.0 |  95.50% |      8692 | 5.05M | 849 |        68 | 66.26K | 1809 |   77.0 | 12.0 |  17.7 | 202.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:33 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  96.59% |      8153 | 5.01M |  912 |        59 | 143.05K | 2983 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:37 |
| Q0L0X40P001   |   40.0 |  97.00% |      8547 | 5.03M |  875 |        58 | 145.16K | 3026 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:38 |
| Q0L0X40P002   |   40.0 |  97.16% |      8359 | 5.04M |  855 |        57 | 137.69K | 2940 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:38 |
| Q0L0X60P000   |   60.0 |  92.19% |      4615 | 4.81M | 1344 |        79 | 201.42K | 3336 |   56.0 | 12.0 |  10.7 | 160.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:34 |
| Q0L0X60P001   |   60.0 |  93.32% |      5033 | 4.87M | 1271 |        74 | 186.31K | 3181 |   56.0 | 12.0 |  10.7 | 160.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:34 |
| Q0L0X60P002   |   60.0 |  92.81% |      4846 | 4.84M | 1311 |        73 | 193.24K | 3296 |   56.0 | 12.0 |  10.7 | 160.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  85.57% |      3045 |  4.5M | 1717 |       341 | 246.85K | 3865 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:32 |
| Q0L0X80P001   |   80.0 |  86.08% |      3269 | 4.53M | 1677 |        87 | 222.36K | 3735 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:30 |
| Q20L60X40P000 |   40.0 |  97.28% |      8387 | 5.04M |  852 |        55 | 132.74K | 2953 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:37 |
| Q20L60X40P001 |   40.0 |  96.99% |      8099 | 5.03M |  895 |        64 | 151.74K | 2929 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:38 |
| Q20L60X40P002 |   40.0 |  97.00% |      8489 | 5.03M |  864 |        61 | 142.45K | 2985 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:39 |
| Q20L60X60P000 |   60.0 |  93.68% |      4986 |  4.9M | 1298 |        61 | 165.15K | 3288 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:32 |
| Q20L60X60P001 |   60.0 |  92.77% |      4713 | 4.85M | 1314 |        66 | 170.76K | 3283 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:32 |
| Q20L60X60P002 |   60.0 |  93.19% |      4971 | 4.87M | 1276 |        65 | 167.89K | 3189 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:32 |
| Q20L60X80P000 |   80.0 |  87.02% |      3107 | 4.57M | 1708 |       114 | 237.31K | 3854 |   74.0 | 17.0 |  13.3 | 216.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:31 |
| Q20L60X80P001 |   80.0 |  86.48% |      3344 | 4.54M | 1626 |       160 |  228.2K | 3657 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:32 |
| Q25L60X40P000 |   40.0 |  98.65% |     13941 | 5.07M |  579 |        53 |  94.83K | 2180 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:36 |
| Q25L60X40P001 |   40.0 |  98.78% |     13010 | 5.09M |  589 |        50 |  96.98K | 2321 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:38 |
| Q25L60X40P002 |   40.0 |  98.66% |     12785 | 5.08M |  587 |        54 |  97.63K | 2262 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:38 |
| Q25L60X60P000 |   60.0 |  97.92% |      9439 | 5.02M |  784 |        66 | 136.95K | 2380 |   58.0 | 12.0 |  11.3 | 164.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:38 |
| Q25L60X60P001 |   60.0 |  97.80% |      9018 | 5.06M |  791 |        60 | 119.84K | 2411 |   57.0 | 12.0 |  11.0 | 162.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:36 |
| Q25L60X60P002 |   60.0 |  97.96% |     10044 | 5.04M |  776 |        70 |    134K | 2323 |   58.0 | 12.0 |  11.3 | 164.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:37 |
| Q25L60X80P000 |   80.0 |  96.09% |      6581 | 4.96M | 1027 |        75 | 153.09K | 2601 |   76.0 | 16.0 |  14.7 | 216.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:35 |
| Q25L60X80P001 |   80.0 |  96.52% |      6794 |    5M | 1009 |        61 | 134.96K | 2642 |   77.0 | 16.0 |  15.0 | 218.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:35 |
| Q30L60X40P000 |   40.0 |  98.97% |     14720 | 5.08M |  549 |        49 |  84.84K | 2165 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:37 |
| Q30L60X40P001 |   40.0 |  98.86% |     13034 | 5.08M |  574 |        42 |  79.26K | 2299 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:39 |
| Q30L60X40P002 |   40.0 |  98.91% |     14005 | 5.08M |  553 |        50 |  89.95K | 2180 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:37 |
| Q30L60X60P000 |   60.0 |  98.76% |     11761 | 5.07M |  645 |        56 | 106.34K | 2307 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:40 |
| Q30L60X60P001 |   60.0 |  98.60% |     12503 | 5.07M |  656 |        50 | 101.19K | 2395 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:41 |
| Q30L60X80P000 |   80.0 |  97.89% |      8899 | 5.02M |  821 |        63 | 129.43K | 2531 |   77.0 | 16.0 |  15.0 | 218.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:39 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  99.21% |     89371 |  5.1M | 122 |        75 | 19.77K | 375 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:34 |
| MRX40P001 |   40.0 |  99.18% |     78413 |  5.1M | 130 |       109 | 27.66K | 415 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:35 |
| MRX40P002 |   40.0 |  99.15% |     70142 |  5.1M | 141 |        95 | 25.54K | 433 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:35 |
| MRX60P000 |   60.0 |  99.04% |     58216 |  5.1M | 157 |        92 | 22.01K | 415 |   59.0 |  9.0 |  13.7 | 154.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:37 |
| MRX60P001 |   60.0 |  98.94% |     49777 | 5.11M | 179 |       106 | 24.88K | 472 |   59.0 |  9.0 |  13.7 | 154.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:36 |
| MRX80P000 |   80.0 |  98.74% |     41892 |  5.1M | 212 |        96 | 24.14K | 500 |   78.0 | 11.0 |  18.7 | 200.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:35 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.89% |     12695 | 5.07M |  610 |        63 | 101.11K | 1996 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |
| Q0L0X40P001   |   40.0 |  98.11% |     14771 | 5.07M |  563 |        69 | 104.35K | 1916 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:39 |
| Q0L0X40P002   |   40.0 |  98.00% |     13813 | 5.07M |  560 |        67 |  97.98K | 1846 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:37 |
| Q0L0X60P000   |   60.0 |  96.25% |      7848 | 5.02M |  881 |        92 | 148.92K | 2308 |   57.0 | 12.0 |  11.0 | 162.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:33 |
| Q0L0X60P001   |   60.0 |  96.47% |      8715 | 5.01M |  829 |        71 | 126.16K | 2160 |   57.0 | 12.0 |  11.0 | 162.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:33 |
| Q0L0X60P002   |   60.0 |  96.28% |      7985 | 5.02M |  879 |        65 | 123.72K | 2252 |   57.0 | 12.0 |  11.0 | 162.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  95.24% |      5864 | 4.95M | 1147 |        75 | 190.02K | 3075 |   76.0 | 16.0 |  14.7 | 216.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:37 |
| Q0L0X80P001   |   80.0 |  95.51% |      6017 | 4.99M | 1133 |        60 | 163.69K | 3034 |   76.0 | 16.0 |  14.7 | 216.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:36 |
| Q20L60X40P000 |   40.0 |  98.11% |     13542 | 5.08M |  582 |        66 |  94.38K | 1863 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:36 |
| Q20L60X40P001 |   40.0 |  97.99% |     13346 | 5.06M |  597 |        66 | 100.34K | 1920 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:36 |
| Q20L60X40P002 |   40.0 |  97.97% |     12885 | 5.07M |  592 |        68 | 100.71K | 1931 |   38.0 |  9.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:37 |
| Q20L60X60P000 |   60.0 |  96.72% |      8119 | 5.03M |  887 |        62 | 123.47K | 2321 |   57.0 | 12.0 |  11.0 | 162.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:33 |
| Q20L60X60P001 |   60.0 |  96.07% |      7612 | 4.99M |  906 |        84 | 138.32K | 2297 |   57.0 | 12.0 |  11.0 | 162.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:33 |
| Q20L60X60P002 |   60.0 |  96.45% |      8293 | 5.01M |  863 |        68 | 129.68K | 2289 |   57.0 | 12.0 |  11.0 | 162.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:35 |
| Q20L60X80P000 |   80.0 |  95.51% |      5658 | 4.97M | 1169 |        65 | 175.28K | 3094 |   76.0 | 16.0 |  14.7 | 216.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:36 |
| Q20L60X80P001 |   80.0 |  95.11% |      6005 | 4.95M | 1113 |        67 | 175.79K | 3021 |   76.0 | 16.0 |  14.7 | 216.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:36 |
| Q25L60X40P000 |   40.0 |  99.04% |     18972 | 5.06M |  424 |        61 |  72.96K | 1481 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |
| Q25L60X40P001 |   40.0 |  99.02% |     19955 | 5.09M |  427 |        55 |  68.69K | 1548 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:36 |
| Q25L60X40P002 |   40.0 |  99.04% |     18933 | 5.09M |  417 |        57 |  66.47K | 1531 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |
| Q25L60X60P000 |   60.0 |  98.72% |     14172 | 5.06M |  552 |        88 | 108.93K | 1665 |   58.0 | 12.0 |  11.3 | 164.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:38 |
| Q25L60X60P001 |   60.0 |  98.55% |     12943 | 5.08M |  569 |        72 |  94.52K | 1642 |   58.0 | 12.0 |  11.3 | 164.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:37 |
| Q25L60X60P002 |   60.0 |  98.78% |     15968 | 5.08M |  526 |        85 | 101.67K | 1579 |   58.0 | 12.0 |  11.3 | 164.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:38 |
| Q25L60X80P000 |   80.0 |  98.13% |     11514 | 5.05M |  654 |        77 | 115.73K | 1917 |   77.0 | 16.0 |  15.0 | 218.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:36 |
| Q25L60X80P001 |   80.0 |  98.38% |     12867 | 5.07M |  611 |        58 |  91.82K | 1852 |   77.0 | 16.0 |  15.0 | 218.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:38 |
| Q30L60X40P000 |   40.0 |  99.25% |     21317 | 5.12M |  387 |        51 |  61.51K | 1490 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q30L60X40P001 |   40.0 |  99.22% |     19912 | 5.13M |  421 |        59 |  70.05K | 1579 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:36 |
| Q30L60X40P002 |   40.0 |  99.26% |     20684 | 5.08M |  387 |        49 |  60.42K | 1520 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |
| Q30L60X60P000 |   60.0 |  99.15% |     17428 | 5.08M |  452 |        64 |  82.41K | 1619 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:39 |
| Q30L60X60P001 |   60.0 |  99.12% |     16375 | 5.07M |  480 |        72 |  94.89K | 1747 |   58.0 | 12.0 |  11.3 | 164.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:40 |
| Q30L60X80P000 |   80.0 |  98.96% |     13099 | 5.06M |  559 |        97 | 115.76K | 1887 |   78.0 | 16.0 |  15.3 | 220.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:42 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  99.18% |    102710 |  5.1M | 108 |        94 | 18.82K | 296 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:33 |
| MRX40P001 |   40.0 |  99.19% |     89375 |  5.1M | 112 |       138 | 25.26K | 320 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:33 |
| MRX40P002 |   40.0 |  99.27% |    107057 | 4.95M | 107 |       112 | 21.28K | 333 |   39.0 |  6.0 |   9.0 | 102.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| MRX60P000 |   60.0 |  99.08% |     82757 |  5.1M | 117 |        85 | 14.14K | 288 |   58.0 |  9.0 |  13.3 | 152.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:34 |
| MRX60P001 |   60.0 |  99.18% |     89403 |  5.1M | 112 |        98 | 16.28K | 299 |   59.0 |  9.0 |  13.7 | 154.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:37 |
| MRX80P000 |   80.0 |  99.05% |     72132 |  5.1M | 122 |       152 | 19.16K | 300 |   78.0 | 11.0 |  18.7 | 200.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:36 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  98.75% |     53353 |  5.1M | 189 |      2179 |  739.8K | 368 |  193.0 | 33.0 |  42.3 | 518.0 |   0:00:51 |
| 7_merge_mr_unitigs_bcalm      |  99.14% |     51900 | 4.44M | 167 |       778 |  19.04K |  22 |  192.0 | 33.0 |  42.0 | 516.0 |   0:01:01 |
| 7_merge_mr_unitigs_superreads |  99.12% |     46500 | 5.09M | 195 |       982 |  38.48K |  40 |  192.0 | 33.0 |  42.0 | 516.0 |   0:01:01 |
| 7_merge_mr_unitigs_tadpole    |  99.19% |     47324 | 4.05M | 159 |       738 |  18.77K |  22 |  193.0 | 33.0 |  42.3 | 518.0 |   0:01:05 |
| 7_merge_unitigs_bcalm         |  99.42% |     57042 |  5.1M | 175 |      1458 | 339.88K | 225 |  190.0 | 35.0 |  40.0 | 520.0 |   0:01:24 |
| 7_merge_unitigs_superreads    |  99.30% |     43244 | 5.13M | 216 |      1474 | 433.88K | 294 |  189.0 | 35.0 |  39.7 | 518.0 |   0:01:21 |
| 7_merge_unitigs_tadpole       |  99.41% |     59925 |  5.1M | 171 |      1988 |  372.7K | 204 |  190.0 | 35.0 |  40.0 | 520.0 |   0:01:20 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  84.42% |      3488 | 4.45M | 1555 |       437 | 227.04K | 3134 |  185.0 | 36.0 |  37.7 | 514.0 |   0:00:37 |
| 8_mr_spades  |  99.20% |    121402 | 4.45M |   74 |       150 |  12.62K |  125 |  139.0 | 18.0 |  34.3 | 350.0 |   0:00:37 |
| 8_megahit    |  98.24% |     35436 | 1.96M |   94 |       165 |   12.1K |  133 |  193.0 | 33.0 |  42.3 | 518.0 |   0:00:39 |
| 8_mr_megahit |  99.38% |    128025 | 4.11M |   63 |       150 |   9.34K |  107 |  139.0 | 18.0 |  34.3 | 350.0 |   0:00:36 |
| 8_platanus   |  98.21% |     31927 | 4.68M |  237 |        74 |  28.91K |  417 |  193.0 | 33.0 |  42.3 | 518.0 |   0:00:40 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 5067172 | 5090491 |    2 |
| Paralogs                 |    1693 |   83291 |   53 |
| Repetitives              |     192 |   15322 |   82 |
| 7_merge_anchors.anchors  |   53353 | 5098515 |  189 |
| 7_merge_anchors.others   |    2179 |  739802 |  368 |
| glue_anchors             |   54193 | 5098343 |  182 |
| fill_anchors             |  125196 | 5186424 |   71 |
| spades.contig            |    2938 | 5836762 | 6077 |
| spades.scaffold          |    2942 | 5837300 | 6071 |
| spades.non-contained     |    3684 | 4676177 | 1583 |
| mr_spades.contig         |  166845 | 5137989 |   77 |
| mr_spades.scaffold       |  220390 | 5138416 |   70 |
| mr_spades.non-contained  |  166845 | 5126753 |   53 |
| megahit.contig           |  149440 | 5139631 |  120 |
| megahit.non-contained    |  149440 | 5119717 |   69 |
| mr_megahit.contig        |  215357 | 5146799 |   71 |
| mr_megahit.non-contained |  215357 | 5136352 |   48 |
| platanus.contig          |   20365 | 5199284 |  627 |
| platanus.scaffold        |   46922 | 5133169 |  223 |
| platanus.non-contained   |   46922 | 5116849 |  183 |


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
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 60 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 250 \
    --uscale 2 \
    --lscale 3 \
    --redo \
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

Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 459.0 |    422 | 1253.7 |                         34.48% |
| R.tadpole.bbtools | 406.5 |    420 |   63.4 |                         66.23% |
| R.genome.picard   | 413.0 |    422 |   39.3 |                             FR |
| R.tadpole.picard  | 407.6 |    421 |   47.5 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 383       | 3826134         | 0.6238       | 67.02   |
| R.31 | 336       | 3842214         | 0.6620       | 68.04   |
| R.41 | 302       | 3868647         | 0.5713       | 68.18   |
| R.51 | 273       | 3883371         | 0.5103       | 68.25   |
| R.61 | 247       | 3925348         | 0.4011       | 68.27   |
| R.71 | 217       | 3836390         | 0.3371       | 68.29   |
| R.81 | 205       | 4277776         | 0.2241       | 68.29   |


Table: statReads

| Name        |     N50 |     Sum |        # |
|:------------|--------:|--------:|---------:|
| Genome      | 3188524 | 4602977 |        7 |
| Paralogs    |    2337 |  146789 |       66 |
| Repetitives |     572 |   57281 |      165 |
| Illumina.R  |     251 |   4.24G | 16881336 |
| trim.R      |     150 |   1.64G | 11594652 |
| Q20L60      |     150 |    1.6G | 11350230 |
| Q25L60      |     141 |   1.42G | 10770536 |
| Q30L60      |     120 |   1.06G |  9475633 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 251 |    4.2G | 16724610 |
| highpass | 251 |   3.38G | 13485806 |
| trim     | 150 |   1.64G | 11594652 |
| filter   | 150 |   1.64G | 11594652 |
| R1       | 166 | 905.39M |  5797326 |
| R2       | 134 | 731.39M |  5797326 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	605450	4.48954%
#Name	Reads	ReadsPct
Reverse_adapter	350589	2.59969%
pcr_dimer	116327	0.86259%
PCR_Primers	63705	0.47239%
TruSeq_Universal_Adapter	45368	0.33641%
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	25532192
#error_kmers	21184881
#genomic_kmers	4347311
#main_peak	270
#genome_size_in_peaks	4820900
#genome_size	4954771
#haploid_genome_size	4954771
#fold_coverage	270
#haploid_fold_coverage	270
#ploidy	1
#percent_repeat_in_peaks	9.852
#percent_repeat	11.089
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |        # |
|:--------------|----:|--------:|---------:|
| clumped       | 150 |   1.64G | 11594502 |
| ecco          | 150 |   1.64G | 11594502 |
| eccc          | 150 |   1.64G | 11594502 |
| ecct          | 150 |   1.63G | 11528402 |
| extended      | 188 |   2.09G | 11528402 |
| merged.raw    | 460 |   2.04G |  4802348 |
| unmerged.raw  | 163 | 296.36M |  1923706 |
| unmerged.trim | 163 | 296.34M |  1923574 |
| M1            | 460 |   2.02G |  4770772 |
| U1            | 175 | 160.46M |   961787 |
| U2            | 148 | 135.88M |   961787 |
| Us            |   0 |       0 |        0 |
| M.cor         | 456 |   2.32G | 11465118 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 193.1 |    189 |  63.7 |         10.84% |
| M.ihist.merge.txt  | 423.8 |    457 |  84.5 |         83.31% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 355.6 |  321.1 |    9.70% | "41" |  4.6M | 4.83M |     1.05 | 0:03'04'' |
| Q20L60.R | 347.9 |  316.4 |    9.05% | "41" |  4.6M | 4.76M |     1.03 | 0:03'02'' |
| Q25L60.R | 308.2 |  293.6 |    4.73% | "37" |  4.6M | 4.58M |     0.99 | 0:02'45'' |
| Q30L60.R | 231.3 |  226.5 |    2.06% | "31" |  4.6M | 4.55M |     0.99 | 0:02'13'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  91.87% |      6810 | 4.21M |  887 |      1008 | 291.44K | 2816 |   35.0 |  8.0 |   6.3 | 102.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:36 |
| Q0L0X40P001   |   40.0 |  92.14% |      6681 | 4.22M |  899 |      1042 | 293.59K | 2777 |   35.0 |  8.0 |   6.3 | 102.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  91.72% |      6686 |  4.2M |  898 |      1025 | 323.35K | 2839 |   35.0 |  8.0 |   6.3 | 102.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q0L0X60P000   |   60.0 |  86.53% |      4822 | 4.12M | 1131 |       752 | 231.85K | 3115 |   53.0 | 12.0 |   9.7 | 154.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:33 |
| Q0L0X60P001   |   60.0 |  86.56% |      4820 | 4.09M | 1139 |       953 | 276.33K | 3176 |   52.0 | 12.0 |   9.3 | 152.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:32 |
| Q0L0X60P002   |   60.0 |  85.74% |      4697 | 4.08M | 1148 |       843 | 281.25K | 3188 |   52.0 | 11.0 |  10.0 | 148.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:32 |
| Q0L0X80P000   |   80.0 |  81.99% |      3707 | 3.98M | 1331 |       806 | 249.18K | 3419 |   69.0 | 15.0 |  13.0 | 198.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:31 |
| Q0L0X80P001   |   80.0 |  80.65% |      3666 | 3.92M | 1311 |       872 | 264.92K | 3362 |   69.0 | 15.0 |  13.0 | 198.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:33 |
| Q0L0X80P002   |   80.0 |  82.60% |      3742 | 3.95M | 1301 |      1001 | 299.31K | 3486 |   69.0 | 15.0 |  13.0 | 198.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:33 |
| Q20L60X40P000 |   40.0 |  91.93% |      7399 | 4.21M |  858 |      1016 | 285.17K | 2597 |   36.0 |  8.0 |   6.7 | 104.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:34 |
| Q20L60X40P001 |   40.0 |  92.29% |      7267 | 4.18M |  846 |      1064 | 279.54K | 2604 |   36.0 |  8.0 |   6.7 | 104.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:33 |
| Q20L60X40P002 |   40.0 |  92.86% |      7326 |  4.2M |  828 |      1069 | 320.41K | 2569 |   36.0 |  8.0 |   6.7 | 104.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q20L60X60P000 |   60.0 |  87.99% |      5169 | 4.14M | 1058 |       937 | 260.28K | 3003 |   53.0 | 12.0 |   9.7 | 154.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:33 |
| Q20L60X60P001 |   60.0 |  88.59% |      5339 | 4.14M | 1009 |      1004 | 256.35K | 2904 |   53.0 | 12.0 |   9.7 | 154.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:33 |
| Q20L60X60P002 |   60.0 |  88.42% |      5419 | 4.16M | 1029 |       857 | 237.86K | 2955 |   53.0 | 12.0 |   9.7 | 154.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:33 |
| Q20L60X80P000 |   80.0 |  84.23% |      4113 | 4.01M | 1248 |       837 | 289.93K | 3307 |   70.0 | 15.0 |  13.3 | 200.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:33 |
| Q20L60X80P001 |   80.0 |  84.41% |      4213 |    4M | 1219 |       852 |    282K | 3263 |   70.0 | 15.0 |  13.3 | 200.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:34 |
| Q20L60X80P002 |   80.0 |  84.38% |      4043 | 4.01M | 1249 |       831 | 302.94K | 3255 |   70.0 | 15.0 |  13.3 | 200.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:33 |
| Q25L60X40P000 |   40.0 |  97.12% |     10466 | 4.23M |  615 |      1795 | 343.77K | 1723 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  97.24% |     12962 | 4.29M |  554 |      2160 | 283.35K | 1600 |   36.0 |  9.0 |   6.0 | 108.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:30 |
| Q25L60X40P002 |   40.0 |  97.42% |     10686 | 4.22M |  608 |      1647 |  327.9K | 1676 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:31 |
| Q25L60X60P000 |   60.0 |  96.82% |     10283 | 4.23M |  646 |      1374 | 358.64K | 1867 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:34 |
| Q25L60X60P001 |   60.0 |  96.50% |     10675 | 4.21M |  608 |      1376 | 328.63K | 1733 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:34 |
| Q25L60X60P002 |   60.0 |  96.98% |     10610 | 4.25M |  606 |      1462 | 323.31K | 1789 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:35 |
| Q25L60X80P000 |   80.0 |  96.27% |     10449 | 4.25M |  642 |      1382 | 303.14K | 1832 |   73.0 | 16.0 |  13.7 | 210.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:36 |
| Q25L60X80P001 |   80.0 |  96.35% |     10331 | 4.21M |  613 |      1388 | 320.13K | 1844 |   73.0 | 16.0 |  13.7 | 210.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:38 |
| Q25L60X80P002 |   80.0 |  95.94% |     10901 | 4.19M |  615 |      1179 | 270.92K | 1846 |   73.0 | 16.0 |  13.7 | 210.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:36 |
| Q30L60X40P000 |   40.0 |  97.99% |     10104 | 4.18M |  652 |      2956 | 373.93K | 1738 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| Q30L60X40P001 |   40.0 |  97.88% |      9821 |  4.2M |  648 |      2884 | 343.18K | 1764 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| Q30L60X40P002 |   40.0 |  97.97% |      9869 | 4.18M |  645 |      2568 | 331.51K | 1722 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:31 |
| Q30L60X60P000 |   60.0 |  98.11% |     10076 | 4.15M |  632 |      2943 | 387.91K | 1609 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:33 |
| Q30L60X60P001 |   60.0 |  98.03% |     10083 | 4.17M |  632 |      2453 | 355.73K | 1625 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:33 |
| Q30L60X60P002 |   60.0 |  98.07% |     10597 | 4.15M |  599 |      3213 | 379.73K | 1590 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  98.23% |     11330 | 4.15M |  581 |      3253 | 378.45K | 1519 |   75.0 | 18.0 |  13.0 | 222.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:36 |
| Q30L60X80P001 |   80.0 |  98.11% |     11179 | 4.16M |  574 |      2895 | 357.64K | 1514 |   74.0 | 18.0 |  12.7 | 220.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:35 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.82% |     26985 | 3.89M | 250 |      7028 | 238.48K | 500 |   36.0 |  7.0 |   7.3 | 100.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:31 |
| MRX40P001 |   40.0 |  97.81% |     29867 | 3.93M | 248 |      5406 | 220.08K | 498 |   36.0 |  7.0 |   7.3 | 100.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:33 |
| MRX40P002 |   40.0 |  97.89% |     31281 | 4.15M | 260 |      2846 | 178.08K | 507 |   36.0 |  7.0 |   7.3 | 100.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:31 |
| MRX60P000 |   60.0 |  97.71% |     30480 | 3.95M | 238 |      6077 | 231.73K | 469 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:31 |
| MRX60P001 |   60.0 |  97.78% |     32673 |  3.8M | 231 |      4263 |  187.4K | 482 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:32 |
| MRX60P002 |   60.0 |  97.92% |     33888 | 4.18M | 239 |      4254 | 181.93K | 497 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:33 |
| MRX80P000 |   80.0 |  97.66% |     32600 |  3.9M | 235 |      6037 | 239.35K | 473 |   73.0 | 13.0 |  15.7 | 198.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:33 |
| MRX80P001 |   80.0 |  97.88% |     30783 |  3.9M | 243 |      4191 | 210.24K | 540 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:33 |
| MRX80P002 |   80.0 |  97.70% |     29605 | 3.92M | 257 |      5399 | 179.05K | 492 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:32 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  96.51% |     10311 | 4.27M | 633 |      1674 | 321.87K | 2290 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:32 |
| Q0L0X40P001   |   40.0 |  96.12% |     10343 | 4.28M | 638 |      1848 | 300.11K | 2285 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:33 |
| Q0L0X40P002   |   40.0 |  96.36% |     10545 | 4.27M | 637 |      1506 | 294.19K | 2410 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:34 |
| Q0L0X60P000   |   60.0 |  94.90% |      8651 | 4.27M | 701 |      1097 | 284.97K | 2637 |   54.0 | 11.0 |  10.7 | 152.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:36 |
| Q0L0X60P001   |   60.0 |  95.02% |      9704 | 4.27M | 690 |      1183 |  310.3K | 2594 |   54.0 | 11.0 |  10.7 | 152.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  94.59% |      8887 | 4.28M | 730 |      1136 |  291.9K | 2755 |   54.0 | 11.0 |  10.7 | 152.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:37 |
| Q0L0X80P000   |   80.0 |  92.61% |      7503 | 4.25M | 818 |      1002 | 255.67K | 2975 |   71.0 | 15.0 |  13.7 | 202.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:39 |
| Q0L0X80P001   |   80.0 |  92.26% |      6989 | 4.24M | 854 |      1012 | 292.47K | 3144 |   71.0 | 15.0 |  13.7 | 202.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:39 |
| Q0L0X80P002   |   80.0 |  92.59% |      7876 | 4.25M | 825 |      1038 | 288.46K | 3046 |   72.0 | 15.0 |  14.0 | 204.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:37 |
| Q20L60X40P000 |   40.0 |  96.00% |      9396 | 4.26M | 679 |      1689 | 270.82K | 2222 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:31 |
| Q20L60X40P001 |   40.0 |  96.49% |      9670 | 4.25M | 665 |      1682 | 284.57K | 2194 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:31 |
| Q20L60X40P002 |   40.0 |  96.53% |     10508 | 4.26M | 630 |      1569 | 297.71K | 2156 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:32 |
| Q20L60X60P000 |   60.0 |  94.87% |      9294 | 4.28M | 687 |      1127 | 281.93K | 2422 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:35 |
| Q20L60X60P001 |   60.0 |  95.25% |      9916 | 4.28M | 652 |      1244 | 279.44K | 2454 |   54.0 | 12.0 |  10.0 | 156.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:35 |
| Q20L60X60P002 |   60.0 |  95.15% |      9693 | 4.29M | 677 |      1188 | 253.81K | 2589 |   54.0 | 12.0 |  10.0 | 156.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:35 |
| Q20L60X80P000 |   80.0 |  92.52% |      8146 | 4.25M | 772 |      1014 | 255.36K | 2720 |   72.0 | 15.0 |  14.0 | 204.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:36 |
| Q20L60X80P001 |   80.0 |  93.53% |      7985 | 4.26M | 774 |      1040 | 259.25K | 2780 |   72.0 | 15.0 |  14.0 | 204.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:37 |
| Q20L60X80P002 |   80.0 |  92.81% |      7968 | 4.25M | 810 |      1007 | 268.84K | 2858 |   72.0 | 15.0 |  14.0 | 204.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  97.20% |      8071 | 4.18M | 723 |      2753 |  355.8K | 2032 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:30 |
| Q25L60X40P001 |   40.0 |  97.34% |      9049 | 4.19M | 701 |      2672 | 334.18K | 2027 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:31 |
| Q25L60X40P002 |   40.0 |  97.14% |      8948 | 4.16M | 692 |      2483 | 336.05K | 1996 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:31 |
| Q25L60X60P000 |   60.0 |  97.63% |     11111 | 4.27M | 607 |      2672 | 340.39K | 1891 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:33 |
| Q25L60X60P001 |   60.0 |  97.22% |     11947 | 4.27M | 578 |      2225 | 318.83K | 1792 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:34 |
| Q25L60X60P002 |   60.0 |  97.62% |     11679 | 4.25M | 571 |      2286 | 316.73K | 1805 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:34 |
| Q25L60X80P000 |   80.0 |  97.34% |     11584 | 4.28M | 594 |      2125 | 322.57K | 1885 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:36 |
| Q25L60X80P001 |   80.0 |  97.42% |     11592 | 4.27M | 579 |      1906 |  338.6K | 1885 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:36 |
| Q25L60X80P002 |   80.0 |  97.26% |     12037 | 4.27M | 573 |      1968 | 307.34K | 1827 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:35 |
| Q30L60X40P000 |   40.0 |  96.28% |      5979 |  4.1M | 944 |      2872 | 356.28K | 2418 |   38.0 | 10.0 |   6.0 | 116.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:29 |
| Q30L60X40P001 |   40.0 |  96.01% |      6292 |  4.1M | 916 |      2713 | 301.72K | 2353 |   38.0 | 10.0 |   6.0 | 116.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  96.29% |      6393 | 4.11M | 916 |      2672 | 346.87K | 2382 |   38.0 | 10.0 |   6.0 | 116.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:29 |
| Q30L60X60P000 |   60.0 |  97.25% |      8070 | 4.16M | 764 |      2882 | 367.84K | 2097 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:32 |
| Q30L60X60P001 |   60.0 |  97.53% |      8457 | 4.18M | 753 |      2672 |  373.4K | 2097 |   57.0 | 14.0 |   9.7 | 170.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:33 |
| Q30L60X60P002 |   60.0 |  97.27% |      8316 | 4.17M | 737 |      3106 | 374.55K | 2052 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:32 |
| Q30L60X80P000 |   80.0 |  97.67% |      9588 | 4.16M | 664 |      2672 | 361.25K | 1936 |   75.0 | 18.0 |  13.0 | 222.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:34 |
| Q30L60X80P001 |   80.0 |  97.86% |      9198 | 4.15M | 669 |      2821 | 370.34K | 1879 |   75.0 | 18.0 |  13.0 | 222.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:34 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.91% |     31592 | 4.19M | 254 |      8856 | 240.75K | 600 |   36.0 |  7.0 |   7.3 | 100.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:30 |
| MRX40P001 |   40.0 |  97.82% |     32800 | 4.34M | 263 |      7359 | 237.43K | 620 |   37.0 |  7.0 |   7.7 | 102.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:30 |
| MRX40P002 |   40.0 |  97.85% |     29115 | 4.29M | 264 |      3200 | 193.59K | 624 |   36.0 |  7.0 |   7.3 | 100.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:31 |
| MRX60P000 |   60.0 |  97.77% |     34809 | 4.27M | 246 |      6077 | 228.04K | 540 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:32 |
| MRX60P001 |   60.0 |  97.82% |     32690 | 4.31M | 253 |      4950 | 195.51K | 566 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:31 |
| MRX60P002 |   60.0 |  97.83% |     35912 | 4.34M | 241 |      5247 | 200.08K | 538 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:33 |
| MRX80P000 |   80.0 |  97.78% |     33415 | 4.27M | 249 |      6059 | 240.09K | 516 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:33 |
| MRX80P001 |   80.0 |  97.78% |     32825 | 4.26M | 247 |      4263 | 195.76K | 525 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:33 |
| MRX80P002 |   80.0 |  97.91% |     33899 | 4.14M | 245 |      8066 | 199.44K | 526 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:34 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.82% |     14126 | 4.28M | 486 |      2752 | 282.43K | 1474 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:32 |
| Q0L0X40P001   |   40.0 |  97.97% |     13978 | 4.29M | 477 |      3439 | 362.83K | 1462 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:32 |
| Q0L0X40P002   |   40.0 |  97.90% |     14763 | 4.28M | 478 |      4364 | 352.74K | 1518 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:32 |
| Q0L0X60P000   |   60.0 |  98.07% |     14644 | 4.31M | 486 |      2470 | 345.57K | 1617 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:36 |
| Q0L0X60P001   |   60.0 |  98.03% |     14640 |  4.3M | 472 |      2679 | 322.61K | 1592 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  97.92% |     14230 | 4.26M | 497 |      2600 | 338.14K | 1624 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:35 |
| Q0L0X80P000   |   80.0 |  98.20% |     13002 | 4.27M | 506 |      2785 |  384.1K | 1749 |   74.0 | 15.0 |  14.7 | 208.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:38 |
| Q0L0X80P001   |   80.0 |  98.14% |     13474 | 4.29M | 527 |      2473 | 390.27K | 1935 |   73.0 | 15.0 |  14.3 | 206.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:39 |
| Q0L0X80P002   |   80.0 |  98.00% |     13586 | 4.28M | 518 |      2670 | 381.98K | 1869 |   73.0 | 15.0 |  14.3 | 206.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| Q20L60X40P000 |   40.0 |  97.95% |     12879 | 4.26M | 526 |      3630 | 320.82K | 1549 |   37.0 |  8.0 |   7.0 | 106.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:32 |
| Q20L60X40P001 |   40.0 |  97.84% |     13149 | 4.29M | 507 |      2952 |  296.4K | 1464 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:31 |
| Q20L60X40P002 |   40.0 |  97.93% |     14160 | 4.27M | 499 |      4227 | 350.18K | 1490 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:32 |
| Q20L60X60P000 |   60.0 |  98.11% |     14146 | 4.29M | 489 |      2823 |  347.9K | 1563 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:36 |
| Q20L60X60P001 |   60.0 |  98.02% |     14290 | 4.24M | 489 |      4411 |  403.1K | 1487 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| Q20L60X60P002 |   60.0 |  98.13% |     13565 |  4.3M | 499 |      2813 | 333.85K | 1589 |   55.0 | 12.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:34 |
| Q20L60X80P000 |   80.0 |  98.22% |     13081 | 4.27M | 521 |      2241 | 366.52K | 1726 |   74.0 | 15.0 |  14.7 | 208.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:36 |
| Q20L60X80P001 |   80.0 |  98.17% |     13415 | 4.22M | 505 |      3320 | 362.73K | 1710 |   74.0 | 15.0 |  14.7 | 208.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:37 |
| Q20L60X80P002 |   80.0 |  98.11% |     12760 | 4.26M | 521 |      2170 | 373.12K | 1668 |   74.0 | 15.0 |  14.7 | 208.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  98.12% |     10437 | 4.22M | 614 |      4504 | 374.95K | 1673 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  97.96% |     11274 | 4.23M | 599 |      2803 | 320.53K | 1664 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q25L60X40P002 |   40.0 |  98.05% |     10809 | 4.19M | 604 |      2914 | 348.68K | 1664 |   37.0 |  9.0 |   6.3 | 110.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:31 |
| Q25L60X60P000 |   60.0 |  98.28% |     12808 | 4.25M | 525 |      3394 | 319.64K | 1480 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:32 |
| Q25L60X60P001 |   60.0 |  98.38% |     14069 | 4.25M | 520 |      3014 |  333.4K | 1479 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:32 |
| Q25L60X60P002 |   60.0 |  98.40% |     13404 | 4.25M | 511 |      3106 | 318.51K | 1493 |   56.0 | 13.0 |  10.0 | 164.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:33 |
| Q25L60X80P000 |   80.0 |  98.39% |     13913 | 4.21M | 503 |      2920 | 310.92K | 1413 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:35 |
| Q25L60X80P001 |   80.0 |  98.51% |     13677 | 4.26M | 499 |      2792 | 320.73K | 1461 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:36 |
| Q25L60X80P002 |   80.0 |  98.50% |     13672 | 4.24M | 505 |      3106 | 329.97K | 1458 |   74.0 | 16.0 |  14.0 | 212.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:36 |
| Q30L60X40P000 |   40.0 |  97.34% |      7704 | 4.14M | 777 |      2881 | 367.82K | 2096 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:31 |
| Q30L60X40P001 |   40.0 |  97.34% |      8081 | 4.19M | 770 |      3082 |  341.6K | 2084 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:32 |
| Q30L60X40P002 |   40.0 |  97.33% |      8213 | 4.18M | 755 |      2801 |  340.8K | 2028 |   38.0 | 10.0 |   6.0 | 116.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:31 |
| Q30L60X60P000 |   60.0 |  97.99% |      9019 | 4.17M | 692 |      2827 | 373.92K | 1842 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:33 |
| Q30L60X60P001 |   60.0 |  97.90% |      9262 | 4.18M | 693 |      2818 | 362.52K | 1849 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:33 |
| Q30L60X60P002 |   60.0 |  97.98% |      9494 | 4.16M | 665 |      2816 |  367.9K | 1849 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  98.26% |     10292 | 4.18M | 626 |      2884 | 373.48K | 1731 |   75.0 | 18.0 |  13.0 | 222.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:35 |
| Q30L60X80P001 |   80.0 |  98.16% |      9873 | 4.17M | 616 |      2895 | 367.84K | 1714 |   75.0 | 18.0 |  13.0 | 222.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:34 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.94% |     32756 | 4.13M | 242 |     12569 | 218.26K | 518 |   36.0 |  7.0 |   7.3 | 100.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| MRX40P001 |   40.0 |  97.83% |     32762 | 4.34M | 261 |      6157 | 223.09K | 551 |   36.0 |  7.0 |   7.3 | 100.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:30 |
| MRX40P002 |   40.0 |  97.90% |     30773 | 4.18M | 249 |      4985 | 167.99K | 522 |   37.0 |  7.0 |   7.7 | 102.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:30 |
| MRX60P000 |   60.0 |  97.89% |     37927 | 4.11M | 225 |      8043 | 229.45K | 463 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:31 |
| MRX60P001 |   60.0 |  97.83% |     34812 | 4.12M | 229 |      5382 | 176.99K | 459 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:30 |
| MRX60P002 |   60.0 |  97.88% |     37969 | 4.28M | 228 |      5530 | 163.04K | 479 |   55.0 | 10.0 |  11.7 | 150.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| MRX80P000 |   80.0 |  97.94% |     36092 | 4.08M | 222 |      8076 | 239.63K | 449 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:33 |
| MRX80P001 |   80.0 |  97.91% |     33684 | 4.18M | 232 |      3381 | 194.82K | 455 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:32 |
| MRX80P002 |   80.0 |  98.00% |     34861 | 4.06M | 225 |     11873 | 218.11K | 446 |   73.0 | 12.0 |  16.3 | 194.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:34 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  92.06% |     19593 | 4.32M | 387 |      8700 |   2.65M | 835 |  293.0 | 52.0 |  63.0 | 794.0 |   0:00:58 |
| 7_merge_mr_unitigs_bcalm      |  96.32% |     20578 | 3.74M | 326 |      9644 | 472.24K | 128 |  293.0 | 54.0 |  61.7 | 802.0 |   0:01:29 |
| 7_merge_mr_unitigs_superreads |  96.32% |     20081 | 3.79M | 332 |     12569 | 487.32K | 143 |  294.0 | 54.0 |  62.0 | 804.0 |   0:01:30 |
| 7_merge_mr_unitigs_tadpole    |  96.19% |     19892 | 3.57M | 310 |     22261 | 739.72K | 127 |  290.0 | 56.0 |  59.3 | 804.0 |   0:01:30 |
| 7_merge_unitigs_bcalm         |  96.30% |     21056 | 4.33M | 371 |      3808 |   1.35M | 587 |  288.0 | 57.0 |  58.0 | 804.0 |   0:01:35 |
| 7_merge_unitigs_superreads    |  96.17% |     20403 | 4.33M | 368 |      2267 |   1.47M | 755 |  287.0 | 56.0 |  58.3 | 798.0 |   0:01:36 |
| 7_merge_unitigs_tadpole       |  96.19% |     21077 | 4.33M | 361 |      7834 |   1.38M | 503 |  288.0 | 58.0 |  57.3 | 808.0 |   0:01:38 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |     Sum |   # | N50Others |     Sum |   # | median |  MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|--------:|----:|----------:|--------:|----:|-------:|-----:|------:|-------:|----------:|
| 8_spades     |  98.83% |     18466 | 326.14K |  29 |     19602 |  89.33K |  48 |  298.0 | 55.0 |  62.7 |  816.0 |   0:00:44 |
| 8_mr_spades  |  99.23% |     58662 |   1.05M |  44 |      8183 |  79.89K |  55 |  470.0 | 68.0 | 111.3 | 1212.0 |   0:00:46 |
| 8_megahit    |  98.10% |     13708 |   1.02M | 120 |     12571 | 184.61K | 184 |  297.0 | 55.0 |  62.3 |  814.0 |   0:00:42 |
| 8_mr_megahit |  99.41% |     30453 |   1.17M |  66 |     19264 | 141.16K |  91 |  469.0 | 68.0 | 111.0 | 1210.0 |   0:00:44 |
| 8_platanus   |  97.70% |     20415 |   1.61M | 135 |      3619 |  183.1K | 213 |  298.0 | 54.0 |  63.3 |  812.0 |   0:00:41 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 3188524 | 4602977 |    7 |
| Paralogs                 |    2337 |  146789 |   66 |
| Repetitives              |     572 |   57281 |  165 |
| 7_merge_anchors.anchors  |   19593 | 4323423 |  387 |
| 7_merge_anchors.others   |    8700 | 2653527 |  835 |
| glue_anchors             |   21018 | 4319668 |  364 |
| fill_anchors             |  110005 | 4341858 |   90 |
| spades.contig            |  250095 | 4577425 |   86 |
| spades.scaffold          |  333463 | 4577753 |   82 |
| spades.non-contained     |  250095 | 4567805 |   44 |
| mr_spades.contig         |  172630 | 4586098 |   73 |
| mr_spades.scaffold       |  204122 | 4586240 |   71 |
| mr_spades.non-contained  |  172630 | 4576305 |   47 |
| megahit.contig           |  151747 | 4573779 |  166 |
| megahit.non-contained    |  151747 | 4544284 |  110 |
| mr_megahit.contig        |  156942 | 4593401 |   95 |
| mr_megahit.non-contained |  156942 | 4578298 |   62 |
| platanus.contig          |    4196 | 4934729 | 3533 |
| platanus.scaffold        |   73363 | 4795432 | 1916 |
| platanus.non-contained   |   97244 | 4529145 |  138 |

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
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 60 80" \
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
| trim.R      |     188 |   1.18G | 6551416 |
| Q20L60      |     188 |   1.16G | 6412990 |
| Q25L60      |     186 |   1.09G | 6164391 |
| Q30L60      |     180 | 964.91M | 5680864 |


Table: statTrimReads

| Name     | N50 |     Sum |       # |
|:---------|----:|--------:|--------:|
| clumpify | 251 |   1.73G | 6883006 |
| highpass | 251 |    1.7G | 6792780 |
| trim     | 188 |   1.19G | 6588208 |
| filter   | 188 |   1.18G | 6551416 |
| R1       | 192 | 607.54M | 3275708 |
| R2       | 183 | 573.18M | 3275708 |
| Rs       |   0 |       0 |       0 |


```text
#R.trim
#Matched	5546294	81.64984%
#Name	Reads	ReadsPct
Reverse_adapter	2691008	39.61571%
pcr_dimer	1546135	22.76145%
PCR_Primers	793401	11.68006%
TruSeq_Universal_Adapter	214101	3.15189%
TruSeq_Adapter_Index_1_6	201895	2.97220%
Nextera_LMP_Read2_External_Adapter	83194	1.22474%
```

```text
#R.filter
#Matched	36770	0.55812%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	36770	0.55812%
```

```text
#R.peaks
#k	31
#unique_kmers	31019710
#error_kmers	26998993
#genomic_kmers	4020717
#main_peak	243
#genome_size_in_peaks	156208614
#genome_size	159153153
#haploid_genome_size	3978828
#fold_coverage	6
#haploid_fold_coverage	241
#ploidy	40
#het_rate	0.00088
#percent_repeat_in_peaks	52.203
#percent_repeat	97.469
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 190 |   1.09G | 5973700 |
| ecco          | 190 |   1.09G | 5973700 |
| eccc          | 190 |   1.09G | 5973700 |
| ecct          | 190 |   1.08G | 5930182 |
| extended      | 229 |   1.31G | 5930182 |
| merged.raw    | 240 | 710.01M | 2943098 |
| unmerged.raw  | 226 |   9.09M |   43986 |
| unmerged.trim | 226 |   9.09M |   43978 |
| M1            | 240 | 524.68M | 2180741 |
| U1            | 239 |   4.96M |   21989 |
| U2            | 213 |   4.14M |   21989 |
| Us            |   0 |       0 |       0 |
| M.cor         | 239 | 535.95M | 4405460 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 198.7 |    193 |  44.4 |         94.73% |
| M.ihist.merge.txt  | 241.2 |    233 |  52.1 |         99.26% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   | 292.7 |  245.1 |   16.26% | "109" | 4.03M | 4.21M |     1.04 | 0:02'07'' |
| Q20L60.R | 287.4 |  243.8 |   15.18% | "109" | 4.03M | 4.17M |     1.03 | 0:02'04'' |
| Q25L60.R | 271.4 |  236.3 |   12.92% | "109" | 4.03M | 4.02M |     1.00 | 0:01'58'' |
| Q30L60.R | 239.3 |  215.2 |   10.06% | "105" | 4.03M | 3.99M |     0.99 | 0:01'44'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  79.39% |      3003 | 3.26M | 1298 |       585 | 190.42K | 2762 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:27 |
| Q0L0X40P001   |   40.0 |  78.97% |      2822 | 3.24M | 1318 |       802 | 208.33K | 2793 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:27 |
| Q0L0X40P002   |   40.0 |  78.93% |      2838 | 3.22M | 1303 |       895 |  213.2K | 2785 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:24 |
| Q0L0X60P000   |   60.0 |  69.44% |      2255 | 2.88M | 1389 |      1001 | 217.54K | 2947 |   55.0 | 14.0 |   9.0 | 166.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:24 |
| Q0L0X60P001   |   60.0 |  68.32% |      2210 | 2.81M | 1390 |      1008 | 250.72K | 2978 |   55.0 | 14.0 |   9.0 | 166.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:24 |
| Q0L0X60P002   |   60.0 |  71.04% |      2177 |  2.9M | 1438 |      1012 | 274.22K | 3076 |   55.0 | 14.0 |   9.0 | 166.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  60.09% |      1918 | 2.46M | 1333 |      1015 | 285.34K | 2889 |   73.0 | 18.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  60.76% |      1918 | 2.47M | 1337 |      1013 | 288.99K | 2886 |   73.0 | 18.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  59.94% |      1889 | 2.47M | 1362 |      1014 | 281.63K | 2900 |   73.0 | 18.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:24 |
| Q20L60X40P000 |   40.0 |  80.38% |      2876 |  3.3M | 1313 |       566 | 183.93K | 2783 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:25 |
| Q20L60X40P001 |   40.0 |  79.83% |      2894 | 3.26M | 1290 |       814 | 202.29K | 2743 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:25 |
| Q20L60X40P002 |   40.0 |  80.70% |      2960 |  3.3M | 1309 |       541 |  189.9K | 2773 |   37.0 | 10.0 |   5.7 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:25 |
| Q20L60X60P000 |   60.0 |  70.42% |      2258 |  2.9M | 1398 |      1006 | 237.31K | 2982 |   55.0 | 14.0 |   9.0 | 166.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:23 |
| Q20L60X60P001 |   60.0 |  71.98% |      2242 | 2.96M | 1434 |      1002 | 247.97K | 3006 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q20L60X60P002 |   60.0 |  71.27% |      2202 |  2.9M | 1400 |      1004 | 272.11K | 3007 |   56.0 | 14.0 |   9.3 | 168.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q20L60X80P000 |   80.0 |  61.88% |      1974 | 2.53M | 1364 |      1012 | 277.23K | 2927 |   73.0 | 18.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:25 |
| Q20L60X80P001 |   80.0 |  61.99% |      2002 | 2.53M | 1339 |      1020 |  285.9K | 2871 |   73.0 | 18.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:26 |
| Q20L60X80P002 |   80.0 |  62.15% |      1960 | 2.55M | 1374 |      1013 |  275.6K | 2935 |   73.0 | 18.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:25 |
| Q25L60X40P000 |   40.0 |  91.12% |      5420 | 3.63M |  885 |       685 | 156.37K | 1887 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:26 |
| Q25L60X40P001 |   40.0 |  91.20% |      5703 | 3.63M |  840 |       727 | 159.64K | 1796 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  91.28% |      5328 | 3.63M |  891 |       759 | 156.55K | 1884 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:27 |
| Q25L60X60P000 |   60.0 |  87.98% |      4205 | 3.52M | 1060 |       701 | 186.19K | 2209 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q25L60X60P001 |   60.0 |  88.21% |      4225 | 3.53M | 1050 |       828 | 187.32K | 2177 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:29 |
| Q25L60X60P002 |   60.0 |  87.60% |      4215 | 3.51M | 1064 |       930 | 199.77K | 2187 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:33 |
| Q25L60X80P000 |   80.0 |  84.75% |      3456 | 3.42M | 1182 |       764 | 208.01K | 2467 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q25L60X80P001 |   80.0 |  84.78% |      3441 | 3.43M | 1200 |       579 | 196.82K | 2467 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:27 |
| Q30L60X40P000 |   40.0 |  93.31% |      6974 | 3.68M |  740 |       822 | 151.11K | 1616 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  92.84% |      7001 | 3.66M |  755 |       608 | 134.07K | 1646 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| Q30L60X40P002 |   40.0 |  93.06% |      6855 | 3.69M |  749 |       496 |  127.3K | 1630 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| Q30L60X60P000 |   60.0 |  90.97% |      5471 | 3.62M |  899 |       834 | 166.89K | 1883 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q30L60X60P001 |   60.0 |  91.10% |      5355 | 3.64M |  897 |       399 | 142.24K | 1910 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q30L60X60P002 |   60.0 |  90.07% |      5338 |  3.6M |  928 |       740 | 153.17K | 1922 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:28 |
| Q30L60X80P000 |   80.0 |  88.47% |      4450 | 3.53M | 1032 |       836 | 182.45K | 2147 |   78.0 | 17.0 |  14.7 | 224.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:28 |
| Q30L60X80P001 |   80.0 |  87.83% |      4306 | 3.54M | 1070 |       597 | 164.88K | 2203 |   78.0 | 17.0 |  14.7 | 224.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.34% |     27834 | 3.69M | 261 |       968 | 82.48K | 486 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:27 |
| MRX40P001 |   40.0 |  96.08% |     26906 | 3.57M | 256 |       792 | 79.47K | 500 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:27 |
| MRX40P002 |   40.0 |  96.07% |     27017 | 3.75M | 267 |       931 | 90.93K | 526 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:28 |
| MRX60P000 |   60.0 |  95.62% |     21737 | 3.68M | 310 |       803 | 85.09K | 593 |   60.0 | 10.0 |  13.3 | 160.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:28 |
| MRX60P001 |   60.0 |  95.68% |     24426 | 3.68M | 292 |       908 | 94.05K | 576 |   60.0 | 10.0 |  13.3 | 160.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:29 |
| MRX80P000 |   80.0 |  95.15% |     18909 | 3.67M | 339 |       871 | 83.97K | 670 |   80.0 | 13.0 |  18.0 | 212.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:28 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.04% |     10784 | 3.74M |  527 |       367 | 125.48K | 1514 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:30 |
| Q0L0X40P001   |   40.0 |  94.92% |     10825 | 3.75M |  534 |       117 | 101.63K | 1449 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:28 |
| Q0L0X40P002   |   40.0 |  94.82% |      9186 | 3.74M |  572 |       218 | 120.31K | 1567 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:29 |
| Q0L0X60P000   |   60.0 |  92.31% |      6587 | 3.68M |  765 |       790 | 163.75K | 1717 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:28 |
| Q0L0X60P001   |   60.0 |  92.05% |      6550 | 3.69M |  785 |       458 | 146.15K | 1749 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:27 |
| Q0L0X60P002   |   60.0 |  92.09% |      6642 | 3.68M |  777 |       691 | 159.04K | 1709 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:28 |
| Q0L0X80P000   |   80.0 |  87.69% |      4446 | 3.53M | 1022 |       567 | 179.19K | 2156 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:28 |
| Q0L0X80P001   |   80.0 |  87.25% |      4278 | 3.52M | 1029 |       744 | 191.01K | 2174 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:28 |
| Q0L0X80P002   |   80.0 |  87.14% |      4250 | 3.51M | 1059 |       731 | 203.36K | 2239 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:28 |
| Q20L60X40P000 |   40.0 |  94.91% |     10247 | 3.73M |  555 |       363 | 119.68K | 1493 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:28 |
| Q20L60X40P001 |   40.0 |  94.94% |     10804 | 3.75M |  513 |       135 | 106.82K | 1419 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:27 |
| Q20L60X40P002 |   40.0 |  95.28% |     10805 | 3.75M |  534 |       190 |  114.9K | 1444 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:27 |
| Q20L60X60P000 |   60.0 |  91.78% |      6425 | 3.68M |  818 |       473 | 149.93K | 1817 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:28 |
| Q20L60X60P001 |   60.0 |  92.19% |      6765 | 3.69M |  777 |       485 | 145.77K | 1707 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:28 |
| Q20L60X60P002 |   60.0 |  91.64% |      6989 | 3.67M |  783 |       350 |  142.4K | 1764 |   58.0 | 13.0 |  10.7 | 168.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:28 |
| Q20L60X80P000 |   80.0 |  88.33% |      4399 | 3.56M | 1040 |       517 |  175.5K | 2205 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:29 |
| Q20L60X80P001 |   80.0 |  88.15% |      4488 | 3.57M | 1017 |       501 | 174.05K | 2127 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:28 |
| Q20L60X80P002 |   80.0 |  88.15% |      4217 | 3.56M | 1063 |       604 | 187.92K | 2245 |   77.0 | 17.0 |  14.3 | 222.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:27 |
| Q25L60X40P000 |   40.0 |  96.27% |     14121 | 3.77M |  431 |       373 | 100.69K | 1252 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:29 |
| Q25L60X40P001 |   40.0 |  96.52% |     14157 | 3.75M |  425 |       573 | 123.42K | 1229 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:30 |
| Q25L60X40P002 |   40.0 |  96.59% |     15853 | 3.77M |  414 |       185 |  94.13K | 1206 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:30 |
| Q25L60X60P000 |   60.0 |  94.80% |     10005 | 3.75M |  569 |       267 | 106.03K | 1312 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:27 |
| Q25L60X60P001 |   60.0 |  94.82% |     10220 | 3.75M |  549 |       669 | 122.59K | 1288 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:28 |
| Q25L60X60P002 |   60.0 |  94.64% |     10053 | 3.75M |  548 |       411 |  108.6K | 1229 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:27 |
| Q25L60X80P000 |   80.0 |  92.95% |      7336 |  3.7M |  730 |       533 | 138.89K | 1543 |   79.0 | 16.0 |  15.7 | 222.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:30 |
| Q25L60X80P001 |   80.0 |  93.07% |      7150 |  3.7M |  732 |       504 |  138.1K | 1538 |   79.0 | 16.0 |  15.7 | 222.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:30 |
| Q30L60X40P000 |   40.0 |  96.80% |     16174 | 3.73M |  397 |       328 | 114.14K | 1273 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:31 |
| Q30L60X40P001 |   40.0 |  96.83% |     14203 |  3.7M |  411 |       432 | 112.27K | 1204 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:31 |
| Q30L60X40P002 |   40.0 |  96.68% |     15073 | 3.76M |  412 |       552 | 128.01K | 1287 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:30 |
| Q30L60X60P000 |   60.0 |  95.41% |     12158 | 3.73M |  484 |       618 | 110.29K | 1167 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:28 |
| Q30L60X60P001 |   60.0 |  95.15% |     10536 | 3.75M |  522 |       325 | 106.24K | 1212 |   60.0 | 13.0 |  11.3 | 172.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:27 |
| Q30L60X60P002 |   60.0 |  95.47% |     10936 | 3.71M |  516 |       614 | 112.36K | 1229 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:27 |
| Q30L60X80P000 |   80.0 |  93.84% |      8382 | 3.72M |  643 |       594 | 122.99K | 1348 |   79.0 | 16.0 |  15.7 | 222.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:29 |
| Q30L60X80P001 |   80.0 |  93.99% |      8074 | 3.71M |  683 |       458 | 126.45K | 1479 |   79.0 | 16.0 |  15.7 | 222.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:29 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.63% |     36908 | 3.79M | 222 |       843 | 77.71K | 465 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:28 |
| MRX40P001 |   40.0 |  96.64% |     42288 |  3.8M | 214 |       627 | 69.53K | 461 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:28 |
| MRX40P002 |   40.0 |  96.78% |     35644 | 3.79M | 227 |       797 | 72.52K | 463 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:28 |
| MRX60P000 |   60.0 |  96.42% |     35781 | 3.73M | 232 |       793 | 69.78K | 439 |   60.0 | 10.0 |  13.3 | 160.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:30 |
| MRX60P001 |   60.0 |  96.43% |     35592 |  3.6M | 211 |      1007 | 76.35K | 405 |   60.0 | 10.0 |  13.3 | 160.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:29 |
| MRX80P000 |   80.0 |  96.40% |     33044 | 3.47M | 228 |       827 | 78.76K | 453 |   79.0 | 12.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:30 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.84% |     14765 | 3.76M | 425 |       882 |  129.2K | 1077 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:27 |
| Q0L0X40P001   |   40.0 |  95.92% |     14744 | 3.77M | 420 |       456 |  95.43K | 1025 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q0L0X40P002   |   40.0 |  95.96% |     14383 | 3.71M | 418 |       705 | 118.59K | 1074 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:28 |
| Q0L0X60P000   |   60.0 |  94.75% |     10016 | 3.75M | 570 |       790 | 129.09K | 1288 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:27 |
| Q0L0X60P001   |   60.0 |  94.38% |     10256 | 3.73M | 562 |       816 | 132.13K | 1230 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| Q0L0X60P002   |   60.0 |  94.59% |      9861 | 3.74M | 571 |       764 | 136.72K | 1249 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  94.54% |      7995 | 3.69M | 673 |       589 | 175.04K | 1751 |   80.0 | 17.0 |  15.3 | 228.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| Q0L0X80P001   |   80.0 |  94.03% |      7601 | 3.68M | 690 |       579 | 170.42K | 1713 |   80.0 | 17.0 |  15.3 | 228.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:29 |
| Q0L0X80P002   |   80.0 |  94.36% |      7411 |  3.7M | 719 |       594 |    174K | 1775 |   80.0 | 16.0 |  16.0 | 224.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:29 |
| Q20L60X40P000 |   40.0 |  95.70% |     13139 | 3.75M | 451 |       757 | 116.78K | 1106 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:27 |
| Q20L60X40P001 |   40.0 |  95.99% |     15279 | 3.73M | 397 |       574 | 113.38K | 1030 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q20L60X40P002 |   40.0 |  95.95% |     13993 | 3.77M | 422 |       665 | 113.95K | 1031 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:27 |
| Q20L60X60P000 |   60.0 |  94.20% |      9521 | 3.74M | 611 |       615 | 126.99K | 1323 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q20L60X60P001 |   60.0 |  94.74% |      9674 | 3.75M | 583 |       553 | 127.78K | 1283 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:27 |
| Q20L60X60P002 |   60.0 |  94.58% |      9881 | 3.71M | 574 |       545 |  115.4K | 1298 |   59.0 | 13.0 |  11.0 | 170.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:27 |
| Q20L60X80P000 |   80.0 |  94.69% |      7940 |  3.7M | 694 |       655 | 176.56K | 1776 |   80.0 | 17.0 |  15.3 | 228.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q20L60X80P001 |   80.0 |  94.81% |      7651 |  3.7M | 692 |       552 | 169.53K | 1751 |   81.0 | 17.0 |  15.7 | 230.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:31 |
| Q20L60X80P002 |   80.0 |  94.72% |      7477 | 3.71M | 698 |       764 | 182.77K | 1774 |   80.0 | 17.0 |  15.3 | 228.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q25L60X40P000 |   40.0 |  96.88% |     16677 | 3.73M | 367 |       856 | 122.31K |  922 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q25L60X40P001 |   40.0 |  96.87% |     19936 | 3.77M | 343 |       700 | 107.57K |  922 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q25L60X40P002 |   40.0 |  96.78% |     20431 | 3.77M | 349 |       729 |  92.55K |  868 |   39.0 |  9.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:28 |
| Q25L60X60P000 |   60.0 |  96.01% |     13432 | 3.75M | 444 |       781 | 123.07K |  979 |   60.0 | 13.0 |  11.3 | 172.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q25L60X60P001 |   60.0 |  96.04% |     14214 | 3.65M | 432 |       929 | 123.48K |  943 |   60.0 | 13.0 |  11.3 | 172.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:29 |
| Q25L60X60P002 |   60.0 |  95.98% |     12188 | 3.76M | 456 |       724 | 120.05K |  951 |   60.0 | 13.0 |  11.3 | 172.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:28 |
| Q25L60X80P000 |   80.0 |  96.30% |     11523 | 3.74M | 497 |       686 | 156.78K | 1286 |   81.0 | 16.0 |  16.3 | 226.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:30 |
| Q25L60X80P001 |   80.0 |  95.87% |     10403 | 3.72M | 536 |       785 | 147.44K | 1275 |   81.0 | 16.0 |  16.3 | 226.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| Q30L60X40P000 |   40.0 |  97.25% |     19322 | 3.66M | 344 |       614 | 108.45K |  896 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:30 |
| Q30L60X40P001 |   40.0 |  96.96% |     17323 | 3.64M | 347 |       735 | 113.76K |  919 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q30L60X40P002 |   40.0 |  97.06% |     18566 | 3.76M | 357 |       622 | 118.36K |  978 |   40.0 |  9.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q30L60X60P000 |   60.0 |  96.58% |     14560 | 3.69M | 411 |       899 | 125.06K |  912 |   60.0 | 13.0 |  11.3 | 172.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q30L60X60P001 |   60.0 |  96.19% |     13505 | 3.77M | 430 |       468 | 100.19K |  960 |   60.0 | 12.0 |  12.0 | 168.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:29 |
| Q30L60X60P002 |   60.0 |  96.27% |     15117 | 3.76M | 419 |       786 | 112.31K |  913 |   60.0 | 13.0 |  11.3 | 172.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:29 |
| Q30L60X80P000 |   80.0 |  96.38% |     12218 | 3.69M | 463 |       891 | 141.47K | 1098 |   81.0 | 16.0 |  16.3 | 226.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:31 |
| Q30L60X80P001 |   80.0 |  96.26% |     11571 | 3.73M | 497 |       718 | 135.97K | 1199 |   81.0 | 16.0 |  16.3 | 226.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.72% |     38698 | 3.73M | 211 |       902 | 76.59K | 394 |   39.0 |  7.0 |   8.3 | 106.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:29 |
| MRX40P001 |   40.0 |  96.88% |     41002 | 3.51M | 198 |       792 | 69.54K | 379 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:29 |
| MRX40P002 |   40.0 |  96.67% |     37568 | 3.64M | 207 |       924 | 70.63K | 387 |   40.0 |  7.0 |   8.7 | 108.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:27 |
| MRX60P000 |   60.0 |  96.72% |     36912 | 3.54M | 208 |       891 | 76.52K | 395 |   60.0 | 10.0 |  13.3 | 160.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:30 |
| MRX60P001 |   60.0 |  96.54% |     37542 | 3.46M | 191 |      1025 | 75.91K | 364 |   60.0 | 10.0 |  13.3 | 160.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:28 |
| MRX80P000 |   80.0 |  96.73% |     34203 | 3.33M | 201 |      1001 |  79.9K | 398 |   81.0 | 13.0 |  18.3 | 214.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:29 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  96.95% |     36914 | 3.79M | 227 |      2012 |   1.06M | 547 |  246.0 | 40.0 |  55.3 | 652.0 |   0:00:48 |
| 7_merge_mr_unitigs_bcalm      |  97.74% |     33384 | 3.78M | 236 |      1179 | 142.85K | 131 |  249.0 | 41.0 |  55.7 | 662.0 |   0:00:58 |
| 7_merge_mr_unitigs_superreads |  97.59% |     31656 | 3.77M | 236 |      1144 | 172.68K | 155 |  248.0 | 41.0 |  55.3 | 660.0 |   0:00:55 |
| 7_merge_mr_unitigs_tadpole    |  97.89% |     33403 | 3.79M | 236 |      1184 | 149.19K | 133 |  249.0 | 41.0 |  55.7 | 662.0 |   0:00:59 |
| 7_merge_unitigs_bcalm         |  98.59% |     37064 |  3.8M | 216 |      1497 | 609.14K | 418 |  245.0 | 43.0 |  53.0 | 662.0 |   0:01:23 |
| 7_merge_unitigs_superreads    |  98.42% |     35887 | 3.82M | 218 |      1608 | 578.62K | 390 |  241.0 | 44.0 |  51.0 | 658.0 |   0:01:13 |
| 7_merge_unitigs_tadpole       |  98.60% |     39312 | 3.81M | 206 |      2118 | 752.97K | 395 |  242.0 | 44.0 |  51.3 | 660.0 |   0:01:21 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  94.17% |     11125 | 3.65M | 520 |       987 | 120.15K | 949 |  245.0 | 43.0 |  53.0 | 662.0 |   0:00:35 |
| 8_mr_spades  |  98.55% |     67938 | 2.29M | 113 |      1041 |  56.32K | 196 |  134.0 | 20.0 |  31.3 | 348.0 |   0:00:31 |
| 8_megahit    |  97.13% |     42358 | 2.74M | 150 |      1091 |  65.36K | 221 |  247.0 | 42.0 |  54.3 | 662.0 |   0:00:36 |
| 8_mr_megahit |  98.97% |     59521 | 1.99M | 101 |      1084 |  55.29K | 169 |  134.0 | 20.0 |  31.3 | 348.0 |   0:00:30 |
| 8_platanus   |  96.56% |     39227 | 3.65M | 210 |       879 |  76.87K | 358 |  248.0 | 41.0 |  55.3 | 660.0 |   0:00:36 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 2961149 | 4033464 |    2 |
| Paralogs                 |    3424 |  119270 |   49 |
| Repetitives              |    1070 |  120471 |  244 |
| 7_merge_anchors.anchors  |   36914 | 3788140 |  227 |
| 7_merge_anchors.others   |    2012 | 1056551 |  547 |
| glue_anchors             |   44007 | 3787733 |  215 |
| fill_anchors             |  104815 | 3816258 |  107 |
| spades.contig            |   12898 | 4165902 | 2176 |
| spades.scaffold          |   12986 | 4166848 | 2166 |
| spades.non-contained     |   13389 | 3870353 |  453 |
| mr_spades.contig         |  112695 | 3954569 |  204 |
| mr_spades.scaffold       |  124669 | 3955015 |  198 |
| mr_spades.non-contained  |  112695 | 3920939 |  106 |
| megahit.contig           |   87668 | 3941731 |  182 |
| megahit.non-contained    |   87668 | 3903364 |  107 |
| mr_megahit.contig        |  203715 | 3978813 |  176 |
| mr_megahit.non-contained |  203715 | 3940173 |   89 |
| platanus.contig          |   45170 | 3985546 |  406 |
| platanus.scaffold        |   59030 | 3934384 |  284 |
| platanus.non-contained   |   59043 | 3891546 |  165 |

