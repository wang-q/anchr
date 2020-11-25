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
| merged.raw    | 236 | 835.96M | 3577134 |
| unmerged.raw  | 210 |  15.17M |   81814 |
| unmerged.trim | 210 |  15.17M |   81792 |
| M1            | 236 | 705.17M | 3020305 |
| U1            | 230 |   8.66M |   40896 |
| U2            | 186 |   6.51M |   40896 |
| Us            |   0 |       0 |       0 |
| M.cor         | 235 | 723.36M | 6122402 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 191.2 |    187 |  46.4 |         92.44% |
| M.ihist.merge.txt  | 233.7 |    227 |  51.7 |         98.87% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 245.2 |  197.7 |   19.34% | "45" | 5.09M | 5.57M |     1.09 | 0:02'11'' |
| Q20L60.R | 237.2 |  194.6 |   18.00% | "47" | 5.09M |  5.5M |     1.08 | 0:02'09'' |
| Q25L60.R | 215.5 |  182.3 |   15.41% | "45" | 5.09M | 5.23M |     1.03 | 0:01'59'' |
| Q30L60.R | 176.2 |  154.8 |   12.14% | "39" | 5.09M | 5.18M |     1.02 | 0:01'43'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  75.99% |      2365 | 4.08M | 1916 |       644 | 206.84K | 4169 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:35 |
| Q0L0X40P001   |   40.0 |  77.19% |      2407 | 4.14M | 1906 |       363 | 199.14K | 4166 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:36 |
| Q0L0X40P002   |   40.0 |  75.55% |      2316 | 4.05M | 1915 |       958 | 216.27K | 4179 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:29 |
| Q0L0X60P000   |   60.0 |  61.19% |      1821 | 3.37M | 1907 |      1002 | 164.21K | 4017 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:27 |
| Q0L0X60P001   |   60.0 |  60.66% |      1859 | 3.34M | 1868 |      1007 | 165.32K | 3938 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:26 |
| Q0L0X60P002   |   60.0 |  61.19% |      1782 | 3.38M | 1950 |      1006 | 161.38K | 4124 |   53.0 | 13.0 |   5.0 | 184.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q0L0X80P000   |   80.0 |  48.78% |      1574 | 2.71M | 1716 |      1004 | 151.04K | 3596 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:27 |
| Q0L0X80P001   |   80.0 |  48.53% |      1568 |  2.7M | 1730 |      1006 | 144.04K | 3631 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:25 |
| Q20L60X40P000 |   40.0 |  78.15% |      2446 | 4.18M | 1876 |       907 | 214.25K | 4112 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:28 |
| Q20L60X40P001 |   40.0 |  77.56% |      2395 | 4.14M | 1910 |       777 | 213.36K | 4185 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:28 |
| Q20L60X40P002 |   40.0 |  77.94% |      2446 | 4.18M | 1894 |       916 | 209.87K | 4092 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:29 |
| Q20L60X60P000 |   60.0 |  65.03% |      1946 | 3.58M | 1935 |      1002 |  150.8K | 4073 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:27 |
| Q20L60X60P001 |   60.0 |  64.10% |      1889 | 3.53M | 1954 |      1001 | 145.08K | 4111 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:26 |
| Q20L60X60P002 |   60.0 |  63.87% |      1918 | 3.49M | 1915 |      1006 | 166.88K | 4060 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:25 |
| Q20L60X80P000 |   80.0 |  52.23% |      1653 | 2.91M | 1778 |      1005 | 133.47K | 3745 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:26 |
| Q20L60X80P001 |   80.0 |  52.42% |      1700 | 2.92M | 1763 |      1003 | 138.45K | 3693 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  95.64% |      6466 | 5.01M | 1062 |        47 |  92.07K | 2657 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:31 |
| Q25L60X40P001 |   40.0 |  95.38% |      6113 | 4.98M | 1083 |        67 | 108.19K | 2675 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:31 |
| Q25L60X40P002 |   40.0 |  95.30% |      6115 | 4.98M | 1105 |        48 |  96.67K | 2704 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:32 |
| Q25L60X60P000 |   60.0 |  92.71% |      4658 | 4.91M | 1336 |        46 |  93.84K | 3019 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:33 |
| Q25L60X60P001 |   60.0 |  91.95% |      4331 | 4.88M | 1420 |        35 |  87.65K | 3177 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:33 |
| Q25L60X60P002 |   60.0 |  92.63% |      4689 |  4.9M | 1347 |        39 |  92.91K | 3118 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:34 |
| Q25L60X80P000 |   80.0 |  89.30% |      3725 | 4.75M | 1573 |        42 | 105.46K | 3444 |   75.0 | 17.0 |   8.0 | 252.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:32 |
| Q25L60X80P001 |   80.0 |  88.88% |      3692 | 4.73M | 1587 |        44 | 108.29K | 3475 |   75.0 | 17.0 |   8.0 | 252.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:33 |
| Q30L60X40P000 |   40.0 |  97.73% |      9710 | 5.06M |  801 |        46 |  75.52K | 2346 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:35 |
| Q30L60X40P001 |   40.0 |  97.61% |      8570 | 5.06M |  870 |        41 |  74.98K | 2478 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:34 |
| Q30L60X40P002 |   40.0 |  97.82% |      8646 | 5.08M |  831 |        39 |   71.6K | 2408 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:34 |
| Q30L60X60P000 |   60.0 |  96.18% |      6687 | 5.03M | 1046 |        30 |  63.84K | 2676 |   57.0 | 13.0 |   6.0 | 192.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:36 |
| Q30L60X60P001 |   60.0 |  96.38% |      6763 | 5.05M | 1052 |        34 |  68.37K | 2682 |   57.0 | 13.0 |   6.0 | 192.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  94.46% |      5084 | 4.96M | 1277 |        39 |  93.57K | 3041 |   76.0 | 17.0 |   8.3 | 254.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:36 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.07% |     23828 |  5.1M | 358 |        50 | 35.33K |  855 |   39.0 |  7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:34 |
| MRX40P001 |   40.0 |  98.24% |     22087 |  5.1M | 357 |        50 | 36.05K |  877 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:33 |
| MRX40P002 |   40.0 |  98.18% |     22293 |  5.1M | 357 |        54 |  37.5K |  855 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:33 |
| MRX60P000 |   60.0 |  96.86% |     13248 | 5.08M | 624 |        53 | 52.04K | 1390 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:34 |
| MRX60P001 |   60.0 |  97.07% |     12925 | 5.08M | 574 |        53 | 49.08K | 1297 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:34 |
| MRX80P000 |   80.0 |  95.31% |      8846 | 5.06M | 855 |        50 | 54.86K | 1833 |   77.0 | 12.0 |  13.7 | 226.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:35 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  96.82% |      8898 | 5.06M |  849 |        38 |  84.26K | 2837 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:36 |
| Q0L0X40P001   |   40.0 |  97.00% |      9150 | 5.07M |  824 |        36 |   78.5K | 2810 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  96.83% |      8263 | 5.06M |  856 |        35 |  81.04K | 2926 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:36 |
| Q0L0X60P000   |   60.0 |  93.07% |      5050 | 4.92M | 1284 |        43 | 108.52K | 3172 |   56.0 | 12.0 |   6.7 | 184.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:33 |
| Q0L0X60P001   |   60.0 |  92.97% |      4750 | 4.94M | 1322 |        33 |  87.63K | 3286 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  92.62% |      4816 | 4.91M | 1312 |        40 | 103.21K | 3209 |   56.0 | 12.0 |   6.7 | 184.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  85.79% |      3185 | 4.61M | 1697 |        51 | 139.15K | 3750 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:32 |
| Q0L0X80P001   |   80.0 |  85.08% |      3132 | 4.57M | 1710 |        47 | 137.53K | 3798 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:41 |   0:00:34 |
| Q20L60X40P000 |   40.0 |  97.36% |      9010 | 5.08M |  784 |        40 |   81.8K | 2701 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:36 |
| Q20L60X40P001 |   40.0 |  97.09% |      8675 | 5.09M |  852 |        36 |   79.1K | 2845 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:34 |
| Q20L60X40P002 |   40.0 |  97.01% |      9035 | 5.05M |  817 |        41 |  90.07K | 2798 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:36 |
| Q20L60X60P000 |   60.0 |  94.06% |      5310 | 4.99M | 1223 |        30 |  73.94K | 3081 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:36 |
| Q20L60X60P001 |   60.0 |  93.04% |      5073 | 4.91M | 1284 |        37 |  95.89K | 3202 |   56.0 | 12.0 |   6.7 | 184.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:35 |
| Q20L60X60P002 |   60.0 |  93.03% |      4900 | 4.92M | 1284 |        33 |  85.96K | 3192 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:34 |
| Q20L60X80P000 |   80.0 |  87.54% |      3388 |  4.7M | 1639 |        44 | 126.29K | 3623 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:33 |
| Q20L60X80P001 |   80.0 |  86.64% |      3294 | 4.63M | 1664 |        44 | 129.94K | 3724 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:41 |   0:00:32 |
| Q25L60X40P000 |   40.0 |  98.75% |     14881 |  5.1M |  535 |        41 |  67.74K | 2165 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:38 |
| Q25L60X40P001 |   40.0 |  98.58% |     13363 | 5.11M |  564 |        34 |  55.93K | 2173 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:37 |
| Q25L60X40P002 |   40.0 |  98.61% |     13849 | 5.11M |  540 |        33 |  54.03K | 2175 |   38.5 |  8.5 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:38 |
| Q25L60X60P000 |   60.0 |  97.81% |      9923 | 5.07M |  734 |        39 |  73.38K | 2303 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:39 |
| Q25L60X60P001 |   60.0 |  97.58% |      9123 | 5.06M |  765 |        39 |  73.25K | 2323 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:37 |
| Q25L60X60P002 |   60.0 |  97.96% |     10109 | 5.08M |  730 |        36 |  67.65K | 2300 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:38 |
| Q25L60X80P000 |   80.0 |  96.07% |      7051 | 5.02M | 1011 |        41 |  88.32K | 2617 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:36 |
| Q25L60X80P001 |   80.0 |  96.16% |      6932 | 5.02M | 1000 |        36 |  78.48K | 2609 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:37 |
| Q30L60X40P000 |   40.0 |  98.97% |     15154 |  5.1M |  509 |        32 |  52.49K | 2084 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:39 |
| Q30L60X40P001 |   40.0 |  98.86% |     14136 |  5.1M |  551 |        33 |  51.33K | 2132 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:37 |
| Q30L60X40P002 |   40.0 |  99.01% |     15066 |  5.1M |  515 |        32 |  50.54K | 2057 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:37 |
| Q30L60X60P000 |   60.0 |  98.65% |     12884 |  5.1M |  612 |        27 |  49.43K | 2268 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:43 |
| Q30L60X60P001 |   60.0 |  98.72% |     11622 | 5.12M |  639 |        29 |  53.65K | 2256 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:40 |
| Q30L60X80P000 |   80.0 |  97.94% |     10037 | 5.08M |  783 |        30 |  61.71K | 2464 |   77.0 | 16.0 |   9.7 | 250.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:40 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  99.15% |     90017 | 5.11M | 105 |        58 | 16.29K | 353 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:35 |
| MRX40P001 |   40.0 |  99.19% |     92752 | 5.11M | 111 |        51 |    14K | 368 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:36 |
| MRX40P002 |   40.0 |  99.06% |     81614 | 5.11M | 125 |        54 | 14.37K | 363 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:34 |
| MRX60P000 |   60.0 |  98.91% |     62225 | 5.11M | 154 |        56 | 15.75K | 425 |   59.0 |  9.0 |  10.7 | 172.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:35 |
| MRX60P001 |   60.0 |  98.94% |     65031 | 5.11M | 150 |        54 | 14.15K | 382 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:36 |
| MRX80P000 |   80.0 |  98.77% |     44973 | 5.11M | 200 |        67 | 18.48K | 511 |   78.0 | 11.0 |  15.0 | 222.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:35 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.96% |     14805 | 5.11M |  529 |        39 |  53.53K | 1828 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:35 |
| Q0L0X40P001   |   40.0 |  97.95% |     15631 | 5.11M |  504 |        39 |  52.87K | 1772 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  97.90% |     13737 | 5.15M |  551 |        40 |  58.63K | 1878 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:37 |
| Q0L0X60P000   |   60.0 |  96.45% |      9043 | 5.07M |  800 |        39 |  66.24K | 2092 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:35 |
| Q0L0X60P001   |   60.0 |  96.32% |      8352 | 5.07M |  831 |        34 |  62.16K | 2184 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  96.41% |      8767 | 5.05M |  827 |        43 |  76.13K | 2225 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:35 |
| Q0L0X80P000   |   80.0 |  95.38% |      6164 | 5.04M | 1087 |        40 | 101.47K | 2935 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:37 |
| Q0L0X80P001   |   80.0 |  94.79% |      5962 | 4.99M | 1140 |        42 | 109.31K | 3019 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:36 |
| Q20L60X40P000 |   40.0 |  98.14% |     16329 | 5.09M |  479 |        43 |  55.26K | 1727 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:37 |
| Q20L60X40P001 |   40.0 |  98.10% |     15623 |  5.1M |  519 |        39 |  53.85K | 1836 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:38 |
| Q20L60X40P002 |   40.0 |  98.15% |     14870 | 5.13M |  524 |        36 |  51.81K | 1862 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:36 |
| Q20L60X60P000 |   60.0 |  96.71% |      9541 | 5.06M |  760 |        38 |  62.84K | 2058 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:34 |
| Q20L60X60P001 |   60.0 |  96.54% |      8925 | 5.08M |  836 |        36 |  64.58K | 2194 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:35 |
| Q20L60X60P002 |   60.0 |  96.45% |      8652 | 5.07M |  849 |        38 |  66.15K | 2162 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:35 |
| Q20L60X80P000 |   80.0 |  95.58% |      6699 | 5.02M | 1033 |        43 | 104.12K | 2860 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:37 |
| Q20L60X80P001 |   80.0 |  95.52% |      6175 | 5.02M | 1095 |        40 |  105.1K | 3045 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  99.12% |     24495 | 5.11M |  340 |        36 |  37.75K | 1392 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:38 |
| Q25L60X40P001 |   40.0 |  99.06% |     21287 | 5.11M |  377 |        39 |  42.74K | 1477 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q25L60X40P002 |   40.0 |  99.03% |     22227 | 5.11M |  350 |        44 |  47.71K | 1398 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |
| Q25L60X60P000 |   60.0 |  98.64% |     16201 | 5.09M |  470 |        41 |  49.68K | 1488 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:39 |
| Q25L60X60P001 |   60.0 |  98.69% |     15171 | 5.12M |  488 |        41 |  53.56K | 1603 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:40 |
| Q25L60X60P002 |   60.0 |  98.77% |     17502 | 5.12M |  478 |        38 |  48.12K | 1545 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:39 |
| Q25L60X80P000 |   80.0 |  98.27% |     12597 | 5.14M |  607 |        39 |  61.08K | 1852 |   77.0 | 16.0 |   9.7 | 250.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:38 |
| Q25L60X80P001 |   80.0 |  98.31% |     13073 | 5.09M |  595 |        34 |  52.61K | 1855 |   78.0 | 16.0 |  10.0 | 252.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:38 |
| Q30L60X40P000 |   40.0 |  99.33% |     25680 | 5.11M |  329 |        52 |  46.34K | 1366 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:38 |
| Q30L60X40P001 |   40.0 |  99.22% |     22736 | 5.13M |  363 |        34 |  35.81K | 1432 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:39 |
| Q30L60X40P002 |   40.0 |  99.28% |     27358 |  5.1M |  334 |        33 |  34.81K | 1364 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:37 |
| Q30L60X60P000 |   60.0 |  99.15% |     20859 | 5.11M |  408 |        33 |  40.71K | 1571 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:41 |
| Q30L60X60P001 |   60.0 |  99.19% |     22133 | 5.12M |  401 |        31 |  37.67K | 1550 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:41 |
| Q30L60X80P000 |   80.0 |  98.95% |     17397 |  5.1M |  477 |        32 |  45.69K | 1726 |   78.0 | 16.0 |  10.0 | 252.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:41 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  99.14% |    111865 | 5.11M |  90 |        74 | 14.32K | 272 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |
| MRX40P001 |   40.0 |  99.22% |    110473 | 5.11M |  95 |        70 | 13.86K | 292 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:35 |
| MRX40P002 |   40.0 |  99.15% |    110503 | 5.11M |  98 |        62 | 12.88K | 284 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:34 |
| MRX60P000 |   60.0 |  99.18% |    111918 | 5.11M |  90 |        78 |  12.3K | 274 |   59.0 |  9.0 |  10.7 | 172.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:38 |
| MRX60P001 |   60.0 |  99.12% |    107284 | 5.11M | 104 |        69 | 11.79K | 283 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:37 |
| MRX80P000 |   80.0 |  99.15% |    107259 | 5.11M |  96 |        95 | 13.06K | 281 |   78.0 | 11.0 |  15.0 | 222.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:37 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  99.03% |     79689 | 5.11M | 121 |      4773 | 500.69K | 186 |  193.0 | 33.0 |  31.3 | 584.0 |   0:00:58 |
| 7_merge_mr_unitigs_bcalm      |  99.16% |     87968 | 5.11M | 117 |      1135 |  10.01K |  10 |  192.0 | 33.0 |  31.0 | 582.0 |   0:01:04 |
| 7_merge_mr_unitigs_superreads |  99.14% |     79691 | 5.11M | 126 |      1051 |  18.25K |  18 |  192.0 | 33.0 |  31.0 | 582.0 |   0:01:03 |
| 7_merge_mr_unitigs_tadpole    |  99.12% |     87968 | 5.08M | 116 |      1287 |   9.08K |   8 |  192.0 | 33.0 |  31.0 | 582.0 |   0:01:01 |
| 7_merge_unitigs_bcalm         |  99.40% |     87976 | 5.11M | 118 |      3188 | 151.32K |  74 |  190.0 | 35.0 |  28.3 | 590.0 |   0:01:23 |
| 7_merge_unitigs_superreads    |  99.35% |     63056 | 5.16M | 161 |      1624 | 227.52K | 140 |  189.0 | 35.0 |  28.0 | 588.0 |   0:01:20 |
| 7_merge_unitigs_tadpole       |  99.45% |     88081 | 5.11M | 114 |     19265 | 271.08K |  77 |  190.0 | 35.0 |  28.3 | 590.0 |   0:01:27 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  84.42% |      3582 | 4.52M | 1556 |       156 | 153.65K | 3131 |  185.0 | 36.0 |  25.7 | 586.0 |   0:00:38 |
| 8_mr_spades  |  99.20% |    127450 | 5.12M |   71 |        92 |   8.24K |  119 |  139.0 | 18.0 |  28.3 | 386.0 |   0:00:37 |
| 8_megahit    |  98.24% |     73343 | 3.17M |   91 |       116 |   9.28K |  148 |  193.0 | 33.0 |  31.3 | 584.0 |   0:00:40 |
| 8_mr_megahit |  99.38% |    142386 | 4.94M |   61 |       109 |   7.45K |  108 |  139.0 | 18.0 |  28.3 | 386.0 |   0:00:39 |
| 8_platanus   |  98.21% |     43191 |  5.1M |  208 |        55 |  18.47K |  386 |  193.0 | 33.0 |  31.3 | 584.0 |   0:00:40 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 5067172 | 5090491 |    2 |
| Paralogs                 |    1693 |   83291 |   53 |
| Repetitives              |     192 |   15322 |   82 |
| 7_merge_anchors.anchors  |   79689 | 5109916 |  121 |
| 7_merge_anchors.others   |    4773 |  500691 |  186 |
| glue_anchors             |   82785 | 5109688 |  115 |
| fill_anchors             |  144405 | 5123381 |   62 |
| spades.contig            |    2938 | 5836762 | 6077 |
| spades.scaffold          |    2942 | 5837300 | 6071 |
| spades.non-contained     |    3684 | 4676177 | 1583 |
| mr_spades.contig         |  166845 | 5137989 |   77 |
| mr_spades.scaffold       |  220390 | 5138416 |   70 |
| mr_spades.non-contained  |  166845 | 5126753 |   53 |
| megahit.contig           |  149440 | 5139746 |  121 |
| megahit.non-contained    |  149440 | 5119832 |   70 |
| mr_megahit.contig        |  215357 | 5146452 |   70 |
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
| merged.raw    | 460 |   2.04G |  4802359 |
| unmerged.raw  | 163 | 296.36M |  1923684 |
| unmerged.trim | 163 | 296.34M |  1923552 |
| M1            | 460 |   2.02G |  4770790 |
| U1            | 175 | 160.46M |   961776 |
| U2            | 148 | 135.88M |   961776 |
| Us            |   0 |       0 |        0 |
| M.cor         | 456 |   2.32G | 11465132 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 193.1 |    189 |  63.7 |         10.84% |
| M.ihist.merge.txt  | 423.8 |    457 |  84.5 |         83.31% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 355.6 |  321.1 |    9.70% | "39" |  4.6M | 4.83M |     1.05 | 0:03'08'' |
| Q20L60.R | 347.9 |  316.4 |    9.05% | "39" |  4.6M | 4.76M |     1.03 | 0:03'04'' |
| Q25L60.R | 308.2 |  293.6 |    4.73% | "35" |  4.6M | 4.58M |     0.99 | 0:02'51'' |
| Q30L60.R | 231.3 |  226.5 |    2.06% | "31" |  4.6M | 4.55M |     0.99 | 0:02'15'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  91.84% |      7318 |  4.3M |  856 |      1010 | 172.75K | 2758 |   35.0 |  8.0 |   5.0 | 118.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:37 |
| Q0L0X40P001   |   40.0 |  91.68% |      7135 |  4.3M |  856 |      1029 | 181.89K | 2695 |   35.0 |  8.0 |   5.0 | 118.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:36 |
| Q0L0X40P002   |   40.0 |  91.21% |      7326 | 4.29M |  846 |      1081 | 160.94K | 2654 |   35.0 |  8.0 |   5.0 | 118.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q0L0X60P000   |   60.0 |  86.63% |      5279 | 4.21M | 1087 |      1012 |  130.6K | 3013 |   52.0 | 12.0 |   5.3 | 176.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:34 |
| Q0L0X60P001   |   60.0 |  86.27% |      5268 |  4.2M | 1108 |       910 | 118.47K | 3098 |   52.0 | 12.0 |   5.3 | 176.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  86.79% |      5239 | 4.21M | 1092 |      1005 | 127.36K | 3101 |   52.0 | 12.0 |   5.3 | 176.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:33 |
| Q0L0X80P000   |   80.0 |  81.63% |      3880 | 4.06M | 1317 |      1000 | 134.62K | 3385 |   69.0 | 15.0 |   8.0 | 228.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:32 |
| Q0L0X80P001   |   80.0 |  81.50% |      3987 | 4.06M | 1288 |       713 | 120.39K | 3296 |   69.0 | 15.0 |   8.0 | 228.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:31 |
| Q0L0X80P002   |   80.0 |  81.26% |      3901 | 4.04M | 1313 |       817 |  134.7K | 3405 |   69.0 | 15.0 |   8.0 | 228.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:33 |
| Q20L60X40P000 |   40.0 |  92.68% |      8215 | 4.31M |  763 |      1156 | 167.01K | 2456 |   35.0 |  8.0 |   5.0 | 118.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:33 |
| Q20L60X40P001 |   40.0 |  92.85% |      8123 | 4.32M |  783 |      1108 | 183.46K | 2480 |   36.0 |  8.0 |   5.0 | 120.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:34 |
| Q20L60X40P002 |   40.0 |  93.18% |      8398 | 4.31M |  768 |      1184 | 191.47K | 2483 |   36.0 |  8.0 |   5.0 | 120.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:33 |
| Q20L60X60P000 |   60.0 |  88.40% |      5998 | 4.25M |  974 |      1026 |  126.1K | 2785 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:35 |
| Q20L60X60P001 |   60.0 |  88.63% |      6531 | 4.25M |  957 |      1017 | 121.87K | 2762 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:34 |
| Q20L60X60P002 |   60.0 |  87.92% |      5874 | 4.24M |  985 |      1002 | 116.91K | 2796 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:33 |
| Q20L60X80P000 |   80.0 |  84.66% |      4562 | 4.14M | 1171 |       757 | 140.56K | 3078 |   70.0 | 15.0 |   8.3 | 230.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:32 |
| Q20L60X80P001 |   80.0 |  83.72% |      4559 | 4.11M | 1171 |       682 | 136.73K | 3114 |   70.0 | 15.0 |   8.3 | 230.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:34 |
| Q20L60X80P002 |   80.0 |  84.22% |      4293 | 4.13M | 1218 |       750 | 135.76K | 3173 |   70.0 | 15.0 |   8.3 | 230.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:36 |
| Q25L60X40P000 |   40.0 |  97.58% |     15251 | 4.37M |  502 |      1992 | 205.44K | 1521 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  97.48% |     13919 | 4.34M |  515 |      1893 | 237.24K | 1592 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q25L60X40P002 |   40.0 |  97.25% |     14080 | 4.35M |  500 |      1861 | 214.84K | 1576 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:33 |
| Q25L60X60P000 |   60.0 |  96.87% |     16285 | 4.34M |  461 |      1883 | 201.86K | 1587 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:35 |
| Q25L60X60P001 |   60.0 |  96.74% |     14717 | 4.35M |  494 |      1862 | 204.41K | 1650 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:36 |
| Q25L60X60P002 |   60.0 |  96.86% |     14914 | 4.35M |  482 |      1710 | 196.46K | 1622 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:37 |
| Q25L60X80P000 |   80.0 |  96.39% |     14592 | 4.34M |  503 |      1484 | 200.73K | 1726 |   73.0 | 16.0 |   8.3 | 242.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:39 |
| Q25L60X80P001 |   80.0 |  96.43% |     14673 | 4.36M |  499 |      1508 | 176.37K | 1756 |   73.0 | 16.0 |   8.3 | 242.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:40 |
| Q25L60X80P002 |   80.0 |  96.33% |     14334 | 4.35M |  534 |      1464 | 185.79K | 1695 |   73.0 | 16.0 |   8.3 | 242.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:37 |
| Q30L60X40P000 |   40.0 |  97.79% |     11483 | 4.32M |  618 |      2825 | 236.76K | 1689 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:30 |
| Q30L60X40P001 |   40.0 |  97.95% |     11816 | 4.31M |  587 |      4161 | 240.75K | 1635 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:31 |
| Q30L60X40P002 |   40.0 |  97.76% |     11164 |  4.3M |  604 |      2441 |  236.1K | 1660 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| Q30L60X60P000 |   60.0 |  98.06% |     14481 | 4.35M |  482 |      5492 | 259.31K | 1433 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:33 |
| Q30L60X60P001 |   60.0 |  98.05% |     15476 | 4.35M |  480 |      6175 | 217.01K | 1429 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:34 |
| Q30L60X60P002 |   60.0 |  98.10% |     14701 | 4.35M |  496 |      6835 | 268.67K | 1431 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  98.18% |     16315 | 4.36M |  446 |      5541 | 246.22K | 1364 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:36 |
| Q30L60X80P001 |   80.0 |  98.21% |     17321 | 4.36M |  446 |      5668 | 251.35K | 1356 |   74.0 | 18.0 |   6.7 | 256.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:36 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.89% |     46872 | 4.37M | 213 |      5393 | 160.29K | 486 |   37.0 |  7.0 |   5.3 | 116.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:44 |
| MRX40P001 |   40.0 |  97.73% |     44623 | 4.37M | 207 |      4798 | 141.34K | 446 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:42 |
| MRX40P002 |   40.0 |  97.82% |     45204 | 4.38M | 205 |      4406 | 142.93K | 481 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:31 |
| MRX60P000 |   60.0 |  97.79% |     45068 | 4.37M | 221 |      4263 | 146.41K | 483 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:37 |
| MRX60P001 |   60.0 |  97.76% |     36732 | 4.39M | 229 |      3022 | 123.13K | 466 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:38 |
| MRX60P002 |   60.0 |  97.86% |     36744 | 4.32M | 232 |      4263 | 160.34K | 514 |   55.0 |  9.0 |   9.3 | 164.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:35 |
| MRX80P000 |   80.0 |  97.67% |     39144 | 4.27M | 228 |      4263 | 141.06K | 464 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:33 |
| MRX80P001 |   80.0 |  97.76% |     34863 | 4.31M | 240 |      2834 | 153.72K | 527 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:34 |
| MRX80P002 |   80.0 |  97.76% |     36745 | 4.23M | 229 |      4937 | 166.79K | 492 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:32 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  96.09% |     10779 | 4.33M | 629 |      1664 | 211.25K | 2218 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:47 |
| Q0L0X40P001   |   40.0 |  96.13% |     10674 | 4.33M | 623 |      1960 | 216.82K | 2283 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:47 |
| Q0L0X40P002   |   40.0 |  96.21% |     10833 | 4.33M | 617 |      1641 | 217.86K | 2229 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:33 |
| Q0L0X60P000   |   60.0 |  94.75% |     10310 | 4.34M | 625 |      1378 | 192.07K | 2432 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:36 |
| Q0L0X60P001   |   60.0 |  94.55% |     10332 | 4.35M | 664 |      1330 | 189.35K | 2614 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  94.74% |     11047 | 4.35M | 624 |      1387 |  204.4K | 2572 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:36 |
| Q0L0X80P000   |   80.0 |  92.20% |      8110 | 4.34M | 791 |      1015 | 141.97K | 2890 |   71.0 | 15.0 |   8.7 | 232.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:38 |
| Q0L0X80P001   |   80.0 |  92.56% |      8226 | 4.36M | 794 |      1028 | 145.13K | 2910 |   71.0 | 15.0 |   8.7 | 232.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:40 |
| Q0L0X80P002   |   80.0 |  92.72% |      7966 | 4.35M | 814 |      1040 | 159.59K | 2961 |   71.0 | 15.0 |   8.7 | 232.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:40 |
| Q20L60X40P000 |   40.0 |  96.55% |     11300 | 4.32M | 612 |      1721 | 226.86K | 2145 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:32 |
| Q20L60X40P001 |   40.0 |  96.30% |     10687 | 4.32M | 623 |      1776 | 238.75K | 2126 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:31 |
| Q20L60X40P002 |   40.0 |  96.31% |     11028 | 4.31M | 611 |      2396 | 227.91K | 2091 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:32 |
| Q20L60X60P000 |   60.0 |  95.16% |     10466 | 4.35M | 633 |      1432 | 170.77K | 2307 |   54.0 | 12.0 |   6.0 | 180.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:35 |
| Q20L60X60P001 |   60.0 |  95.25% |     11360 | 4.35M | 625 |      1277 | 172.58K | 2335 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:37 |
| Q20L60X60P002 |   60.0 |  95.03% |     10728 | 4.34M | 618 |      1380 | 175.17K | 2293 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:37 |
| Q20L60X80P000 |   80.0 |  93.43% |      8625 | 4.34M | 754 |      1232 | 161.92K | 2703 |   72.0 | 15.0 |   9.0 | 234.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:38 |
| Q20L60X80P001 |   80.0 |  93.09% |      9076 | 4.33M | 739 |      1214 |  177.9K | 2658 |   72.0 | 15.0 |   9.0 | 234.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:37 |
| Q20L60X80P002 |   80.0 |  92.90% |      8141 | 4.33M | 754 |      1163 | 151.99K | 2688 |   72.0 | 15.0 |   9.0 | 234.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  97.29% |      9524 |  4.3M | 703 |      3218 | 231.75K | 1984 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:31 |
| Q25L60X40P001 |   40.0 |  97.01% |      9731 | 4.28M | 673 |      2700 | 263.32K | 1917 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:29 |
| Q25L60X40P002 |   40.0 |  97.07% |      9629 | 4.28M | 685 |      3166 | 242.36K | 1977 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:32 |
| Q25L60X60P000 |   60.0 |  97.50% |     13074 | 4.34M | 542 |      2250 |  210.5K | 1716 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:35 |
| Q25L60X60P001 |   60.0 |  97.50% |     13124 | 4.34M | 556 |      2176 | 200.42K | 1770 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:34 |
| Q25L60X60P002 |   60.0 |  97.36% |     13125 | 4.34M | 524 |      2409 | 227.69K | 1727 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:34 |
| Q25L60X80P000 |   80.0 |  97.48% |     14174 | 4.35M | 506 |      1976 | 226.38K | 1732 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:37 |
| Q25L60X80P001 |   80.0 |  97.42% |     13652 | 4.35M | 506 |      2602 | 266.64K | 1766 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:39 |
| Q25L60X80P002 |   80.0 |  97.40% |     14636 | 4.35M | 506 |      2109 | 254.61K | 1741 |   75.0 | 16.0 |   9.0 | 246.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:38 |
| Q30L60X40P000 |   40.0 |  96.42% |      6474 |  4.2M | 923 |      2825 | 256.95K | 2375 |   38.0 | 10.0 |   5.0 | 136.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:30 |
| Q30L60X40P001 |   40.0 |  96.41% |      6933 |  4.2M | 873 |      2820 | 242.41K | 2328 |   38.0 | 10.0 |   5.0 | 136.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  96.15% |      6324 | 4.19M | 925 |      2650 | 239.81K | 2321 |   38.0 | 10.0 |   5.0 | 136.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:30 |
| Q30L60X60P000 |   60.0 |  97.52% |      9339 | 4.28M | 703 |      4983 | 246.83K | 1970 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:33 |
| Q30L60X60P001 |   60.0 |  97.31% |      8883 | 4.28M | 717 |      3130 |  246.7K | 2035 |   57.0 | 14.0 |   5.0 | 198.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:32 |
| Q30L60X60P002 |   60.0 |  97.54% |      9442 | 4.28M | 699 |      3593 | 257.36K | 1957 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:32 |
| Q30L60X80P000 |   80.0 |  97.66% |     11671 | 4.32M | 601 |      5792 | 251.42K | 1761 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:35 |
| Q30L60X80P001 |   80.0 |  97.62% |     10341 | 4.31M | 624 |      2965 | 247.03K | 1822 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:36 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.73% |     39640 | 4.37M | 225 |      4263 | 144.22K | 547 |   37.0 |  7.0 |   5.3 | 116.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:43 |
| MRX40P001 |   40.0 |  97.95% |     45130 | 4.37M | 214 |      4796 |  146.2K | 590 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:44 |
| MRX40P002 |   40.0 |  97.77% |     42013 | 4.36M | 215 |      4263 | 115.78K | 554 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:30 |
| MRX60P000 |   60.0 |  97.78% |     44335 | 4.34M | 214 |      2964 | 126.39K | 514 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:32 |
| MRX60P001 |   60.0 |  97.82% |     42038 | 4.36M | 212 |      4802 | 132.44K | 527 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:32 |
| MRX60P002 |   60.0 |  97.90% |     42185 | 4.37M | 215 |      2975 | 113.96K | 521 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:34 |
| MRX80P000 |   80.0 |  97.71% |     39639 | 4.36M | 223 |      5247 | 177.02K | 509 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:33 |
| MRX80P001 |   80.0 |  97.88% |     38210 | 4.34M | 222 |      4263 | 167.28K | 512 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:34 |
| MRX80P002 |   80.0 |  97.76% |     42024 | 4.36M | 212 |      5406 | 171.37K | 505 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:32 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.90% |     16206 | 4.36M | 440 |      3827 | 253.01K | 1462 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:38 |
| Q0L0X40P001   |   40.0 |  97.77% |     16643 | 4.35M | 437 |      3423 | 262.25K | 1389 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:38 |
| Q0L0X40P002   |   40.0 |  97.82% |     17170 | 4.35M | 428 |      2740 | 248.66K | 1410 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:32 |
| Q0L0X60P000   |   60.0 |  98.00% |     21036 | 4.37M | 382 |      2965 | 228.53K | 1474 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:37 |
| Q0L0X60P001   |   60.0 |  97.95% |     19166 | 4.37M | 389 |      3164 | 246.44K | 1453 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  97.97% |     19114 | 4.35M | 402 |      3384 | 271.02K | 1540 |   55.0 | 11.0 |   7.3 | 176.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:35 |
| Q0L0X80P000   |   80.0 |  97.91% |     17636 | 4.37M | 432 |      2529 | 246.77K | 1650 |   73.0 | 15.0 |   9.3 | 236.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:37 |
| Q0L0X80P001   |   80.0 |  98.01% |     17767 | 4.37M | 430 |      3298 | 291.26K | 1675 |   73.0 | 15.0 |   9.3 | 236.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:39 |
| Q0L0X80P002   |   80.0 |  97.99% |     16273 | 4.37M | 430 |      3592 | 279.52K | 1694 |   73.0 | 15.0 |   9.3 | 236.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:39 |
| Q20L60X40P000 |   40.0 |  97.86% |     16464 | 4.37M | 461 |      2576 | 243.43K | 1486 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:33 |
| Q20L60X40P001 |   40.0 |  97.90% |     17070 | 4.34M | 447 |      2489 | 246.67K | 1457 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:31 |
| Q20L60X40P002 |   40.0 |  97.97% |     16836 | 4.34M | 453 |      3547 | 217.33K | 1449 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:32 |
| Q20L60X60P000 |   60.0 |  98.10% |     19361 | 4.36M | 392 |      3448 | 255.95K | 1452 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:36 |
| Q20L60X60P001 |   60.0 |  98.15% |     20569 | 4.37M | 403 |      2894 |  231.3K | 1438 |   56.0 | 12.0 |   6.7 | 184.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:35 |
| Q20L60X60P002 |   60.0 |  98.20% |     20575 | 4.38M | 396 |      2670 | 233.89K | 1415 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:36 |
| Q20L60X80P000 |   80.0 |  98.11% |     18501 | 4.36M | 416 |      3606 | 277.16K | 1577 |   74.0 | 15.0 |   9.7 | 238.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:37 |
| Q20L60X80P001 |   80.0 |  98.20% |     19299 | 4.36M | 410 |      3387 | 249.57K | 1564 |   74.0 | 15.0 |   9.7 | 238.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:38 |
| Q20L60X80P002 |   80.0 |  98.13% |     16812 | 4.36M | 432 |      3189 | 259.33K | 1617 |   73.0 | 15.0 |   9.3 | 236.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:38 |
| Q25L60X40P000 |   40.0 |  97.95% |     12112 | 4.33M | 571 |      2905 | 192.17K | 1577 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:30 |
| Q25L60X40P001 |   40.0 |  97.98% |     12392 | 4.32M | 537 |      2815 | 222.36K | 1589 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q25L60X40P002 |   40.0 |  98.01% |     12748 | 4.32M | 552 |      4783 | 247.35K | 1614 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:32 |
| Q25L60X60P000 |   60.0 |  98.38% |     17779 | 4.36M | 433 |      3502 | 207.84K | 1375 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:36 |
| Q25L60X60P001 |   60.0 |  98.29% |     17952 | 4.35M | 441 |      4194 | 213.48K | 1377 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:34 |
| Q25L60X60P002 |   60.0 |  98.32% |     18251 | 4.36M | 428 |      5112 | 216.24K | 1331 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:33 |
| Q25L60X80P000 |   80.0 |  98.47% |     19592 | 4.36M | 396 |      3436 |  205.7K | 1342 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:36 |
| Q25L60X80P001 |   80.0 |  98.49% |     18275 | 4.35M | 408 |      5245 | 265.01K | 1363 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:37 |
| Q25L60X80P002 |   80.0 |  98.47% |     18698 | 4.36M | 416 |      5539 | 250.33K | 1354 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| Q30L60X40P000 |   40.0 |  97.43% |      8118 | 4.27M | 778 |      2824 | 258.82K | 2094 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q30L60X40P001 |   40.0 |  97.41% |      8682 | 4.27M | 731 |      2953 | 243.46K | 1996 |   38.0 | 10.0 |   5.0 | 136.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  97.14% |      8036 | 4.26M | 773 |      2440 | 234.32K | 2074 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:30 |
| Q30L60X60P000 |   60.0 |  97.90% |     11572 | 4.32M | 601 |      4987 | 248.76K | 1706 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| Q30L60X60P001 |   60.0 |  97.93% |     12508 | 4.32M | 586 |      6107 | 225.43K | 1725 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:33 |
| Q30L60X60P002 |   60.0 |  98.01% |     10957 | 4.32M | 602 |      3593 | 246.83K | 1736 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:33 |
| Q30L60X80P000 |   80.0 |  98.17% |     13564 | 4.34M | 528 |      6107 | 232.74K | 1592 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:36 |
| Q30L60X80P001 |   80.0 |  98.20% |     13758 | 4.34M | 519 |      2965 | 237.04K | 1596 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:35 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.82% |     42597 | 4.37M | 216 |      5392 | 157.34K | 498 |   37.0 |  7.0 |   5.3 | 116.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:36 |
| MRX40P001 |   40.0 |  97.91% |     46786 | 4.37M | 199 |      5245 | 142.99K | 504 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:36 |
| MRX40P002 |   40.0 |  97.92% |     45737 | 4.36M | 193 |      5245 | 111.27K | 484 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:32 |
| MRX60P000 |   60.0 |  97.92% |     45165 | 4.37M | 208 |      2670 | 124.24K | 473 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:32 |
| MRX60P001 |   60.0 |  97.93% |     50069 | 4.37M | 188 |      4261 |  111.9K | 448 |   54.0 | 10.0 |   8.0 | 168.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:32 |
| MRX60P002 |   60.0 |  97.99% |     45182 | 4.37M | 206 |      2975 | 117.02K | 477 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:35 |
| MRX80P000 |   80.0 |  97.87% |     45045 | 4.29M | 206 |      8198 | 160.56K | 445 |   74.0 | 12.0 |  12.7 | 220.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:35 |
| MRX80P001 |   80.0 |  97.87% |     42465 | 4.34M | 208 |      5245 | 147.37K | 448 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:35 |
| MRX80P002 |   80.0 |  97.87% |     45205 | 4.35M | 197 |      8395 | 166.57K | 423 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:33 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  92.93% |     31409 | 4.37M | 275 |     11586 |   1.56M | 433 |  296.0 | 53.0 |  45.7 | 910.0 |   0:01:07 |
| 7_merge_mr_unitigs_bcalm      |  96.78% |     33046 | 4.32M | 253 |     14502 | 532.37K | 122 |  290.0 | 54.0 |  42.7 | 904.0 |   0:01:28 |
| 7_merge_mr_unitigs_superreads |  97.13% |     39190 | 4.32M | 233 |     19279 | 592.36K | 114 |  290.0 | 56.0 |  40.7 | 916.0 |   0:01:52 |
| 7_merge_mr_unitigs_tadpole    |  96.77% |     34582 | 4.28M | 238 |     14202 | 419.09K |  89 |  290.0 | 55.0 |  41.7 | 910.0 |   0:01:28 |
| 7_merge_unitigs_bcalm         |  96.90% |     31352 | 4.37M | 278 |      4774 | 795.15K | 275 |  290.0 | 56.0 |  40.7 | 916.0 |   0:01:41 |
| 7_merge_unitigs_superreads    |  96.71% |     33081 | 4.37M | 252 |      2712 | 863.28K | 396 |  287.0 | 56.0 |  39.7 | 910.0 |   0:01:40 |
| 7_merge_unitigs_tadpole       |  96.90% |     32802 | 4.37M | 264 |      5762 | 776.02K | 242 |  289.0 | 56.0 |  40.3 | 914.0 |   0:01:35 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |     Sum |   # | N50Others |     Sum |   # | median |  MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|--------:|----:|----------:|--------:|----:|-------:|-----:|------:|-------:|----------:|
| 8_spades     |  98.83% |     29628 | 530.91K |  38 |      8152 |  80.54K |  60 |  298.0 | 55.0 |  44.3 |  926.0 |   0:00:49 |
| 8_mr_spades  |  99.23% |     71110 |   1.61M |  54 |      9161 |  56.47K |  64 |  470.0 | 68.0 |  88.7 | 1348.0 |   0:00:52 |
| 8_megahit    |  98.14% |     26501 |   1.81M | 132 |      8309 | 121.13K | 214 |  297.0 | 55.0 |  44.0 |  924.0 |   0:00:52 |
| 8_mr_megahit |  99.45% |     45159 |    1.8M |  79 |      8353 |  65.14K | 104 |  469.0 | 68.0 |  88.3 | 1346.0 |   0:01:00 |
| 8_platanus   |  97.70% |     30534 |   2.77M | 161 |      3068 |  168.2K | 254 |  298.0 | 54.0 |  45.3 |  920.0 |   0:00:47 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 3188524 | 4602977 |    7 |
| Paralogs                 |    2337 |  146789 |   66 |
| Repetitives              |     572 |   57281 |  165 |
| 7_merge_anchors.anchors  |   31409 | 4369712 |  275 |
| 7_merge_anchors.others   |   11586 | 1561718 |  433 |
| glue_anchors             |   32654 | 4367264 |  257 |
| fill_anchors             |  131501 | 4431543 |   84 |
| spades.contig            |  250095 | 4577425 |   86 |
| spades.scaffold          |  333463 | 4577753 |   82 |
| spades.non-contained     |  250095 | 4567805 |   44 |
| mr_spades.contig         |  172630 | 4586343 |   75 |
| mr_spades.scaffold       |  204122 | 4586485 |   73 |
| mr_spades.non-contained  |  172630 | 4576422 |   48 |
| megahit.contig           |  151747 | 4572637 |  158 |
| megahit.non-contained    |  151747 | 4546933 |  110 |
| mr_megahit.contig        |  156942 | 4590788 |   88 |
| mr_megahit.non-contained |  156942 | 4579130 |   63 |
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
| merged.raw    | 240 | 710.01M | 2943100 |
| unmerged.raw  | 226 |   9.09M |   43982 |
| unmerged.trim | 226 |   9.09M |   43972 |
| M1            | 240 | 524.67M | 2180698 |
| U1            | 239 |   4.96M |   21986 |
| U2            | 213 |   4.14M |   21986 |
| Us            |   0 |       0 |       0 |
| M.cor         | 239 | 535.94M | 4405368 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 198.7 |    193 |  44.4 |         94.73% |
| M.ihist.merge.txt  | 241.2 |    233 |  52.1 |         99.26% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   | 292.7 |  245.1 |   16.26% | "109" | 4.03M | 4.21M |     1.04 | 0:02'04'' |
| Q20L60.R | 287.4 |  243.8 |   15.18% | "109" | 4.03M | 4.17M |     1.03 | 0:02'02'' |
| Q25L60.R | 271.4 |  236.3 |   12.92% | "107" | 4.03M | 4.02M |     1.00 | 0:01'59'' |
| Q30L60.R | 239.3 |  215.2 |   10.06% | "103" | 4.03M | 3.99M |     0.99 | 0:01'48'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  79.96% |      2911 |  3.3M | 1315 |       986 | 171.68K | 2805 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:31 |
| Q0L0X40P001   |   40.0 |  80.07% |      2954 | 3.31M | 1318 |       785 | 163.76K | 2850 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q0L0X40P002   |   40.0 |  78.82% |      2737 | 3.28M | 1358 |       552 | 161.49K | 2884 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:25 |
| Q0L0X60P000   |   60.0 |  70.05% |      2275 |    3M | 1444 |      1001 |  125.7K | 3015 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:25 |
| Q0L0X60P001   |   60.0 |  69.31% |      2228 | 2.98M | 1457 |       923 | 120.21K | 3049 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:25 |
| Q0L0X60P002   |   60.0 |  69.39% |      2235 | 2.97M | 1450 |      1001 | 121.78K | 3025 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  60.27% |      1950 |  2.6M | 1405 |      1007 | 150.12K | 2947 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  60.00% |      1903 | 2.59M | 1425 |      1003 | 136.23K | 2995 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  59.86% |      1902 |  2.6M | 1421 |      1009 | 144.46K | 2978 |   72.0 | 18.0 |   6.0 | 252.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:24 |
| Q20L60X40P000 |   40.0 |  80.30% |      2987 | 3.32M | 1278 |       829 | 163.64K | 2716 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:24 |
| Q20L60X40P001 |   40.0 |  80.24% |      3001 | 3.33M | 1298 |       665 | 150.11K | 2727 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:25 |
| Q20L60X40P002 |   40.0 |  81.00% |      2903 | 3.34M | 1301 |       955 | 171.83K | 2774 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:25 |
| Q20L60X60P000 |   60.0 |  71.53% |      2356 | 3.06M | 1427 |       529 | 105.29K | 2960 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:31 |
| Q20L60X60P001 |   60.0 |  71.40% |      2338 | 3.06M | 1444 |       672 | 113.97K | 3023 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:25 |
| Q20L60X60P002 |   60.0 |  70.46% |      2321 | 3.01M | 1433 |       587 | 111.68K | 2978 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:25 |
| Q20L60X80P000 |   80.0 |  63.29% |      2012 | 2.73M | 1448 |      1007 | 136.78K | 3018 |   74.0 | 18.0 |   6.7 | 256.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:26 |
| Q20L60X80P001 |   80.0 |  63.44% |      1980 | 2.72M | 1454 |      1008 | 153.48K | 3069 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:25 |
| Q20L60X80P002 |   80.0 |  62.48% |      1952 |  2.7M | 1455 |      1003 | 145.61K | 3041 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:24 |
| Q25L60X40P000 |   40.0 |  91.13% |      6006 | 3.69M |  843 |       104 |  89.99K | 1840 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:27 |
| Q25L60X40P001 |   40.0 |  90.91% |      5788 | 3.67M |  867 |       663 | 108.02K | 1870 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| Q25L60X40P002 |   40.0 |  91.31% |      5850 | 3.68M |  833 |       133 |  91.52K | 1823 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:27 |
| Q25L60X60P000 |   60.0 |  87.70% |      4387 | 3.61M | 1058 |       230 | 103.32K | 2200 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q25L60X60P001 |   60.0 |  87.77% |      4636 | 3.59M | 1027 |       382 |  105.8K | 2152 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q25L60X60P002 |   60.0 |  87.83% |      4596 |  3.6M | 1020 |       575 | 106.01K | 2120 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q25L60X80P000 |   80.0 |  84.67% |      3637 |  3.5M | 1187 |       443 | 113.84K | 2462 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q25L60X80P001 |   80.0 |  84.55% |      3736 | 3.48M | 1157 |       788 | 121.82K | 2406 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q30L60X40P000 |   40.0 |  92.77% |      7601 | 3.73M |  709 |       151 |  77.77K | 1562 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| Q30L60X40P001 |   40.0 |  93.33% |      7633 | 3.74M |  705 |       228 |  76.58K | 1559 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:27 |
| Q30L60X40P002 |   40.0 |  92.92% |      7352 | 3.73M |  726 |       384 |  85.45K | 1614 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| Q30L60X60P000 |   60.0 |  90.40% |      5503 | 3.68M |  876 |        80 |  74.75K | 1828 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:28 |
| Q30L60X60P001 |   60.0 |  90.47% |      5762 | 3.69M |  885 |        72 |  73.21K | 1867 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:28 |
| Q30L60X60P002 |   60.0 |  90.00% |      5653 | 3.66M |  883 |       500 |  89.69K | 1857 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:28 |
| Q30L60X80P000 |   80.0 |  88.44% |      4583 | 3.62M | 1015 |        80 |  82.54K | 2126 |   78.0 | 17.0 |   9.0 | 258.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:29 |
| Q30L60X80P001 |   80.0 |  87.91% |      4497 |  3.6M | 1035 |       361 |  96.06K | 2185 |   78.0 | 17.0 |   9.0 | 258.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:27 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  95.83% |     30198 | 3.77M | 240 |       539 | 48.15K | 502 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRX40P001 |   40.0 |  96.11% |     27775 | 3.76M | 241 |       521 | 51.65K | 512 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRX40P002 |   40.0 |  96.30% |     33376 | 3.75M | 232 |       799 | 57.29K | 475 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:32 |
| MRX60P000 |   60.0 |  95.37% |     25321 | 3.81M | 275 |       514 | 43.67K | 572 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:31 |
| MRX60P001 |   60.0 |  95.85% |     27500 |  3.8M | 272 |       883 | 54.98K | 545 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:28 |
| MRX80P000 |   80.0 |  95.02% |     22127 |  3.8M | 316 |       578 | 51.01K | 656 |   80.0 | 13.0 |  13.7 | 238.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:30 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.06% |     11255 | 3.79M |  501 |        53 | 58.37K | 1466 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:30 |
| Q0L0X40P001   |   40.0 |  95.02% |     11597 | 3.79M |  529 |        73 |  69.4K | 1520 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:30 |
| Q0L0X40P002   |   40.0 |  95.00% |     10683 | 3.79M |  526 |        64 |  64.6K | 1479 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:28 |
| Q0L0X60P000   |   60.0 |  91.98% |      7069 | 3.74M |  758 |        70 | 77.05K | 1740 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:29 |
| Q0L0X60P001   |   60.0 |  91.83% |      6506 | 3.76M |  807 |        57 | 71.63K | 1796 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:28 |
| Q0L0X60P002   |   60.0 |  91.62% |      6359 | 3.73M |  780 |        55 | 67.66K | 1735 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:29 |
| Q0L0X80P000   |   80.0 |  87.77% |      4697 | 3.62M | 1003 |       116 | 98.06K | 2151 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:30 |
| Q0L0X80P001   |   80.0 |  87.41% |      4384 | 3.62M | 1078 |        79 | 96.79K | 2302 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:29 |
| Q0L0X80P002   |   80.0 |  87.66% |      4454 | 3.64M | 1062 |        54 | 83.86K | 2217 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:28 |
| Q20L60X40P000 |   40.0 |  95.25% |     12626 | 3.78M |  487 |        87 |  68.8K | 1389 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:29 |
| Q20L60X40P001 |   40.0 |  95.21% |     11259 | 3.78M |  522 |        64 | 66.07K | 1489 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:28 |
| Q20L60X40P002 |   40.0 |  94.97% |     11837 | 3.79M |  494 |        60 | 58.67K | 1404 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:28 |
| Q20L60X60P000 |   60.0 |  91.86% |      7319 | 3.74M |  736 |        57 | 66.89K | 1658 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:28 |
| Q20L60X60P001 |   60.0 |  92.37% |      7296 | 3.77M |  763 |        44 | 59.24K | 1708 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:29 |
| Q20L60X60P002 |   60.0 |  92.38% |      7062 | 3.76M |  773 |        54 | 67.12K | 1742 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:30 |
| Q20L60X80P000 |   80.0 |  87.88% |      4684 | 3.64M | 1017 |        49 | 79.96K | 2132 |   78.0 | 17.0 |   9.0 | 258.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:27 |
| Q20L60X80P001 |   80.0 |  87.68% |      4486 | 3.63M | 1039 |        64 | 88.21K | 2190 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:30 |
| Q20L60X80P002 |   80.0 |  88.43% |      4530 | 3.65M | 1039 |        56 | 86.04K | 2196 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:28 |
| Q25L60X40P000 |   40.0 |  96.57% |     17592 |  3.8M |  370 |        68 | 53.53K | 1168 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:30 |
| Q25L60X40P001 |   40.0 |  96.04% |     15713 |  3.8M |  415 |        60 | 48.03K | 1181 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:28 |
| Q25L60X40P002 |   40.0 |  96.16% |     16381 | 3.79M |  385 |        69 | 53.33K | 1193 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:30 |
| Q25L60X60P000 |   60.0 |  94.16% |     11402 | 3.79M |  523 |       111 | 56.91K | 1170 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:28 |
| Q25L60X60P001 |   60.0 |  94.70% |     11214 |  3.8M |  539 |        68 | 55.06K | 1293 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:30 |
| Q25L60X60P002 |   60.0 |  94.82% |     11448 |  3.8M |  515 |        53 |    47K | 1232 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:30 |
| Q25L60X80P000 |   80.0 |  92.58% |      7725 | 3.75M |  716 |        93 | 76.32K | 1514 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:29 |
| Q25L60X80P001 |   80.0 |  92.62% |      7944 | 3.73M |  690 |       157 | 85.62K | 1485 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  96.59% |     17930 |  3.8M |  345 |        66 | 49.81K | 1156 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:29 |
| Q30L60X40P001 |   40.0 |  96.79% |     19267 | 3.81M |  350 |        58 | 48.56K | 1172 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:32 |
| Q30L60X40P002 |   40.0 |  96.48% |     18063 |  3.8M |  361 |        62 | 51.31K | 1181 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:31 |
| Q30L60X60P000 |   60.0 |  95.42% |     13524 | 3.79M |  472 |       335 | 70.99K | 1184 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:30 |
| Q30L60X60P001 |   60.0 |  95.43% |     12277 | 3.82M |  478 |        49 | 43.51K | 1153 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:31 |
| Q30L60X60P002 |   60.0 |  95.62% |     14021 | 3.81M |  447 |        73 | 46.27K | 1085 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:29 |
| Q30L60X80P000 |   80.0 |  93.88% |      9080 | 3.77M |  608 |        78 | 64.64K | 1336 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:30 |
| Q30L60X80P001 |   80.0 |  93.84% |      9579 | 3.78M |  616 |        67 | 62.92K | 1366 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:29 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.01% |     44087 | 3.82M | 196 |       174 | 41.11K | 468 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:31 |
| MRX40P001 |   40.0 |  96.64% |     42392 | 3.82M | 191 |       404 | 40.35K | 454 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:29 |
| MRX40P002 |   40.0 |  96.89% |     41140 | 3.82M | 192 |       452 | 46.81K | 450 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:30 |
| MRX60P000 |   60.0 |  96.46% |     39324 | 3.83M | 200 |       485 |  37.1K | 429 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:28 |
| MRX60P001 |   60.0 |  96.53% |     39303 | 3.83M | 197 |       653 | 44.66K | 431 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:28 |
| MRX80P000 |   80.0 |  96.01% |     32117 | 3.82M | 221 |       795 | 45.61K | 455 |   80.0 | 13.0 |  13.7 | 238.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:29 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.87% |     16453 | 3.81M | 375 |       130 | 58.16K |  982 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:30 |
| Q0L0X40P001   |   40.0 |  95.98% |     16687 | 3.82M | 388 |       116 |  61.6K | 1058 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:29 |
| Q0L0X40P002   |   40.0 |  95.98% |     16186 | 3.81M | 393 |       156 | 64.38K | 1055 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q0L0X60P000   |   60.0 |  94.50% |     10332 | 3.79M | 541 |       454 | 72.19K | 1248 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q0L0X60P001   |   60.0 |  94.29% |     10166 |  3.8M | 585 |       127 | 70.63K | 1331 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q0L0X60P002   |   60.0 |  94.46% |     10016 |  3.8M | 570 |       231 | 67.42K | 1288 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:29 |
| Q0L0X80P000   |   80.0 |  94.37% |      8240 | 3.77M | 650 |       119 | 90.87K | 1699 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:32 |
| Q0L0X80P001   |   80.0 |  94.55% |      7978 | 3.78M | 668 |       109 | 89.34K | 1759 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:31 |
| Q0L0X80P002   |   80.0 |  94.41% |      8574 | 3.78M | 664 |       101 | 86.94K | 1711 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q20L60X40P000 |   40.0 |  96.03% |     17422 | 3.81M | 361 |       475 | 63.67K |  951 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:30 |
| Q20L60X40P001 |   40.0 |  95.82% |     16030 |  3.8M | 392 |       552 | 66.85K | 1006 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q20L60X40P002 |   40.0 |  95.97% |     16451 | 3.81M | 382 |       851 | 74.34K | 1013 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q20L60X60P000 |   60.0 |  94.39% |     11254 | 3.79M | 535 |       156 | 65.23K | 1210 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:29 |
| Q20L60X60P001 |   60.0 |  94.64% |     10229 | 3.81M | 547 |       104 | 60.83K | 1241 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:29 |
| Q20L60X60P002 |   60.0 |  94.50% |     10237 |  3.8M | 552 |        94 | 60.53K | 1248 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q20L60X80P000 |   80.0 |  94.08% |      8568 | 3.77M | 680 |       536 | 98.68K | 1709 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| Q20L60X80P001 |   80.0 |  94.27% |      8060 | 3.78M | 684 |       106 | 88.24K | 1743 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| Q20L60X80P002 |   80.0 |  94.55% |      8476 | 3.77M | 651 |       197 | 94.63K | 1712 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:31 |
| Q25L60X40P000 |   40.0 |  96.83% |     23379 | 3.81M | 313 |       227 | 55.65K |  861 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q25L60X40P001 |   40.0 |  96.60% |     22236 | 3.81M | 329 |       100 | 48.25K |  872 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q25L60X40P002 |   40.0 |  96.75% |     20290 | 3.81M | 298 |       311 | 54.39K |  839 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:30 |
| Q25L60X60P000 |   60.0 |  95.93% |     15910 | 3.82M | 395 |       486 | 58.84K |  925 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:31 |
| Q25L60X60P001 |   60.0 |  95.83% |     15511 | 3.82M | 403 |       119 | 50.17K |  932 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q25L60X60P002 |   60.0 |  95.89% |     15891 | 3.81M | 379 |       407 | 54.43K |  912 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:28 |
| Q25L60X80P000 |   80.0 |  95.96% |     13293 | 3.79M | 474 |       429 | 78.66K | 1201 |   81.0 | 16.0 |  11.0 | 258.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q25L60X80P001 |   80.0 |  95.90% |     13188 | 3.79M | 468 |       557 | 84.15K | 1194 |   81.0 | 16.0 |  11.0 | 258.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  97.07% |     25919 | 3.82M | 277 |       393 | 54.19K |  856 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:29 |
| Q30L60X40P001 |   40.0 |  97.06% |     24132 | 3.82M | 289 |        84 | 44.68K |  836 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  96.88% |     22225 | 3.81M | 299 |       203 | 50.72K |  833 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q30L60X60P000 |   60.0 |  96.31% |     18645 | 3.81M | 366 |       648 | 68.39K |  886 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:30 |
| Q30L60X60P001 |   60.0 |  96.50% |     16391 | 3.83M | 375 |        94 | 48.04K |  911 |   60.0 | 12.0 |   8.0 | 192.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q30L60X60P002 |   60.0 |  96.37% |     18043 | 3.82M | 355 |       614 | 50.59K |  832 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:29 |
| Q30L60X80P000 |   80.0 |  96.16% |     15336 | 3.81M | 422 |       260 | 65.34K | 1034 |   81.0 | 16.0 |  11.0 | 258.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q30L60X80P001 |   80.0 |  96.24% |     14813 | 3.82M | 425 |        96 | 57.75K | 1049 |   81.0 | 16.0 |  11.0 | 258.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:32 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.57% |     46231 | 3.83M | 191 |       765 | 47.99K | 402 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:27 |
| MRX40P001 |   40.0 |  96.73% |     44195 | 3.83M | 186 |       765 | 46.72K | 393 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:29 |
| MRX40P002 |   40.0 |  96.84% |     44832 | 3.82M | 187 |       825 | 51.97K | 390 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:29 |
| MRX60P000 |   60.0 |  96.66% |     39474 | 3.79M | 190 |       863 | 51.49K | 392 |   61.0 | 10.0 |  10.3 | 182.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:29 |
| MRX60P001 |   60.0 |  96.74% |     44856 | 3.61M | 173 |       980 |  48.7K | 362 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:30 |
| MRX80P000 |   80.0 |  96.56% |     40400 | 3.83M | 188 |       800 | 45.86K | 410 |   81.0 | 13.0 |  14.0 | 240.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:30 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.36% |     46540 | 3.83M | 171 |      1903 |  678.8K | 353 |  246.0 | 40.0 |  42.0 | 732.0 |   0:00:49 |
| 7_merge_mr_unitigs_bcalm      |  97.86% |     44132 | 3.83M | 185 |      1088 |  67.52K |  66 |  247.0 | 40.0 |  42.3 | 734.0 |   0:00:59 |
| 7_merge_mr_unitigs_superreads |  97.75% |     46522 | 3.83M | 179 |      1066 |  76.85K |  77 |  246.0 | 41.0 |  41.0 | 738.0 |   0:00:57 |
| 7_merge_mr_unitigs_tadpole    |  97.70% |     46522 | 3.83M | 176 |      1046 |  68.93K |  68 |  248.0 | 42.0 |  40.7 | 748.0 |   0:00:56 |
| 7_merge_unitigs_bcalm         |  98.62% |     58738 | 3.84M | 162 |      1530 | 355.69K | 227 |  242.0 | 44.0 |  36.7 | 748.0 |   0:01:20 |
| 7_merge_unitigs_superreads    |  98.40% |     45086 | 3.84M | 177 |      1293 |    344K | 275 |  241.0 | 45.0 |  35.3 | 752.0 |   0:01:17 |
| 7_merge_unitigs_tadpole       |  98.66% |     58812 | 3.84M | 162 |      1840 | 432.71K | 241 |  244.0 | 44.0 |  37.3 | 752.0 |   0:01:22 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  94.17% |     12021 | 3.78M | 491 |       840 | 91.35K | 939 |  245.0 | 43.0 |  38.7 | 748.0 |   0:00:38 |
| 8_mr_spades  |  98.55% |     71096 | 3.12M | 126 |       911 | 49.67K | 227 |  134.0 | 20.0 |  24.7 | 388.0 |   0:00:32 |
| 8_megahit    |  97.18% |     47647 | 3.29M | 146 |      1034 |  46.7K | 225 |  247.0 | 42.0 |  40.3 | 746.0 |   0:00:36 |
| 8_mr_megahit |  98.92% |     92324 | 2.88M | 104 |      1066 | 50.95K | 188 |  134.0 | 20.0 |  24.7 | 388.0 |   0:00:30 |
| 8_platanus   |  96.56% |     52703 | 3.85M | 166 |       162 | 40.15K | 326 |  248.0 | 41.0 |  41.7 | 742.0 |   0:00:36 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 2961149 | 4033464 |    2 |
| Paralogs                 |    3424 |  119270 |   49 |
| Repetitives              |    1070 |  120471 |  244 |
| 7_merge_anchors.anchors  |   46540 | 3834580 |  171 |
| 7_merge_anchors.others   |    1903 |  678804 |  353 |
| glue_anchors             |   58984 | 3833948 |  152 |
| fill_anchors             |  114960 | 3866303 |   95 |
| spades.contig            |   12898 | 4165902 | 2176 |
| spades.scaffold          |   12986 | 4166848 | 2166 |
| spades.non-contained     |   13389 | 3870353 |  453 |
| mr_spades.contig         |  112695 | 3954692 |  204 |
| mr_spades.scaffold       |  124669 | 3955138 |  198 |
| mr_spades.non-contained  |  112695 | 3920861 |  105 |
| megahit.contig           |   87114 | 3941693 |  182 |
| megahit.non-contained    |   87114 | 3906692 |  108 |
| mr_megahit.contig        |  203715 | 3976930 |  171 |
| mr_megahit.non-contained |  203715 | 3938723 |   89 |
| platanus.contig          |   45170 | 3985546 |  406 |
| platanus.scaffold        |   59030 | 3934384 |  284 |
| platanus.non-contained   |   59043 | 3891546 |  165 |

