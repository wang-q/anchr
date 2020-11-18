# Assemble genomes from GAGE-B data sets

[TOC levels=1-3]: # ""

- [Assemble genomes from GAGE-B data sets](#assemble-genomes-from-gage-b-data-sets)
- [*Bacillus cereus* ATCC 10987](#bacillus-cereus-atcc-10987)
  - [Bcer_mi: reference](#bcer_mi-reference)
  - [Bcer_mi: download](#bcer_mi-download)
  - [Bcer_mi: template](#bcer_mi-template)
  - [Bcer_mi: run](#bcer_mi-run)
- [*Mycobacterium abscessus* 6G-0125-R](#mycobacterium-abscessus-6g-0125-r)
  - [Mabs_mi: reference](#mabs_mi-reference)
  - [Mabs_mi: download](#mabs_mi-download)
  - [Mabs_mi: template](#mabs_mi-template)
  - [Mabs_mi: run](#mabs_mi-run)
- [*Rhodobacter sphaeroides* 2.4.1](#rhodobacter-sphaeroides-241)
  - [Rsph_mi: reference](#rsph_mi-reference)
  - [Rsph_mi: download](#rsph_mi-download)
  - [Rsph_mi: template](#rsph_mi-template)
  - [Rsph_mi: run](#rsph_mi-run)
- [*Vibrio cholerae* CP1032(5)](#vibrio-cholerae-cp10325)
  - [Vcho_mi: reference](#vcho_mi-reference)
  - [Vcho_mi: download](#vcho_mi-download)
  - [Vcho_mi: template](#vcho_mi-template)
  - [Vcho_mi: run](#vcho_mi-run)

* Rsync to hpcc

```shell script
for D in Bcer_mi Mabs_mi Rsph_mi Vcho_mi; do
    rsync -avP \
        ~/data/anchr/${D}/ \
        wangq@202.119.37.251:data/anchr/${D}
done

```


# *Bacillus cereus* ATCC 10987

## Bcer_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Bcer_mi/1_genome
cd ~/data/anchr/Bcer_mi/1_genome

cp ~/data/anchr/ref/Bcer/genome.fa .
cp ~/data/anchr/ref/Bcer/paralogs.fa .
cp ~/data/anchr/ref/Bcer/repetitives.fa .

```

## Bcer_mi: download

* Illumina

```shell script
cd ~/data/anchr/Bcer_mi

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/datasets/B_cereus_MiSeq.tar.gz

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
cd ~/data/anchr/Bcer_mi

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/B_cereus_MiSeq.tar.gz

tar xvfz B_cereus_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz mira_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz sga_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz soap_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz spades_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz velvet_ctg.fasta

```

## Bcer_mi: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_mi

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

## Bcer_mi: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_mi

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Bcer_mi bash 0_bsub.sh
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

## Mabs_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Mabs_mi/1_genome
cd ~/data/anchr/Mabs_mi/1_genome

cp ~/data/anchr/ref/Mabs/genome.fa .
cp ~/data/anchr/ref/Mabs/paralogs.fa .
cp ~/data/anchr/ref/Mabs/repetitives.fa .

```

## Mabs_mi: download

* Illumina

```shell script
cd ~/data/anchr/Mabs_mi

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/datasets/M_abscessus_MiSeq.tar.gz

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
cd ~/data/anchr/Mabs_mi

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/M_abscessus_MiSeq.tar.gz

tar xvfz M_abscessus_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz mira_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz sga_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz soap_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz spades_ctg.fasta
tar xvfz M_abscessus_MiSeq.tar.gz velvet_ctg.fasta

```

## Mabs_mi: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Mabs_mi

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

## Mabs_mi: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Mabs_mi

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Mabs_mi bash 0_bsub.sh
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

## Rsph_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Rsph_mi/1_genome
cd ~/data/anchr/Rsph_mi/1_genome

cp ~/data/anchr/ref/Rsph/genome.fa .
cp ~/data/anchr/ref/Rsph/paralogs.fa .
cp ~/data/anchr/ref/Rsph/repetitives.fa .

```

## Rsph_mi: download

* Illumina

```shell script
cd ~/data/anchr/Rsph_mi

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/datasets/R_sphaeroides_MiSeq.tar.gz

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
cd ~/data/anchr/Rsph_mi

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/R_sphaeroides_MiSeq.tar.gz

tar xvfz R_sphaeroides_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz mira_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz sga_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz soap_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz spades_ctg.fasta
tar xvfz R_sphaeroides_MiSeq.tar.gz velvet_ctg.fasta

```


## Rsph_mi: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_mi

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

## Rsph_mi: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_mi

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Rsph_mi bash 0_bsub.sh
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

## Vcho_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Vcho_mi/1_genome
cd ~/data/anchr/Vcho_mi/1_genome

cp ~/data/anchr/ref/Vcho/genome.fa .
cp ~/data/anchr/ref/Vcho/paralogs.fa .
cp ~/data/anchr/ref/Vcho/repetitives.fa .

```

## Vcho_mi: download

* Illumina

```shell script
cd ~/data/anchr/Vcho_mi

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/datasets/V_cholerae_MiSeq.tar.gz

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
cd ~/data/anchr/Vcho_mi

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/V_cholerae_MiSeq.tar.gz

tar xvfz V_cholerae_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz mira_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz sga_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz soap_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz spades_ctg.fasta
tar xvfz V_cholerae_MiSeq.tar.gz velvet_ctg.fasta

```


## Vcho_mi: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Vcho_mi

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

## Vcho_mi: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Vcho_mi

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

