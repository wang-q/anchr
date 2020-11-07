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

anchr template \
    --genome 5432652 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --tile" \
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
| R.genome.bbtools  | 589.6 |    583 | 726.0 |                         96.28% |
| R.tadpole.bbtools | 568.3 |    575 | 154.4 |                         87.24% |
| R.genome.picard   | 582.1 |    585 | 146.5 |                             FR |
| R.tadpole.picard  | 574.2 |    578 | 147.3 |                             FR |


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

| Name       |     N50 |     Sum |       # |
|:-----------|--------:|--------:|--------:|
| Genome     | 5224283 | 5432652 |       2 |
| Paralogs   |    2295 |  220468 |     101 |
| Illumina.R |     251 | 481.02M | 2080000 |
| trim.R     |     250 | 403.34M | 1804252 |
| Q20L60     |     250 | 395.62M | 1753944 |
| Q25L60     |     250 | 378.46M | 1701397 |
| Q30L60     |     250 | 343.34M | 1606801 |


Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 480.99M | 2079856 |
| filteredbytile | 251 | 461.81M | 1998582 |
| trim           | 250 | 403.38M | 1804444 |
| filter         | 250 | 403.34M | 1804252 |
| R1             | 250 | 208.73M |  902126 |
| R2             | 247 | 194.61M |  902126 |
| Rs             |   0 |       0 |       0 |


```text
#R.trim
#Matched	5538	0.27710%
#Name	Reads	ReadsPct
Reverse_adapter	4568	0.22856%
```

```text
#R.filter
#Matched	99	0.00549%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	14060984
#error_kmers	8846626
#genomic_kmers	5214358
#main_peak	68
#genome_size_in_peaks	5295678
#genome_size	5414375
#haploid_genome_size	5414375
#fold_coverage	68
#haploid_fold_coverage	68
#ploidy	1
#percent_repeat_in_peaks	1.549
#percent_repeat	1.999
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 250 | 403.34M | 1804252 |
| ecco          | 250 | 403.33M | 1804252 |
| eccc          | 250 | 403.33M | 1804252 |
| ecct          | 250 | 399.27M | 1782610 |
| extended      | 290 | 469.67M | 1782610 |
| merged.raw    | 586 | 315.63M |  582457 |
| unmerged.raw  | 285 | 149.35M |  617696 |
| unmerged.trim | 285 | 149.34M |  617674 |
| M1            | 586 |  315.6M |  582402 |
| U1            | 290 |  79.75M |  308837 |
| U2            | 270 |  69.59M |  308837 |
| Us            |   0 |       0 |       0 |
| M.cor         | 518 | 465.52M | 1782478 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 362.6 |    389 |  97.8 |         19.13% |
| M.ihist.merge.txt  | 541.9 |    564 | 120.1 |         65.35% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   |  74.2 |   64.6 |   13.01% | "127" | 5.43M | 5.35M |     0.98 | 0:00'47'' |
| Q20L60.R |  72.8 |   64.5 |   11.47% | "127" | 5.43M | 5.35M |     0.98 | 0:00'46'' |
| Q25L60.R |  69.7 |   63.6 |    8.70% | "127" | 5.43M | 5.34M |     0.98 | 0:00'45'' |
| Q30L60.R |  63.2 |   59.6 |    5.75% | "127" | 5.43M | 5.34M |     0.98 | 0:00'41'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  97.27% |     16489 | 5.07M | 528 |       363 | 109.31K | 1148 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:01:07 |   0:00:38 |
| Q0L0X50P000    |   50.0 |  97.30% |     17564 | 5.14M | 482 |       161 |  92.65K | 1128 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:40 |
| Q0L0X60P000    |   60.0 |  97.30% |     17715 |  5.2M | 474 |       123 |   86.1K | 1121 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:39 |
| Q0L0XallP000   |   64.6 |  97.30% |     17089 | 5.22M | 496 |       122 |  90.88K | 1133 |   64.0 | 7.0 |  14.3 | 127.5 | "31,41,51,61,71,81" |   0:01:30 |   0:00:39 |
| Q20L60X40P000  |   40.0 |  97.52% |     20610 | 5.16M | 450 |       305 |  85.51K | 1038 |   39.0 | 5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:37 |
| Q20L60X50P000  |   50.0 |  97.50% |     19000 | 5.24M | 467 |       159 |  87.42K | 1067 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:39 |
| Q20L60X60P000  |   60.0 |  97.52% |     20202 | 5.24M | 443 |       122 |  81.11K | 1058 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:41 |
| Q20L60XallP000 |   64.5 |  97.51% |     20693 | 5.23M | 424 |       133 |  78.27K | 1038 |   63.0 | 8.0 |  13.0 | 126.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:42 |
| Q25L60X40P000  |   40.0 |  97.85% |     17419 | 4.66M | 443 |       416 |  98.22K | 1008 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:01:03 |   0:00:39 |
| Q25L60X50P000  |   50.0 |  97.88% |     20603 | 5.02M | 402 |       244 |  84.04K |  988 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:42 |
| Q25L60X60P000  |   60.0 |  97.93% |     21736 |  5.1M | 386 |       156 |  79.43K |  983 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:43 |
| Q25L60XallP000 |   63.6 |  97.97% |     22251 | 5.08M | 368 |       113 |  68.93K |  959 |   63.0 | 8.0 |  13.0 | 126.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:42 |
| Q30L60X40P000  |   40.0 |  98.20% |     19616 |  4.9M | 395 |       422 |   84.9K |  958 |   40.0 | 5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:40 |
| Q30L60X50P000  |   50.0 |  98.27% |     20721 | 4.89M | 381 |       205 |  73.35K |  957 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:40 |
| Q30L60XallP000 |   59.6 |  98.34% |     22476 | 4.96M | 360 |       153 |   73.1K |  952 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:41 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.40% |     31851 | 5.31M | 284 |       166 | 46.15K | 568 |   39.0 |  5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:35 |
| MRX40P001  |   40.0 |  97.47% |     29651 | 5.21M | 281 |       123 | 43.34K | 586 |   39.0 |  5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:35 |
| MRX50P000  |   50.0 |  97.36% |     31689 | 5.29M | 281 |       142 | 46.28K | 590 |   49.0 |  6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:38 |
| MRX60P000  |   60.0 |  97.38% |     30363 | 5.31M | 275 |       124 | 41.73K | 587 |   59.0 |  8.0 |  11.7 | 118.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:37 |
| MRXallP000 |   85.7 |  97.31% |     32332 | 5.25M | 271 |       100 | 38.87K | 602 |   84.0 | 11.0 |  17.0 | 168.0 | "31,41,51,61,71,81" |   0:02:05 |   0:00:35 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  97.61% |     20482 |  5.2M | 420 |       119 | 78.17K | 1179 |   39.0 | 5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:39 |
| Q0L0X50P000    |   50.0 |  97.64% |     20711 | 5.21M | 408 |       124 | 74.05K | 1078 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:40 |
| Q0L0X60P000    |   60.0 |  97.64% |     21665 | 5.27M | 395 |       113 | 67.74K | 1022 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:42 |
| Q0L0XallP000   |   64.6 |  97.66% |     20697 | 5.24M | 420 |       118 | 77.23K | 1076 |   64.0 | 7.0 |  14.3 | 127.5 | "31,41,51,61,71,81" |   0:01:22 |   0:00:41 |
| Q20L60X40P000  |   40.0 |  97.59% |     20681 |  5.2M | 421 |       101 | 73.02K | 1190 |   39.0 | 5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:39 |
| Q20L60X50P000  |   50.0 |  97.58% |     18990 | 5.21M | 446 |       204 | 93.82K | 1124 |   49.0 | 5.0 |  11.3 |  96.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:39 |
| Q20L60X60P000  |   60.0 |  97.67% |     22525 | 5.29M | 388 |       150 | 72.84K | 1032 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:41 |
| Q20L60XallP000 |   64.5 |  97.69% |     20720 | 5.29M | 419 |       120 | 75.86K | 1061 |   64.0 | 7.0 |  14.3 | 127.5 | "31,41,51,61,71,81" |   0:01:21 |   0:00:42 |
| Q25L60X40P000  |   40.0 |  97.66% |     20215 |  5.2M | 423 |       158 | 81.78K | 1202 |   39.0 | 5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:39 |
| Q25L60X50P000  |   50.0 |  97.73% |     21430 | 5.14M | 388 |       132 | 71.74K | 1063 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:39 |
| Q25L60X60P000  |   60.0 |  97.75% |     22350 | 5.26M | 382 |       107 | 65.51K | 1004 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:40 |
| Q25L60XallP000 |   63.6 |  97.80% |     22379 | 5.23M | 381 |       100 | 64.88K | 1011 |   63.0 | 7.0 |  14.0 | 126.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:40 |
| Q30L60X40P000  |   40.0 |  98.06% |     20224 | 5.23M | 425 |       160 | 92.71K | 1267 |   40.0 | 5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:39 |
| Q30L60X50P000  |   50.0 |  98.16% |     21453 | 5.24M | 387 |        96 | 67.65K | 1136 |   50.0 | 6.0 |  10.7 | 100.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:40 |
| Q30L60XallP000 |   59.6 |  98.21% |     23704 | 5.27M | 370 |       109 | 67.43K | 1048 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:41 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.32% |     33067 | 5.25M | 260 |       162 | 39.19K | 511 |   39.0 |  5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:33 |
| MRX40P001  |   40.0 |  97.56% |     32733 | 5.25M | 253 |       133 | 34.75K | 522 |   39.0 |  5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:34 |
| MRX50P000  |   50.0 |  97.35% |     37668 |  5.3M | 249 |       161 | 36.57K | 495 |   49.0 |  6.5 |   9.8 |  98.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:36 |
| MRX60P000  |   60.0 |  97.41% |     34247 |  5.3M | 250 |       142 | 36.73K | 515 |   58.0 |  7.0 |  12.3 | 116.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:36 |
| MRXallP000 |   85.7 |  97.23% |     37916 | 5.25M | 246 |       125 | 32.52K | 494 |   84.0 | 10.5 |  17.5 | 168.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:35 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  98.26% |     19303 | 5.02M | 443 |       251 | 92.78K | 1013 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:40 |
| Q0L0X50P000    |   50.0 |  98.25% |     19083 | 5.07M | 451 |       202 |  89.9K | 1030 |   49.0 | 5.0 |  11.3 |  96.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:41 |
| Q0L0X60P000    |   60.0 |  98.17% |     20726 | 5.17M | 410 |       117 | 71.93K |  981 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:41 |
| Q0L0XallP000   |   64.6 |  98.17% |     20301 | 5.19M | 437 |       120 | 78.55K | 1014 |   64.0 | 7.0 |  14.3 | 127.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:41 |
| Q20L60X40P000  |   40.0 |  98.21% |     21307 |  5.2M | 411 |       435 | 80.66K |  944 |   39.0 | 5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:38 |
| Q20L60X50P000  |   50.0 |  98.25% |     20715 | 5.26M | 424 |       188 | 79.92K |  974 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:41 |
| Q20L60X60P000  |   60.0 |  98.25% |     21660 | 5.29M | 409 |       127 | 76.17K |  993 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:42 |
| Q20L60XallP000 |   64.5 |  98.23% |     21847 | 5.22M | 392 |       145 | 73.05K |  968 |   63.0 | 8.0 |  13.0 | 126.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:40 |
| Q25L60X40P000  |   40.0 |  98.48% |     20692 | 5.11M | 401 |       381 | 79.61K |  947 |   39.0 | 5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:40 |
| Q25L60X50P000  |   50.0 |  98.47% |     21734 | 5.05M | 382 |       226 | 78.48K |  939 |   49.0 | 6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:42 |
| Q25L60X60P000  |   60.0 |  98.47% |     22577 |  5.2M | 376 |       158 | 75.12K |  946 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:43 |
| Q25L60XallP000 |   63.6 |  98.48% |     24097 | 5.16M | 356 |       114 | 65.77K |  919 |   63.0 | 8.0 |  13.0 | 126.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:41 |
| Q30L60X40P000  |   40.0 |  98.53% |     19313 | 5.13M | 448 |       578 | 97.81K | 1001 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:39 |   0:00:40 |
| Q30L60X50P000  |   50.0 |  98.54% |     21932 | 5.24M | 378 |       158 | 69.51K |  948 |   50.0 | 6.0 |  10.7 | 100.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:40 |
| Q30L60XallP000 |   59.6 |  98.56% |     23752 | 5.13M | 359 |       153 | 68.37K |  915 |   59.0 | 7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:43 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.35% |     33069 |  5.3M | 263 |       218 | 38.14K | 486 |   39.0 |  5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:42 |
| MRX40P001  |   40.0 |  97.74% |     33199 | 5.25M | 252 |       218 | 35.98K | 475 |   39.0 |  5.0 |   8.0 |  78.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:42 |
| MRX50P000  |   50.0 |  97.33% |     36661 |  5.3M | 256 |       195 | 39.22K | 492 |   49.0 |  6.0 |  10.3 |  98.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:37 |
| MRX60P000  |   60.0 |  97.28% |     34247 |  5.3M | 256 |       146 | 35.27K | 502 |   59.0 |  7.0 |  12.7 | 118.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:37 |
| MRXallP000 |   85.7 |  97.23% |     34494 | 5.25M | 253 |       122 | 32.81K | 510 |   84.0 | 11.0 |  17.0 | 168.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:36 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  97.06% |     26074 | 5.26M | 326 |      1276 | 191.72K | 131 |   64.0 | 8.0 |  13.3 | 128.0 |   0:00:49 |
| 7_merge_mr_unitigs_bcalm      |  97.93% |     25836 | 5.14M | 321 |      1049 |  34.75K |  34 |   64.0 | 8.0 |  13.3 | 128.0 |   0:00:45 |
| 7_merge_mr_unitigs_superreads |  97.72% |     24492 | 5.08M | 328 |      1398 |  51.76K |  35 |   64.0 | 9.0 |  12.3 | 128.0 |   0:00:45 |
| 7_merge_mr_unitigs_tadpole    |  97.84% |     25284 | 4.98M | 316 |      1131 |  39.16K |  36 |   64.0 | 8.0 |  13.3 | 128.0 |   0:00:55 |
| 7_merge_unitigs_bcalm         |  97.86% |     24483 | 5.27M | 352 |      7272 | 159.01K |  74 |   64.0 | 8.0 |  13.3 | 128.0 |   0:00:50 |
| 7_merge_unitigs_superreads    |  97.84% |     23742 | 5.27M | 366 |      1031 |  86.01K |  91 |   64.0 | 8.0 |  13.3 | 128.0 |   0:00:51 |
| 7_merge_unitigs_tadpole       |  97.82% |     24492 | 5.27M | 352 |      1287 | 131.71K |  91 |   64.0 | 8.0 |  13.3 | 128.0 |   0:00:52 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |    Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|-------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.96% |     36330 |  1.32M |  61 |      1624 | 19.97K |  86 |   64.0 |  9.5 |  11.8 | 128.0 |   0:00:36 |
| 8_mr_spades  |  98.77% |     73928 |  5.28M | 139 |      1186 | 17.81K | 234 |   84.0 | 12.0 |  16.0 | 168.0 |   0:00:34 |
| 8_megahit    |  98.68% |     27031 |  3.94M | 244 |       559 | 40.75K | 329 |   64.0 |  8.0 |  13.3 | 128.0 |   0:00:35 |
| 8_mr_megahit |  98.80% |     64033 |  5.29M | 161 |       973 | 22.28K | 292 |   84.0 | 11.0 |  17.0 | 168.0 |   0:00:34 |
| 8_platanus   |  97.36% |     34226 | 677.2K |  36 |      1283 | 10.63K |  48 |   64.0 |  9.0 |  12.3 | 128.0 |   0:00:33 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 5224283 | 5432652 |   2 |
| Paralogs                 |    2295 |  220468 | 101 |
| 7_merge_anchors.anchors  |   26074 | 5261113 | 326 |
| 7_merge_anchors.others   |    1276 |  191715 | 131 |
| glue_anchors             |   26074 | 5260914 | 323 |
| fill_anchors             |   72010 | 5300408 | 124 |
| spades.contig            |  209173 | 5367611 | 155 |
| spades.scaffold          |  285140 | 5368010 | 140 |
| spades.non-contained     |  246457 | 5350249 |  58 |
| mr_spades.contig         |  100015 | 5367784 | 124 |
| mr_spades.scaffold       |  284296 | 5374450 |  62 |
| mr_spades.non-contained  |  100015 | 5361229 | 103 |
| megahit.contig           |   59600 | 5361716 | 209 |
| megahit.non-contained    |   59600 | 5340827 | 156 |
| mr_megahit.contig        |   75019 | 5382321 | 165 |
| mr_megahit.non-contained |   75019 | 5369915 | 139 |
| platanus.contig          |   18988 | 5423306 | 677 |
| platanus.scaffold        |  283954 | 5400050 | 274 |
| platanus.non-contained   |  283954 | 5346127 |  39 |


# *Mycobacterium abscessus* 6G-0125-R

## Mabs_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Mabs_mi/1_genome
cd ~/data/anchr/Mabs_mi/1_genome

cp ~/data/anchr/ref/Mabs/genome.fa .
cp ~/data/anchr/ref/Mabs/paralogs.fa .

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

anchr template \
    --genome 5090491 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --tile" \
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

# *Rhodobacter sphaeroides* 2.4.1

## Rsph_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Rsph_mi/1_genome
cd ~/data/anchr/Rsph_mi/1_genome

cp ~/data/anchr/ref/Rsph/genome.fa .
cp ~/data/anchr/ref/Rsph/paralogs.fa .

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
    --unitigger "superreads bcalm" \
    --statp 2 \
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

anchr template \
    --genome 4033464 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --tile" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 50 all" \
    --unitigger "superreads bcalm" \
    --statp 2 \
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

