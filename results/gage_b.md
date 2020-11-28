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
    --trim "--dedupe --tile --cutoff 5 --cutk 31" \
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
| trim.R      |     250 | 404.36M | 1807384 |
| Q20L60      |     250 | 396.79M | 1758525 |
| Q25L60      |     250 | 379.57M | 1706208 |
| Q30L60      |     250 | 344.25M | 1611233 |


Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 480.99M | 2079856 |
| filteredbytile | 251 | 463.51M | 2005704 |
| highpass       | 251 | 459.84M | 1989878 |
| trim           | 250 | 404.41M | 1807576 |
| filter         | 250 | 404.36M | 1807384 |
| R1             | 250 | 209.26M |  903692 |
| R2             | 247 |  195.1M |  903692 |
| Rs             |   0 |       0 |       0 |


```text
#R.trim
#Matched	1247	0.06267%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	99	0.00548%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	13975806
#error_kmers	8743688
#genomic_kmers	5232118
#main_peak	64
#genome_size_in_peaks	5313826
#genome_size	5551203
#haploid_genome_size	5551203
#fold_coverage	64
#haploid_fold_coverage	64
#ploidy	1
#percent_repeat_in_peaks	1.571
#percent_repeat	2.617
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 250 | 404.36M | 1807382 |
| ecco          | 250 | 404.36M | 1807382 |
| eccc          | 250 | 404.36M | 1807382 |
| ecct          | 250 | 400.39M | 1786476 |
| extended      | 290 | 470.96M | 1786476 |
| merged.raw    | 586 | 316.57M |  584078 |
| unmerged.raw  | 285 | 149.75M |  618320 |
| unmerged.trim | 285 | 149.74M |  618292 |
| M1            | 586 | 316.54M |  584024 |
| U1            | 290 |  79.96M |  309146 |
| U2            | 270 |  69.78M |  309146 |
| Us            |   0 |       0 |       0 |
| M.cor         | 518 | 466.86M | 1786340 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 362.0 |    388 |  97.6 |         19.50% |
| M.ihist.merge.txt  | 542.0 |    564 | 120.0 |         65.39% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   |  74.4 |   64.8 |   12.98% | "127" | 5.43M | 5.35M |     0.98 | 0:01'00'' |
| Q20L60.R |  73.0 |   64.7 |   11.48% | "127" | 5.43M | 5.35M |     0.98 | 0:00'52'' |
| Q25L60.R |  69.9 |   63.8 |    8.71% | "127" | 5.43M | 5.34M |     0.98 | 0:00'50'' |
| Q30L60.R |  63.4 |   59.7 |    5.75% | "127" | 5.43M | 5.34M |     0.98 | 0:00'49'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  97.42% |     22565 | 5.31M | 364 |        47 | 30.01K |  959 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:38 |
| Q0L0X50P000    |   50.0 |  97.39% |     21975 |  5.3M | 392 |        60 | 42.38K | 1027 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:39 |
| Q0L0X60P000    |   60.0 |  97.37% |     20912 |  5.3M | 390 |        58 |  41.6K | 1024 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:39 |
| Q0L0XallP000   |   64.8 |  97.36% |     20488 | 5.29M | 414 |        63 |  49.8K | 1069 |   64.0 | 12.0 |   9.3 | 200.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:40 |
| Q20L60X40P000  |   40.0 |  97.66% |     24120 | 5.31M | 339 |        50 | 31.24K |  926 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:37 |
| Q20L60X50P000  |   50.0 |  97.55% |     23574 |  5.3M | 359 |        58 | 38.12K |  964 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:37 |
| Q20L60X60P000  |   60.0 |  97.57% |     23574 |  5.3M | 361 |        56 | 38.87K |  983 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:39 |
| Q20L60XallP000 |   64.7 |  97.54% |     22352 |  5.3M | 376 |        65 | 48.05K | 1001 |   64.0 | 12.0 |   9.3 | 200.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:38 |
| Q25L60X40P000  |   40.0 |  97.92% |     29564 | 5.33M | 302 |        55 | 31.35K |  896 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:37 |
| Q25L60X50P000  |   50.0 |  97.94% |     26549 | 5.33M | 317 |        60 | 38.15K |  933 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:39 |
| Q25L60X60P000  |   60.0 |  98.04% |     28543 | 5.33M | 309 |        54 |  36.2K |  932 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:44 |
| Q25L60XallP000 |   63.8 |  98.01% |     25870 | 5.31M | 317 |        58 | 38.98K |  948 |   63.0 | 12.0 |   9.0 | 198.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:41 |
| Q30L60X40P000  |   40.0 |  98.29% |     28356 | 5.31M | 305 |        59 | 35.63K |  878 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:37 |
| Q30L60X50P000  |   50.0 |  98.33% |     29701 | 5.33M | 297 |        57 | 34.69K |  898 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:40 |
| Q30L60XallP000 |   59.7 |  98.40% |     31700 | 5.33M | 286 |        57 | 36.04K |  904 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:44 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.62% |     37925 | 5.31M | 240 |        81 | 22.62K | 519 |   39.0 |  7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:34 |
| MRX40P001  |   40.0 |  97.54% |     40827 | 5.31M | 237 |        64 | 18.73K | 454 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:32 |
| MRX50P000  |   50.0 |  97.58% |     37545 | 5.31M | 249 |        90 | 22.64K | 547 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:36 |
| MRX60P000  |   60.0 |  97.53% |     37952 | 5.31M | 250 |        92 | 20.41K | 527 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:37 |
| MRXallP000 |   85.9 |  97.46% |     36340 | 5.33M | 253 |        90 |  18.9K | 537 |   84.0 | 15.0 |  13.0 | 258.0 | "31,41,51,61,71,81" |   0:02:00 |   0:00:35 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  97.59% |     24324 | 5.34M | 344 |        48 | 30.22K |  984 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:38 |
| Q0L0X50P000    |   50.0 |  97.45% |     22709 | 5.31M | 362 |        53 | 35.77K | 1003 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:37 |
| Q0L0X60P000    |   60.0 |  97.40% |     22353 | 5.31M | 375 |        48 | 33.76K | 1019 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:40 |
| Q0L0XallP000   |   64.8 |  97.38% |     21304 |  5.3M | 393 |        58 | 43.22K | 1063 |   64.0 | 12.0 |   9.3 | 200.0 | "31,41,51,61,71,81" |   0:01:44 |   0:00:39 |
| Q20L60X40P000  |   40.0 |  97.61% |     25791 | 5.31M | 330 |        44 | 26.42K |  899 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:36 |
| Q20L60X50P000  |   50.0 |  97.55% |     25614 |  5.3M | 344 |        47 | 30.05K |  955 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:38 |
| Q20L60X60P000  |   60.0 |  97.56% |     24098 | 5.31M | 355 |        50 | 33.34K |  986 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:39 |
| Q20L60XallP000 |   64.7 |  97.55% |     23567 |  5.3M | 367 |        58 | 41.72K | 1012 |   64.0 | 12.0 |   9.3 | 200.0 | "31,41,51,61,71,81" |   0:01:42 |   0:00:39 |
| Q25L60X40P000  |   40.0 |  97.77% |     29303 | 5.32M | 310 |        43 | 25.31K |  891 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:36 |
| Q25L60X50P000  |   50.0 |  97.78% |     28611 | 5.32M | 312 |        50 | 30.19K |  898 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:37 |
| Q25L60X60P000  |   60.0 |  97.80% |     28634 | 5.31M | 309 |        48 | 29.74K |  879 |   60.0 | 12.0 |   8.0 | 192.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:39 |
| Q25L60XallP000 |   63.8 |  97.86% |     26733 | 5.33M | 312 |        51 | 32.26K |  920 |   63.0 | 12.0 |   9.0 | 198.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:41 |
| Q30L60X40P000  |   40.0 |  98.09% |     26940 | 5.35M | 319 |        45 | 29.25K |  927 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:36 |
| Q30L60X50P000  |   50.0 |  98.18% |     28611 | 5.32M | 304 |        47 | 30.69K |  936 |   50.0 | 10.0 |   6.7 | 160.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:38 |
| Q30L60XallP000 |   59.7 |  98.20% |     30381 | 5.32M | 291 |        45 | 28.57K |  891 |   60.0 | 12.0 |   8.0 | 192.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:37 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.81% |     39857 | 5.34M | 232 |        87 | 18.83K | 442 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:32 |
| MRX40P001  |   40.0 |  97.77% |     42781 | 5.31M | 229 |        86 | 18.97K | 412 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:31 |
| MRX50P000  |   50.0 |  97.44% |     39857 | 5.34M | 236 |        80 | 16.34K | 466 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:33 |
| MRX60P000  |   60.0 |  97.42% |     41648 | 5.38M | 238 |        98 | 17.86K | 463 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:35 |
| MRXallP000 |   85.9 |  97.38% |     39857 | 5.31M | 246 |        97 | 16.29K | 489 |   85.0 | 15.0 |  13.3 | 260.0 | "31,41,51,61,71,81" |   0:01:54 |   0:00:35 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  98.33% |     29507 | 5.35M | 294 |        49 | 28.88K | 863 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:36 |
| Q0L0X50P000    |   50.0 |  98.25% |     27886 | 5.31M | 313 |        59 | 37.42K | 906 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:37 |
| Q0L0X60P000    |   60.0 |  98.18% |     27094 | 5.31M | 328 |        57 | 37.48K | 933 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:39 |
| Q0L0XallP000   |   64.8 |  98.17% |     25414 | 5.31M | 348 |        59 | 41.06K | 977 |   64.0 | 13.0 |   8.3 | 206.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:38 |
| Q20L60X40P000  |   40.0 |  98.29% |     30477 | 5.33M | 297 |        51 | 29.62K | 854 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:36 |
| Q20L60X50P000  |   50.0 |  97.90% |     29506 | 5.31M | 312 |        50 | 30.97K | 883 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |
| Q20L60X60P000  |   60.0 |  98.27% |     28025 | 5.31M | 318 |        55 | 36.07K | 917 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:41 |
| Q20L60XallP000 |   64.7 |  98.25% |     26737 | 5.31M | 330 |        56 | 37.95K | 947 |   64.0 | 13.0 |   8.3 | 206.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:41 |
| Q25L60X40P000  |   40.0 |  98.46% |     29701 | 5.33M | 300 |        64 | 36.28K | 862 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:36 |
| Q25L60X50P000  |   50.0 |  98.49% |     31133 | 5.33M | 294 |        56 | 34.34K | 870 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| Q25L60X60P000  |   60.0 |  98.50% |     32273 | 5.31M | 293 |        57 | 36.36K | 890 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:39 |
| Q25L60XallP000 |   63.8 |  98.52% |     29510 | 5.31M | 301 |        59 |  38.8K | 920 |   63.0 | 12.0 |   9.0 | 198.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:38 |
| Q30L60X40P000  |   40.0 |  98.81% |     28611 | 5.31M | 298 |        59 | 36.59K | 902 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:37 |
| Q30L60X50P000  |   50.0 |  98.61% |     31192 | 5.31M | 288 |        51 | 32.31K | 889 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:37 |
| Q30L60XallP000 |   59.7 |  98.61% |     32356 | 5.31M | 282 |        48 | 30.22K | 870 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:39 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.83% |     44440 | 5.32M | 217 |       107 |  18.4K | 404 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:32 |
| MRX40P001  |   40.0 |  97.77% |     44376 | 5.36M | 219 |        83 | 16.18K | 378 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:33 |
| MRX50P000  |   50.0 |  97.76% |     42799 | 5.34M | 224 |       131 | 19.98K | 444 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:32 |
| MRX60P000  |   60.0 |  97.77% |     42132 | 5.34M | 228 |       131 | 19.13K | 443 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:34 |
| MRXallP000 |   85.9 |  97.70% |     41648 | 5.32M | 239 |       120 | 18.37K | 484 |   85.0 | 15.0 |  13.3 | 260.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:35 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.76% |     36737 |  5.3M | 254 |     41729 | 291.14K | 43 |   64.0 | 12.0 |   9.3 | 200.0 |   0:00:42 |
| 7_merge_mr_unitigs_bcalm      |  98.15% |     37087 | 5.21M | 252 |     41793 |  70.11K |  9 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:49 |
| 7_merge_mr_unitigs_superreads |  97.83% |     34002 | 5.15M | 259 |      1091 |    6.8K |  7 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:43 |
| 7_merge_mr_unitigs_tadpole    |  98.14% |     42769 | 5.21M | 232 |     41634 |  53.78K | 12 |   63.0 | 13.0 |   8.0 | 204.0 |   0:00:43 |
| 7_merge_unitigs_bcalm         |  98.04% |     32860 | 5.31M | 273 |     36517 |  63.09K | 19 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:49 |
| 7_merge_unitigs_superreads    |  98.10% |     32344 | 5.31M | 273 |     41729 |  64.41K | 24 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:46 |
| 7_merge_unitigs_tadpole       |  98.02% |     32710 | 5.31M | 277 |     62295 |  81.41K | 21 |   64.0 | 12.0 |   9.3 | 200.0 |   0:00:45 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.96% |     49586 | 2.14M |  80 |       899 |  9.35K | 106 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:33 |
| 8_mr_spades  |  98.81% |     78041 | 5.35M | 130 |       217 | 10.11K | 212 |   85.0 | 15.0 |  13.3 | 260.0 |   0:00:36 |
| 8_megahit    |  98.69% |     37648 | 4.85M | 218 |       125 | 21.73K | 356 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:40 |
| 8_mr_megahit |  98.85% |     65403 | 5.36M | 153 |       451 | 13.53K | 265 |   85.0 | 15.0 |  13.3 | 260.0 |   0:00:36 |
| 8_platanus   |  97.38% |     43046 | 1.03M |  49 |        61 |  2.73K |  66 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:35 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 5224283 | 5432652 |   2 |
| Paralogs                 |    2295 |  220468 | 101 |
| Repetitives              |    2461 |  113050 | 173 |
| 7_merge_anchors.anchors  |   36737 | 5304696 | 254 |
| 7_merge_anchors.others   |   41729 |  291138 |  43 |
| glue_anchors             |   36737 | 5304175 | 253 |
| fill_anchors             |   61075 | 5317582 | 159 |
| spades.contig            |  207470 | 5366804 | 153 |
| spades.scaffold          |  285416 | 5367163 | 139 |
| spades.non-contained     |  207470 | 5349666 |  58 |
| mr_spades.contig         |  100015 | 5367895 | 128 |
| mr_spades.scaffold       |  284294 | 5374592 |  66 |
| mr_spades.non-contained  |  100015 | 5361433 | 105 |
| megahit.contig           |   59732 | 5360219 | 204 |
| megahit.non-contained    |   59732 | 5341409 | 158 |
| mr_megahit.contig        |   75019 | 5388594 | 186 |
| mr_megahit.non-contained |   75019 | 5369027 | 141 |
| platanus.contig          |   18988 | 5420878 | 665 |
| platanus.scaffold        |  284635 | 5398027 | 263 |
| platanus.non-contained   |  284635 | 5346865 |  39 |


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
    --redo \
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
| trim.R      |     177 | 280.98M | 1707706 |
| Q20L60      |     177 | 272.37M | 1643445 |
| Q25L60      |     174 | 248.49M | 1547604 |
| Q30L60      |     165 | 204.66M | 1371769 |


Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 511.87M | 2039328 |
| filteredbytile | 251 | 487.98M | 1944160 |
| highpass       | 251 | 473.94M | 1888212 |
| trim           | 177 | 281.94M | 1712272 |
| filter         | 177 | 280.98M | 1707706 |
| R1             | 187 | 150.12M |  853853 |
| R2             | 167 | 130.85M |  853853 |
| Rs             |   0 |       0 |       0 |


```text
#R.trim
#Matched	1394914	73.87486%
#Name	Reads	ReadsPct
Reverse_adapter	722797	38.27944%
pcr_dimer	390154	20.66262%
TruSeq_Universal_Adapter	112820	5.97496%
PCR_Primers	98601	5.22192%
TruSeq_Adapter_Index_1_6	46377	2.45613%
Nextera_LMP_Read2_External_Adapter	14169	0.75039%
TruSeq_Adapter_Index_11	4992	0.26438%
```

```text
#R.filter
#Matched	4566	0.26666%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	4566	0.26666%
```

```text
#R.peaks
#k	31
#unique_kmers	12872099
#error_kmers	7867805
#genomic_kmers	5004294
#main_peak	41
#genome_size_in_peaks	5000925
#genome_size	5549817
#haploid_genome_size	5549817
#fold_coverage	41
#haploid_fold_coverage	41
#ploidy	1
#percent_repeat_in_peaks	0.034
#percent_repeat	4.197
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 177 | 279.37M | 1697058 |
| ecco          | 177 | 279.31M | 1697058 |
| eccc          | 177 | 279.31M | 1697058 |
| ecct          | 176 | 270.46M | 1646896 |
| extended      | 214 | 335.98M | 1646896 |
| merged.raw    | 235 | 189.59M |  814700 |
| unmerged.raw  | 210 |   3.25M |   17496 |
| unmerged.trim | 210 |   3.25M |   17488 |
| M1            | 235 |  181.4M |  779644 |
| U1            | 230 |   1.85M |    8744 |
| U2            | 188 |    1.4M |    8744 |
| Us            |   0 |       0 |       0 |
| M.cor         | 234 | 185.43M | 1576776 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 190.5 |    186 |  46.4 |         92.94% |
| M.ihist.merge.txt  | 232.7 |    226 |  51.5 |         98.94% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   |  55.2 |   44.6 |   19.26% | "45" | 5.09M | 5.15M |     1.01 | 0:00'44'' |
| Q20L60.R |  53.5 |   43.9 |   17.89% | "45" | 5.09M | 5.15M |     1.01 | 0:00'40'' |
| Q25L60.R |  48.8 |   41.6 |   14.87% | "43" | 5.09M | 5.14M |     1.01 | 0:00'37'' |
| Q30L60.R |  40.2 |   35.5 |   11.70% | "39" | 5.09M | 5.13M |     1.01 | 0:00'34'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  96.15% |      7471 | 5.03M | 927 |        47 | 78.97K | 2241 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:41 |
| Q0L0XallP000   |   44.6 |  95.85% |      7168 | 5.02M | 969 |        38 | 65.84K | 2272 |   42.0 | 10.0 |   5.0 | 144.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:42 |
| Q20L60X40P000  |   40.0 |  96.50% |      7813 | 5.04M | 888 |        57 | 84.27K | 2152 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:34 |
| Q20L60XallP000 |   43.9 |  96.20% |      7434 | 5.03M | 927 |        40 | 66.51K | 2192 |   42.0 | 10.0 |   5.0 | 144.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:33 |
| Q25L60X40P000  |   40.0 |  97.52% |      9409 | 5.06M | 766 |        45 | 67.32K | 2003 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:33 |
| Q25L60XallP000 |   41.6 |  97.41% |      9360 | 5.06M | 773 |        45 | 65.86K | 2007 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:33 |
| Q30L60XallP000 |   35.5 |  98.65% |     13699 |  5.1M | 598 |        50 | 65.29K | 1839 |   34.0 |  8.0 |   5.0 | 116.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:36 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   36.4 |  99.02% |     38945 | 5.11M | 227 |        48 | 20.56K | 580 |   35.0 | 6.0 |   5.7 | 106.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:35 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  97.75% |     10064 | 5.07M | 729 |        38 | 71.79K | 2544 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:36 |
| Q0L0XallP000   |   44.6 |  97.38% |      8805 | 5.07M | 809 |        35 | 66.07K | 2520 |   43.0 | 10.0 |   5.0 | 146.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:36 |
| Q20L60X40P000  |   40.0 |  97.93% |     10279 | 5.08M | 723 |        43 | 80.96K | 2498 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:37 |
| Q20L60XallP000 |   43.9 |  97.54% |      9310 | 5.06M | 784 |        38 | 70.32K | 2509 |   42.0 |  9.0 |   5.0 | 138.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:37 |
| Q25L60X40P000  |   40.0 |  98.28% |     11749 | 5.09M | 675 |        38 | 71.29K | 2489 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:39 |
| Q25L60XallP000 |   41.6 |  98.14% |     11128 | 5.09M | 691 |        39 | 71.27K | 2464 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:37 |
| Q30L60XallP000 |   35.5 |  98.73% |     12172 | 5.09M | 647 |        37 | 69.56K | 2530 |   34.0 |  8.0 |   5.0 | 116.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:38 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   36.4 |  99.32% |     67083 | 5.14M | 139 |        47 | 12.72K | 404 |   35.0 | 6.0 |   5.7 | 106.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:36 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  98.46% |     15592 | 5.09M | 518 |        43 |  56.3K | 1705 |   38.0 | 9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:36 |
| Q0L0XallP000   |   44.6 |  98.24% |     13962 | 5.08M | 576 |        50 | 67.57K | 1759 |   43.0 | 9.0 |   5.3 | 140.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:37 |
| Q20L60X40P000  |   40.0 |  98.52% |     15444 |  5.1M | 529 |        46 | 58.31K | 1736 |   39.0 | 9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:38 |
| Q20L60XallP000 |   43.9 |  98.27% |     14178 | 5.09M | 573 |        41 | 52.84K | 1764 |   42.0 | 9.0 |   5.0 | 138.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:38 |
| Q25L60X40P000  |   40.0 |  98.80% |     16491 | 5.11M | 489 |        37 | 50.89K | 1788 |   38.0 | 9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:37 |
| Q25L60XallP000 |   41.6 |  98.73% |     16356 | 5.13M | 505 |        36 | 48.81K | 1771 |   40.0 | 9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:37 |
| Q30L60XallP000 |   35.5 |  99.14% |     17511 | 5.12M | 458 |        43 |  61.9K | 1901 |   34.0 | 8.0 |   5.0 | 116.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   36.4 |  99.35% |     79695 | 5.13M | 125 |        51 | 11.64K | 346 |   35.0 | 6.0 |   5.7 | 106.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:36 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |  # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|-------:|---:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  99.44% |     77345 | 5.11M | 133 |      6262 | 72.52K | 36 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:43 |
| 7_merge_mr_unitigs_bcalm      |  99.02% |     62462 |  5.1M | 158 |       696 |  2.29K |  3 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:32 |
| 7_merge_mr_unitigs_superreads |  98.69% |     34328 | 5.09M | 246 |      1093 |  6.79K |  7 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:33 |
| 7_merge_mr_unitigs_tadpole    |  99.06% |     66766 |  5.1M | 146 |      1092 |  3.55K |  4 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:32 |
| 7_merge_unitigs_bcalm         |  98.99% |     20670 | 5.11M | 380 |      6262 | 62.97K | 24 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:43 |
| 7_merge_unitigs_superreads    |  98.59% |     16997 |  5.1M | 492 |      1012 | 41.18K | 45 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:42 |
| 7_merge_unitigs_tadpole       |  99.19% |     29682 | 5.11M | 274 |      1007 | 23.36K | 25 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:44 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|----------:|
| 8_spades     |  99.32% |    100016 | 3.78M |  87 |       115 |  8.39K | 136 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:34 |
| 8_mr_spades  |  99.36% |    103525 | 5.11M |  93 |        64 |  7.19K | 178 |   36.0 | 6.0 |   6.0 | 108.0 |   0:00:33 |
| 8_megahit    |  99.20% |     76765 | 4.61M | 123 |        90 | 10.58K | 201 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:34 |
| 8_mr_megahit |  99.62% |     95731 | 5.13M |  96 |        62 |  5.85K | 185 |   36.0 | 6.0 |   6.0 | 108.0 |   0:00:30 |
| 8_platanus   |  98.28% |     23229 | 5.07M | 328 |        46 | 24.23K | 632 |   43.0 | 9.0 |   5.3 | 140.0 |   0:00:32 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 5067172 | 5090491 |   2 |
| Paralogs                 |    1693 |   83291 |  53 |
| Repetitives              |     192 |   15322 |  82 |
| 7_merge_anchors.anchors  |   77345 | 5113714 | 133 |
| 7_merge_anchors.others   |    6262 |   72521 |  36 |
| glue_anchors             |   87935 | 5112668 | 131 |
| fill_anchors             |  104681 | 5114246 | 101 |
| spades.contig            |  145002 | 5134933 |  89 |
| spades.scaffold          |  145002 | 5135123 |  85 |
| spades.non-contained     |  145002 | 5123061 |  62 |
| mr_spades.contig         |  110924 | 5135670 | 116 |
| mr_spades.scaffold       |  113271 | 5135970 | 113 |
| mr_spades.non-contained  |  110924 | 5121841 |  87 |
| megahit.contig           |  107304 | 5138499 | 142 |
| megahit.non-contained    |  107304 | 5118130 |  90 |
| mr_megahit.contig        |   95831 | 5139653 | 105 |
| mr_megahit.non-contained |   95831 | 5133034 |  92 |
| platanus.contig          |   16541 | 5191584 | 890 |
| platanus.scaffold        |   25324 | 5129772 | 407 |
| platanus.non-contained   |   25593 | 5097323 | 304 |


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
| trim.R      |     148 |  200.1M | 1452706 |
| Q20L60      |     148 | 193.66M | 1401466 |
| Q25L60      |     139 | 169.12M | 1304628 |
| Q30L60      |     119 | 125.02M | 1123194 |


Table: statTrimReads

| Name     | N50 |     Sum |       # |
|:---------|----:|--------:|--------:|
| clumpify | 251 | 447.53M | 1782994 |
| trim     | 148 |  200.1M | 1452706 |
| filter   | 148 |  200.1M | 1452706 |
| R1       | 164 | 100.23M |  655190 |
| R2       | 133 |  81.52M |  655190 |
| Rs       | 141 |  18.34M |  142326 |


```text
#R.trim
#Matched	113970	6.39206%
#Name	Reads	ReadsPct
Reverse_adapter	81598	4.57646%
pcr_dimer	14481	0.81217%
PCR_Primers	8081	0.45323%
TruSeq_Universal_Adapter	5665	0.31772%
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
#error_kmers	3605043
#genomic_kmers	4414304
#main_peak	30
#genome_size_in_peaks	4263517
#genome_size	5264941
#haploid_genome_size	5264941
#fold_coverage	30
#haploid_fold_coverage	30
#ploidy	1
#percent_repeat_in_peaks	0.376
#percent_repeat	12.772
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 148 | 200.09M | 1452579 |
| ecco          | 148 | 199.84M | 1452578 |
| ecct          | 148 | 198.72M | 1444000 |
| extended      | 186 | 255.79M | 1444000 |
| merged.raw    | 455 | 197.38M |  475527 |
| unmerged.raw  | 172 |  80.09M |  492946 |
| unmerged.trim | 172 |  80.07M |  492605 |
| M1            | 455 |  197.2M |  475127 |
| U1            | 172 |  19.67M |  121605 |
| U2            | 151 |  17.53M |  121605 |
| Us            | 182 |  42.86M |  249395 |
| M.cor         | 443 | 277.99M | 1692254 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 184.3 |    179 |  66.0 |         10.54% |
| M.ihist.merge.txt  | 415.1 |    452 |  89.0 |         65.86% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   |  43.5 |   38.7 |   11.10% | "39" |  4.6M | 4.55M |     0.99 | 0:00'37'' |
| Q20L60.R |  42.1 |   37.9 |    9.98% | "39" |  4.6M | 4.55M |     0.99 | 0:00'32'' |
| Q25L60.R |  36.8 |   34.9 |    5.03% | "35" |  4.6M | 4.54M |     0.99 | 0:00'29'' |
| Q30L60.R |  27.2 |   26.6 |    2.20% | "31" |  4.6M | 4.52M |     0.98 | 0:00'26'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  97.65% |     17835 | 4.34M | 400 |      2369 | 212.26K | 1303 |   27.0 | 6.0 |   5.0 |  90.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:40 |
| Q0L0XallP000   |   38.7 |  97.71% |     23038 | 4.36M | 339 |      2672 | 213.91K | 1214 |   35.0 | 7.0 |   5.0 | 112.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:41 |
| Q20L60X30P000  |   30.0 |  97.75% |     18402 | 4.34M | 402 |      2673 | 218.89K | 1250 |   27.0 | 6.0 |   5.0 |  90.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| Q20L60XallP000 |   37.9 |  97.78% |     24784 | 4.35M | 329 |      2877 | 201.89K | 1126 |   35.0 | 7.0 |   5.0 | 112.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q25L60X30P000  |   30.0 |  98.26% |     14880 | 4.32M | 497 |      4716 | 280.68K | 1406 |   28.0 | 6.0 |   5.0 |  92.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:31 |
| Q25L60XallP000 |   34.9 |  98.36% |     17779 | 4.34M | 426 |      5833 | 260.62K | 1284 |   32.0 | 7.0 |   5.0 | 106.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:33 |
| Q30L60XallP000 |   26.6 |  97.60% |      8381 | 4.23M | 753 |      2734 | 340.04K | 1845 |   25.0 | 6.0 |   5.0 |  86.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:29 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.39% |     20508 | 4.35M | 364 |      4332 |  140.5K | 739 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRX30P001  |   30.0 |  97.40% |     21605 | 4.34M | 357 |      3817 | 154.53K | 743 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRXallP000 |   60.4 |  97.36% |     22279 | 4.34M | 345 |      5800 | 167.85K | 749 |   55.0 | 9.0 |   9.3 | 164.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:29 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  97.17% |     11350 | 4.31M |  585 |      2815 | 256.99K | 1916 |   28.0 | 6.0 |   5.0 |  92.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:30 |
| Q0L0XallP000   |   38.7 |  97.53% |     15324 | 4.34M |  476 |      2795 | 233.05K | 1748 |   36.0 | 7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:33 |
| Q20L60X30P000  |   30.0 |  97.05% |     10392 | 4.29M |  627 |      2798 | 289.47K | 1937 |   28.0 | 6.0 |   5.0 |  92.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:30 |
| Q20L60XallP000 |   37.9 |  97.55% |     14639 | 4.33M |  495 |      3417 | 266.97K | 1666 |   35.0 | 7.0 |   5.0 | 112.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:31 |
| Q25L60X30P000  |   30.0 |  96.95% |      8398 | 4.25M |  772 |      4300 | 283.37K | 2137 |   28.0 | 7.0 |   5.0 |  98.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:30 |
| Q25L60XallP000 |   34.9 |  97.60% |      9972 |  4.3M |  661 |      5747 | 334.13K | 1935 |   33.0 | 7.0 |   5.0 | 108.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:30 |
| Q30L60XallP000 |   26.6 |  93.57% |      4494 | 4.01M | 1150 |      2438 | 285.52K | 2746 |   25.0 | 6.0 |   5.0 |  86.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:27 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.41% |     18774 | 4.35M | 404 |      5004 | 173.16K | 899 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:30 |
| MRX30P001  |   30.0 |  97.35% |     19071 | 4.34M | 389 |      6101 | 196.85K | 857 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:28 |
| MRXallP000 |   60.4 |  97.44% |     20747 | 4.33M | 353 |      6101 |  199.4K | 792 |   55.0 | 8.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:30 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X30P000    |   30.0 |  98.24% |     16897 | 4.35M | 447 |      3996 | 248.69K | 1499 |   28.0 | 6.0 |   5.0 |  92.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:31 |
| Q0L0XallP000   |   38.7 |  98.43% |     22283 | 4.36M | 350 |      4356 | 254.58K | 1328 |   36.0 | 7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:33 |
| Q20L60X30P000  |   30.0 |  98.14% |     15870 | 4.32M | 468 |      4349 | 292.05K | 1521 |   28.0 | 6.0 |   5.0 |  92.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:31 |
| Q20L60XallP000 |   37.9 |  98.38% |     19964 | 4.35M | 378 |      4490 | 282.47K | 1335 |   35.0 | 7.0 |   5.0 | 112.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:32 |
| Q25L60X30P000  |   30.0 |  97.83% |     11151 | 4.32M | 613 |      2792 |  248.1K | 1739 |   28.0 | 7.0 |   5.0 |  98.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:31 |
| Q25L60XallP000 |   34.9 |  98.14% |     13362 | 4.33M | 528 |      5824 | 305.44K | 1580 |   33.0 | 7.0 |   5.0 | 108.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:31 |
| Q30L60XallP000 |   26.6 |  96.22% |      6627 | 4.15M | 908 |      2431 | 317.11K | 2305 |   25.0 | 6.0 |   5.0 |  86.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:29 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.50% |     19948 | 4.34M | 379 |      5224 | 161.26K | 801 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:28 |
| MRX30P001  |   30.0 |  97.55% |     19862 | 4.34M | 372 |      5475 | 173.65K | 799 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:29 |
| MRXallP000 |   60.4 |  97.54% |     21368 | 4.33M | 344 |      6100 | 177.53K | 759 |   55.0 | 8.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:29 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  91.17% |     36740 | 4.41M | 246 |      4689 | 360.25K | 133 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:33 |
| 7_merge_mr_unitigs_bcalm      |  92.62% |     21353 | 4.33M | 349 |      5224 | 203.81K |  65 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:35 |
| 7_merge_mr_unitigs_superreads |  92.48% |     22986 | 4.39M | 328 |      5487 | 188.38K |  62 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:34 |
| 7_merge_mr_unitigs_tadpole    |  92.68% |     22194 | 4.33M | 335 |      5819 | 199.44K |  57 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:35 |
| 7_merge_unitigs_bcalm         |  90.95% |     17117 | 4.33M | 424 |      5224 | 311.57K | 102 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:30 |
| 7_merge_unitigs_superreads    |  91.58% |     30174 | 4.35M | 271 |      4689 | 308.94K | 114 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:31 |
| 7_merge_unitigs_tadpole       |  91.26% |     22998 | 4.34M | 319 |      5018 | 335.15K | 115 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:30 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 8_spades     |  99.20% |     45176 | 1.94M |  83 |      7794 |  86.43K | 125 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:29 |
| 8_mr_spades  |  99.13% |     30528 | 3.93M | 216 |     16123 |  151.9K | 304 |   55.0 | 8.0 |  10.3 | 158.0 |   0:00:31 |
| 8_megahit    |  98.54% |     37964 |  4.1M | 230 |      5525 | 153.99K | 387 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:28 |
| 8_mr_megahit |  99.34% |     26589 | 4.38M | 287 |     16123 | 164.89K | 523 |   55.0 | 8.0 |  10.3 | 158.0 |   0:00:31 |
| 8_platanus   |  96.34% |     28637 |  2.7M | 170 |      4785 | 157.09K | 241 |   35.0 | 7.0 |   5.0 | 112.0 |   0:00:29 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 3188524 | 4602977 |    7 |
| Paralogs                 |    2337 |  146789 |   66 |
| Repetitives              |     572 |   57281 |  165 |
| 7_merge_anchors.anchors  |   36740 | 4407623 |  246 |
| 7_merge_anchors.others   |    4689 |  360250 |  133 |
| glue_anchors             |   37591 | 4404044 |  240 |
| fill_anchors             |   48535 | 4406444 |  183 |
| spades.contig            |  150729 | 4576779 |  136 |
| spades.scaffold          |  172916 | 4577123 |  131 |
| spades.non-contained     |  150729 | 4562257 |   71 |
| mr_spades.contig         |   55603 | 4566224 |  170 |
| mr_spades.scaffold       |   89512 | 4567395 |  121 |
| mr_spades.non-contained  |   55603 | 4555133 |  149 |
| megahit.contig           |   52830 | 4572904 |  245 |
| megahit.non-contained    |   52830 | 4541309 |  182 |
| mr_megahit.contig        |   31157 | 4576803 |  282 |
| mr_megahit.non-contained |   31157 | 4563775 |  255 |
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
| trim.R      |     189 | 258.49M | 1426392 |
| Q20L60      |     189 | 254.11M | 1398505 |
| Q25L60      |     187 |  240.9M | 1348208 |
| Q30L60      |     181 |    214M | 1250494 |


Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 397.98M | 1585566 |
| filteredbytile | 251 | 374.64M | 1492590 |
| highpass       | 251 | 370.39M | 1475660 |
| trim           | 189 |  260.1M | 1433970 |
| filter         | 189 | 258.49M | 1426392 |
| R1             | 193 | 132.85M |  713196 |
| R2             | 185 | 125.64M |  713196 |
| Rs             |   0 |       0 |       0 |


```text
#R.trim
#Matched	1204585	81.63025%
#Name	Reads	ReadsPct
Reverse_adapter	583287	39.52719%
pcr_dimer	337494	22.87072%
PCR_Primers	173218	11.73834%
TruSeq_Adapter_Index_1_6	44566	3.02007%
TruSeq_Universal_Adapter	44294	3.00164%
Nextera_LMP_Read2_External_Adapter	18244	1.23633%
```

```text
#R.filter
#Matched	7576	0.52832%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	7576	0.52832%
```

```text
#R.peaks
#k	31
#unique_kmers	9654640
#error_kmers	5899734
#genomic_kmers	3754906
#main_peak	52
#genome_size_in_peaks	3822682
#genome_size	4125197
#haploid_genome_size	4125197
#fold_coverage	52
#haploid_fold_coverage	52
#ploidy	1
#percent_repeat_in_peaks	1.924
#percent_repeat	4.220
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 189 | 252.85M | 1391904 |
| ecco          | 189 | 252.82M | 1391904 |
| eccc          | 189 | 252.82M | 1391904 |
| ecct          | 189 | 250.78M | 1380570 |
| extended      | 228 | 305.77M | 1380570 |
| merged.raw    | 239 | 164.16M |  685565 |
| unmerged.raw  | 227 |   1.97M |    9440 |
| unmerged.trim | 227 |   1.97M |    9438 |
| M1            | 239 | 152.62M |  637789 |
| U1            | 240 |   1.07M |    4719 |
| U2            | 214 | 901.01K |    4719 |
| Us            |   0 |       0 |       0 |
| M.cor         | 238 | 155.23M | 1285016 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 197.1 |    191 |  44.4 |         95.41% |
| M.ihist.merge.txt  | 239.5 |    232 |  51.5 |         99.32% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   |  64.1 |   54.1 |   15.64% | "109" | 4.03M | 3.94M |     0.98 | 0:00'42'' |
| Q20L60.R |  63.0 |   53.7 |   14.70% | "111" | 4.03M | 3.92M |     0.97 | 0:00'36'' |
| Q25L60.R |  59.7 |   52.4 |   12.25% | "109" | 4.03M | 3.92M |     0.97 | 0:00'36'' |
| Q30L60.R |  53.1 |   48.0 |    9.55% | "103" | 4.03M | 3.92M |     0.97 | 0:00'35'' |


Table: statUnitigsSuperreads.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  93.59% |      8675 | 3.77M | 638 |       426 |  80.3K | 1358 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:34 |
| Q0L0X50P000    |   50.0 |  92.86% |      7819 | 3.77M | 684 |       102 | 68.73K | 1421 |   49.0 | 11.0 |   5.3 | 164.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:35 |
| Q0L0XallP000   |   54.1 |  92.78% |      7789 | 3.78M | 695 |        77 |  60.9K | 1442 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:29 |
| Q20L60X40P000  |   40.0 |  96.07% |     22584 | 3.81M | 317 |       528 | 60.36K |  748 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:28 |
| Q20L60X50P000  |   50.0 |  95.42% |     18102 | 3.82M | 363 |       146 | 46.18K |  806 |   50.0 | 11.0 |   5.7 | 166.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:29 |
| Q20L60XallP000 |   53.7 |  95.23% |     16576 | 3.81M | 391 |       575 | 60.89K |  838 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:29 |
| Q25L60X40P000  |   40.0 |  96.59% |     23667 | 3.74M | 300 |       609 | 54.38K |  718 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:28 |
| Q25L60X50P000  |   50.0 |  96.13% |     19245 | 3.82M | 335 |       157 | 43.29K |  755 |   50.0 | 11.0 |   5.7 | 166.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:29 |
| Q25L60XallP000 |   52.4 |  96.04% |     18128 | 3.81M | 350 |       720 | 58.56K |  775 |   52.0 | 11.0 |   6.3 | 170.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:29 |
| Q30L60X40P000  |   40.0 |  96.80% |     27778 | 3.81M | 275 |       967 | 55.74K |  700 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:29 |
| Q30L60XallP000 |   48.0 |  96.44% |     21003 | 3.82M | 306 |       309 | 44.93K |  711 |   48.0 | 10.0 |   6.0 | 156.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:27 |


Table: statMRUnitigsSuperreads.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   38.5 |  96.68% |     41046 | 3.83M | 188 |       909 | 40.22K | 401 |   39.0 | 7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:28 |


Table: statUnitigsBcalm.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  95.87% |     13169 | 3.81M | 453 |        54 | 50.33K | 1307 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:29 |
| Q0L0X50P000    |   50.0 |  94.70% |      9759 | 3.81M | 577 |        58 | 57.33K | 1429 |   49.0 | 11.0 |   5.3 | 164.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:30 |
| Q0L0XallP000   |   54.1 |  94.36% |      9626 | 3.81M | 608 |        51 | 54.45K | 1486 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:30 |
| Q20L60X40P000  |   40.0 |  97.43% |     31325 | 3.82M | 250 |        90 | 46.21K |  846 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:32 |
| Q20L60X50P000  |   50.0 |  97.02% |     30501 | 3.83M | 251 |       137 | 52.11K |  747 |   50.0 | 10.0 |   6.7 | 160.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:32 |
| Q20L60XallP000 |   53.7 |  96.93% |     32560 | 3.83M | 246 |       109 | 43.88K |  709 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:32 |
| Q25L60X40P000  |   40.0 |  97.66% |     31195 | 3.83M | 240 |        67 | 38.69K |  849 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:32 |
| Q25L60X50P000  |   50.0 |  97.48% |     35825 | 3.84M | 225 |        94 | 38.01K |  713 |   50.0 | 11.0 |   5.7 | 166.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:32 |
| Q25L60XallP000 |   52.4 |  97.42% |     35823 | 3.83M | 229 |       182 | 48.11K |  694 |   53.0 | 11.0 |   6.7 | 172.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:33 |
| Q30L60X40P000  |   40.0 |  97.65% |     34756 | 3.82M | 228 |       101 | 44.59K |  847 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:33 |
| Q30L60XallP000 |   48.0 |  97.68% |     39336 | 3.83M | 222 |       212 | 52.23K |  750 |   49.0 | 10.0 |   6.3 | 158.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:33 |


Table: statMRUnitigsBcalm.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   38.5 |  96.97% |     49347 | 3.84M | 170 |       457 | 38.36K | 401 |   39.0 | 7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:27 |


Table: statUnitigsTadpole.md

| Name           | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:---------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000    |   40.0 |  96.60% |     18508 | 3.83M | 361 |       484 | 60.56K |  917 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:29 |
| Q0L0X50P000    |   50.0 |  95.71% |     13576 | 3.84M | 466 |        79 | 49.42K | 1048 |   50.0 | 11.0 |   5.7 | 166.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:28 |
| Q0L0XallP000   |   54.1 |  95.54% |     13298 | 3.84M | 483 |        71 | 45.68K | 1079 |   54.0 | 12.0 |   6.0 | 180.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:30 |
| Q20L60X40P000  |   40.0 |  97.69% |     37279 | 3.83M | 221 |       903 | 57.65K |  658 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:33 |
| Q20L60X50P000  |   50.0 |  97.34% |     39643 | 3.85M | 215 |       189 | 37.86K |  565 |   50.0 | 11.0 |   5.7 | 166.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:33 |
| Q20L60XallP000 |   53.7 |  97.32% |     35197 | 3.84M | 229 |       603 | 49.46K |  578 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:31 |
| Q25L60X40P000  |   40.0 |  97.58% |     35238 | 3.84M | 219 |       521 | 47.92K |  646 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:31 |
| Q25L60X50P000  |   50.0 |  97.45% |     41077 | 3.85M | 217 |       136 | 36.11K |  572 |   50.0 | 11.0 |   5.7 | 166.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:31 |
| Q25L60XallP000 |   52.4 |  97.42% |     40557 | 3.76M | 219 |       753 | 49.17K |  569 |   53.0 | 11.0 |   6.7 | 172.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:36 |
| Q30L60X40P000  |   40.0 |  97.57% |     40681 | 3.84M | 207 |       278 | 42.96K |  667 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:33 |
| Q30L60XallP000 |   48.0 |  97.61% |     40680 | 3.84M | 203 |       266 | 40.67K |  588 |   48.0 | 10.0 |   6.0 | 156.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:31 |


Table: statMRUnitigsTadpole.md

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRXallP000 |   38.5 |  96.62% |     59534 | 3.84M | 163 |       905 | 34.81K | 330 |   39.0 | 7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:26 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.21% |     47140 | 3.83M | 167 |      1118 | 102.97K | 92 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:34 |
| 7_merge_mr_unitigs_bcalm      |  95.58% |     42378 | 3.82M | 183 |      1076 |  25.95K | 23 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:29 |
| 7_merge_mr_unitigs_superreads |  95.41% |     35955 | 3.82M | 204 |      1046 |  29.27K | 27 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:28 |
| 7_merge_mr_unitigs_tadpole    |  95.66% |     44878 | 3.82M | 178 |      1070 |   26.2K | 24 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:29 |
| 7_merge_unitigs_bcalm         |  97.72% |     42570 | 3.83M | 189 |      1108 |  55.99K | 50 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:38 |
| 7_merge_unitigs_superreads    |  97.45% |     37274 | 3.82M | 220 |      1118 |   80.2K | 73 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:37 |
| 7_merge_unitigs_tadpole       |  97.74% |     44902 |  3.8M | 180 |      1118 |  59.45K | 55 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:36 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.36% |     90125 | 3.25M | 120 |      1082 | 35.82K | 187 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:30 |
| 8_mr_spades  |  98.43% |     90099 | 3.87M | 139 |      1051 | 38.18K | 268 |   39.0 |  7.0 |   6.0 | 120.0 |   0:00:35 |
| 8_megahit    |  97.91% |     47116 | 3.73M | 162 |      1034 | 45.26K | 254 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:28 |
| 8_mr_megahit |  99.06% |     90226 | 3.89M | 123 |      1075 | 38.51K | 232 |   39.0 |  7.0 |   6.0 | 120.0 |   0:00:28 |
| 8_platanus   |  96.53% |     47050 | 3.84M | 177 |      1034 | 36.85K | 347 |   54.0 | 11.0 |   7.0 | 174.0 |   0:00:30 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 2961149 | 4033464 |   2 |
| Paralogs                 |    3424 |  119270 |  49 |
| Repetitives              |    1070 |  120471 | 244 |
| 7_merge_anchors.anchors  |   47140 | 3834259 | 167 |
| 7_merge_anchors.others   |    1118 |  102972 |  92 |
| glue_anchors             |   47606 | 3834197 | 165 |
| fill_anchors             |   82588 | 3836005 | 130 |
| spades.contig            |  124666 | 3945980 | 189 |
| spades.scaffold          |  124666 | 3946280 | 186 |
| spades.non-contained     |  124666 | 3917216 |  91 |
| mr_spades.contig         |   92221 | 3946942 | 228 |
| mr_spades.scaffold       |   94852 | 3947488 | 221 |
| mr_spades.non-contained  |   92221 | 3911904 | 129 |
| megahit.contig           |   82706 | 3941024 | 190 |
| megahit.non-contained    |   82706 | 3907965 | 123 |
| mr_megahit.contig        |   92352 | 3959784 | 166 |
| mr_megahit.non-contained |   92352 | 3931854 | 111 |
| platanus.contig          |   47235 | 3985037 | 565 |
| platanus.scaffold        |   55299 | 3930954 | 332 |
| platanus.non-contained   |   58949 | 3873293 | 175 |


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

