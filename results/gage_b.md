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
| merged.raw    | 236 | 835.96M | 3577133 |
| unmerged.raw  | 210 |  15.17M |   81816 |
| unmerged.trim | 210 |  15.17M |   81794 |
| M1            | 236 | 705.17M | 3020306 |
| U1            | 230 |   8.66M |   40897 |
| U2            | 186 |   6.51M |   40897 |
| Us            |   0 |       0 |       0 |
| M.cor         | 235 | 723.36M | 6122406 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 191.2 |    187 |  46.4 |         92.44% |
| M.ihist.merge.txt  | 233.7 |    227 |  51.7 |         98.87% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 245.2 |  197.7 |   19.34% | "45" | 5.09M | 5.57M |     1.09 | 0:02'16'' |
| Q20L60.R | 237.2 |  194.6 |   18.00% | "47" | 5.09M |  5.5M |     1.08 | 0:02'13'' |
| Q25L60.R | 215.5 |  182.3 |   15.41% | "45" | 5.09M | 5.23M |     1.03 | 0:02'03'' |
| Q30L60.R | 176.2 |  154.8 |   12.14% | "39" | 5.09M | 5.18M |     1.02 | 0:01'43'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  75.86% |      2352 | 4.07M | 1891 |       550 | 204.91K | 4108 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:36 |
| Q0L0X40P001   |   40.0 |  75.32% |      2366 | 4.03M | 1900 |       922 | 223.27K | 4183 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:36 |
| Q0L0X40P002   |   40.0 |  76.40% |      2287 | 4.08M | 1937 |      1001 | 225.92K | 4212 |   36.0 | 10.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:28 |
| Q0L0X60P000   |   60.0 |  60.05% |      1883 |  3.3M | 1847 |      1005 | 165.09K | 3898 |   53.0 | 13.0 |   5.0 | 184.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q0L0X60P001   |   60.0 |  61.18% |      1819 | 3.36M | 1902 |      1002 | 158.78K | 4000 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q0L0X60P002   |   60.0 |  62.31% |      1865 | 3.44M | 1914 |      1001 |  149.4K | 4057 |   53.0 | 13.0 |   5.0 | 184.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q0L0X80P000   |   80.0 |  47.76% |      1589 | 2.66M | 1673 |      1007 | 149.44K | 3510 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:28 |
| Q0L0X80P001   |   80.0 |  48.69% |      1618 | 2.72M | 1702 |      1005 | 130.45K | 3564 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:27 |
| Q20L60X40P000 |   40.0 |  77.78% |      2459 | 4.15M | 1888 |       388 | 204.01K | 4160 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:29 |
| Q20L60X40P001 |   40.0 |  77.89% |      2470 | 4.16M | 1888 |       833 | 220.48K | 4198 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:29 |
| Q20L60X40P002 |   40.0 |  77.78% |      2443 | 4.17M | 1908 |       607 | 208.09K | 4205 |   36.0 |  9.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:28 |
| Q20L60X60P000 |   60.0 |  64.52% |      1941 | 3.54M | 1936 |      1006 | 158.01K | 4134 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:27 |
| Q20L60X60P001 |   60.0 |  64.14% |      1890 | 3.52M | 1941 |      1003 | 161.35K | 4146 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q20L60X60P002 |   60.0 |  63.11% |      1869 | 3.49M | 1941 |      1002 | 144.64K | 4044 |   53.0 | 14.0 |   5.0 | 190.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q20L60X80P000 |   80.0 |  52.51% |      1669 | 2.91M | 1782 |      1005 | 144.61K | 3789 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:27 |
| Q20L60X80P001 |   80.0 |  51.61% |      1647 | 2.86M | 1770 |      1005 | 148.09K | 3691 |   70.0 | 18.0 |   5.3 | 248.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  95.68% |      6437 |    5M | 1048 |        46 |  90.78K | 2607 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:33 |
| Q25L60X40P001 |   40.0 |  95.41% |      6360 | 4.99M | 1068 |        46 |  90.73K | 2638 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:31 |
| Q25L60X40P002 |   40.0 |  95.22% |      6015 | 4.99M | 1123 |        45 |  93.27K | 2725 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:32 |
| Q25L60X60P000 |   60.0 |  92.65% |      4636 |  4.9M | 1339 |        34 |  79.89K | 3022 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:36 |
| Q25L60X60P001 |   60.0 |  92.17% |      4518 | 4.89M | 1388 |        34 |  85.18K | 3147 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:34 |
| Q25L60X60P002 |   60.0 |  92.90% |      4533 | 4.92M | 1396 |        40 |  92.24K | 3136 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:35 |
| Q25L60X80P000 |   80.0 |  89.67% |      3689 | 4.77M | 1551 |        48 | 110.81K | 3453 |   75.0 | 17.0 |   8.0 | 252.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:35 |
| Q25L60X80P001 |   80.0 |  89.23% |      3695 | 4.75M | 1616 |        40 | 109.35K | 3574 |   75.0 | 17.0 |   8.0 | 252.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:33 |
| Q30L60X40P000 |   40.0 |  97.61% |      8821 | 5.06M |  837 |        40 |  71.78K | 2384 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:34 |
| Q30L60X40P001 |   40.0 |  97.62% |      8425 | 5.07M |  829 |        43 |  78.02K | 2456 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:35 |
| Q30L60X40P002 |   40.0 |  97.89% |      9129 | 5.07M |  814 |        39 |  67.87K | 2368 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:35 |
| Q30L60X60P000 |   60.0 |  96.09% |      6395 | 5.03M | 1063 |        35 |  70.47K | 2701 |   57.0 | 13.0 |   6.0 | 192.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:36 |
| Q30L60X60P001 |   60.0 |  96.31% |      6805 | 5.05M | 1028 |        28 |  60.29K | 2639 |   57.0 | 13.0 |   6.0 | 192.0 | "31,41,51,61,71,81" |   0:01:03 |   0:00:36 |
| Q30L60X80P000 |   80.0 |  94.37% |      5183 | 4.95M | 1260 |        40 |  94.19K | 3035 |   76.0 | 17.0 |   8.3 | 254.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:36 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.30% |     23842 |  5.1M | 333 |        54 | 35.34K |  796 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:32 |
| MRX40P001 |   40.0 |  98.24% |     23707 | 5.11M | 357 |        46 | 32.15K |  863 |   39.0 |  7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:33 |
| MRX40P002 |   40.0 |  98.06% |     22558 | 5.12M | 383 |        51 | 36.71K |  884 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:33 |
| MRX60P000 |   60.0 |  97.18% |     13655 | 5.09M | 547 |        52 |  44.7K | 1213 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:34 |
| MRX60P001 |   60.0 |  96.92% |     12163 |  5.1M | 619 |        49 | 46.03K | 1356 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:36 |
| MRX80P000 |   80.0 |  95.57% |      9106 | 5.06M | 830 |        50 | 54.08K | 1790 |   77.0 | 12.0 |  13.7 | 226.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:36 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  96.92% |      8383 | 5.07M |  849 |        41 |  91.34K | 2868 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:37 |
| Q0L0X40P001   |   40.0 |  96.93% |      8402 | 5.07M |  845 |        38 |  85.61K | 2859 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:39 |
| Q0L0X40P002   |   40.0 |  96.96% |      8326 | 5.07M |  848 |        41 |  92.91K | 2934 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:37 |
| Q0L0X60P000   |   60.0 |  92.44% |      4781 | 4.91M | 1311 |        44 | 110.39K | 3210 |   56.0 | 12.0 |   6.7 | 184.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:36 |
| Q0L0X60P001   |   60.0 |  93.16% |      4987 | 4.94M | 1305 |        31 |  82.98K | 3232 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:36 |
| Q0L0X60P002   |   60.0 |  92.94% |      4971 | 4.92M | 1286 |        40 | 103.82K | 3202 |   56.0 | 12.0 |   6.7 | 184.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:37 |
| Q0L0X80P000   |   80.0 |  84.76% |      3098 | 4.57M | 1717 |        50 | 140.57K | 3772 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:43 |   0:00:34 |
| Q0L0X80P001   |   80.0 |  86.02% |      3242 | 4.61M | 1679 |        49 | 138.48K | 3775 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:42 |   0:00:34 |
| Q20L60X40P000 |   40.0 |  97.03% |      8772 | 5.06M |  818 |        46 |  96.57K | 2795 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:38 |
| Q20L60X40P001 |   40.0 |  97.25% |      9139 | 5.08M |  795 |        40 |  87.41K | 2827 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:37 |
| Q20L60X40P002 |   40.0 |  96.88% |      8702 | 5.07M |  872 |        44 | 101.79K | 2995 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:39 |
| Q20L60X60P000 |   60.0 |  93.42% |      5180 | 4.95M | 1251 |        39 |  93.98K | 3154 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:35 |
| Q20L60X60P001 |   60.0 |  93.14% |      5109 | 4.96M | 1259 |        33 |  83.34K | 3174 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:37 |
| Q20L60X60P002 |   60.0 |  93.44% |      4907 | 4.93M | 1286 |        44 | 113.28K | 3261 |   56.0 | 12.0 |   6.7 | 184.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:36 |
| Q20L60X80P000 |   80.0 |  86.96% |      3394 | 4.64M | 1647 |        58 | 133.43K | 3742 |   74.0 | 17.0 |   7.7 | 250.0 | "31,41,51,61,71,81" |   0:01:46 |   0:00:35 |
| Q20L60X80P001 |   80.0 |  86.26% |      3223 | 4.63M | 1718 |        41 | 128.23K | 3823 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:44 |   0:00:34 |
| Q25L60X40P000 |   40.0 |  98.74% |     14085 | 5.15M |  544 |        35 |  58.26K | 2188 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:39 |
| Q25L60X40P001 |   40.0 |  98.73% |     14024 | 5.13M |  547 |        36 |  60.59K | 2192 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:38 |
| Q25L60X40P002 |   40.0 |  98.58% |     13700 | 5.08M |  579 |        35 |  58.86K | 2181 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:37 |
| Q25L60X60P000 |   60.0 |  97.89% |      9511 | 5.12M |  739 |        38 |  70.29K | 2296 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:38 |
| Q25L60X60P001 |   60.0 |  97.67% |      9622 | 5.07M |  786 |        40 |  78.59K | 2370 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:36 |
| Q25L60X60P002 |   60.0 |  98.06% |     10216 | 5.07M |  738 |        36 |     72K | 2354 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:40 |
| Q25L60X80P000 |   80.0 |  96.21% |      7021 | 5.03M |  989 |        33 |  71.53K | 2584 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:01:48 |   0:00:37 |
| Q25L60X80P001 |   80.0 |  96.17% |      6989 | 5.04M | 1018 |        42 |  93.78K | 2687 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:01:51 |   0:00:37 |
| Q30L60X40P000 |   40.0 |  98.79% |     14088 | 5.09M |  525 |        33 |   53.7K | 2097 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:38 |
| Q30L60X40P001 |   40.0 |  98.95% |     15511 | 5.14M |  519 |        31 |  48.48K | 2077 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:38 |
| Q30L60X40P002 |   40.0 |  98.90% |     15105 |  5.1M |  527 |        29 |  45.33K | 2061 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:37 |
| Q30L60X60P000 |   60.0 |  98.59% |     12069 |  5.1M |  621 |        30 |  54.07K | 2337 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:44 |
| Q30L60X60P001 |   60.0 |  98.73% |     13247 | 5.11M |  611 |        28 |  49.78K | 2249 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:42 |
| Q30L60X80P000 |   80.0 |  97.82% |      9145 | 5.07M |  799 |        31 |  63.94K | 2486 |   77.0 | 16.0 |   9.7 | 250.0 | "31,41,51,61,71,81" |   0:01:46 |   0:00:39 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  99.13% |     81424 | 5.11M | 118 |        55 | 16.57K | 373 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:34 |
| MRX40P001 |   40.0 |  99.22% |     90130 | 5.13M | 111 |        51 |  15.5K | 383 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:36 |
| MRX40P002 |   40.0 |  99.21% |     82486 | 5.11M | 120 |        57 | 17.17K | 411 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:39 |
| MRX60P000 |   60.0 |  98.97% |     57844 | 5.11M | 148 |        64 | 16.71K | 404 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:37 |
| MRX60P001 |   60.0 |  98.93% |     60744 | 5.13M | 157 |        54 | 15.03K | 432 |   58.0 |  9.0 |  10.3 | 170.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:37 |
| MRX80P000 |   80.0 |  98.74% |     45415 | 5.11M | 199 |        63 | 15.76K | 488 |   78.0 | 11.0 |  15.0 | 222.0 | "31,41,51,61,71,81" |   0:01:42 |   0:00:37 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.08% |     16289 | 5.13M |  510 |        46 |  63.02K | 1814 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:38 |
| Q0L0X40P001   |   40.0 |  97.93% |     14675 |  5.1M |  528 |        39 |  56.14K | 1828 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:37 |
| Q0L0X40P002   |   40.0 |  98.00% |     14522 |  5.1M |  518 |        42 |  57.51K | 1814 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:36 |
| Q0L0X60P000   |   60.0 |  96.06% |      8665 | 5.07M |  817 |        40 |  68.69K | 2137 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:35 |
| Q0L0X60P001   |   60.0 |  96.54% |      9061 | 5.08M |  799 |        37 |  64.64K | 2119 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  96.61% |      9086 | 5.07M |  806 |        36 |  62.77K | 2190 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:37 |
| Q0L0X80P000   |   80.0 |  94.85% |      6091 | 5.01M | 1110 |        41 | 104.67K | 2998 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:39 |
| Q0L0X80P001   |   80.0 |  95.24% |      6624 | 5.01M | 1059 |        43 | 104.62K | 2895 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:38 |
| Q20L60X40P000 |   40.0 |  98.13% |     15457 |  5.1M |  504 |        43 |  56.58K | 1747 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |
| Q20L60X40P001 |   40.0 |  98.16% |     15672 | 5.09M |  481 |        43 |  59.06K | 1730 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:36 |
| Q20L60X40P002 |   40.0 |  98.03% |     15372 | 5.12M |  540 |        49 |  69.67K | 1891 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:37 |
| Q20L60X60P000 |   60.0 |  96.63% |      9124 | 5.08M |  819 |        40 |  70.76K | 2183 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:34 |
| Q20L60X60P001 |   60.0 |  96.59% |      9437 | 5.08M |  803 |        35 |   63.5K | 2194 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:34 |
| Q20L60X60P002 |   60.0 |  96.57% |      8292 | 5.07M |  847 |        40 |  72.33K | 2245 |   57.0 | 12.0 |   7.0 | 186.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:35 |
| Q20L60X80P000 |   80.0 |  95.54% |      6376 | 5.02M | 1072 |        45 |  109.3K | 2900 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:38 |
| Q20L60X80P001 |   80.0 |  95.47% |      6091 | 5.02M | 1101 |        39 |  99.48K | 2998 |   76.0 | 16.0 |   9.3 | 248.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  99.06% |     23343 | 5.14M |  355 |        39 |  40.44K | 1376 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:39 |
| Q25L60X40P001 |   40.0 |  99.13% |     21635 | 5.14M |  354 |        40 |  43.23K | 1442 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:39 |
| Q25L60X40P002 |   40.0 |  99.01% |     23574 |  5.1M |  367 |        44 |  46.82K | 1439 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:37 |
| Q25L60X60P000 |   60.0 |  98.57% |     15856 | 5.12M |  490 |        42 |  51.83K | 1562 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:39 |
| Q25L60X60P001 |   60.0 |  98.57% |     14924 | 5.11M |  515 |        40 |  52.66K | 1632 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:39 |
| Q25L60X60P002 |   60.0 |  98.76% |     16112 | 5.11M |  499 |        40 |  52.62K | 1575 |   58.0 | 12.0 |   7.3 | 188.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:41 |
| Q25L60X80P000 |   80.0 |  98.23% |     12940 | 5.11M |  591 |        36 |   55.1K | 1818 |   77.0 | 16.0 |   9.7 | 250.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:40 |
| Q25L60X80P001 |   80.0 |  98.31% |     12472 | 5.13M |  625 |        38 |  59.92K | 1899 |   77.0 | 16.0 |   9.7 | 250.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:41 |
| Q30L60X40P000 |   40.0 |  99.21% |     22654 | 5.12M |  351 |        35 |  39.17K | 1486 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:41 |
| Q30L60X40P001 |   40.0 |  99.20% |     22753 | 5.13M |  363 |        32 |  35.21K | 1429 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:38 |
| Q30L60X40P002 |   40.0 |  99.28% |     23120 | 5.12M |  365 |        34 |  37.99K | 1437 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:37 |
| Q30L60X60P000 |   60.0 |  99.09% |     20087 | 5.11M |  411 |        31 |  36.65K | 1573 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:42 |
| Q30L60X60P001 |   60.0 |  99.14% |     20567 | 5.11M |  402 |        30 |  37.15K | 1546 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:41 |
| Q30L60X80P000 |   80.0 |  98.93% |     16369 |  5.1M |  486 |        31 |  45.75K | 1797 |   78.0 | 16.0 |  10.0 | 252.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:42 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  99.16% |    116480 | 5.13M |  96 |        63 | 14.05K | 288 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:35 |
| MRX40P001 |   40.0 |  99.27% |    117645 | 5.18M |  95 |        61 | 14.84K | 317 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:39 |
| MRX40P002 |   40.0 |  99.14% |    102717 | 5.11M | 102 |        62 | 15.05K | 325 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| MRX60P000 |   60.0 |  99.15% |     98730 | 5.13M | 104 |        75 | 13.64K | 295 |   59.0 |  9.0 |  10.7 | 172.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:39 |
| MRX60P001 |   60.0 |  99.15% |     98853 | 5.11M | 102 |        60 | 10.74K | 293 |   59.0 |  9.0 |  10.7 | 172.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:38 |
| MRX80P000 |   80.0 |  99.13% |     94964 | 5.13M | 105 |        82 | 12.33K | 305 |   78.0 | 11.0 |  15.0 | 222.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:39 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  99.01% |     79692 | 5.11M | 119 |     14396 | 938.02K | 214 |  193.0 | 33.0 |  31.3 | 584.0 |   0:01:00 |
| 7_merge_mr_unitigs_bcalm      |  99.23% |     87966 | 5.11M | 117 |     16484 |  21.41K |   6 |  192.0 | 33.0 |  31.0 | 582.0 |   0:01:09 |
| 7_merge_mr_unitigs_superreads |  99.16% |     75129 | 5.11M | 125 |      1060 |  14.63K |  15 |  191.0 | 33.0 |  30.7 | 580.0 |   0:01:07 |
| 7_merge_mr_unitigs_tadpole    |  99.12% |     82774 | 4.77M | 111 |     44227 |  53.83K |  10 |  192.0 | 33.0 |  31.0 | 582.0 |   0:01:06 |
| 7_merge_unitigs_bcalm         |  99.38% |     92832 | 5.11M | 117 |      9290 | 402.18K | 101 |  188.0 | 36.0 |  26.7 | 592.0 |   0:01:36 |
| 7_merge_unitigs_superreads    |  99.33% |     68842 | 5.13M | 153 |      1941 | 299.95K | 167 |  189.0 | 35.0 |  28.0 | 588.0 |   0:01:25 |
| 7_merge_unitigs_tadpole       |  99.44% |     93581 | 5.11M | 111 |     27430 | 434.91K |  74 |  189.0 | 36.0 |  27.0 | 594.0 |   0:01:35 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  84.42% |      3582 | 4.52M | 1556 |       156 | 153.65K | 3131 |  185.0 | 36.0 |  25.7 | 586.0 |   0:00:39 |
| 8_mr_spades  |  99.20% |    127450 | 5.12M |   71 |        92 |   8.24K |  119 |  139.0 | 18.0 |  28.3 | 386.0 |   0:00:39 |
| 8_megahit    |  98.24% |     73343 | 3.17M |   90 |       116 |   9.23K |  146 |  193.0 | 33.0 |  31.3 | 584.0 |   0:00:42 |
| 8_mr_megahit |  99.39% |    142388 | 4.94M |   61 |       109 |   7.48K |  108 |  139.0 | 18.0 |  28.3 | 386.0 |   0:00:37 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 5067172 | 5090491 |    2 |
| Paralogs                 |    1693 |   83291 |   53 |
| Repetitives              |     192 |   15322 |   82 |
| 7_merge_anchors.anchors  |   79692 | 5107054 |  119 |
| 7_merge_anchors.others   |   14396 |  938021 |  214 |
| glue_anchors             |   82779 | 5106943 |  116 |
| fill_anchors             |  122243 | 5291097 |   70 |
| spades.contig            |    2938 | 5836762 | 6077 |
| spades.scaffold          |    2942 | 5837300 | 6071 |
| spades.non-contained     |    3684 | 4676177 | 1583 |
| mr_spades.contig         |  166845 | 5137989 |   77 |
| mr_spades.scaffold       |  220390 | 5138416 |   70 |
| mr_spades.non-contained  |  166845 | 5126753 |   53 |
| megahit.contig           |  149440 | 5139604 |  120 |
| megahit.non-contained    |  149440 | 5119438 |   68 |
| mr_megahit.contig        |  215357 | 5146723 |   70 |
| mr_megahit.non-contained |  215357 | 5137091 |   48 |


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
| merged.raw    | 460 |   2.04G |  4802337 |
| unmerged.raw  | 163 | 296.36M |  1923728 |
| unmerged.trim | 163 | 296.34M |  1923596 |
| M1            | 460 |   2.02G |  4770795 |
| U1            | 175 | 160.46M |   961798 |
| U2            | 148 | 135.88M |   961798 |
| Us            |   0 |       0 |        0 |
| M.cor         | 456 |   2.32G | 11465186 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 193.1 |    189 |  63.7 |         10.84% |
| M.ihist.merge.txt  | 423.8 |    457 |  84.5 |         83.31% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 355.6 |  321.1 |    9.70% | "39" |  4.6M | 4.83M |     1.05 | 0:03'05'' |
| Q20L60.R | 347.9 |  316.4 |    9.05% | "39" |  4.6M | 4.76M |     1.03 | 0:03'06'' |
| Q25L60.R | 308.2 |  293.6 |    4.73% | "35" |  4.6M | 4.58M |     0.99 | 0:02'51'' |
| Q30L60.R | 231.3 |  226.5 |    2.06% | "31" |  4.6M | 4.55M |     0.99 | 0:02'22'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  91.50% |      7371 |  4.3M |  873 |      1000 | 160.56K | 2812 |   35.0 |  8.0 |   5.0 | 118.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:47 |
| Q0L0X40P001   |   40.0 |  91.72% |      7338 |  4.3M |  817 |      1034 | 159.62K | 2620 |   36.0 |  8.0 |   5.0 | 120.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:47 |
| Q0L0X40P002   |   40.0 |  91.79% |      7053 |  4.3M |  870 |      1037 | 193.26K | 2755 |   35.0 |  8.0 |   5.0 | 118.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:34 |
| Q0L0X60P000   |   60.0 |  86.75% |      4990 | 4.22M | 1115 |      1003 | 127.22K | 3122 |   52.0 | 12.0 |   5.3 | 176.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:35 |
| Q0L0X60P001   |   60.0 |  86.36% |      5324 | 4.18M | 1078 |      1000 | 129.18K | 3074 |   52.0 | 12.0 |   5.3 | 176.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:35 |
| Q0L0X60P002   |   60.0 |  86.52% |      5143 |  4.2M | 1111 |       694 | 132.93K | 3118 |   52.0 | 11.0 |   6.3 | 170.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:37 |
| Q0L0X80P000   |   80.0 |  81.72% |      3802 | 4.06M | 1331 |       850 | 128.01K | 3361 |   69.0 | 15.0 |   8.0 | 228.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:34 |
| Q0L0X80P001   |   80.0 |  80.90% |      3839 | 4.05M | 1298 |       538 |  116.3K | 3304 |   69.0 | 15.0 |   8.0 | 228.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:33 |
| Q0L0X80P002   |   80.0 |  81.73% |      3741 | 4.06M | 1343 |       583 | 128.18K | 3400 |   69.0 | 15.0 |   8.0 | 228.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:32 |
| Q20L60X40P000 |   40.0 |  92.92% |      8049 |  4.3M |  802 |      1050 | 182.28K | 2644 |   36.0 |  8.0 |   5.0 | 120.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:35 |
| Q20L60X40P001 |   40.0 |  92.97% |      7994 | 4.32M |  801 |      1108 | 189.03K | 2598 |   36.0 |  8.0 |   5.0 | 120.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:33 |
| Q20L60X40P002 |   40.0 |  92.10% |      7832 | 4.29M |  791 |      1043 | 170.94K | 2514 |   36.0 |  8.0 |   5.0 | 120.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:31 |
| Q20L60X60P000 |   60.0 |  88.02% |      5819 | 4.25M | 1022 |       875 | 114.31K | 2897 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:34 |
| Q20L60X60P001 |   60.0 |  88.07% |      5667 | 4.25M | 1014 |      1003 | 133.37K | 2892 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:34 |
| Q20L60X60P002 |   60.0 |  88.40% |      6239 | 4.25M |  961 |      1027 | 119.35K | 2777 |   53.0 | 12.0 |   5.7 | 178.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:34 |
| Q20L60X80P000 |   80.0 |  83.56% |      4483 | 4.11M | 1203 |       746 |  133.4K | 3153 |   70.0 | 15.0 |   8.3 | 230.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:34 |
| Q20L60X80P001 |   80.0 |  83.91% |      4688 | 4.09M | 1150 |       853 | 134.39K | 3099 |   70.0 | 15.0 |   8.3 | 230.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:34 |
| Q20L60X80P002 |   80.0 |  83.99% |      4561 | 4.14M | 1201 |       399 | 123.44K | 3097 |   70.0 | 15.0 |   8.3 | 230.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:35 |
| Q25L60X40P000 |   40.0 |  97.38% |     13710 | 4.34M |  508 |      1623 | 214.24K | 1577 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  97.47% |     13545 | 4.33M |  514 |      2110 | 232.09K | 1561 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:32 |
| Q25L60X40P002 |   40.0 |  97.62% |     14201 | 4.34M |  512 |      2161 | 230.21K | 1578 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:33 |
| Q25L60X60P000 |   60.0 |  96.79% |     14681 | 4.34M |  500 |      1703 | 208.75K | 1631 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:35 |
| Q25L60X60P001 |   60.0 |  97.04% |     14322 | 4.34M |  499 |      1839 |  222.4K | 1644 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:36 |
| Q25L60X60P002 |   60.0 |  96.78% |     14553 | 4.35M |  494 |      1672 |  184.9K | 1577 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:37 |
| Q25L60X80P000 |   80.0 |  96.28% |     14335 | 4.36M |  534 |      1458 | 187.32K | 1714 |   73.0 | 16.0 |   8.3 | 242.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:40 |
| Q25L60X80P001 |   80.0 |  96.36% |     13514 | 4.35M |  522 |      1535 |    190K | 1641 |   73.0 | 16.0 |   8.3 | 242.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:38 |
| Q25L60X80P002 |   80.0 |  96.22% |     14669 | 4.36M |  510 |      1484 | 183.11K | 1692 |   73.0 | 16.0 |   8.3 | 242.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:39 |
| Q30L60X40P000 |   40.0 |  97.96% |     11209 | 4.31M |  588 |      2445 | 235.48K | 1608 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:31 |
| Q30L60X40P001 |   40.0 |  97.73% |     10949 |  4.3M |  607 |      2730 | 248.86K | 1662 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:31 |
| Q30L60X40P002 |   40.0 |  98.01% |     11377 | 4.31M |  618 |      3601 | 275.61K | 1655 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:31 |
| Q30L60X60P000 |   60.0 |  98.26% |     15316 | 4.36M |  475 |      4603 | 225.67K | 1398 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:32 |
| Q30L60X60P001 |   60.0 |  98.04% |     14819 | 4.35M |  492 |      4781 | 264.01K | 1386 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:36 |
| Q30L60X60P002 |   60.0 |  98.12% |     15314 | 4.36M |  482 |      3176 | 215.75K | 1423 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  98.18% |     17385 | 4.36M |  438 |      6108 | 264.34K | 1357 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:36 |
| Q30L60X80P001 |   80.0 |  98.17% |     16828 | 4.36M |  455 |      4697 | 258.74K | 1347 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:35 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.96% |     44621 | 4.27M | 216 |      4795 |  154.9K | 476 |   37.0 |  7.0 |   5.3 | 116.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:33 |
| MRX40P001 |   40.0 |  97.78% |     44876 | 4.38M | 210 |      2510 | 135.49K | 487 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:31 |
| MRX40P002 |   40.0 |  97.85% |     44734 | 4.37M | 202 |      4820 | 166.91K | 481 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:31 |
| MRX60P000 |   60.0 |  97.89% |     44626 | 4.38M | 230 |      4369 | 148.57K | 492 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:35 |
| MRX60P001 |   60.0 |  97.69% |     40502 | 4.16M | 210 |      4814 | 155.73K | 471 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:32 |
| MRX60P002 |   60.0 |  97.88% |     44626 | 4.26M | 210 |      3304 | 137.93K | 471 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:35 |
| MRX80P000 |   80.0 |  97.76% |     39782 | 4.37M | 245 |      2672 |  152.6K | 528 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:36 |
| MRX80P001 |   80.0 |  97.77% |     34852 | 4.32M | 233 |      4263 | 158.79K | 518 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:34 |
| MRX80P002 |   80.0 |  97.81% |     36031 | 4.28M | 237 |      4263 | 155.19K | 491 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:35 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.94% |     11329 | 4.34M | 611 |      1696 | 217.97K | 2183 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:42 |
| Q0L0X40P001   |   40.0 |  96.18% |     11207 | 4.33M | 610 |      1605 | 202.81K | 2168 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:42 |
| Q0L0X40P002   |   40.0 |  96.27% |     11316 | 4.34M | 619 |      1725 | 222.51K | 2240 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:36 |
| Q0L0X60P000   |   60.0 |  94.77% |     10164 | 4.36M | 680 |      1402 | 195.47K | 2669 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:41 |
| Q0L0X60P001   |   60.0 |  94.73% |     10811 | 4.34M | 645 |      1438 | 214.66K | 2536 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:42 |
| Q0L0X60P002   |   60.0 |  94.99% |      9915 | 4.35M | 667 |      1257 | 200.19K | 2638 |   54.0 | 11.0 |   7.0 | 174.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:38 |
| Q0L0X80P000   |   80.0 |  92.67% |      8385 | 4.34M | 790 |      1075 | 159.67K | 2956 |   72.0 | 15.0 |   9.0 | 234.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:42 |
| Q0L0X80P001   |   80.0 |  92.66% |      7976 | 4.34M | 795 |      1188 | 160.87K | 2902 |   71.0 | 15.0 |   8.7 | 232.0 | "31,41,51,61,71,81" |   0:01:31 |   0:00:41 |
| Q0L0X80P002   |   80.0 |  92.21% |      7753 | 4.34M | 833 |      1063 | 160.63K | 2968 |   71.0 | 15.0 |   8.7 | 232.0 | "31,41,51,61,71,81" |   0:01:43 |   0:00:41 |
| Q20L60X40P000 |   40.0 |  96.40% |     11386 | 4.32M | 602 |      1985 | 251.31K | 2086 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:34 |
| Q20L60X40P001 |   40.0 |  96.33% |     10651 | 4.35M | 650 |      1896 | 218.72K | 2235 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:34 |
| Q20L60X40P002 |   40.0 |  96.35% |     10043 | 4.31M | 642 |      1639 | 218.07K | 2137 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:33 |
| Q20L60X60P000 |   60.0 |  94.88% |     11141 | 4.34M | 620 |      1361 | 160.93K | 2370 |   54.0 | 12.0 |   6.0 | 180.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:37 |
| Q20L60X60P001 |   60.0 |  94.95% |     10735 | 4.35M | 646 |      1382 | 181.33K | 2433 |   54.0 | 12.0 |   6.0 | 180.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:37 |
| Q20L60X60P002 |   60.0 |  94.98% |     10925 | 4.36M | 647 |      1798 | 185.41K | 2392 |   54.0 | 12.0 |   6.0 | 180.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:43 |
| Q20L60X80P000 |   80.0 |  93.09% |      8926 | 4.32M | 752 |      1142 | 169.07K | 2776 |   72.0 | 15.0 |   9.0 | 234.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:39 |
| Q20L60X80P001 |   80.0 |  93.07% |      9165 | 4.33M | 728 |      1203 | 173.33K | 2660 |   72.0 | 15.0 |   9.0 | 234.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:40 |
| Q20L60X80P002 |   80.0 |  93.21% |      9318 | 4.32M | 733 |      1091 | 156.92K | 2614 |   72.0 | 15.0 |   9.0 | 234.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:39 |
| Q25L60X40P000 |   40.0 |  97.08% |      9518 | 4.29M | 684 |      2552 | 259.76K | 1932 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  97.29% |      9392 | 4.28M | 670 |      3619 | 268.04K | 1941 |   38.0 |  9.0 |   5.0 | 130.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:34 |
| Q25L60X40P002 |   40.0 |  97.01% |      9288 | 4.29M | 696 |      2467 | 238.75K | 1958 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:32 |
| Q25L60X60P000 |   60.0 |  97.48% |     13318 | 4.34M | 523 |      2195 | 223.33K | 1645 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:37 |
| Q25L60X60P001 |   60.0 |  97.73% |     12738 | 4.35M | 553 |      2997 | 237.59K | 1789 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:37 |
| Q25L60X60P002 |   60.0 |  97.55% |     12834 | 4.34M | 535 |      2906 | 224.92K | 1741 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:36 |
| Q25L60X80P000 |   80.0 |  97.31% |     15448 | 4.35M | 487 |      2108 | 212.79K | 1636 |   75.0 | 16.0 |   9.0 | 246.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:37 |
| Q25L60X80P001 |   80.0 |  97.51% |     14298 | 4.34M | 496 |      2170 | 235.61K | 1741 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:41 |
| Q25L60X80P002 |   80.0 |  97.36% |     14165 | 4.35M | 497 |      1912 | 224.04K | 1730 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:40 |
| Q30L60X40P000 |   40.0 |  96.39% |      6338 | 4.19M | 928 |      2938 | 269.42K | 2369 |   38.0 | 10.0 |   5.0 | 136.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:32 |
| Q30L60X40P001 |   40.0 |  96.33% |      6788 | 4.19M | 894 |      2963 | 273.54K | 2281 |   38.0 | 10.0 |   5.0 | 136.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:31 |
| Q30L60X40P002 |   40.0 |  96.19% |      6724 | 4.18M | 911 |      2731 | 249.79K | 2335 |   38.0 | 10.0 |   5.0 | 136.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:32 |
| Q30L60X60P000 |   60.0 |  97.60% |      8859 | 4.29M | 713 |      4602 | 255.54K | 1958 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:35 |
| Q30L60X60P001 |   60.0 |  97.49% |      8873 |  4.3M | 718 |      2991 | 230.88K | 1935 |   57.0 | 14.0 |   5.0 | 198.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:34 |
| Q30L60X60P002 |   60.0 |  97.28% |      9310 | 4.29M | 704 |      3120 | 244.15K | 1990 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:35 |
| Q30L60X80P000 |   80.0 |  97.64% |     10766 | 4.31M | 611 |      3612 | 243.99K | 1830 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:38 |
| Q30L60X80P001 |   80.0 |  97.91% |     11007 | 4.32M | 621 |      4583 | 243.72K | 1781 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:35 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.88% |     46626 | 4.37M | 213 |      2672 | 116.13K | 551 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:33 |
| MRX40P001 |   40.0 |  97.88% |     44334 | 4.37M | 212 |      5367 | 152.61K | 557 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:33 |
| MRX40P002 |   40.0 |  97.76% |     42710 | 4.36M | 213 |      4813 | 144.29K | 570 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:31 |
| MRX60P000 |   60.0 |  97.75% |     44626 | 4.38M | 221 |      4369 | 150.32K | 523 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:32 |
| MRX60P001 |   60.0 |  97.77% |     44622 | 4.36M | 206 |      4807 | 151.26K | 516 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:34 |
| MRX60P002 |   60.0 |  97.82% |     44631 | 4.38M | 213 |      3581 | 145.89K | 512 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:36 |
| MRX80P000 |   80.0 |  97.91% |     45188 | 4.37M | 225 |      4407 | 163.05K | 534 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:38 |   0:00:37 |
| MRX80P001 |   80.0 |  97.90% |     36403 | 4.35M | 218 |      5247 |    171K | 533 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:37 |
| MRX80P002 |   80.0 |  97.88% |     42108 | 4.38M | 220 |      4918 | 158.83K | 506 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:35 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.76% |     18347 | 4.36M | 430 |      4454 | 255.32K | 1441 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:32 |
| Q0L0X40P001   |   40.0 |  97.87% |     17245 | 4.35M | 427 |      4570 | 234.65K | 1430 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:33 |
| Q0L0X40P002   |   40.0 |  97.80% |     16471 | 4.35M | 444 |      2992 |  248.1K | 1474 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:34 |
| Q0L0X60P000   |   60.0 |  97.94% |     20566 | 4.38M | 393 |      2731 | 236.76K | 1508 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:43 |
| Q0L0X60P001   |   60.0 |  98.06% |     19953 | 4.37M | 394 |      3505 | 241.37K | 1442 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:40 |
| Q0L0X60P002   |   60.0 |  98.04% |     20306 | 4.37M | 389 |      2964 | 228.94K | 1517 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:37 |
| Q0L0X80P000   |   80.0 |  98.00% |     19132 | 4.38M | 412 |      3125 | 262.24K | 1672 |   74.0 | 15.0 |   9.7 | 238.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:41 |
| Q0L0X80P001   |   80.0 |  98.12% |     18205 | 4.37M | 426 |      2897 | 245.58K | 1724 |   73.0 | 15.0 |   9.3 | 236.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:41 |
| Q0L0X80P002   |   80.0 |  97.90% |     17632 | 4.37M | 432 |      2995 | 271.86K | 1726 |   73.0 | 15.0 |   9.3 | 236.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:41 |
| Q20L60X40P000 |   40.0 |  97.94% |     17924 | 4.35M | 448 |      3326 | 239.42K | 1469 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:34 |
| Q20L60X40P001 |   40.0 |  97.83% |     16462 | 4.35M | 459 |      2948 |  214.1K | 1483 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:34 |
| Q20L60X40P002 |   40.0 |  97.91% |     15882 | 4.34M | 459 |      5214 | 277.83K | 1470 |   37.0 |  8.0 |   5.0 | 122.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:34 |
| Q20L60X60P000 |   60.0 |  98.08% |     20241 | 4.36M | 401 |      3414 | 230.43K | 1462 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:39 |
| Q20L60X60P001 |   60.0 |  98.16% |     18859 | 4.37M | 409 |      3086 | 272.41K | 1521 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:39 |
| Q20L60X60P002 |   60.0 |  98.11% |     17941 | 4.37M | 419 |      3828 | 261.05K | 1483 |   55.0 | 12.0 |   6.3 | 182.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:38 |
| Q20L60X80P000 |   80.0 |  97.97% |     17422 | 4.36M | 437 |      2754 | 247.73K | 1645 |   73.0 | 15.0 |   9.3 | 236.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:38 |
| Q20L60X80P001 |   80.0 |  98.17% |     18026 | 4.37M | 409 |      4597 | 297.67K | 1619 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:41 |
| Q20L60X80P002 |   80.0 |  98.27% |     17627 | 4.36M | 427 |      3630 | 280.24K | 1627 |   74.0 | 15.0 |   9.7 | 238.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:41 |
| Q25L60X40P000 |   40.0 |  98.07% |     12173 | 4.32M | 550 |      3248 | 248.96K | 1562 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  98.05% |     13102 | 4.32M | 531 |      5208 | 232.63K | 1572 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:33 |
| Q25L60X40P002 |   40.0 |  98.00% |     12870 | 4.32M | 552 |      3164 | 209.24K | 1578 |   37.0 |  9.0 |   5.0 | 128.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:32 |
| Q25L60X60P000 |   60.0 |  98.36% |     16760 | 4.35M | 434 |      4536 | 205.14K | 1317 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:38 |
| Q25L60X60P001 |   60.0 |  98.30% |     18325 | 4.35M | 426 |      5241 | 227.76K | 1347 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:36 |
| Q25L60X60P002 |   60.0 |  98.45% |     17636 | 4.36M | 428 |      4783 | 216.38K | 1381 |   56.0 | 13.0 |   5.7 | 190.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:36 |
| Q25L60X80P000 |   80.0 |  98.44% |     18820 | 4.35M | 402 |      3659 | 219.81K | 1295 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:38 |
| Q25L60X80P001 |   80.0 |  98.43% |     18776 | 4.36M | 409 |      4715 | 223.16K | 1329 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:38 |
| Q25L60X80P002 |   80.0 |  98.44% |     18780 | 4.37M | 410 |      2670 | 205.81K | 1347 |   74.0 | 16.0 |   8.7 | 244.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:40 |
| Q30L60X40P000 |   40.0 |  97.35% |      8322 | 4.27M | 750 |      2927 | 265.08K | 2041 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:31 |
| Q30L60X40P001 |   40.0 |  97.35% |      8177 | 4.26M | 743 |      2775 | 260.82K | 1965 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  97.34% |      8299 | 4.27M | 768 |      2940 | 234.88K | 2007 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:32 |
| Q30L60X60P000 |   60.0 |  97.89% |     11189 | 4.34M | 602 |      4603 | 233.36K | 1702 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:36 |
| Q30L60X60P001 |   60.0 |  97.96% |     11512 | 4.33M | 598 |      3851 | 224.09K | 1678 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:35 |
| Q30L60X60P002 |   60.0 |  97.91% |     11669 | 4.32M | 586 |      3155 | 229.79K | 1696 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  98.13% |     13331 | 4.33M | 531 |      4606 | 258.68K | 1596 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:37 |
| Q30L60X80P001 |   80.0 |  98.22% |     13105 | 4.34M | 536 |      4697 | 247.44K | 1589 |   75.0 | 18.0 |   7.0 | 258.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:35 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.84% |     48572 | 4.37M | 198 |      3732 | 112.99K | 471 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:31 |
| MRX40P001 |   40.0 |  97.89% |     42389 | 4.37M | 214 |      4409 | 111.09K | 467 |   37.0 |  7.0 |   5.3 | 116.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:31 |
| MRX40P002 |   40.0 |  97.89% |     45175 | 4.37M | 200 |      5836 | 185.51K | 487 |   36.0 |  7.0 |   5.0 | 114.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:30 |
| MRX60P000 |   60.0 |  97.87% |     46921 | 4.37M | 205 |      5408 | 162.25K | 457 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:32 |
| MRX60P001 |   60.0 |  97.89% |     47121 | 4.37M | 196 |      5084 | 178.01K | 476 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:33 |
| MRX60P002 |   60.0 |  98.06% |     45190 | 4.37M | 199 |      2923 | 116.27K | 457 |   55.0 | 10.0 |   8.3 | 170.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:34 |
| MRX80P000 |   80.0 |  97.85% |     46334 | 4.25M | 200 |      5245 | 159.12K | 445 |   73.0 | 13.0 |  11.3 | 224.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:32 |
| MRX80P001 |   80.0 |  97.93% |     42029 | 4.32M | 203 |      4828 | 153.56K | 449 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:35 |
| MRX80P002 |   80.0 |  97.87% |     44622 | 4.37M | 207 |      8003 | 157.12K | 438 |   73.0 | 12.0 |  12.3 | 218.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:34 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  92.74% |     31611 | 4.37M | 274 |      6869 |    1.4M | 423 |  295.0 | 53.0 |  45.3 | 908.0 |   0:01:02 |
| 7_merge_mr_unitigs_bcalm      |  97.04% |     32839 | 4.25M | 253 |      5227 | 367.98K | 112 |  294.0 | 54.0 |  44.0 | 912.0 |   0:01:31 |
| 7_merge_mr_unitigs_superreads |  97.12% |     39187 | 4.31M | 230 |     21322 | 659.86K | 119 |  291.0 | 56.0 |  41.0 | 918.0 |   0:01:42 |
| 7_merge_mr_unitigs_tadpole    |  96.92% |     41754 | 4.32M | 229 |     25996 | 561.39K |  99 |  291.0 | 57.0 |  40.0 | 924.0 |   0:01:32 |
| 7_merge_unitigs_bcalm         |  97.09% |     30564 | 4.37M | 282 |      4470 | 756.15K | 269 |  289.0 | 56.0 |  40.3 | 914.0 |   0:01:40 |
| 7_merge_unitigs_superreads    |  96.87% |     33060 | 4.37M | 255 |      2653 | 894.37K | 406 |  289.0 | 56.0 |  40.3 | 914.0 |   0:01:49 |
| 7_merge_unitigs_tadpole       |  96.91% |     30193 | 4.36M | 284 |      5256 |  709.2K | 242 |  292.0 | 54.0 |  43.3 | 908.0 |   0:01:39 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |     Sum |   # | N50Others |     Sum |   # | median |  MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|--------:|----:|----------:|--------:|----:|-------:|-----:|------:|-------:|----------:|
| 8_spades     |  98.83% |     29628 | 530.91K |  38 |      8152 |  80.54K |  60 |  298.0 | 55.0 |  44.3 |  926.0 |   0:00:45 |
| 8_mr_spades  |  99.23% |     71110 |    1.6M |  53 |      9161 |  56.39K |  62 |  470.0 | 68.0 |  88.7 | 1348.0 |   0:00:48 |
| 8_megahit    |  98.11% |     26501 |   1.81M | 130 |      8309 | 122.33K | 212 |  297.0 | 55.0 |  44.0 |  924.0 |   0:00:44 |
| 8_mr_megahit |  99.45% |     45159 |    1.8M |  79 |      8353 |  65.14K | 104 |  469.0 | 68.0 |  88.3 | 1346.0 |   0:00:49 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 3188524 | 4602977 |    7 |
| Paralogs                 |    2337 |  146789 |   66 |
| Repetitives              |     572 |   57281 |  165 |
| 7_merge_anchors.anchors  |   31611 | 4370419 |  274 |
| 7_merge_anchors.others   |    6869 | 1402167 |  423 |
| glue_anchors             |   31612 | 4367142 |  263 |
| fill_anchors             |   92767 | 4421405 |  122 |
| spades.contig            |  250095 | 4577425 |   86 |
| spades.scaffold          |  333463 | 4577753 |   82 |
| spades.non-contained     |  250095 | 4567805 |   44 |
| mr_spades.contig         |  172630 | 4586098 |   73 |
| mr_spades.scaffold       |  204122 | 4586240 |   71 |
| mr_spades.non-contained  |  172630 | 4576305 |   47 |
| megahit.contig           |  151747 | 4573489 |  164 |
| megahit.non-contained    |  151747 | 4544657 |  110 |
| mr_megahit.contig        |  156942 | 4590155 |   86 |
| mr_megahit.non-contained |  156942 | 4579130 |   63 |


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
| merged.raw    | 240 | 710.02M | 2943106 |
| unmerged.raw  | 226 |   9.09M |   43970 |
| unmerged.trim | 226 |   9.09M |   43962 |
| M1            | 240 | 524.67M | 2180703 |
| U1            | 239 |   4.95M |   21981 |
| U2            | 213 |   4.14M |   21981 |
| Us            |   0 |       0 |       0 |
| M.cor         | 239 | 535.94M | 4405368 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 198.7 |    193 |  44.4 |         94.73% |
| M.ihist.merge.txt  | 241.2 |    233 |  52.1 |         99.26% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   | 292.7 |  245.1 |   16.26% | "111" | 4.03M | 4.21M |     1.04 | 0:02'08'' |
| Q20L60.R | 287.4 |  243.8 |   15.18% | "113" | 4.03M | 4.17M |     1.03 | 0:02'07'' |
| Q25L60.R | 271.4 |  236.3 |   12.92% | "109" | 4.03M | 4.02M |     1.00 | 0:02'00'' |
| Q30L60.R | 239.3 |  215.2 |   10.06% | "105" | 4.03M | 3.99M |     0.99 | 0:01'48'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  79.71% |      2856 | 3.28M | 1313 |      1007 | 190.17K | 2790 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:27 |
| Q0L0X40P001   |   40.0 |  78.85% |      2812 | 3.27M | 1342 |       770 | 169.63K | 2849 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:27 |
| Q0L0X40P002   |   40.0 |  79.71% |      2956 | 3.29M | 1301 |       786 | 170.35K | 2779 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:25 |
| Q0L0X60P000   |   60.0 |  69.68% |      2254 | 2.98M | 1451 |      1006 | 140.06K | 3037 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:25 |
| Q0L0X60P001   |   60.0 |  69.21% |      2273 | 2.96M | 1445 |      1001 | 125.15K | 3029 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:24 |
| Q0L0X60P002   |   60.0 |  69.51% |      2348 | 2.97M | 1394 |      1001 | 122.69K | 2912 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  60.03% |      1893 | 2.59M | 1417 |      1012 | 152.71K | 2989 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  60.40% |      1895 | 2.61M | 1431 |      1003 | 133.98K | 2991 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  60.83% |      1913 | 2.63M | 1434 |      1004 |  138.9K | 3021 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:25 |
| Q20L60X40P000 |   40.0 |  80.27% |      2878 | 3.33M | 1338 |       692 | 160.79K | 2822 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:24 |
| Q20L60X40P001 |   40.0 |  80.26% |      2896 | 3.33M | 1297 |       436 | 146.12K | 2763 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:27 |
| Q20L60X40P002 |   40.0 |  79.95% |      2933 | 3.31M | 1298 |       538 | 151.96K | 2732 |   37.0 | 10.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:25 |
| Q20L60X60P000 |   60.0 |  71.12% |      2276 | 3.07M | 1471 |       628 | 110.89K | 3055 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:24 |
| Q20L60X60P001 |   60.0 |  70.87% |      2296 | 3.04M | 1441 |       740 | 112.05K | 3007 |   55.0 | 14.0 |   5.0 | 194.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:25 |
| Q20L60X60P002 |   60.0 |  70.24% |      2263 |    3M | 1435 |      1002 | 121.83K | 3002 |   56.0 | 14.0 |   5.0 | 196.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:26 |
| Q20L60X80P000 |   80.0 |  62.47% |      1997 | 2.69M | 1448 |      1008 | 149.16K | 3028 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:24 |
| Q20L60X80P001 |   80.0 |  61.62% |      2010 | 2.65M | 1388 |      1006 | 141.62K | 2913 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:25 |
| Q20L60X80P002 |   80.0 |  62.65% |      1988 |  2.7M | 1437 |      1005 | 139.74K | 3021 |   73.0 | 18.0 |   6.3 | 254.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:25 |
| Q25L60X40P000 |   40.0 |  91.40% |      6344 |  3.7M |  842 |       154 |  91.34K | 1821 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:27 |
| Q25L60X40P001 |   40.0 |  91.08% |      5800 | 3.69M |  859 |       135 |  93.69K | 1843 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  91.21% |      6158 | 3.68M |  852 |       614 | 113.78K | 1874 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:28 |
| Q25L60X60P000 |   60.0 |  87.87% |      4481 | 3.61M | 1042 |        70 |  90.14K | 2168 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:28 |
| Q25L60X60P001 |   60.0 |  87.62% |      4417 |  3.6M | 1051 |       534 | 110.84K | 2194 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:28 |
| Q25L60X60P002 |   60.0 |  87.90% |      4709 |  3.6M | 1004 |       587 | 108.86K | 2106 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:28 |
| Q25L60X80P000 |   80.0 |  84.61% |      3518 |  3.5M | 1188 |       195 |  110.4K | 2475 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:28 |
| Q25L60X80P001 |   80.0 |  84.82% |      3682 | 3.51M | 1167 |       415 | 113.18K | 2435 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  92.74% |      7468 | 3.73M |  701 |       296 |  75.48K | 1525 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:25 |
| Q30L60X40P001 |   40.0 |  92.72% |      7146 | 3.73M |  727 |       159 |  80.32K | 1596 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| Q30L60X40P002 |   40.0 |  92.74% |      7170 | 3.73M |  734 |       257 |  82.38K | 1604 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:28 |
| Q30L60X60P000 |   60.0 |  90.56% |      5521 | 3.68M |  875 |       156 |  82.31K | 1856 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:30 |
| Q30L60X60P001 |   60.0 |  90.38% |      5494 | 3.68M |  895 |        71 |  75.22K | 1884 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:28 |
| Q30L60X60P002 |   60.0 |  90.45% |      5539 | 3.68M |  922 |       221 |  87.05K | 1936 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:29 |
| Q30L60X80P000 |   80.0 |  87.95% |      4395 |  3.6M | 1025 |       443 |  95.88K | 2148 |   78.0 | 17.0 |   9.0 | 258.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:27 |
| Q30L60X80P001 |   80.0 |  88.44% |      4340 | 3.63M | 1059 |        76 |  87.64K | 2216 |   78.0 | 17.0 |   9.0 | 258.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:30 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.11% |     28285 |  3.8M | 250 |       852 | 64.21K | 522 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| MRX40P001 |   40.0 |  96.02% |     33419 | 3.74M | 237 |       789 | 60.08K | 475 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRX40P002 |   40.0 |  95.92% |     30749 | 3.78M | 243 |       635 | 53.07K | 503 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:28 |
| MRX60P000 |   60.0 |  95.76% |     23783 | 3.81M | 281 |       848 | 59.86K | 586 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:30 |
| MRX60P001 |   60.0 |  95.60% |     24989 | 3.81M | 278 |       574 | 46.02K | 567 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:30 |
| MRX80P000 |   80.0 |  95.35% |     21539 |  3.8M | 307 |       856 | 58.85K | 628 |   80.0 | 13.0 |  13.7 | 238.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:30 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |    # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|-----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  94.85% |     10931 | 3.77M |  525 |        64 |  64.34K | 1459 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:31 |
| Q0L0X40P001   |   40.0 |  94.89% |     11347 | 3.78M |  529 |        92 |  74.07K | 1473 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:29 |
| Q0L0X40P002   |   40.0 |  94.32% |     11106 | 3.78M |  527 |        72 |  66.82K | 1436 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:29 |
| Q0L0X60P000   |   60.0 |  91.94% |      6967 | 3.74M |  770 |        73 |  77.12K | 1742 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:29 |
| Q0L0X60P001   |   60.0 |  91.52% |      6868 | 3.73M |  785 |        75 |  78.57K | 1771 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:31 |
| Q0L0X60P002   |   60.0 |  91.82% |      6834 | 3.75M |  758 |        47 |   61.1K | 1697 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:29 |
| Q0L0X80P000   |   80.0 |  87.44% |      4393 | 3.62M | 1065 |        83 |  96.65K | 2251 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:29 |
| Q0L0X80P001   |   80.0 |  87.84% |      4326 | 3.64M | 1072 |        53 |  84.96K | 2238 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:28 |
| Q0L0X80P002   |   80.0 |  87.81% |      4170 | 3.65M | 1082 |        44 |  77.66K | 2269 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:29 |
| Q20L60X40P000 |   40.0 |  94.95% |     11324 | 3.77M |  514 |       104 |  74.78K | 1430 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:30 |
| Q20L60X40P001 |   40.0 |  95.00% |     11784 | 3.78M |  486 |        58 |  60.61K | 1401 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:29 |
| Q20L60X40P002 |   40.0 |  95.09% |     12450 | 3.78M |  489 |        80 |  67.87K | 1390 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:30 |
| Q20L60X60P000 |   60.0 |  92.33% |      7099 | 3.75M |  751 |        70 |  73.09K | 1676 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:30 |
| Q20L60X60P001 |   60.0 |  91.90% |      7216 | 3.73M |  740 |       106 |  84.15K | 1714 |   58.0 | 13.0 |   6.3 | 194.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:29 |
| Q20L60X60P002 |   60.0 |  92.20% |      6934 | 3.75M |  758 |        67 |  72.54K | 1692 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:29 |
| Q20L60X80P000 |   80.0 |  88.15% |      4543 | 3.64M | 1042 |        58 |  86.48K | 2192 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:30 |
| Q20L60X80P001 |   80.0 |  88.36% |      4785 | 3.64M |  988 |        85 |  92.76K | 2111 |   77.0 | 17.0 |   8.7 | 256.0 | "31,41,51,61,71,81" |   0:01:34 |   0:00:28 |
| Q20L60X80P002 |   80.0 |  88.14% |      4888 | 3.62M |  971 |       345 | 105.18K | 2033 |   77.0 | 16.0 |   9.7 | 250.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:29 |
| Q25L60X40P000 |   40.0 |  96.09% |     16510 | 3.81M |  395 |        53 |  45.74K | 1207 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:30 |
| Q25L60X40P001 |   40.0 |  96.29% |     16325 | 3.82M |  379 |        73 |  52.96K | 1148 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:29 |
| Q25L60X40P002 |   40.0 |  96.31% |     16437 | 3.81M |  381 |        71 |   55.6K | 1198 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:29 |
| Q25L60X60P000 |   60.0 |  94.75% |     11136 |  3.8M |  525 |        51 |  47.23K | 1242 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:31 |
| Q25L60X60P001 |   60.0 |  94.72% |     11292 | 3.81M |  525 |        65 |  50.81K | 1207 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:30 |
| Q25L60X60P002 |   60.0 |  94.93% |     10806 | 3.79M |  513 |       100 |  67.58K | 1233 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:29 |
| Q25L60X80P000 |   80.0 |  93.15% |      7888 | 3.76M |  699 |       102 |     81K | 1500 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:31 |
| Q25L60X80P001 |   80.0 |  92.93% |      7612 | 3.76M |  690 |       103 |  79.79K | 1472 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:30 |
| Q30L60X40P000 |   40.0 |  96.69% |     19522 |  3.8M |  355 |        50 |  43.75K | 1155 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:31 |
| Q30L60X40P001 |   40.0 |  96.86% |     17183 | 3.81M |  364 |        62 |  52.37K | 1232 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  96.61% |     18323 |  3.8M |  362 |        60 |  47.87K | 1122 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:30 |
| Q30L60X60P000 |   60.0 |  95.42% |     12262 |  3.8M |  475 |        75 |  56.16K | 1151 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:31 |
| Q30L60X60P001 |   60.0 |  95.22% |     11780 |  3.8M |  479 |        80 |   57.4K | 1175 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:31 |
| Q30L60X60P002 |   60.0 |  95.58% |     13225 | 3.81M |  471 |        71 |  55.64K | 1168 |   59.0 | 12.0 |   7.7 | 190.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:30 |
| Q30L60X80P000 |   80.0 |  93.98% |      9309 | 3.78M |  607 |        77 |  65.25K | 1349 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:29 |
| Q30L60X80P001 |   80.0 |  93.97% |      8358 | 3.78M |  636 |       100 |  74.42K | 1415 |   79.0 | 16.0 |  10.3 | 254.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:30 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.95% |     41036 |  3.8M | 183 |       867 | 51.31K | 443 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:31 |
| MRX40P001 |   40.0 |  96.70% |     46098 | 3.82M | 194 |       634 | 50.93K | 459 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:30 |
| MRX40P002 |   40.0 |  96.87% |     39348 | 3.82M | 193 |       226 | 46.01K | 475 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:30 |
| MRX60P000 |   60.0 |  96.60% |     43064 | 3.77M | 180 |       774 | 48.66K | 397 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:28 |
| MRX60P001 |   60.0 |  96.48% |     39653 | 3.83M | 198 |       557 | 41.25K | 433 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:30 |
| MRX80P000 |   80.0 |  96.36% |     35515 | 3.81M | 213 |       876 | 54.93K | 439 |   80.0 | 12.0 |  14.7 | 232.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:33 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.79% |     16406 | 3.79M | 406 |       568 |  75.17K | 1043 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q0L0X40P001   |   40.0 |  95.86% |     15754 |  3.8M | 391 |       451 |  73.43K | 1057 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:29 |
| Q0L0X40P002   |   40.0 |  95.72% |     15270 | 3.81M | 403 |       161 |  61.45K | 1022 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:27 |
| Q0L0X60P000   |   60.0 |  94.40% |     10214 |  3.8M | 557 |       269 |  69.67K | 1249 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:29 |
| Q0L0X60P001   |   60.0 |  94.21% |     10561 | 3.79M | 567 |       123 |  68.44K | 1301 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:29 |
| Q0L0X60P002   |   60.0 |  94.36% |     10422 |  3.8M | 547 |       112 |  63.04K | 1257 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q0L0X80P000   |   80.0 |  94.31% |      8642 | 3.77M | 674 |       475 | 100.39K | 1747 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:33 |
| Q0L0X80P001   |   80.0 |  93.99% |      7933 | 3.76M | 683 |       118 |  96.64K | 1718 |   80.0 | 16.0 |  10.7 | 256.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:34 |
| Q0L0X80P002   |   80.0 |  94.43% |      7931 | 3.79M | 691 |        91 |  85.25K | 1734 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:34 |
| Q20L60X40P000 |   40.0 |  95.83% |     16164 | 3.79M | 395 |       681 |  76.22K | 1017 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q20L60X40P001 |   40.0 |  95.96% |     17288 |  3.8M | 364 |       141 |  60.63K |  982 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:30 |
| Q20L60X40P002 |   40.0 |  96.10% |     17278 | 3.81M | 368 |       499 |  69.99K | 1014 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:32 |
| Q20L60X60P000 |   60.0 |  94.45% |     11005 | 3.79M | 537 |       380 |  67.59K | 1215 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q20L60X60P001 |   60.0 |  94.48% |     10236 |  3.8M | 541 |       548 |  74.23K | 1252 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:31 |
| Q20L60X60P002 |   60.0 |  94.85% |     11929 | 3.81M | 531 |       118 |  63.78K | 1227 |   59.0 | 13.0 |   6.7 | 196.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:30 |
| Q20L60X80P000 |   80.0 |  94.42% |      8165 | 3.77M | 663 |        90 |  83.57K | 1705 |   80.0 | 17.0 |   9.7 | 262.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:32 |
| Q20L60X80P001 |   80.0 |  95.16% |      9642 | 3.79M | 610 |       281 |  93.86K | 1654 |   81.0 | 17.0 |  10.0 | 264.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| Q20L60X80P002 |   80.0 |  94.11% |      8756 | 3.75M | 644 |       454 | 106.47K | 1670 |   80.0 | 16.0 |  10.7 | 256.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:32 |
| Q25L60X40P000 |   40.0 |  96.78% |     21952 |  3.8M | 320 |        99 |  47.02K |  894 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:29 |
| Q25L60X40P001 |   40.0 |  96.70% |     21095 | 3.82M | 311 |       126 |   47.9K |  808 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:31 |
| Q25L60X40P002 |   40.0 |  96.84% |     22002 | 3.81M | 299 |       381 |  56.21K |  855 |   39.0 |  9.0 |   5.0 | 132.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:30 |
| Q25L60X60P000 |   60.0 |  95.92% |     15129 | 3.82M | 406 |        87 |  44.34K |  945 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:28 |
| Q25L60X60P001 |   60.0 |  96.10% |     15451 | 3.82M | 393 |       459 |  54.87K |  936 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:29 |
| Q25L60X60P002 |   60.0 |  96.01% |     16924 | 3.78M | 377 |       534 |   61.8K |  904 |   60.0 | 13.0 |   7.0 | 198.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q25L60X80P000 |   80.0 |  95.90% |     13944 |  3.8M | 461 |       215 |  72.52K | 1163 |   81.0 | 16.0 |  11.0 | 258.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:32 |
| Q25L60X80P001 |   80.0 |  96.09% |     11592 | 3.82M | 480 |       435 |   73.9K | 1244 |   81.0 | 17.0 |  10.0 | 264.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:33 |
| Q30L60X40P000 |   40.0 |  97.22% |     25125 | 3.82M | 285 |       123 |  46.28K |  830 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q30L60X40P001 |   40.0 |  97.14% |     22083 | 3.82M | 291 |       119 |  48.82K |  842 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  97.10% |     25393 | 3.82M | 292 |       127 |  50.29K |  855 |   40.0 |  9.0 |   5.0 | 134.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:32 |
| Q30L60X60P000 |   60.0 |  96.29% |     18114 | 3.82M | 354 |       505 |  52.84K |  823 |   60.0 | 12.0 |   8.0 | 192.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:30 |
| Q30L60X60P001 |   60.0 |  96.36% |     18401 | 3.81M | 360 |       426 |  58.89K |  843 |   60.0 | 12.0 |   8.0 | 192.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:30 |
| Q30L60X60P002 |   60.0 |  96.34% |     17541 | 3.82M | 358 |       171 |  54.81K |  866 |   60.0 | 12.0 |   8.0 | 192.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:32 |
| Q30L60X80P000 |   80.0 |  96.48% |     14463 | 3.81M | 423 |       526 |   72.3K | 1048 |   82.0 | 16.0 |  11.3 | 260.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:33 |
| Q30L60X80P001 |   80.0 |  96.29% |     14422 | 3.82M | 444 |       486 |  73.07K | 1125 |   81.0 | 16.0 |  11.0 | 258.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:32 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.84% |     42385 | 3.55M | 171 |      1001 | 54.02K | 379 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:28 |
| MRX40P001 |   40.0 |  96.56% |     47072 | 3.69M | 184 |       864 | 53.93K | 385 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:28 |
| MRX40P002 |   40.0 |  96.67% |     41709 | 3.61M | 180 |       717 | 48.63K | 376 |   40.0 |  7.0 |   6.3 | 122.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:28 |
| MRX60P000 |   60.0 |  96.60% |     55001 | 3.56M | 164 |      1001 | 48.37K | 345 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:30 |
| MRX60P001 |   60.0 |  96.51% |     41041 | 3.62M | 178 |       775 | 46.47K | 381 |   60.0 | 10.0 |  10.0 | 180.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:29 |
| MRX80P000 |   80.0 |  96.54% |     41069 | 3.65M | 187 |       864 | 50.02K | 387 |   81.0 | 13.0 |  14.0 | 240.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:30 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.06% |     46530 | 3.83M | 170 |      1929 | 648.63K | 344 |  245.0 | 40.0 |  41.7 | 730.0 |   0:00:50 |
| 7_merge_mr_unitigs_bcalm      |  97.64% |     44124 | 3.83M | 179 |      1138 |  73.11K |  70 |  247.0 | 42.0 |  40.3 | 746.0 |   0:00:58 |
| 7_merge_mr_unitigs_superreads |  97.69% |     42436 | 3.82M | 189 |      1081 |  89.18K |  84 |  246.0 | 40.0 |  42.0 | 732.0 |   0:00:57 |
| 7_merge_mr_unitigs_tadpole    |  97.79% |     46519 | 3.83M | 182 |      1148 |  74.88K |  70 |  248.0 | 41.0 |  41.7 | 742.0 |   0:00:58 |
| 7_merge_unitigs_bcalm         |  98.56% |     58736 | 3.84M | 165 |      1179 | 293.18K | 226 |  244.0 | 44.0 |  37.3 | 752.0 |   0:01:23 |
| 7_merge_unitigs_superreads    |  98.41% |     43165 | 3.85M | 185 |      1358 | 359.18K | 273 |  242.0 | 44.0 |  36.7 | 748.0 |   0:01:21 |
| 7_merge_unitigs_tadpole       |  98.58% |     58777 | 3.84M | 163 |      1814 | 428.91K | 236 |  242.0 | 44.0 |  36.7 | 748.0 |   0:01:21 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  94.17% |     12021 | 3.78M | 491 |       840 | 91.35K | 939 |  245.0 | 43.0 |  38.7 | 748.0 |   0:00:36 |
| 8_mr_spades  |  98.52% |     71096 | 3.12M | 125 |       911 |  49.7K | 226 |  134.0 | 20.0 |  24.7 | 388.0 |   0:00:33 |
| 8_megahit    |  97.19% |     47647 | 3.29M | 151 |       958 | 45.72K | 233 |  247.0 | 42.0 |  40.3 | 746.0 |   0:00:38 |
| 8_mr_megahit |  98.92% |     92324 | 2.88M | 103 |      1084 | 52.68K | 189 |  134.0 | 20.0 |  24.7 | 388.0 |   0:00:32 |


Table: statFinal

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 2961149 | 4033464 |    2 |
| Paralogs                 |    3424 |  119270 |   49 |
| Repetitives              |    1070 |  120471 |  244 |
| 7_merge_anchors.anchors  |   46530 | 3834139 |  170 |
| 7_merge_anchors.others   |    1929 |  648626 |  344 |
| glue_anchors             |   53312 | 3833735 |  161 |
| fill_anchors             |   84202 | 3834961 |  119 |
| spades.contig            |   12898 | 4165902 | 2176 |
| spades.scaffold          |   12986 | 4166848 | 2166 |
| spades.non-contained     |   13389 | 3870353 |  453 |
| mr_spades.contig         |  112695 | 3954942 |  206 |
| mr_spades.scaffold       |  124669 | 3955388 |  200 |
| mr_spades.non-contained  |  112695 | 3919714 |  105 |
| megahit.contig           |   87668 | 3941712 |  181 |
| megahit.non-contained    |   87668 | 3907337 |  111 |
| mr_megahit.contig        |  203715 | 3978805 |  178 |
| mr_megahit.non-contained |  203715 | 3939555 |   91 |

