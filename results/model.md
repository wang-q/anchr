# Assemble genomes of model organisms by `anchr`

<!-- TOC -->
* [Assemble genomes of model organisms by `anchr`](#assemble-genomes-of-model-organisms-by-anchr)
  * [*Mycoplasma genitalium* G37](#mycoplasma-genitalium-g37)
    * [g37: reference](#g37-reference)
    * [g37: download](#g37-download)
    * [g37: template](#g37-template)
    * [g37: run](#g37-run)
  * [*E. coli* str. K-12 substr. MG1655](#e-coli-str-k-12-substr-mg1655)
    * [mg1655: reference](#mg1655-reference)
    * [mg1655: download](#mg1655-download)
    * [mg1655: template](#mg1655-template)
    * [mg1655: run](#mg1655-run)
  * [*E. coli* str. K-12 substr. DH5alpha](#e-coli-str-k-12-substr-dh5alpha)
    * [dh5alpha: reference](#dh5alpha-reference)
    * [dh5alpha: download](#dh5alpha-download)
    * [dh5alpha: template](#dh5alpha-template)
    * [dh5alpha: run](#dh5alpha-run)
<!-- TOC -->

## *Mycoplasma genitalium* G37

### g37: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/g37/1_genome
cd ~/data/anchr/g37/1_genome

cp ~/data/anchr/ref/g37/genome.fa .
cp ~/data/anchr/ref/g37/paralogs.fa .

```

### g37: download

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

```

| name | srx       | platform | layout | ilength | srr       |   spots | bases  |
|------|-----------|----------|--------|--------:|-----------|--------:|--------|
| G37  | ERX452667 | ILLUMINA | PAIRED |     447 | ERR486835 | 680,644 | 97.37M |

* Illumina

```shell
cd ~/data/anchr/g37

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/ERR486835_1.fastq.gz R1.fq.gz
ln -s ../ena/ERR486835_2.fastq.gz R2.fq.gz

```

### g37: template

* template

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=g37

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 580076 \
    --parallel 8 \
    --xmx 12g \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "bcalm bifrost superreads tadpole" \
    --statp 2 \
    --readl 125 \
    --uscale 2 \
    --lscale 3 \
    --redo

```

### g37: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=g37

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge 9_markdown/statReads.md
# rm -fr 4_down_sampling 6_down_sampling

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

| Group             |   Mean | Median |   STDev | Pairs%/Orientation |
|-------------------|-------:|-------:|--------:|-------------------:|
| R.genome.bbtools  | 3491.1 |    466 | 14779.2 |             96.05% |
| R.tadpole.bbtools |  468.8 |    452 |   136.0 |             83.54% |
| R.genome.picard   |  475.2 |    460 |   133.4 |                 FR |
| R.tadpole.picard  |  468.6 |    452 |   132.0 |                 FR |

Table: statInsertSize

| K    | property              |        min |        max |
|------|-----------------------|-----------:|-----------:|
| R.21 | Homozygous (a)        |            |       100% |
|      | Genome Haploid Length |            | 577,872 bp |
|      | Genome Repeat Length  |   3,638 bp |   3,643 bp |
|      | Genome Unique Length  | 573,811 bp | 574,653 bp |
|      | Model Fit             |   92.6091% |   93.3366% |
|      | Read Error Rate       |            |  0.137192% |
|      | Kmer Cov              |            |      148.8 |
| R.51 | Homozygous (a)        |            |       100% |
|      | Genome Haploid Length |            | 578,025 bp |
|      | Genome Repeat Length  |            |       0 bp |
|      | Genome Unique Length  |            | 578,025 bp |
|      | Model Fit             |   95.6385% |   95.7255% |
|      | Read Error Rate       |            | 0.0942288% |
|      | Kmer Cov              |            |      112.2 |
| R.81 | Homozygous (a)        |            |       100% |
|      | Genome Haploid Length |            | 578,301 bp |
|      | Genome Repeat Length  |            |       0 bp |
|      | Genome Unique Length  |            | 578,301 bp |
|      | Model Fit             |   97.0931% |   97.1094% |
|      | Read Error Rate       |            | 0.0778153% |
|      | Kmer Cov              |            |       77.3 |

Table: statFastK

| Name       |    N50 |     Sum |      # |
|------------|-------:|--------:|-------:|
| genome     | 580076 | 580.08K |      1 |
| repetitive |    184 |   6.91K |     41 |
| Illumina.R |    150 |  102.1M | 680644 |
| trim.R     |    150 | 101.56M | 677478 |
| Q0L0       |    150 | 101.56M | 677478 |
| Q25L60     |    150 |  98.27M | 658139 |
| Q30L60     |    150 |  95.05M | 638376 |

Table: statReads

| Name     | N50 |     Sum |      # |
|----------|----:|--------:|-------:|
| clumpify | 150 | 102.01M | 680076 |
| highpass | 150 | 101.67M | 677786 |
| trim     | 150 | 101.56M | 677478 |
| filter   | 150 | 101.56M | 677478 |
| R1       | 150 |   50.8M | 338739 |
| R2       | 150 |  50.76M | 338739 |
| Rs       |   0 |       0 |      0 |

Table: statTrimReads

```text
#R.trim
#Matched        602     0.08882%
#Name   Reads   ReadsPct
```

```text
#R.filter
#Matched        0       0.00000%
#Name   Reads   ReadsPct
```

| Name          | N50 |     Sum |      # |
|---------------|----:|--------:|-------:|
| clumped       | 150 | 101.56M | 677478 |
| ecco          | 150 | 101.56M | 677478 |
| eccc          | 150 | 101.56M | 677478 |
| ecct          | 150 |  98.41M | 656430 |
| extended      | 190 | 124.33M | 656430 |
| merged.raw    | 443 |  77.84M | 183173 |
| unmerged.raw  | 190 |  54.78M | 290084 |
| unmerged.trim | 190 |  54.78M | 290084 |
| M1            | 443 |   77.8M | 183085 |
| U1            | 190 |   27.5M | 145042 |
| U2            | 190 |  27.28M | 145042 |
| Us            |   0 |       0 |      0 |
| M.cor         | 363 | 132.77M | 656254 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
|--------------------|------:|-------:|------:|-------:|
| M.ihist.merge1.txt | 248.0 |    254 |  23.8 |  2.80% |
| M.ihist.merge.txt  | 425.0 |    431 |  64.8 | 55.81% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% |  Kmer |   RealG |    EstG | Est/Real | RunTime |
|----------|------:|-------:|---------:|------:|--------:|--------:|---------:|--------:|
| Q0L0.R   | 175.1 |  160.8 |    8.16% | "105" | 580.08K | 584.05K |     1.01 | 0:00:18 |
| Q25L60.R | 169.4 |  160.4 |    5.31% | "105" | 580.08K | 582.89K |     1.00 | 0:00:17 |
| Q30L60.R | 163.9 |  157.3 |    4.05% | "105" | 580.08K | 582.28K |     1.00 | 0:00:17 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |     Sum |   # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|--------:|----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.962 |     16645 |  561.5K |  54 |     3.62K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P001   |   40.0 |  0.963 |     15870 | 560.35K |  54 |     4.22K |     40 |   5 |   8.3 |   110 |
| Q0L0X40P002   |   40.0 |  0.961 |     16544 | 551.23K |  50 |     3.39K |     39 |   6 |     7 |   114 |
| Q0L0X80P000   |   80.0 |  0.934 |      7710 | 546.38K | 100 |      8.3K |     79 |   9 |  17.3 |   212 |
| Q0L0X80P001   |   80.0 |  0.938 |      7387 |  549.5K |  99 |     7.18K |     79 |   9 |  17.3 |   212 |
| Q25L60X40P000 |   40.0 |  0.963 |     14828 |  561.2K |  50 |     3.15K |     39 |   5 |     8 |   108 |
| Q25L60X40P001 |   40.0 |  0.965 |     22403 | 560.41K |  45 |     3.32K |     40 |   5 |   8.3 |   110 |
| Q25L60X40P002 |   40.0 |  0.963 |     25424 | 559.68K |  44 |     3.16K |     39 |   5 |     8 |   108 |
| Q25L60X80P000 |   80.0 |  0.939 |      7531 | 549.84K | 100 |     6.78K |     79 |   9 |  17.3 |   212 |
| Q25L60X80P001 |   80.0 |  0.935 |      7523 | 547.71K | 103 |     8.53K |     79 |   9 |  17.3 |   212 |
| Q30L60X40P000 |   40.0 |  0.962 |     20001 |    555K |  48 |     2.93K |     39 |   6 |     7 |   114 |
| Q30L60X40P001 |   40.0 |  0.964 |     17542 | 561.74K |  50 |     3.26K |     39 |   6 |     7 |   114 |
| Q30L60X40P002 |   40.0 |  0.966 |     22885 | 559.47K |  39 |     2.96K |     40 |   6 |   7.3 |   116 |
| Q30L60X80P000 |   80.0 |  0.937 |      7324 | 548.48K | 102 |     7.56K |     79 |  10 |  16.3 |   218 |

Table: statUnitigsBcalm.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.951 |     31168 | 561.41K | 32 |     2.27K |     39 |   6 |     7 |   114 |
| MRX40P001 |   40.0 |   0.95 |     20278 | 560.84K | 40 |     3.79K |     39 |   5 |     8 |   108 |
| MRX40P002 |   40.0 |  0.953 |     34102 | 561.41K | 31 |     3.05K |     39 |   5 |     8 |   108 |
| MRX80P000 |   80.0 |  0.935 |     13951 |  554.4K | 62 |     5.36K |     78 |   9 |    17 |   210 |
| MRX80P001 |   80.0 |  0.937 |     16095 | 555.84K | 59 |     5.06K |     78 |   9 |    17 |   210 |

Table: statMRUnitigsBcalm.md

| Name          | CovCor | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.967 |     55036 | 561.57K | 24 |     2.04K |     40 |   5 |   8.3 |   110 |
| Q0L0X40P001   |   40.0 |  0.969 |     31700 | 561.88K | 27 |     2.26K |     40 |   5 |   8.3 |   110 |
| Q0L0X40P002   |   40.0 |  0.969 |     41526 | 556.84K | 22 |     2.05K |     40 |   5 |   8.3 |   110 |
| Q0L0X80P000   |   80.0 |  0.954 |     17281 | 554.89K | 53 |        5K |     79 |   9 |  17.3 |   212 |
| Q0L0X80P001   |   80.0 |  0.957 |     16521 | 557.12K | 51 |     3.94K |     79 |   9 |  17.3 |   212 |
| Q25L60X40P000 |   40.0 |  0.968 |     55036 | 561.56K | 19 |      1.6K |     40 |   5 |   8.3 |   110 |
| Q25L60X40P001 |   40.0 |  0.967 |     37688 | 560.76K | 25 |     2.05K |     40 |   5 |   8.3 |   110 |
| Q25L60X40P002 |   40.0 |  0.969 |     39432 | 561.66K | 26 |     2.04K |     39 |   5 |     8 |   108 |
| Q25L60X80P000 |   80.0 |  0.955 |     22617 | 554.97K | 48 |     4.59K |     79 |   9 |  17.3 |   212 |
| Q25L60X80P001 |   80.0 |  0.953 |     19834 | 556.14K | 52 |     3.68K |     79 |   9 |  17.3 |   212 |
| Q30L60X40P000 |   40.0 |  0.967 |     36712 | 561.71K | 23 |     1.67K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P001 |   40.0 |  0.967 |     33379 | 561.15K | 27 |        2K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P002 |   40.0 |   0.97 |     47384 |  562.3K | 21 |     1.76K |     40 |   6 |   7.3 |   116 |
| Q30L60X80P000 |   80.0 |  0.956 |     16042 | 557.01K | 55 |     4.08K |     79 |   9 |  17.3 |   212 |

Table: statUnitigsBifrost.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.949 |     30386 | 560.29K | 32 |     2.44K |     39 |   6 |     7 |   114 |
| MRX40P001 |   40.0 |  0.948 |     24091 | 559.18K | 40 |     3.45K |     39 |   5 |     8 |   108 |
| MRX40P002 |   40.0 |   0.95 |     27042 | 560.48K | 36 |     3.38K |     39 |   5 |     8 |   108 |
| MRX80P000 |   80.0 |  0.934 |     13950 | 554.78K | 58 |     4.41K |     78 |   9 |    17 |   210 |
| MRX80P001 |   80.0 |  0.937 |     19853 | 553.93K | 53 |     5.93K |     78 |   9 |    17 |   210 |

Table: statMRUnitigsBifrost.md

| Name          | CovCor | Mapped | N50Anchor |     Sum |   # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|--------:|----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.939 |      8290 | 550.49K |  95 |     6.46K |     39 |   6 |     7 |   114 |
| Q0L0X40P001   |   40.0 |  0.941 |      6927 | 551.38K | 100 |     6.42K |     39 |   5 |     8 |   108 |
| Q0L0X40P002   |   40.0 |  0.941 |      7850 |    552K |  95 |     5.45K |     39 |   6 |     7 |   114 |
| Q0L0X80P000   |   80.0 |  0.911 |      5170 |  533.4K | 127 |     9.06K |     79 |   9 |  17.3 |   212 |
| Q0L0X80P001   |   80.0 |  0.918 |      5083 | 531.24K | 127 |     9.27K |     79 |   9 |  17.3 |   212 |
| Q25L60X40P000 |   40.0 |  0.947 |      8503 | 554.39K |  83 |     4.65K |     39 |   6 |     7 |   114 |
| Q25L60X40P001 |   40.0 |  0.944 |      8444 | 551.68K |  91 |     7.36K |     40 |   5 |   8.3 |   110 |
| Q25L60X40P002 |   40.0 |  0.943 |      8051 |    545K |  87 |     6.72K |     39 |   5 |     8 |   108 |
| Q25L60X80P000 |   80.0 |  0.914 |      5004 | 529.22K | 127 |     8.78K |     79 |   9 |  17.3 |   212 |
| Q25L60X80P001 |   80.0 |  0.914 |      5086 | 528.78K | 127 |    10.09K |     79 |   9 |  17.3 |   212 |
| Q30L60X40P000 |   40.0 |  0.943 |      7561 | 553.55K |  93 |     5.01K |     39 |   6 |     7 |   114 |
| Q30L60X40P001 |   40.0 |  0.947 |      7570 | 549.16K |  96 |     5.28K |     39 |   6 |     7 |   114 |
| Q30L60X40P002 |   40.0 |  0.947 |      8776 | 554.33K |  82 |     4.75K |     39 |   6 |     7 |   114 |
| Q30L60X80P000 |   80.0 |  0.916 |      4977 | 531.43K | 125 |     8.11K |     79 |  10 |  16.3 |   218 |

Table: statUnitigsSuperreads.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.943 |     20012 | 559.49K | 48 |     3.57K |     39 |   6 |     7 |   114 |
| MRX40P001 |   40.0 |  0.941 |     16105 | 555.28K | 56 |     5.97K |     39 |   5 |     8 |   108 |
| MRX40P002 |   40.0 |  0.943 |     17813 | 555.05K | 49 |     5.77K |     39 |   5 |     8 |   108 |
| MRX80P000 |   80.0 |  0.931 |     12105 | 551.71K | 71 |     7.11K |     79 |   9 |  17.3 |   212 |
| MRX80P001 |   80.0 |  0.929 |     12558 | 543.61K | 71 |     8.42K |     78 |   9 |    17 |   210 |

Table: statMRUnitigsSuperreads.md

| Name          | CovCor | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.971 |     27112 |    563K | 32 |    12.06K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P001   |   40.0 |   0.97 |     34506 | 561.96K | 26 |     2.47K |     40 |   5 |   8.3 |   110 |
| Q0L0X40P002   |   40.0 |  0.972 |     27099 | 555.38K | 29 |     3.34K |     39 |   6 |     7 |   114 |
| Q0L0X80P000   |   80.0 |  0.952 |     11428 | 555.58K | 72 |     9.94K |     79 |   9 |  17.3 |   212 |
| Q0L0X80P001   |   80.0 |  0.957 |     10227 |  558.6K | 70 |     5.57K |     79 |   9 |  17.3 |   212 |
| Q25L60X40P000 |   40.0 |  0.975 |     30330 | 562.75K | 27 |     2.93K |     40 |   6 |   7.3 |   116 |
| Q25L60X40P001 |   40.0 |  0.971 |     31181 | 560.38K | 26 |      2.5K |     40 |   5 |   8.3 |   110 |
| Q25L60X40P002 |   40.0 |  0.974 |     30360 | 555.71K | 24 |     2.63K |     39 |   5 |     8 |   108 |
| Q25L60X80P000 |   80.0 |  0.958 |     13232 | 555.44K | 64 |     6.44K |     79 |   9 |  17.3 |   212 |
| Q25L60X80P001 |   80.0 |  0.955 |     10691 | 557.93K | 70 |     5.41K |     79 |   9 |  17.3 |   212 |
| Q30L60X40P000 |   40.0 |  0.969 |     39424 | 559.56K | 23 |     2.36K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P001 |   40.0 |  0.973 |     36733 | 562.13K | 27 |     2.84K |     40 |   5 |   8.3 |   110 |
| Q30L60X40P002 |   40.0 |  0.974 |     40818 | 562.62K | 24 |     2.53K |     40 |   5 |   8.3 |   110 |
| Q30L60X80P000 |   80.0 |  0.958 |     10286 | 559.02K | 71 |      5.2K |     79 |  10 |  16.3 |   218 |

Table: statUnitigsTadpole.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.956 |     41644 | 561.88K | 23 |     2.53K |     39 |   5 |     8 |   108 |
| MRX40P001 |   40.0 |  0.963 |     34437 | 561.39K | 25 |     3.48K |     40 |   5 |   8.3 |   110 |
| MRX40P002 |   40.0 |   0.96 |     39348 | 561.54K | 24 |     3.03K |     40 |   5 |   8.3 |   110 |
| MRX80P000 |   80.0 |  0.955 |     20286 | 559.04K | 42 |    23.35K |     80 |   9 |  17.7 |   214 |
| MRX80P001 |   80.0 |  0.953 |     24000 | 558.97K | 41 |     4.67K |     80 |   9 |  17.7 |   214 |

Table: statMRUnitigsTadpole.md

| Name                          | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|-------------------------------|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| 7_merge_anchors               |   0.99 |     55041 | 562.59K | 17 |     35.8K |    161 |  16 |  37.7 |   418 |
| 7_merge_mr_unitigs_bcalm      |  0.986 |     48771 | 562.08K | 18 |         0 |    162 |  16 |    38 |   420 |
| 7_merge_mr_unitigs_bifrost    |  0.986 |     48778 | 562.27K | 19 |     1.11K |    162 |  16 |    38 |   420 |
| 7_merge_mr_unitigs_superreads |  0.986 |     41943 | 562.55K | 21 |     2.25K |    160 |  17 |  36.3 |   422 |
| 7_merge_mr_unitigs_tadpole    |  0.988 |     48811 | 562.12K | 18 |    18.98K |    159 |  17 |    36 |   420 |
| 7_merge_unitigs_bcalm         |  0.989 |     55039 | 562.99K | 18 |     5.68K |    157 |  16 |  36.3 |   410 |
| 7_merge_unitigs_bifrost       |  0.988 |     55037 | 562.43K | 17 |         0 |    158 |  16 |  36.7 |   412 |
| 7_merge_unitigs_superreads    |  0.988 |     34510 | 563.84K | 28 |     2.19K |    161 |  17 |  36.7 |   424 |
| 7_merge_unitigs_tadpole       |   0.99 |     55057 |  562.5K | 17 |     9.05K |    158 |  16 |  36.7 |   412 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |     Sum |  # | SumOthers | median | MAD | lower | upper |
|--------------|-------:|----------:|--------:|---:|----------:|-------:|----:|------:|------:|
| 8_spades     |  0.984 |     83396 | 572.71K |  8 |       667 |    160 |  16 |  37.3 |   416 |
| 8_mr_spades  |  0.982 |    395772 | 577.43K |  4 |       515 |    226 |  22 |  53.3 |   584 |
| 8_megahit    |  0.972 |     54910 | 562.14K | 15 |     4.62K |    160 |  16 |  37.3 |   416 |
| 8_mr_megahit |  0.984 |    357625 | 579.08K |  5 |      3.4K |    226 |  22 |  53.3 |   584 |

Table: statOtherAnchors.md

| Name                     |    N50 |     Sum |  # |
|--------------------------|-------:|--------:|---:|
| Genome                   | 580076 | 580.08K |  1 |
| repetitive               |    184 |   6.91K | 41 |
| 7_merge_anchors.anchors  |  55041 | 562.59K | 17 |
| 7_merge_anchors.others   |   9048 |   35.8K |  6 |
| spades.contig            | 163847 | 581.72K | 41 |
| spades.scaffold          | 163847 | 581.72K | 41 |
| spades.non-contained     | 163847 | 573.38K |  6 |
| mr_spades.contig         | 395878 | 579.84K |  7 |
| mr_spades.scaffold       | 395878 | 579.84K |  7 |
| mr_spades.non-contained  | 395878 | 577.95K |  4 |
| megahit.contig           |  54990 | 584.37K | 65 |
| megahit.non-contained    |  54990 | 566.76K | 17 |
| mr_megahit.contig        | 357750 | 598.09K | 41 |
| mr_megahit.non-contained | 357750 | 582.48K |  7 |

Table: statFinal

## *E. coli* str. K-12 substr. MG1655

### mg1655: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/mg1655/1_genome
cd ~/data/anchr/mg1655/1_genome

cp ~/data/anchr/ref/mg1655/genome.fa .
cp ~/data/anchr/ref/mg1655/paralogs.fa .

```

### mg1655: download

```shell
mkdir -p ~/data/anchr/mg1655/ena
cd ~/data/anchr/mg1655/ena

aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz
aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz

```

* Illumina

```shell
cd ~/data/anchr/mg1655

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz R1.fq.gz
ln -s ../ena/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz R2.fq.gz

```

### mg1655: template

* Rsync to hpcc

```shell
rsync -avP \
    ~/data/anchr/mg1655/ \
    wangq@202.119.37.251:data/anchr/mg1655

# rsync -avP wangq@202.119.37.251:data/anchr/mg1655/ ~/data/anchr/mg1655

```

* template

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4641652 \
    --parallel 8 \
    --xmx 24g \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --tile --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --bwa "Q25L60" \
    --gatk \
    \
    --cov "40 80" \
    --unitigger "bcalm bifrost superreads" \
    --statp 2 \
    --readl 151 \
    --uscale 2 \
    --lscale 3 \
    \
    --extend \
    \
    --busco

```

### mg1655: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge 9_markdown/statReads.md
# rm -fr 4_down_sampling 6_down_sampling

bash 0_script/1_repetitive.sh

# BASE_NAME=mg1655 bash 0_script/0_bsub.sh
#bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_script/0_master.sh"
# bkill -J "${BASE_NAME}-*"

bash 0_script/0_master.sh
# bash 0_script/0_cleanup.sh

#bash 0_script/9_busco.sh

```

| Group             |  Mean | Median |  STDev | Pairs%/Orientation |
|-------------------|------:|-------:|-------:|-------------------:|
| R.genome.bbtools  | 345.2 |    298 | 1928.5 |             96.08% |
| R.tadpole.bbtools | 294.9 |    296 |   22.0 |             81.44% |
| R.genome.picard   | 298.2 |    298 |   18.0 |                 FR |
| R.tadpole.picard  | 295.0 |    296 |   21.6 |                 FR |

Table: statInsertSize

| K    | property              |          min |          max |
|------|-----------------------|-------------:|-------------:|
| R.21 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,478,083 bp |
|      | Genome Repeat Length  |   136,755 bp |   136,883 bp |
|      | Genome Unique Length  | 4,339,242 bp | 4,343,288 bp |
|      | Model Fit             |     97.2749% |     97.3934% |
|      | Read Error Rate       |              |    0.531821% |
|      | Kmer Cov              |              |        299.7 |
| R.51 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,385,263 bp |
|      | Genome Repeat Length  |    91,444 bp |    91,569 bp |
|      | Genome Unique Length  | 4,290,813 bp | 4,296,704 bp |
|      | Model Fit             |     97.3222% |     97.6265% |
|      | Read Error Rate       |              |    0.326106% |
|      | Kmer Cov              |              |        223.4 |
| R.81 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,304,378 bp |
|      | Genome Repeat Length  |    60,194 bp |    60,313 bp |
|      | Genome Unique Length  | 4,239,931 bp | 4,248,326 bp |
|      | Model Fit             |     96.7883% |     97.4791% |
|      | Read Error Rate       |              |    0.272491% |
|      | Kmer Cov              |              |        151.5 |

Table: statFastK

| Name       |     N50 |     Sum |        # |
|------------|--------:|--------:|---------:|
| genome     | 4641652 |   4.64M |        1 |
| paralogs   |    2003 | 260.35K |      131 |
| repetitive |    1235 |  91.99K |      169 |
| Illumina.R |     151 |   1.73G | 11458940 |
| trim.R     |     149 |   1.44G | 10457850 |
| Q0L0       |     149 |   1.44G | 10457850 |
| Q25L60     |     148 |   1.33G | 10022456 |
| Q30L60     |     128 |   1.11G |  9372881 |

Table: statReads

| Name           | N50 |     Sum |        # |
|----------------|----:|--------:|---------:|
| clumpify       | 151 |   1.73G | 11439000 |
| filteredbytile | 151 |   1.69G | 11172376 |
| highpass       | 151 |   1.68G | 11096920 |
| trim           | 149 |   1.44G | 10457850 |
| filter         | 149 |   1.44G | 10457850 |
| R1             | 150 | 741.67M |  5228925 |
| R2             | 144 | 695.43M |  5228925 |
| Rs             |   0 |       0 |        0 |

Table: statTrimReads

```text
#R.trim
#Matched        17628   0.15885%
#Name   Reads   ReadsPct
```

```text
#R.filter
#Matched        0       0.00000%
#Name   Reads   ReadsPct
```

| Name          | N50 |    Sum |        # |
|---------------|----:|-------:|---------:|
| clumped       | 149 |  1.44G | 10456678 |
| ecco          | 149 |  1.44G | 10456678 |
| eccc          | 149 |  1.44G | 10456678 |
| ecct          | 149 |  1.43G | 10410416 |
| extended      | 189 |  1.84G | 10410416 |
| merged.raw    | 339 |  1.73G |  5136391 |
| unmerged.raw  | 174 | 20.15M |   137634 |
| unmerged.trim | 174 | 20.14M |   137580 |
| M1            | 339 |  1.72G |  5108318 |
| U1            | 181 | 10.58M |    68790 |
| U2            | 168 |  9.55M |    68790 |
| Us            |   0 |      0 |        0 |
| M.cor         | 338 |  1.75G | 10354216 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
|--------------------|------:|-------:|------:|-------:|
| M.ihist.merge1.txt | 248.2 |    263 |  40.4 |  2.42% |
| M.ihist.merge.txt  | 337.7 |    338 |  19.3 | 98.68% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real | RunTime |
|----------|------:|-------:|---------:|-----:|------:|------:|---------:|--------:|
| Q0L0.R   | 309.6 |  287.5 |    7.13% | "95" | 4.64M | 4.66M |     1.00 | 0:02:51 |
| Q25L60.R | 285.8 |  274.7 |    3.90% | "89" | 4.64M | 4.57M |     0.98 | 0:02:41 |
| Q30L60.R | 239.4 |  234.4 |    2.09% | "73" | 4.64M | 4.55M |     0.98 | 0:02:16 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |   0.98 |     19489 | 4.53M | 384 |    56.77K |     40 |   7 |   6.3 |   122 |
| Q0L0X40P001   |   40.0 |  0.979 |     18219 | 4.53M | 389 |    60.95K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P002   |   40.0 |  0.978 |     17755 | 4.53M | 401 |    40.78K |     39 |   7 |     6 |   120 |
| Q0L0X80P000   |   80.0 |  0.958 |      9214 | 4.46M | 708 |    58.92K |     79 |  11 |  15.3 |   224 |
| Q0L0X80P001   |   80.0 |  0.959 |      8933 | 4.46M | 718 |    62.87K |     79 |  12 |  14.3 |   230 |
| Q0L0X80P002   |   80.0 |   0.96 |      9180 | 4.48M | 705 |    51.64K |     79 |  12 |  14.3 |   230 |
| Q25L60X40P000 |   40.0 |  0.984 |     23623 | 4.54M | 311 |    53.15K |     40 |   7 |   6.3 |   122 |
| Q25L60X40P001 |   40.0 |  0.983 |     22120 | 4.54M | 331 |    52.09K |     40 |   7 |   6.3 |   122 |
| Q25L60X40P002 |   40.0 |  0.984 |     23822 | 4.54M | 322 |    49.33K |     40 |   7 |   6.3 |   122 |
| Q25L60X80P000 |   80.0 |  0.977 |     15808 | 4.51M | 449 |    49.69K |     79 |  12 |  14.3 |   230 |
| Q25L60X80P001 |   80.0 |  0.976 |     15319 | 4.52M | 455 |    46.37K |     79 |  12 |  14.3 |   230 |
| Q25L60X80P002 |   80.0 |  0.978 |     16666 | 4.51M | 447 |    39.95K |     80 |  12 |  14.7 |   232 |
| Q30L60X40P000 |   40.0 |  0.985 |     23577 | 4.53M | 326 |    49.16K |     40 |   8 |   5.3 |   128 |
| Q30L60X40P001 |   40.0 |  0.984 |     23046 | 4.53M | 325 |    45.81K |     40 |   8 |   5.3 |   128 |
| Q30L60X40P002 |   40.0 |  0.984 |     22191 | 4.53M | 346 |     48.1K |     40 |   8 |   5.3 |   128 |
| Q30L60X80P000 |   80.0 |  0.986 |     31123 | 4.53M | 252 |    44.31K |     80 |  14 |  12.7 |   244 |
| Q30L60X80P001 |   80.0 |  0.985 |     29653 | 4.53M | 247 |    33.34K |     80 |  14 |  12.7 |   244 |

Table: statUnitigsBcalm.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.968 |     65283 | 4.52M | 129 |    18.68K |     39 |   6 |     7 |   114 |
| MRX40P001 |   40.0 |  0.967 |     57785 | 4.52M | 129 |    19.55K |     39 |   6 |     7 |   114 |
| MRX40P002 |   40.0 |  0.968 |     58708 | 4.52M | 138 |    23.23K |     40 |   6 |   7.3 |   116 |
| MRX80P000 |   80.0 |  0.965 |     43693 | 4.51M | 175 |    28.69K |     79 |  10 |  16.3 |   218 |
| MRX80P001 |   80.0 |  0.963 |     46313 | 4.51M | 185 |    30.08K |     79 |  10 |  16.3 |   218 |
| MRX80P002 |   80.0 |  0.965 |     44530 | 4.51M | 178 |     28.4K |     79 |   9 |  17.3 |   212 |

Table: statMRUnitigsBcalm.md

| Name          | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.985 |     56446 | 4.51M | 151 |    32.06K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P001   |   40.0 |  0.984 |     57795 | 4.51M | 149 |    29.33K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P002   |   40.0 |  0.984 |     58627 | 4.51M | 157 |    27.02K |     40 |   7 |   6.3 |   122 |
| Q0L0X80P000   |   80.0 |  0.981 |     59618 | 4.52M | 145 |    22.59K |     80 |  11 |  15.7 |   226 |
| Q0L0X80P001   |   80.0 |  0.981 |     57926 | 4.53M | 144 |    22.46K |     79 |  11 |  15.3 |   224 |
| Q0L0X80P002   |   80.0 |  0.982 |     57910 | 4.53M | 151 |    22.16K |     79 |  11 |  15.3 |   224 |
| Q25L60X40P000 |   40.0 |  0.986 |     44142 | 4.52M | 180 |    30.87K |     40 |   7 |   6.3 |   122 |
| Q25L60X40P001 |   40.0 |  0.986 |     41602 | 4.53M | 178 |    30.85K |     40 |   7 |   6.3 |   122 |
| Q25L60X40P002 |   40.0 |  0.985 |     40710 | 4.53M | 188 |    31.21K |     40 |   7 |   6.3 |   122 |
| Q25L60X80P000 |   80.0 |  0.986 |     58773 | 4.53M | 156 |    26.48K |     80 |  12 |  14.7 |   232 |
| Q25L60X80P001 |   80.0 |  0.985 |     54563 | 4.53M | 156 |    29.18K |     80 |  12 |  14.7 |   232 |
| Q25L60X80P002 |   80.0 |  0.986 |     57786 | 4.51M | 152 |       27K |     80 |  12 |  14.7 |   232 |
| Q30L60X40P000 |   40.0 |  0.986 |     33949 | 4.53M | 243 |    51.15K |     40 |   8 |   5.3 |   128 |
| Q30L60X40P001 |   40.0 |  0.986 |     31239 | 4.54M | 250 |     41.1K |     40 |   8 |   5.3 |   128 |
| Q30L60X40P002 |   40.0 |  0.985 |     31012 | 4.54M | 256 |    41.67K |     40 |   8 |   5.3 |   128 |
| Q30L60X80P000 |   80.0 |  0.987 |     37684 | 4.53M | 198 |    31.01K |     80 |  14 |  12.7 |   244 |
| Q30L60X80P001 |   80.0 |  0.986 |     38607 | 4.53M | 198 |    28.51K |     80 |  14 |  12.7 |   244 |

Table: statUnitigsBifrost.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.968 |     73567 | 4.52M | 112 |    16.74K |     39 |   6 |     7 |   114 |
| MRX40P001 |   40.0 |  0.967 |     65265 | 4.52M | 114 |     17.2K |     39 |   6 |     7 |   114 |
| MRX40P002 |   40.0 |  0.969 |     67291 | 4.51M | 112 |    19.87K |     40 |   6 |   7.3 |   116 |
| MRX80P000 |   80.0 |  0.967 |     63088 | 4.52M | 119 |    19.74K |     79 |  10 |  16.3 |   218 |
| MRX80P001 |   80.0 |  0.967 |     63108 | 4.52M | 120 |    19.57K |     79 |  10 |  16.3 |   218 |
| MRX80P002 |   80.0 |  0.968 |     58677 | 4.52M | 121 |    20.21K |     79 |   9 |  17.3 |   212 |

Table: statMRUnitigsBifrost.md

| Name          | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|------:|-----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |   0.96 |      9749 | 4.47M |  670 |    44.27K |     39 |   7 |     6 |   120 |
| Q0L0X40P001   |   40.0 |  0.961 |      9934 | 4.48M |  661 |    43.37K |     39 |   7 |     6 |   120 |
| Q0L0X40P002   |   40.0 |  0.961 |      9085 |  4.5M |  718 |    47.39K |     39 |   7 |     6 |   120 |
| Q0L0X80P000   |   80.0 |  0.921 |      4994 |  4.3M | 1154 |    90.18K |     79 |  12 |  14.3 |   230 |
| Q0L0X80P001   |   80.0 |  0.921 |      4893 | 4.32M | 1162 |    86.23K |     79 |  12 |  14.3 |   230 |
| Q0L0X80P002   |   80.0 |  0.922 |      4867 | 4.32M | 1172 |    84.83K |     79 |  12 |  14.3 |   230 |
| Q25L60X40P000 |   40.0 |  0.978 |     16832 | 4.53M |  424 |    38.99K |     40 |   7 |   6.3 |   122 |
| Q25L60X40P001 |   40.0 |  0.978 |     15595 | 4.51M |  436 |    41.91K |     40 |   7 |   6.3 |   122 |
| Q25L60X40P002 |   40.0 |  0.979 |     18062 | 4.54M |  401 |    37.37K |     39 |   7 |     6 |   120 |
| Q25L60X80P000 |   80.0 |  0.966 |     10782 | 4.49M |  656 |    45.54K |     79 |  12 |  14.3 |   230 |
| Q25L60X80P001 |   80.0 |  0.964 |     10634 | 4.48M |  638 |    45.79K |     79 |  12 |  14.3 |   230 |
| Q25L60X80P002 |   80.0 |  0.967 |     10595 |  4.5M |  653 |    50.13K |     79 |  12 |  14.3 |   230 |
| Q30L60X40P000 |   40.0 |  0.985 |     31973 | 4.53M |  242 |     31.1K |     40 |   8 |   5.3 |   128 |
| Q30L60X40P001 |   40.0 |  0.986 |     30840 | 4.53M |  264 |    38.34K |     40 |   7 |   6.3 |   122 |
| Q30L60X40P002 |   40.0 |  0.985 |     28735 | 4.52M |  247 |    34.87K |     40 |   7 |   6.3 |   122 |
| Q30L60X80P000 |   80.0 |  0.985 |     31253 | 4.53M |  246 |    28.71K |     79 |  14 |  12.3 |   242 |
| Q30L60X80P001 |   80.0 |  0.985 |     29025 | 4.55M |  254 |    31.86K |     79 |  13 |  13.3 |   236 |

Table: statUnitigsSuperreads.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |   0.97 |     44458 | 4.52M | 175 |     28.7K |     39 |   6 |     7 |   114 |
| MRX40P001 |   40.0 |  0.969 |     43277 | 4.51M | 188 |     31.4K |     39 |   6 |     7 |   114 |
| MRX40P002 |   40.0 |   0.97 |     46393 | 4.51M | 177 |    30.54K |     39 |   6 |     7 |   114 |
| MRX80P000 |   80.0 |  0.963 |     29043 |  4.5M | 271 |     48.3K |     79 |  10 |  16.3 |   218 |
| MRX80P001 |   80.0 |  0.964 |     31261 | 4.49M | 251 |    45.37K |     79 |  10 |  16.3 |   218 |
| MRX80P002 |   80.0 |  0.964 |     30789 |  4.5M | 255 |    44.16K |     79 |  10 |  16.3 |   218 |

Table: statMRUnitigsSuperreads.md

| Name                          | Mapped | N50Anchor | Sum |      # | SumOthers | median | MAD | lower | upper |
|-------------------------------|-------:|----------:|----:|-------:|----------:|-------:|----:|------:|------:|
| 7_merge_anchors               |  97548 |     4.54M | 101 |      0 |
| 7_merge_mr_unitigs_bcalm      |  78592 |     4.52M | 101 |  5.96K |
| 7_merge_mr_unitigs_bifrost    |  78601 |     4.53M | 107 |  7.09K |
| 7_merge_mr_unitigs_superreads |  67346 |     4.52M | 111 |  6.81K |
| 7_merge_unitigs_bcalm         |  63414 |     4.52M | 125 | 61.14K |
| 7_merge_unitigs_bifrost       |  78639 |     4.52M | 111 | 58.34K |
| 7_merge_unitigs_superreads    |  63418 |     4.54M | 124 | 21.15K |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|--------------|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| 8_spades     |  0.986 |     95432 | 4.54M |  96 |    15.58K |    288 |  35 |    61 |   786 |
| 8_mr_spades  |  0.984 |    148467 | 4.55M |  71 |    22.66K |    375 |  36 |    89 |   966 |
| 8_megahit    |  0.982 |     65378 | 4.53M | 124 |    18.97K |    288 |  35 |    61 |   786 |
| 8_mr_megahit |  0.989 |    132785 | 4.57M |  71 |    20.45K |    375 |  36 |    89 |   966 |

Table: statOtherAnchors.md

| Name                     |     N50 |     Sum |   # |
|--------------------------|--------:|--------:|----:|
| Genome                   | 4641652 |   4.64M |   1 |
| Paralogs                 |    2003 | 260.35K | 131 |
| repetitive               |    1235 |  91.99K | 169 |
| 7_merge_anchors.anchors  |   97548 |   4.54M | 101 |
| glue_anchors             |  105670 |   4.53M |  91 |
| fill_anchors             |  110256 |   4.53M |  85 |
| spades.contig            |  125607 |   4.58M | 140 |
| spades.scaffold          |  132608 |   4.58M | 136 |
| spades.non-contained     |  125607 |   4.56M |  75 |
| mr_spades.contig         |  148607 |   4.59M | 146 |
| mr_spades.scaffold       |  148607 |   4.59M | 144 |
| mr_spades.non-contained  |  148607 |   4.57M |  75 |
| megahit.contig           |   89901 |   4.57M | 148 |
| megahit.non-contained    |   89901 |   4.55M | 105 |
| mr_megahit.contig        |  132896 |   4.61M | 120 |
| mr_megahit.non-contained |  132896 |   4.59M |  69 |

Table: statFinal

Table: statBusco run_bacteria_odb10

| NAME                |   C |   S | D | F | M | Total |
|:--------------------|----:|----:|--:|--:|--:|------:|
| Genome              | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_superreads    | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_bcalm         | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_tadpole       | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_mr_superreads | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_mr_bcalm      | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_mr_tadpole    | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_anchors       | 124 | 124 | 0 | 0 | 0 |   124 |
| glue_anchors        | 124 | 124 | 0 | 0 | 0 |   124 |
| fill_anchors        | 124 | 124 | 0 | 0 | 0 |   124 |
| spades              | 124 | 124 | 0 | 0 | 0 |   124 |
| mr_spades           | 124 | 124 | 0 | 0 | 0 |   124 |
| megahit             | 124 | 124 | 0 | 0 | 0 |   124 |
| mr_megahit          | 124 | 124 | 0 | 0 | 0 |   124 |

Table: statBusco run_enterobacterales_odb10

| NAME                |   C |   S | D | F | M | Total |
|:--------------------|----:|----:|--:|--:|--:|------:|
| Genome              | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_superreads    | 437 | 435 | 2 | 3 | 0 |   440 |
| merge_bcalm         | 437 | 435 | 2 | 3 | 0 |   440 |
| merge_tadpole       | 437 | 435 | 2 | 3 | 0 |   440 |
| merge_mr_superreads | 438 | 436 | 2 | 2 | 0 |   440 |
| merge_mr_bcalm      | 438 | 436 | 2 | 2 | 0 |   440 |
| merge_mr_tadpole    | 438 | 436 | 2 | 2 | 0 |   440 |
| merge_anchors       | 438 | 436 | 2 | 2 | 0 |   440 |
| glue_anchors        | 439 | 437 | 2 | 1 | 0 |   440 |
| fill_anchors        | 439 | 437 | 2 | 1 | 0 |   440 |
| spades              | 440 | 438 | 2 | 0 | 0 |   440 |
| mr_spades           | 440 | 438 | 2 | 0 | 0 |   440 |
| megahit             | 440 | 438 | 2 | 0 | 0 |   440 |
| mr_megahit          | 440 | 438 | 2 | 0 | 0 |   440 |

## *E. coli* str. K-12 substr. DH5alpha

### dh5alpha: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/dh5alpha/1_genome
cd ~/data/anchr/dh5alpha/1_genome

cp ~/data/anchr/ref/dh5alpha/genome.fa .
cp ~/data/anchr/ref/dh5alpha/paralogs.fa .

```

### dh5alpha: download

```shell
cd ~/data/anchr/dh5alpha

mkdir -p ena
cd ena

cat << EOF > source.csv
SRP251726,dh5alpha,HiSeq 2500 PE125
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

rgr md ena_info.tsv --fmt

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name     | srx        | platform        | layout | ilength | srr         | spot    | base  |
|:---------|:-----------|:----------------|:-------|:--------|:------------|:--------|:------|
| dh5alpha | SRX7856678 | ILLUMINA        | PAIRED |         | SRR11245239 | 5881654 | 1.37G |
| dh5alpha | SRX7856679 | OXFORD_NANOPORE | SINGLE |         | SRR11245238 | 346489  | 3.35G |

* Illumina

```shell
cd ~/data/anchr/dh5alpha

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR11245239_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR11245239_2.fastq.gz R2.fq.gz

```

### dh5alpha: template

* template

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4583637 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "bcalm bifrost superreads" \
    --statp 2 \
    --readl 125 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend \
    \
    --busco

```

### dh5alpha: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge 9_markdown/statReads.md
# rm -fr 4_down_sampling 6_down_sampling

bash 0_script/1_repetitive.sh

# BASE_NAME=dh5alpha bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
# bkill -J "${BASE_NAME}-*"

# bash 0_master.sh
# bash 0_cleanup.sh

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
| repetitive |    1200 |  115795 |      213 |
| Illumina.R |     125 |   1.47G | 11763308 |
| trim.R     |     125 |   1.37G | 10962262 |
| Q25L60     |     125 |   1.25G | 10280938 |
| Q30L60     |     125 |   1.13G |  9405535 |

Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 125 |   1.37G | 10970530 |
| highpass | 125 |   1.37G | 10966136 |
| trim     | 125 |   1.37G | 10962262 |
| filter   | 125 |   1.37G | 10962262 |
| R1       | 125 |    683M |  5481131 |
| R2       | 125 | 683.63M |  5481131 |
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
#unique_kmers	35779038
#error_kmers	31288377
#genomic_kmers	4490661
#main_peak	221
#genome_size_in_peaks	9162708
#genome_size	10213708
#haploid_genome_size	5106854
#fold_coverage	96
#haploid_fold_coverage	221
#ploidy	2
#het_rate	0.00012
#percent_repeat_in_peaks	2.524
#percent_repeat	55.129
#start	center	stop	max	volume
```

Table: statMergeReads

| Name          | N50 |     Sum |        # |
|:--------------|----:|--------:|---------:|
| clumped       | 125 |   1.37G | 10959444 |
| ecco          | 125 |   1.37G | 10959444 |
| eccc          | 125 |   1.37G | 10959444 |
| ecct          | 125 |   1.37G | 10952912 |
| extended      | 165 |    1.8G | 10952912 |
| merged.raw    | 343 |    1.1G |  3511622 |
| unmerged.raw  | 165 | 645.89M |  3929668 |
| unmerged.trim | 165 | 645.89M |  3929668 |
| M1            | 343 |   1.06G |  3404325 |
| U1            | 165 | 322.85M |  1964834 |
| U2            | 165 | 323.04M |  1964834 |
| Us            |   0 |       0 |        0 |
| M.cor         | 250 |   1.71G | 10738318 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 180.2 |    180 |  32.0 |         24.10% |
| M.ihist.merge.txt  | 312.5 |    310 |  87.7 |         64.12% |

Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 298.2 |  262.2 |   12.07% | "87" | 4.58M |  4.6M |     1.00 | 0:02'44'' |
| Q25L60.R | 273.7 |  253.6 |    7.35% | "87" | 4.58M | 4.53M |     0.99 | 0:02'37'' |
| Q30L60.R | 246.5 |  232.4 |    5.71% | "87" | 4.58M | 4.52M |     0.99 | 0:02'14'' |

Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.18% |     31096 | 4.45M | 262 |        48 | 29.67K |  643 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| Q0L0X40P001   |   40.0 |  97.13% |     27144 | 4.45M | 294 |        38 | 27.09K |  695 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:38 |
| Q0L0X40P002   |   40.0 |  97.15% |     26026 | 4.46M | 264 |        33 |  20.4K |  635 |   39.0 | 6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:29 |
| Q0L0X80P000   |   80.0 |  95.95% |     12946 | 4.42M | 527 |        35 | 44.34K | 1119 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:33 |
| Q0L0X80P001   |   80.0 |  95.85% |     12796 | 4.42M | 546 |        33 |  41.4K | 1165 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:32 |
| Q0L0X80P002   |   80.0 |  95.91% |     13279 | 4.42M | 519 |        33 | 40.39K | 1109 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:33 |
| Q25L60X40P000 |   40.0 |  97.85% |     43277 | 4.46M | 175 |        39 | 18.54K |  490 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:33 |
| Q25L60X40P001 |   40.0 |  97.67% |     42111 | 4.46M | 192 |        38 | 19.84K |  496 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| Q25L60X40P002 |   40.0 |  97.70% |     40217 | 4.46M | 195 |        35 | 18.86K |  540 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:32 |
| Q25L60X80P000 |   80.0 |  97.17% |     27143 | 4.45M | 296 |        32 | 22.85K |  674 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:32 |
| Q25L60X80P001 |   80.0 |  97.03% |     21632 | 4.45M | 329 |        35 | 28.96K |  740 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:33 |
| Q25L60X80P002 |   80.0 |  96.92% |     26330 | 4.45M | 307 |        32 | 24.85K |  683 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:31 |
| Q30L60X40P000 |   40.0 |  97.78% |     46178 | 4.45M | 163 |        36 | 16.43K |  458 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| Q30L60X40P001 |   40.0 |  97.79% |     47349 | 4.47M | 180 |        49 | 21.01K |  504 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:33 |
| Q30L60X40P002 |   40.0 |  97.83% |     44914 | 4.46M | 171 |     24692 | 41.58K |  506 |   39.0 | 6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:31 |
| Q30L60X80P000 |   80.0 |  97.29% |     27104 | 4.45M | 280 |        32 | 23.52K |  640 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:32 |
| Q30L60X80P001 |   80.0 |  97.29% |     27205 | 4.45M | 275 |        34 | 25.26K |  642 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:32 |

Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.63% |     67275 | 4.45M | 122 |        81 | 20.95K | 331 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:34 |
| MRX40P001 |   40.0 |  97.59% |     67328 | 4.45M | 117 |        91 | 23.48K | 323 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:35 |
| MRX40P002 |   40.0 |  97.58% |     67196 | 4.46M | 118 |        76 | 18.09K | 332 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:36 |
| MRX80P000 |   80.0 |  97.48% |     64119 | 4.46M | 121 |        72 | 16.97K | 324 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:38 |
| MRX80P001 |   80.0 |  97.59% |     67250 | 4.46M | 113 |        77 | 17.27K | 325 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:38 |
| MRX80P002 |   80.0 |  97.56% |     64176 | 4.46M | 121 |        84 | 19.81K | 335 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:40 |

Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.29% |     57013 | 4.46M | 141 |        48 | 32.68K | 697 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:39 |
| Q0L0X40P001   |   40.0 |  98.25% |     54641 | 4.46M | 149 |        46 | 34.29K | 710 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:40 |
| Q0L0X40P002   |   40.0 |  98.24% |     58771 | 4.46M | 131 |        91 | 37.62K | 683 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:46 |
| Q0L0X80P000   |   80.0 |  97.00% |     35763 | 4.46M | 220 |        43 | 22.87K | 471 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:34 |
| Q0L0X80P001   |   80.0 |  96.98% |     41931 | 4.46M | 199 |        39 | 18.45K | 414 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:33 |   0:00:33 |
| Q0L0X80P002   |   80.0 |  97.05% |     40711 | 4.45M | 195 |        59 |  22.8K | 427 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:28 |   0:00:32 |
| Q25L60X40P000 |   40.0 |  98.40% |     58790 | 4.46M | 132 |        42 | 28.14K | 692 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:42 |
| Q25L60X40P001 |   40.0 |  98.45% |     63674 | 4.46M | 129 |     39382 | 69.62K | 711 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:40 |
| Q25L60X40P002 |   40.0 |  98.50% |     63682 | 4.46M | 134 |        35 | 24.75K | 735 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:41 |
| Q25L60X80P000 |   80.0 |  97.40% |     51384 | 4.46M | 157 |        40 | 15.07K | 338 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:33 |
| Q25L60X80P001 |   80.0 |  97.25% |     43280 | 4.46M | 171 |        50 | 18.74K | 379 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:33 |
| Q25L60X80P002 |   80.0 |  97.24% |     52376 | 4.36M | 154 |        48 | 16.56K | 338 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:25 |   0:00:33 |
| Q30L60X40P000 |   40.0 |  98.43% |     60844 | 4.46M | 132 |        43 | 30.05K | 713 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:40 |
| Q30L60X40P001 |   40.0 |  98.44% |     60874 | 4.47M | 136 |        47 | 34.16K | 775 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:43 |
| Q30L60X40P002 |   40.0 |  98.47% |     60871 | 4.46M | 132 |        39 |  28.7K | 773 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:40 |
| Q30L60X80P000 |   80.0 |  97.34% |     54822 | 4.46M | 149 |        66 | 16.58K | 327 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:32 |
| Q30L60X80P001 |   80.0 |  97.44% |     53170 | 4.46M | 157 |        43 | 16.11K | 358 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:33 |

Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.16% |     73681 | 4.45M | 118 |       601 | 23.46K | 259 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:30 |
| MRX40P001 |   40.0 |  97.16% |     67315 | 4.45M | 116 |       198 | 21.29K | 254 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:32 |
| MRX40P002 |   40.0 |  97.09% |     67196 | 4.46M | 118 |        87 | 15.22K | 247 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:31 |
| MRX80P000 |   80.0 |  96.85% |     67364 | 4.28M | 112 |        74 | 11.51K | 228 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:01:41 |   0:00:36 |
| MRX80P001 |   80.0 |  96.84% |     73702 | 4.46M | 111 |        81 | 12.61K | 226 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:33 |
| MRX80P002 |   80.0 |  96.84% |     66225 | 4.46M | 118 |        80 | 13.06K | 236 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:34 |

Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.12% |     64138 | 4.46M | 120 |      1017 | 26.83K | 410 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:34 |
| Q0L0X40P001   |   40.0 |  98.09% |     60860 | 4.45M | 123 |       358 |  25.5K | 404 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:36 |
| Q0L0X40P002   |   40.0 |  98.03% |     67329 | 4.46M | 115 |        46 | 17.04K | 384 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:35 |
| Q0L0X80P000   |   80.0 |  97.16% |     60844 | 4.46M | 133 |        61 | 15.94K | 287 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:33 |
| Q0L0X80P001   |   80.0 |  97.18% |     57835 | 4.46M | 133 |        42 | 13.36K | 282 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q0L0X80P002   |   80.0 |  97.25% |     64115 | 4.46M | 125 |       384 | 16.32K | 268 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:31 |
| Q25L60X40P000 |   40.0 |  98.32% |     63680 | 4.46M | 121 |        64 | 22.67K | 468 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q25L60X40P001 |   40.0 |  98.31% |     66049 | 4.46M | 118 |        93 | 25.49K | 443 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:37 |
| Q25L60X40P002 |   40.0 |  98.26% |     63682 | 4.46M | 122 |        41 | 19.75K | 473 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |
| Q25L60X80P000 |   80.0 |  97.34% |     63593 | 4.46M | 123 |        40 | 10.78K | 245 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:32 |
| Q25L60X80P001 |   80.0 |  97.30% |     58770 | 4.46M | 130 |        57 | 14.49K | 284 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:34 |
| Q25L60X80P002 |   80.0 |  97.36% |     63641 | 4.46M | 119 |       541 | 16.04K | 255 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:34 |
| Q30L60X40P000 |   40.0 |  98.25% |     63654 | 4.46M | 123 |        45 | 19.77K | 481 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:38 |
| Q30L60X40P001 |   40.0 |  98.31% |     64116 | 4.46M | 119 |      1343 | 41.65K | 496 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:39 |
| Q30L60X40P002 |   40.0 |  98.36% |     63678 | 4.46M | 122 |        46 | 21.18K | 498 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:37 |
| Q30L60X80P000 |   80.0 |  97.36% |     67334 | 4.46M | 118 |        50 | 11.82K | 247 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:33 |
| Q30L60X80P001 |   80.0 |  97.46% |     63673 | 4.46M | 122 |        51 | 12.96K | 271 |   80.0 | 9.0 |  17.7 | 214.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:33 |

Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.00% |     73666 | 4.45M | 117 |        92 | 15.38K | 230 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| MRX40P001 |   40.0 |  96.98% |     67328 | 4.45M | 115 |       446 | 22.53K | 240 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:30 |
| MRX40P002 |   40.0 |  97.04% |     67196 | 4.46M | 117 |       103 | 16.34K | 236 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:30 |
| MRX80P000 |   80.0 |  96.87% |     73678 | 4.46M | 115 |        86 | 14.15K | 234 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:31 |
| MRX80P001 |   80.0 |  96.86% |     73689 | 4.46M | 111 |       103 | 13.52K | 229 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:33 |
| MRX80P002 |   80.0 |  96.85% |     66221 | 4.46M | 116 |        81 | 12.82K | 234 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:32 |

Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.99% |     67358 | 4.46M | 112 |     24692 | 203.78K | 42 |  262.0 | 22.0 |  65.3 | 656.0 |   0:00:57 |
| 7_merge_mr_unitigs_bcalm      |  98.74% |     67362 | 4.46M | 112 |      1345 |  17.74K | 12 |  263.0 | 22.0 |  65.7 | 658.0 |   0:01:17 |
| 7_merge_mr_unitigs_superreads |  98.75% |     67361 | 4.46M | 111 |      3330 |   14.7K |  8 |  263.0 | 23.0 |  64.7 | 664.0 |   0:01:23 |
| 7_merge_mr_unitigs_tadpole    |  98.85% |     73698 | 4.46M | 110 |      1377 |  16.47K | 10 |  260.0 | 22.0 |  64.7 | 652.0 |   0:01:26 |
| 7_merge_unitigs_bcalm         |  98.93% |     67344 | 4.46M | 111 |     39382 | 154.36K | 40 |  263.0 | 22.0 |  65.7 | 658.0 |   0:01:35 |
| 7_merge_unitigs_superreads    |  98.98% |     67353 | 4.46M | 113 |     23775 | 126.43K | 32 |  262.0 | 22.0 |  65.3 | 656.0 |   0:01:43 |
| 7_merge_unitigs_tadpole       |  98.88% |     67350 | 4.46M | 111 |     19926 |  77.76K | 26 |  262.0 | 22.0 |  65.3 | 656.0 |   0:01:36 |

Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.45% |    114635 | 4.48M |  78 |      1048 | 18.49K | 157 |  263.0 | 22.0 |  65.7 | 658.0 |   0:00:41 |
| 8_mr_spades  |  98.93% |    132590 | 4.22M |  64 |       682 |  18.9K | 123 |  376.0 | 28.0 |  97.3 | 920.0 |   0:00:46 |
| 8_megahit    |  98.47% |     73665 | 4.46M | 116 |      1031 | 25.17K | 231 |  263.0 | 22.0 |  65.7 | 658.0 |   0:00:41 |
| 8_mr_megahit |  99.22% |    132754 |  4.5M |  75 |       673 | 22.14K | 147 |  376.0 | 28.0 |  97.3 | 920.0 |   0:00:46 |

Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 4583637 | 4583637 |   1 |
| Paralogs                 |    1737 |  188158 | 111 |
| repetitive               |    1200 |  115795 | 213 |
| 7_merge_anchors.anchors  |   67358 | 4463149 | 112 |
| 7_merge_anchors.others   |   24692 |  203783 |  42 |
| glue_anchors             |   85589 | 4462785 | 104 |
| fill_anchors             |  102355 | 4462986 |  91 |
| spades.contig            |  125538 | 4514902 | 165 |
| spades.scaffold          |  143522 | 4515652 | 155 |
| spades.non-contained     |  132337 | 4493720 |  79 |
| mr_spades.contig         |  155808 | 4522785 |  87 |
| mr_spades.scaffold       |  203812 | 4522995 |  84 |
| mr_spades.non-contained  |  155808 | 4512983 |  60 |
| megahit.contig           |   82955 | 4508030 | 164 |
| megahit.non-contained    |   82955 | 4488852 | 115 |
| mr_megahit.contig        |  133730 | 4557477 | 156 |
| mr_megahit.non-contained |  133730 | 4522636 |  72 |

Table: statBusco run_bacteria_odb10

| NAME                |   C |   S | D | F | M | Total |
|:--------------------|----:|----:|--:|--:|--:|------:|
| Genome              | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_superreads    | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_bcalm         | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_tadpole       | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_mr_superreads | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_mr_bcalm      | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_mr_tadpole    | 124 | 124 | 0 | 0 | 0 |   124 |
| merge_anchors       | 124 | 124 | 0 | 0 | 0 |   124 |
| glue_anchors        | 124 | 124 | 0 | 0 | 0 |   124 |
| fill_anchors        | 124 | 124 | 0 | 0 | 0 |   124 |
| spades              | 124 | 124 | 0 | 0 | 0 |   124 |
| mr_spades           | 124 | 124 | 0 | 0 | 0 |   124 |
| megahit             | 124 | 124 | 0 | 0 | 0 |   124 |
| mr_megahit          | 124 | 124 | 0 | 0 | 0 |   124 |

Table: statBusco run_enterobacterales_odb10

| NAME                |   C |   S | D | F | M | Total |
|:--------------------|----:|----:|--:|--:|--:|------:|
| Genome              | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_superreads    | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_bcalm         | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_tadpole       | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_mr_superreads | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_mr_bcalm      | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_mr_tadpole    | 440 | 438 | 2 | 0 | 0 |   440 |
| merge_anchors       | 440 | 438 | 2 | 0 | 0 |   440 |
| glue_anchors        | 440 | 438 | 2 | 0 | 0 |   440 |
| fill_anchors        | 440 | 438 | 2 | 0 | 0 |   440 |
| spades              | 440 | 438 | 2 | 0 | 0 |   440 |
| mr_spades           | 440 | 438 | 2 | 0 | 0 |   440 |
| megahit             | 440 | 438 | 2 | 0 | 0 |   440 |
| mr_megahit          | 440 | 438 | 2 | 0 | 0 |   440 |
