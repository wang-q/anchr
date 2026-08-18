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
  * [*Bacillus cereus* ATCC 10987](#bacillus-cereus-atcc-10987)
    * [bcer: reference](#bcer-reference)
    * [bcer: download](#bcer-download)
    * [bcer: template](#bcer-template)
    * [bcer: run](#bcer-run)
  * [*Rhodobacter sphaeroides* 2.4.1](#rhodobacter-sphaeroides-241)
    * [rsph: reference](#rsph-reference)
    * [rsph: download](#rsph-download)
    * [rsph: template](#rsph-template)
    * [rsph: run](#rsph-run)
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

anchr ena meta source.csv > ena_info.json
anchr ena manifest ena_info.json

tva to md ena_info.tsv --fmt

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
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig bcalm" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

### g37: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=g37

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

| Group    |  Mean | Median | STDev | Pairs% | Orientation |
| -------- | ----: | -----: | ----: | -----: | ----------- |
| R.genome | 466.5 |    451 | 127.8 | 97.74% | FR          |
| R.contig | 462.1 |    447 | 126.1 | 92.63% | FR          |

Table: statInsertSize

| K    | property              |        min |        max |
| ---- | --------------------- | ---------: | ---------: |
| R.21 | Homozygous (a)        |            |       100% |
|      | Genome Haploid Length |            | 577,874 bp |
|      | Genome Repeat Length  |   3,639 bp |   3,645 bp |
|      | Genome Unique Length  | 573,825 bp | 574,639 bp |
|      | Model Fit             |   92.5827% |   93.3366% |
|      | Read Error Rate       |            |  0.137189% |
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

| chr       | chrLength |  size | coverage |
| --------- | --------: | ----: | -------: |
| NC_000908 |    580076 | 16875 |   0.0291 |
| all       |    580076 | 16875 |   0.0291 |

Table: statRepetitive

| Name       |    N50 |     Sum |      # |
| ---------- | -----: | ------: | -----: |
| genome     | 580076 | 580.08K |      1 |
| paralogs   |   1567 |  11.53K |      8 |
| repetitive |    499 |  16.88K |     53 |
| Illumina.R |    150 |  102.1M | 680644 |
| trim.R     |    150 | 101.44M | 676644 |
| Q0L0       |    150 | 101.44M | 676644 |
| Q25L60     |    150 |  98.23M | 657747 |
| Q30L60     |    150 |  95.03M | 638160 |

Table: statReads

| Name     |  N50 |     Sum |      # |
| -------- | ---: | ------: | -----: |
| clumpify |  150 | 102.01M | 680076 |
| highpass |  150 | 101.54M | 676954 |
| trim     |  150 | 101.44M | 676644 |
| filter   |  150 | 101.44M | 676644 |
| R1       |  150 |  50.74M | 338322 |
| R2       |  150 |   50.7M | 338322 |
| Rs       |    0 |       0 |      0 |

Table: statTrimReads

```text
#R.trim
#Matched	589	0.08701%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

| Name          |  N50 |     Sum |      # |
| ------------- | ---: | ------: | -----: |
| clumped       |  150 | 101.44M | 676644 |
| ecco          |  150 | 101.44M | 676644 |
| ecct          |  150 |  97.96M | 653378 |
| extended      |  190 | 123.74M | 653378 |
| merged.raw    |  442 |  77.19M | 181762 |
| unmerged.raw  |  190 |  54.74M | 289854 |
| unmerged.trim |  190 |   54.7M | 289808 |
| M1            |  442 |  77.19M | 181762 |
| U1            |  190 |  27.46M | 144904 |
| U2            |  190 |  27.24M | 144904 |
| Us            |    0 |       0 |      0 |
| M.cor         |  362 | 132.08M | 653332 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
| ------------------ | ----: | -----: | ----: | -----: |
| M.ihist.merge1.txt | 258.3 |    264 |  24.7 |  4.23% |
| M.ihist.merge.txt  | 424.6 |    430 |  64.8 | 55.68% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% |   RealG | RunTime |
| -------- | ----: | -----: | -------: | ------: | ------: |
| Q0L0.R   | 174.9 |  160.8 |    8.07% | 580.08K | 0:00:08 |
| Q25L60.R | 169.4 |  160.5 |    5.26% | 580.08K | 0:00:08 |
| Q30L60.R | 163.9 |  157.3 |    4.01% | 580.08K | 0:00:08 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.968 |     34511 | 562.96K |   25 |      1208 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P001   |   40.0 |  0.968 |     37560 |  563.1K |   28 |      1121 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P002   |   40.0 |  0.966 |     34521 | 562.41K |   30 |      1214 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X80P000   |   80.0 |  0.953 |     15170 | 559.05K |   63 |      1264 |     79 |    9 |  17.3 | 212.0 |
| Q0L0X80P001   |   80.0 |  0.954 |     16579 | 562.11K |   65 |      1292 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X40P000 |   40.0 |  0.966 |     32042 | 562.55K |   32 |      1483 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P001 |   40.0 |  0.968 |     33729 | 562.68K |   27 |      1166 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P002 |   40.0 |  0.969 |     39442 | 563.27K |   24 |      1268 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X80P000 |   80.0 |  0.948 |     12537 | 560.87K |   71 |      1589 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X80P001 |   80.0 |  0.955 |     11089 | 561.83K |   76 |      1496 |     79 |    9 |  17.3 | 212.0 |
| Q30L60X40P000 |   40.0 |  0.968 |     37702 | 562.62K |   24 |       989 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P001 |   40.0 |  0.966 |     31414 | 562.74K |   29 |      1324 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P002 |   40.0 |  0.967 |     34603 | 562.49K |   28 |      1159 |     40 |    6 |   7.3 | 116.0 |
| Q30L60X80P000 |   80.0 |  0.957 |     20428 | 560.95K |   60 |      1314 |     79 |    9 |  17.3 | 212.0 |

Table: statUnitigsMultik.md

| Name          | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.963 |     17769 | 562.36K |   48 |      2193 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P001   |   40.0 |  0.965 |     18533 | 562.77K |   45 |      1761 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P002   |   40.0 |  0.962 |     16790 | 561.29K |   47 |      2101 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X80P000   |   80.0 |  0.937 |      7767 | 553.93K |  100 |      2024 |     79 |    9 |  17.3 | 212.0 |
| Q0L0X80P001   |   80.0 |  0.938 |      7787 | 554.33K |   98 |      2044 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X40P000 |   40.0 |  0.962 |     15038 | 561.94K |   53 |      2110 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P001 |   40.0 |  0.962 |     19198 | 560.71K |   44 |      1887 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P002 |   40.0 |  0.965 |     20810 | 562.46K |   45 |      2053 |     40 |    6 |   7.3 | 116.0 |
| Q25L60X80P000 |   80.0 |  0.938 |      7626 | 554.17K |   98 |      1989 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X80P001 |   80.0 |  0.938 |      6947 | 555.97K |  107 |      2040 |     79 |   10 |  16.3 | 218.0 |
| Q30L60X40P000 |   40.0 |  0.964 |     22529 | 562.11K |   42 |      1852 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P001 |   40.0 |  0.963 |     20410 | 561.75K |   49 |      2101 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P002 |   40.0 |  0.962 |     17563 |  561.6K |   49 |      1980 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X80P000 |   80.0 |  0.939 |      7859 | 554.62K |   97 |      2038 |     79 |    9 |  17.3 | 212.0 |

Table: statUnitigsUnitig.md

| Name          | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.963 |     17769 | 562.36K |   48 |      2193 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P001   |   40.0 |  0.965 |     18533 | 562.77K |   45 |      1761 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P002   |   40.0 |  0.962 |     16790 | 561.29K |   47 |      2101 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X80P000   |   80.0 |  0.937 |      7767 | 553.93K |  100 |      2024 |     79 |    9 |  17.3 | 212.0 |
| Q0L0X80P001   |   80.0 |  0.938 |      7787 | 554.33K |   98 |      2044 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X40P000 |   40.0 |  0.962 |     15038 | 561.94K |   53 |      2110 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P001 |   40.0 |  0.962 |     19198 | 560.71K |   44 |      1887 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P002 |   40.0 |  0.965 |     20810 | 562.46K |   45 |      2053 |     40 |    6 |   7.3 | 116.0 |
| Q25L60X80P000 |   80.0 |  0.938 |      7626 | 554.17K |   98 |      1989 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X80P001 |   80.0 |  0.938 |      6947 | 555.97K |  107 |      2040 |     79 |   10 |  16.3 | 218.0 |
| Q30L60X40P000 |   40.0 |  0.964 |     22529 | 562.11K |   42 |      1852 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P001 |   40.0 |  0.963 |     20410 | 561.75K |   49 |      2101 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P002 |   40.0 |  0.962 |     17563 |  561.6K |   49 |      1980 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X80P000 |   80.0 |  0.939 |      7859 | 554.62K |   97 |      2038 |     79 |    9 |  17.3 | 212.0 |

Table: statUnitigsBcalm.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.685 |     31166 | 561.83K |   34 |      2651 |     39 |    5 |   8.0 | 108.0 |
| MRX40P001 |   40.0 |  0.685 |     30478 | 560.63K |   32 |      2298 |     39 |    5 |   8.0 | 108.0 |
| MRX40P002 |   40.0 |  0.688 |     36795 | 562.71K |   33 |      2318 |     39 |    5 |   8.0 | 108.0 |
| MRX80P000 |   80.0 |  0.669 |      8764 | 558.49K |   82 |      3061 |     78 |    9 |  17.0 | 210.0 |
| MRX80P001 |   80.0 |  0.675 |     12905 | 560.18K |   68 |      2190 |     78 |   10 |  16.0 | 216.0 |

Table: statMRUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.681 |     19129 | 561.27K |   48 |      3264 |     39 |    5 |   8.0 | 108.0 |
| MRX40P001 |   40.0 |  0.682 |     22261 | 560.86K |   45 |      3025 |     39 |    5 |   8.0 | 108.0 |
| MRX40P002 |   40.0 |  0.686 |     28917 | 561.25K |   37 |      2672 |     39 |    5 |   8.0 | 108.0 |
| MRX80P000 |   80.0 |  0.665 |      8113 | 556.39K |   87 |      3151 |     78 |    9 |  17.0 | 210.0 |
| MRX80P001 |   80.0 |  0.670 |     10739 | 556.75K |   77 |      2316 |     78 |   10 |  16.0 | 216.0 |

Table: statMRUnitigsUnitig.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.681 |     19129 | 561.27K |   48 |      3264 |     39 |    5 |   8.0 | 108.0 |
| MRX40P001 |   40.0 |  0.682 |     22261 | 560.86K |   45 |      3025 |     39 |    5 |   8.0 | 108.0 |
| MRX40P002 |   40.0 |  0.686 |     28917 | 561.25K |   37 |      2672 |     39 |    5 |   8.0 | 108.0 |
| MRX80P000 |   80.0 |  0.665 |      8113 | 556.39K |   87 |      3151 |     78 |    9 |  17.0 | 210.0 |
| MRX80P001 |   80.0 |  0.670 |     10739 | 556.75K |   77 |      2316 |     78 |   10 |  16.0 | 216.0 |

Table: statMRUnitigsBcalm.md

| Name                      | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.690 |     55098 | 563.79K |   15 |      1082 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_mr_unitigs_bcalm  |  0.689 |     48832 | 573.07K |   21 |      1272 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_mr_unitigs_multik |  0.690 |     55049 | 565.52K |   17 |      1252 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_mr_unitigs_unitig |  0.689 |     48832 | 573.07K |   21 |      1272 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_unitigs_bcalm     |  0.690 |     48853 | 564.67K |   19 |      1409 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_unitigs_multik    |  0.690 |     55098 | 563.71K |   15 |      1161 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_unitigs_unitig    |  0.690 |     48853 | 564.67K |   19 |      1409 |    226 |   22 |  53.3 | 584.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.988 |     83431 | 573.62K |    9 |      7382 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_spades  |  0.710 |    580031 | 580.03K |    1 |       601 |    225 |   22 |  53.0 | 582.0 |
| 8_megahit    |  0.977 |     39288 | 564.03K |   24 |     15271 |    159 |   17 |  36.0 | 420.0 |
| 8_mr_megahit |  0.710 |    319000 | 579.88K |    9 |     19141 |    223 |   23 |  51.3 | 584.0 |

Table: statOtherAnchors.md

| Name                    |    N50 |     Sum |    # |
| ----------------------- | -----: | ------: | ---: |
| Genome                  | 580076 | 580.08K |    1 |
| Paralogs                |   1567 |  11.53K |    8 |
| repetitive              |    499 |  16.88K |   53 |
| 7_merge_anchors.anchors |  55098 | 563.79K |   15 |
| spades.contig           | 163847 |    581K |   38 |
| spades.scaffold         | 163847 |    581K |   38 |
| mr_spades.contig        | 580506 | 580.63K |    2 |
| mr_spades.scaffold      | 580506 | 580.63K |    2 |
| megahit.contig          |  39730 |  579.3K |   56 |
| mr_megahit.contig       | 319186 | 599.03K |   44 |

Table: statFinal

* Assembly quality by QUAST

| Assembly         | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| ---------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik     |        15 |  179712 | 563707 |  55098 |     0 | 96.997 | 1.000 |   0.00 |   76.44 |      19.20 |
| merge_mr_multik  |        17 |  114442 | 565518 |  55049 |     0 | 96.963 | 1.003 |   0.00 |   77.46 |      19.14 |
| merge_unitig     |        19 |  107330 | 564672 |  48853 |     0 | 96.970 | 1.002 |   0.00 |   77.37 |      19.17 |
| merge_mr_unitig  |        21 |   79910 | 573070 |  48832 |     0 | 96.917 | 1.018 |   0.00 |   76.22 |      18.53 |
| merge_bcalm      |        19 |  107330 | 564672 |  48853 |     0 | 96.970 | 1.002 |   0.00 |   77.37 |      19.17 |
| merge_mr_bcalm   |        21 |   79910 | 573070 |  48832 |     0 | 96.917 | 1.018 |   0.00 |   76.22 |      18.53 |
| merge_anchors    |        15 |  179712 | 563791 |  55098 |     0 | 97.004 | 1.000 |   0.00 |   76.78 |      19.55 |
| spades           |        35 |  236302 | 580752 | 163847 |     1 | 99.056 | 1.002 |   0.00 |  222.74 |      66.49 |
| mr_spades        |         2 |  580506 | 580632 | 580506 |     0 | 100.000 | 1.001 |   0.00 |  300.77 |      91.64 |
| megahit          |        56 |   89673 | 579301 |  39730 |     4 | 97.507 | 1.004 |   0.00 |   85.23 |      18.31 |
| mr_megahit       |        44 |  319186 | 599025 | 319186 |     3 | 99.893 | 1.025 |   0.00 |  311.14 |      88.39 |

Table: statQuast

### g37: unitig / bcalm / multik 三链对比（2026-08-16）

同一批 G37 reads 上，唯一变量是 unitigger：`asm multik`（每主 K
31..81 构建骨架、更大 k 验证、跨主 K 合并）vs `asm unitig`（自研 BCALM
语义，每 k 独立 unitigs）vs 外部 bcalm（K = 31..81），下游
anchor / OLC 合并参数完全相同。

#### 最终合并（QUAST，merge 链）

| 指标        | multik | unitig | bcalm |
| :---------- | -----: | -----: | -----: |
| # contigs   |  **15** |    19 |    19 |
| N50         | **55,098** | 48,853 | 48,853 |
| Largest     | **179,712** | 107,330 | 107,330 |
| # misassemblies | 0 | 0 | 0 |
| GF%         | **96.997** | 96.970 | 96.970 |
| mm/100 kbp  |   76.44 |  77.37 |  77.37 |

#### 结论

* **multik 全面占优**：N50 比 unitig/bcalm 高 13%（55.1K vs 48.9K），
  contigs 更少（15 vs 19），三链均 0 mis——与 MG1655 的结论一致
  （N50 +42%，见 `notes/benchmarks/mg1655-unitig-bcalm-multik.md`）；
* **unitig 与 bcalm 等价**：N50 / contigs / GF 完全一致（`asm unitig`
  是 bcalm 的 Rust 复刻，自研可替代外部依赖）；
* **merge 链（MR）同样如此**：multik 55.0K vs unitig/bcalm 48.8K；
* **mm/100k 三链一致（76-77）**：来自 reads 与参考 NC_000908 的系统性
  差异，与组装链无关；spades 系（222-300）明显更高。

### g37: 新链（2026-08-16 追加：K128 + 气泡合并 + extend --min-len 1000 + min-contig-len 200）

模板 multik 分支自 2026-08-16 起自动包含：主 K `31..121 128 160`（6_ MR
链；4_ 链 `31..91`，bcalm 分支 `31..121`）、multik 气泡合并（默认
`--merge-similar 0.95 --merge-len 20`）、olc `--min-contig-len 200`（仅
multik 分支，bcalm/unitig 分支保持 1000）、`asm extend --min-len 1000`
（短于 1000 bp 的 contig 不延伸，避免重复区嵌合）。详细机制与中间实验见
`notes/benchmarks/g37-megahit-spades.md` §7/§8/§10。

7 组 MR（MRX40P000-004 + MRX80P000-001）QUAST（新链 = 31..128+160 +
气泡合并 + extend --min-len 1000 + min200）：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_mr_multik（旧，K31..81） | 17 | 114442 | 565518 |  55049 | 0 | 96.963 | 1.003 | 0.00 | 77.46 | 19.14 |
| merge_mr_multik（新链） | 13 | 236494 | 583630 | 121532 | 0 | 98.849 | 1.000 | 0.00 | 226.03 | 66.37 |
| mr_spades（旧运行） |  2 | 580506 | 580632 | 580506 | 0 | 100.000 | 1.001 | 0.00 | 300.77 | 91.64 |
| mr_megahit（旧运行） | 44 | 319186 | 599025 | 319186 | 3 | 99.893 | 1.025 | 0.00 | 311.14 | 88.39 |

要点：
* 新链 N50 55.0K→**121.5K**（+121%）、GF 96.963→**98.849%**、**0 mis
  保持**，并跨复制起点接出 236.5K 最大 contig；K160 长 unitig 层提升
  N50，配合 indel 容忍的近重复合并（banded identity）Dup 1.000；
* mm/100k 上升（77→226）是 QUAST 对参考的口径问题：新覆盖集中在低复杂
  区，reads 与 NC_000908 在那里差异 ~5%，consensus 忠实反映 reads
  （reads-vs-contigs 一致率 99.96%）；bwa 口径下真缺口仅 ~1.6K bp；
* QUAST minimap 会漏对齐低复杂度区短 contig，报告的未覆盖 bp 偏大。

### g37: multik 单次调用（--all-masters）全流程门禁（2026-08-17 追加）

`asm multik --all-masters`（auto 阶梯：N50 408 → 31..192）单次调用，
7 组 MR anchor 合并口径，guide 与否结果完全相同：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 单次调用 auto（guide 与无 guide 一致） | 58 | 179800 | 580756 | 55170 | 0 | 97.943 | 1.003 | 0.00 | 108.81 | 28.78 |
| 08-16 旧链基线（merge_multik） | 15 | 179712 | 563707 | 55098 | 0 | 96.997 | 1.000 | 0.00 | 76.44 | 19.20 |

要点：
* **0 mis 保持**，N50 与旧 multik 链持平（55.2K vs 55.1K），GF +0.95；
* guide 与无 guide 输出完全一致，与 MG1655 结论互相印证；
* 本流程为快速门禁复现（手写 anchor 合并脚本），contigs 数多于模板
  7_merge 链属预期；单组验证 N50 121K（auto 31..192）。

### g37: multik 性能优化门禁（2026-08-17 晚追加）

walk 滚动窗口 + succ 索引（删除 400 MB HashMap index）、
`remove_unsupported` 按 unitig 并行、pass0 与 rounds 并发（详见
`notes/design/asm-multik.md` §性能）。multik 输出与优化前字节级一致，
全链复跑：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | -------: | ---------: |
| 性能优化后（auto 31..192） | 58 | 179800 | 581098 | 55170 | 0 | 97.979 | 107.51 | 28.59 |
| 08-17 门禁记录（同口径） | 58 | 179800 | 580756 | 55170 | 0 | 97.943 | 108.81 | 28.78 |

* **0 mis 保持**，逐指标与 08-17 门禁一致（GF +0.036 源于其间合入的
  `remove_unsupported` run>=2 过剪修复，非本次性能改动）；
* 单组 multik 4.2–6.5 s（-p 8、2 组并发），峰值 ~2.0 GB/进程。

### mg1655: multik vs bcalm 对照（2026-08-15）

同一批 `6_down_sampling` reads（MRX40P000/P001/P002 + MRX80P000/P001，
5 组），唯一变量是 unitigger：`asm multik` vs `asm unitig`（自研 BCALM
语义）vs 外部 bcalm（每 k 31..81 独立 unitigs + `asm olc --unitigs` 跨 k
合并），下游 anchor/OLC 合并参数完全相同。

#### unitig / anchor 阶段（N50 / 条数 / Sum）

| 组        | multik unitigs | unitig unitigs | bcalm unitigs |
| --------- | -------------: | -------------: | ------------: |
| MRX40P000 |  21238 / 1455  |  61235 / 131   |  56477 / 142  |
| MRX40P001 |  19872 / 1455  |  61235 / 129   |  54908 / 147  |
| MRX40P002 |  21238 / 1455  |  59716 / 130   |  57961 / 143  |
| MRX80P000 |  19331 / 1525  |  53834 / 158   |  42412 / 194  |
| MRX80P001 |  20068 / 1525  |  50795 / 163   |  42234 / 195  |

multik 的 unitigs 比 bcalm/unitig 短约 2.5-3 倍、条数多约 9-11 倍；
`asm unitig`（自研）与 bcalm 等价且略优；两链 anchor Sum ≈ 4.53M，比
multik 的 4.47M 更接近基因组。

#### 最终合并（`asm olc --unitigs`，5 组合并）与 QUAST

| 指标      | unitig 链（自研） | bcalm 链（现代） | multik 链（现代，同输入） | legacy bcalm | legacy merge_anchors |
| --------- | ----------------: | ---------------: | ------------------------: | ------------: | --------------------: |
| # contigs |               102 |              108 |                       317 |           101 |                   103 |
| N50       |           105,719 |           95,478 |                    23,403 |        78,596 |                95,484 |
| Largest   |          246,019 |          202,937 |                    73,030 |       174,107 |               204,605 |
| # mis     |                 1 |                1 |                         4 |             0 |                     0 |
| GF%       |            97.77 |            97.77 |                     96.44 |         97.26 |                 97.67 |
| Dup       |           1.000 |           1.001 |                      1.003 |         1.000 |                  1.000 |
| mm/100k   |            0.27 |            0.04 |                      0.18 |          0.00 |                  0.40 |
| indel/100k|            0.29 |            0.23 |                      0.02 |          0.02 |                  0.13 |

结论：N50 差距基本全部来自 unitig 阶段——`asm unitig` 链（105.7K）与
bcalm 链（95.5K）都追平并超过 legacy（95.5K），multik 只有 23.4K；mis
从 multik 的 4 降到 unitig/bcalm 链的 1（legacy 0）；`asm unitig` 与外部
bcalm 等价（自研可替代外部依赖）。

#### multik k0 修复（2026-08-15 追加）

multik 碎片化的根因是 pass 0 骨架冻结在 k0（21/31）：迭代轮只验证/合并
k0 的图，从不重新组装（单跑 `asm unitig -k 81` 有 53.5K/705 条，multik
迭代却输出 21K/1455 条）。`auto_ks` k0 从 `N50/10` 改为 `N50/3`
（clamp 31..51）后，MG1655 5 组 k0=51 全链：unitig N50 46-60K（原
19-21K）、merge N50 **65.8K**（原 23.4K，128 contigs）、GF 97.36%
（原 96.44）、Dup 1.000；mis 仍 4（重复区/环状 quast 误报 + merge 阶段
重复序列 exact-overlap 错连，机制分析见 `todo.md`，未修）。

#### Dup 修复（2026-08-15 追加）

原两链 Dup 1.07-1.08 偏高：anchor 阶段无近似重复 contig，重复由
`asm olc --unitigs` 跨组合并产生（跨组 anchors 边界不一致，exact overlap
检测连不上，残留同区域不同边界的 contig 对）。`asm olc` 现在在 consensus
后增加近似 overlap 合并（`consensus::merge_overlapping_contigs`：31-mer
定位主导 offset + 头部锚定 + 重叠区 ≥99% 一致才合并，嵌合 contig 的
多块对齐被拒绝）：unitig 链 Dup **1.079 → 1.000**、bcalm 链 **1.068 →
1.001**，GF 97.77% 不变，unitig 链 N50 105.7K → 110.2K（93 contigs），
bcalm 链 N50 95.5K → 88.1K（102 contigs）。mis 仍为 1（contig_26
relocation，嵌合修复属另一任务）。

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

* template

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 4641652 \
    --parallel 8 \
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
    --merge \
    \
    --bwa "Q25L60" \
    \
    --cov "40 80" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

### mg1655: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

| Group    |  Mean | Median | STDev | Pairs% | Orientation |
| -------- | ----: | -----: | ----: | -----: | ----------- |
| R.genome | 298.1 |    298 |  17.7 | 99.40% | FR          |
| R.contig | 297.8 |    297 |  17.7 | 97.75% | FR          |

Table: statInsertSize

| K    | property              |          min |          max |
| ---- | --------------------- | -----------: | -----------: |
| R.21 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,477,162 bp |
|      | Genome Repeat Length  |   135,854 bp |   135,942 bp |
|      | Genome Unique Length  | 4,339,855 bp | 4,342,674 bp |
|      | Model Fit             |     95.3903% |     97.3934% |
|      | Read Error Rate       |              |    0.531924% |
|      | Kmer Cov              |              |        299.7 |
| R.51 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,384,614 bp |
|      | Genome Repeat Length  |    90,815 bp |    90,901 bp |
|      | Genome Unique Length  | 4,291,722 bp | 4,295,793 bp |
|      | Model Fit             |     95.9309% |     97.6265% |
|      | Read Error Rate       |              |    0.326151% |
|      | Kmer Cov              |              |        223.4 |
| R.81 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,303,893 bp |
|      | Genome Repeat Length  |    59,727 bp |    59,810 bp |
|      | Genome Unique Length  | 4,241,157 bp | 4,247,095 bp |
|      | Model Fit             |     95.9005% |     97.4791% |
|      | Read Error Rate       |              |    0.272518% |
|      | Kmer Cov              |              |        151.5 |

Table: statFastK

| chr       | chrLength |   size | coverage |
| --------- | --------: | -----: | -------: |
| NC_000913 |   4641652 | 124422 |   0.0268 |
| all       |   4641652 | 124422 |   0.0268 |

Table: statRepetitive

| Name       |     N50 |     Sum |        # |
| ---------- | ------: | ------: | -------: |
| genome     | 4641652 |   4.64M |        1 |
| paralogs   |    1737 |    193K |      112 |
| repetitive |    1265 | 124.42K |      155 |
| Illumina.R |     151 |   1.73G | 11458940 |
| trim.R     |     149 |   1.46G | 10636412 |
| Q0L0       |     149 |   1.46G | 10636412 |
| Q25L60     |     148 |   1.35G | 10185038 |
| Q30L60     |     128 |   1.12G |  9506942 |

Table: statReads

| Name     |  N50 |     Sum |        # |
| -------- | ---: | ------: | -------: |
| clumpify |  151 |   1.73G | 11439000 |
| highpass |  151 |    1.7G | 11272376 |
| trim     |  149 |   1.46G | 10636412 |
| filter   |  149 |   1.46G | 10636412 |
| R1       |  150 | 753.89M |  5318206 |
| R2       |  144 | 706.08M |  5318206 |
| Rs       |    0 |       0 |        0 |

Table: statTrimReads

```text
#R.trim
#Matched	17107	0.15176%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

| Name          |  N50 |    Sum |        # |
| ------------- | ---: | -----: | -------: |
| clumped       |  149 |  1.46G | 10635230 |
| ecco          |  149 |  1.46G | 10635230 |
| ecct          |  149 |  1.45G | 10566976 |
| extended      |  189 |  1.87G | 10566976 |
| merged.raw    |  339 |  1.76G |  5210725 |
| unmerged.raw  |  172 | 21.12M |   145526 |
| unmerged.trim |  170 | 20.22M |   140094 |
| M1            |  339 |  1.76G |  5210571 |
| U1            |  175 | 10.62M |    70047 |
| U2            |  165 |   9.6M |    70047 |
| Us            |    0 |      0 |        0 |
| M.cor         |  338 |  1.78G | 10561236 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
| ------------------ | ----: | -----: | ----: | -----: |
| M.ihist.merge1.txt | 271.3 |    277 |  24.3 | 10.45% |
| M.ihist.merge.txt  | 337.6 |    338 |  19.3 | 98.77% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real | RunTime |
| -------- | ----: | -----: | -------: | ---: | ----: | ----: | -------: | ------: |
| Q0L0.R   | 314.5 |  291.7 |    7.27% | "24" | 4.64M | 4.47M |     0.96 | 0:02:46 |
| Q25L60.R | 289.9 |  278.4 |    3.95% | "24" | 4.64M | 4.42M |     0.95 | 0:02:31 |
| Q30L60.R | 242.2 |  237.1 |    2.12% | "24" | 4.64M | 4.37M |     0.94 | 0:02:05 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.964 |     13614 | 4.48M |  545 |    125438 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P001   |   40.0 |  0.964 |     13663 | 4.48M |  546 |    127750 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P002   |   40.0 |  0.965 |     14305 | 4.48M |  511 |    119182 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X80P000   |   80.0 |  0.959 |     10077 | 4.48M |  761 |    152824 |     79 |   12 |  14.3 | 230.0 |
| Q0L0X80P001   |   80.0 |  0.959 |      9538 | 4.48M |  756 |    149797 |     79 |   12 |  14.3 | 230.0 |
| Q0L0X80P002   |   80.0 |  0.959 |      9004 | 4.47M |  769 |    156952 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X40P000 |   40.0 |  0.969 |     14007 | 4.47M |  513 |    115498 |     39 |    7 |   6.0 | 120.0 |
| Q25L60X40P001 |   40.0 |  0.970 |     14179 | 4.48M |  508 |    116620 |     39 |    7 |   6.0 | 120.0 |
| Q25L60X40P002 |   40.0 |  0.969 |     13648 | 4.47M |  513 |    118120 |     39 |    7 |   6.0 | 120.0 |
| Q25L60X80P000 |   80.0 |  0.968 |     12815 | 4.48M |  582 |    118105 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X80P001 |   80.0 |  0.968 |     12039 | 4.48M |  596 |    122840 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X80P002 |   80.0 |  0.969 |     12705 | 4.48M |  583 |    121945 |     79 |   12 |  14.3 | 230.0 |
| Q30L60X40P000 |   40.0 |  0.967 |     11397 | 4.44M |  603 |    114763 |     39 |    8 |   5.0 | 126.0 |
| Q30L60X40P001 |   40.0 |  0.968 |     11784 | 4.45M |  590 |    113624 |     39 |    8 |   5.0 | 126.0 |
| Q30L60X40P002 |   40.0 |  0.966 |     11228 | 4.44M |  598 |    114871 |     39 |    8 |   5.0 | 126.0 |
| Q30L60X80P000 |   80.0 |  0.972 |     13557 | 4.48M |  530 |    106155 |     79 |   14 |  12.3 | 242.0 |
| Q30L60X80P001 |   80.0 |  0.972 |     13097 | 4.48M |  529 |    107238 |     79 |   14 |  12.3 | 242.0 |

Table: statUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.481 |     21714 | 4.47M |  351 |    123635 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.481 |     22510 | 4.47M |  352 |    123687 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.481 |     21169 | 4.47M |  356 |    122936 |     39 |    6 |   7.0 | 114.0 |
| MRX80P000 |   80.0 |  0.480 |     20006 | 4.48M |  373 |    115009 |     78 |   10 |  16.0 | 216.0 |
| MRX80P001 |   80.0 |  0.480 |     19961 | 4.48M |  375 |    112524 |     78 |   10 |  16.0 | 216.0 |
| MRX80P002 |   80.0 |  0.480 |     20779 | 4.48M |  373 |    113540 |     78 |   10 |  16.0 | 216.0 |

Table: statMRUnitigsMultik.md

| Name                      | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.479 |     28019 | 4.61M |  278 |      5841 |    381 |   38 |  89.0 | 990.0 |
| 7_merge_mr_unitigs_multik |  0.477 |     23594 | 4.53M |  317 |      4097 |    380 |   38 |  88.7 | 988.0 |
| 7_merge_unitigs_multik    |  0.479 |     23956 | 4.73M |  318 |      6618 |    380 |   38 |  88.7 | 988.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.989 |     82215 | 4.55M |  121 |     23296 |    292 |   35 |  62.3 | 794.0 |
| 8_mr_spades  |  0.500 |    145104 | 4.56M |   84 |     25914 |    382 |   37 |  90.3 | 986.0 |
| 8_megahit    |  0.985 |     35210 | 4.54M |  210 |     46173 |    291 |   36 |  61.0 | 798.0 |
| 8_mr_megahit |  0.502 |    126306 | 4.58M |   87 |     31206 |    382 |   37 |  90.3 | 986.0 |

Table: statOtherAnchors.md

| Name                    |     N50 |     Sum |    # |
| ----------------------- | ------: | ------: | ---: |
| Genome                  | 4641652 |   4.64M |    1 |
| Paralogs                |    1737 |    193K |  112 |
| repetitive              |    1265 | 124.42K |  155 |
| 7_merge_anchors.anchors |   28019 |   4.61M |  278 |
| spades.contig           |  125607 |   4.57M |  135 |
| spades.scaffold         |  132608 |   4.57M |  131 |
| mr_spades.contig        |  148607 |   4.59M |  152 |
| mr_spades.scaffold      |  148607 |   4.59M |  149 |
| megahit.contig          |   43891 |   4.59M |  273 |
| mr_megahit.contig       |  126312 |   4.61M |  138 |

Table: statFinal

* Assembly quality by QUAST

| Assembly       | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik   |       318 |   81358 | 4726538 |  23956 |    12 | 96.78 | 1.052 |   0.00 |    0.63 |       0.19 |
| merge_mr_multik |      317 |   73037 | 4527602 |  23594 |     5 | 96.49 | 1.011 |   0.00 |    0.24 |       0.02 |
| merge_anchors  |       278 |   83824 | 4606970 |  28019 |    12 | 96.84 | 1.025 |   0.00 |    0.63 |       0.20 |
| spades         |       125 |  224028 | 4572988 | 125607 |     0 | 98.47 | 1.000 |   0.00 |    0.74 |       0.17 |
| mr_spades      |       152 |  284843 | 4588125 | 148607 |     0 | 98.64 | 1.002 |   0.00 |    0.92 |       0.17 |
| megahit        |       273 |  175838 | 4588105 |  43891 |    92 | 98.26 | 1.004 |   0.00 |    3.73 |       0.50 |
| mr_megahit     |       138 |  311797 | 4611104 | 126312 |     2 | 98.89 | 1.004 |   0.00 |    2.13 |       0.26 |

Table: statQuast

### mg1655: 新链处理与对比（2026-08-16 追加）

同一批 5 组 MR reads（`6_down_sampling/MRX40P000/P001/P002` +
`MRX80P000/P001`），当前 multik 链（K31..121+128+160 + 气泡合并 + extend
`--min-len 1000` + `--min-contig-len 200`）QUAST：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 新 multik 链（5 组） | 90 | 268333 | 4584494 | 123988 | 0 | 98.523 | 1.002 | 0.00 | 0.22 | 0.09 |
| 旧 multik 链（同 5 组，K31..81） | 107 | — | — | 95478 | 0 | 97.61 | ~1.00 | — | — | — |
| merge_mr_multik（旧运行，全组） | 317 | 73037 | 4527602 | 23594 | 5 | 96.49 | 1.011 | 0.00 | 0.24 | 0.02 |
| spades（旧运行） | 125 | 224028 | 4572988 | 125607 | 0 | 98.47 | 1.000 | 0.00 | 0.74 | 0.17 |
| mr_spades（旧运行） | 152 | 284843 | 4588125 | 148607 | 0 | 98.64 | 1.002 | 0.00 | 0.92 | 0.17 |
| megahit（旧运行） | 273 | 175838 | 4588105 | 43891 | 92 | 98.26 | 1.004 | 0.00 | 3.73 | 0.50 |
| mr_megahit（旧运行） | 138 | 311797 | 4611104 | 126312 | 2 | 98.89 | 1.004 | 0.00 | 2.13 | 0.26 |

要点（详细分析见 `notes/benchmarks/mg1655-process-compare.md`）：
* 新链 N50 95.5K→**124.0K**（+30%）、GF 97.61→**98.523%**、**0 mis**；
  同 5 组输入口径下已超过 megahit（82.8K / 1 mis）并与 spades 持平，
  仅 mr_spades（148.6K，全量输入）更高；
* **extend 必须 `--min-len 1000`**：extend 短碎片会把 1.2 Mb 处重复元件
  拷贝接成嵌合体（238 bp 碎片被长成 1,238 bp relocation，3 mis）；加门槛
  后 0 mis（G37 同步验证）；
* 处理过程中顺带修复 `asm olc` consensus 的 O(n²×L) 性能瓶颈（种子索引
  预筛，输出逐位一致），MG1655 单组 9 主池从数小时降到 3 分钟。

### mg1655: multik 单次调用（--all-masters）门禁（2026-08-17 追加）

`asm multik --all-masters` 单次调用取代模板 per-master 循环（k-major
顺序共享每 k 计数表、rounds 跨 master 并行、计数表即建即弃）。同 5 组
anchor 合并口径，三个变体与 08-16 基线对比：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 08-16 基线（per-master 链） |  90 | 268333 | 4584494 | 123988 | 0 | 98.523 | 1.002 | 0.00 | 0.22 | 0.09 |
| 单次调用 guide 31..192 | 181 | 174306 | 4603907 |  79617 | 0 | 98.880 | 1.003 | 0.00 | 0.83 | 0.07 |
| 单次调用 无guide 31..192 | 182 | 174306 | 4604139 |  79625 | 0 | 98.877 | 1.003 | 0.00 | 0.83 | 0.07 |
| **单次调用 auto（模板默认）** | 178 | 174306 | 4603590 |  79617 | 0 | 98.852 | 1.003 | 0.00 | 1.80 | 0.07 |

要点：
* **三变体全部 0 mis**；guide 与否逐指标几乎相同（GF 差 0.003）→ 5 组
  anchor 投票下 guide 无贡献，模板移除 `--use-guide`（速度 ~2.1×）；
* **单组对照无回归**：新旧同 k（31..128）单组质量完全一致（均 2 mis、
  GF 99.40、N50 112,514），重构本身质量中性；
* **N50 79.6K vs 08-16 链 124.0K**：~~归因为验证密度~~（08-17 下午证伪，
  实为跨组 olc 缺 `--unitigs` 的测量错误，见下节）；
* 计时（release、-p 8、单组）：auto 140–172 s（无 guide）vs guide
  ~360 s vs 旧 per-master 串行 7:23（-p 4、9 次全量计数）；
* 内存：计数表即建即弃后峰值 6.6 GB（40×组）/ 10.8 GB（80×组）；
  曾因缓存全部 K 张表达 26.6 GB 并在 5 组并发时 OOM，已改；
* auto 阶梯改为固定梯 `31..192` 截断于 `clamp(N50/2, 81, 192)`
  （本数据 N50 339 → 31..160）。旧公式 0.8×N50 给出 51..251，高 k
  master 被残余错误打碎（N50 9.4K、5 mis），不可用。

### mg1655: N50 差异归因终局——跨组 olc 缺 `--unitigs`（2026-08-17 下午追加）

上表 79.6K 与 run>=2 修复后的 96.4K，均出自 /tmp 实验脚本（run.sh /
fix_chain.sh）：跨组 anchors 合并用了 `asm olc --list-files`（无
`--unitigs`），anchors 被**当作 reads 重新走 S0 多 k 组装**（切碎重组），
故 N50 塌缩。正式模板 `7_merge_anchors.tera.sh` 一直带 `--unitigs`，
流程本身无此问题。补 `--unitigs` 重跑跨组（multik 输出与 anchors 不变）：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 08-16 基线（per-master 链） |  90 | 268333 | 4584494 | 123988 | 0 | 98.523 | 1.002 | 0.00 | 0.22 | 0.09 |
| **单次调用 auto（修复跨组后）** |  90 | 268842 | 4624071 | 118731 | 0 | **99.072** | 1.005 | 0.00 | 2.06 | 0.06 |

要点：
* 单次调用 multik 达到基线水平：contigs 数持平（90）、N50 118.7K
  （-4.2%）、GF **+0.55 pp**（99.072）、Largest 268,842 持平、**0 mis**；
* **验证密度中性**：step 1（每 k 验证）vs step 3（every-third，复刻
  08-16 链密度）全链输出仅一处 166 kb 重复区序列不同，QUAST 逐指标
  相同——上节"归因为验证密度"证伪；
* **guide / last-k cut 均中性**：`--use-guide`（233 s/组）与移除 last-k
  cut 的全链与默认配置逐指标相同（96,447 → 修复跨组后 118,731 同）；
  5 组 anchor 投票抹平 multik 层的这些差异。代码恢复原状（cut 保留、
  实验 `--validate-step` 开关移除）；
* 单组 anchor 对照：08-16 链 95 条 / N50 112,781 vs 新链 96 条 /
  N50 112,781（相同）——单组层完全等价，差异全部在跨组测量错误；
* 复现：`/tmp/mg1655_fix/`（fix 后 multik）+ `olc --unitigs --list-files
  anchors.list` + extend + quast（quast_fix_u）；变体对照在
  `/tmp/mg1655_vstep`、`/tmp/mg1655_vguide`、`/tmp/mg1655_nocut`。

### mg1655: multik 性能优化门禁（2026-08-17 晚追加）

`--all-masters` 单次调用的三个性能改动（walk 滚动 fw/rc 窗口 + 分类期
预联唯一后继下标、`remove_unsupported` 按 unitig 并行、pass0(k) 与
earlier-masters rounds(k) rayon::join 并发），multik 输出与优化前
**字节级一致**，全链复跑（跨组 olc 带 `--unitigs`，同上节口径）：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | -------: | ---------: |
| 性能优化后（auto 31..160） | 90 | 268842 | 4624071 | 118731 | 0 | 99.072 | 2.06 | 0.06 |
| 上节基线（优化前同口径） | 90 | 268842 | 4624071 | 118731 | 0 | 99.072 | 2.06 | 0.06 |

* 逐指标完全相同（multik 字节级一致 → 下游全同），**0 mis 保持**；
* 单组计时（release、-p 8、单进程）：105.8 s → **48.9 s**（2.2×）；
  5 组链 2 并发下每组 59–65 s；峰值内存不变（~11.7 GB/进程）；
* 分解：DFA walk 7.3→1.9 s（k=160 隔离）、classify 1.5→1.1 s（删
  HashMap 构建）、46 轮 `remove_unsupported` CPU 104.9 s 转入并行。

### mg1655/g37: chain-per-master 调度修复与内存受控（2026-08-17 深夜追加）

上节 chain-per-master 重构引入两个缺陷，本节修复并完成门禁：

1. **builder 越界 panic**：全部表按序到齐时 `buffered[next]` 越界
   （next==len），所有组 multik 直接失败——即此前"门禁链缺 QUAST 报告"
   的根因（multik 失败 → 下游全部静默跳过）；
2. **并发建表无界**：10 张 reads 表同时构建，单进程峰值 15.9 GB。
   改为 `BUILD_WINDOW=3` 的有界前瞻（构建+待投递 ≤3 张表，仍与
   chains 的 round 工作重叠），峰值降到 **12.8 GB**（MG1655）/
   2.1–2.6 GB（G37）；
3. 非顶部 chain 不再保留 last_table（仅 top master 的 `cut` 需要）；
4. 单 k 路径（`--kmer K` 无后续验证轮）`cut` 复用 pass0 表，k=160
   单 k 调用 16.3→11.6 s；
5. 计时开关 `ANCHR_MULTIK_TIMING` 覆盖 pass0/cut/finalize/probe/
   builder（此前只有 round）。

门禁（修复后全链复跑，跨组 olc `--unitigs` 正式口径）：

| Dataset | Assembly | # contigs | N50 | # mis | GF% | 对照 |
| ------- | -------- | --------: | ---: | ----: | --: | ---- |
| MG1655  | 5 组 anchor 合并 | 90 | 118731 | **0** | 99.072 | 与 08-17 基线逐指标相同 |
| G37     | 7 组 anchor 合并（--unitigs） | — | 134724 | **1** | 98.774 | 见下方归因 |
| G37     | 7 组快速口径（无 --unitigs） | 58 | 55170 | **0** | 97.979 | 与 08-17 晚门禁逐指标相同 |

* MG1655 multik 输出与修复前字节级一致；G37 P000/P002 multik 输出与
  昨晚提交版（452ec01）二进制字节级一致——本节全部改动质量中性；
* `-p` 扫描：输出在 p8/p16/p24 下字节级一致（确定性）；但 all-masters
  单组 p16/p24 反而更慢（25.9/24.7 s vs p8 22.8 s）——10 chain 并发下
  `remove_unsupported` 随机访问在 p16 劣化 3×（内存带宽瓶颈），
  **单组 p8 即最优**，32 核靠多组并发利用；
* 计时：MG1655 单组 multik 22.8 s（-p 8；原 105.8 s），聚合分解
  rounds 63 s + 表构建 ~31 s + pass0 22 s + finalize 13 s ≈ 130 s 池
  工作量（8 线程饱和，利用率 ~76%）。

#### g37 `--unitigs` 口径 1 mis 归因（既有问题，非本节引入）

* 嵌合位于 **MRX40P002 单组** `olc --unitigs` 输出 contig_2
  （134,487 bp；relocation 52938↔134693，跨 ~82 kb 重复区），经
  anchor 进入 7 组合并被保留；其余 6 组 anchor 均 0 mis；
* P002 multik 输出与昨晚二进制字节一致 → 嵌合自 08-17 单次调用链
  （auto 31..192）即存在。此前未暴露：快速口径（无 `--unitigs` 的
  S0 重组装）会切碎重组 anchors 掩盖单组嵌合；08-16 per-master 模板
  链（K≤160）记录为 0 mis，是否同样产生该嵌合未验证；
* 遗留：单组嵌合的 master-k 归因与修复（如 anchor 层跨组投票对该
  位点的处理）留待专项实验，按实验纪律需 A/B 验证。

### g37/mg1655: `olc --cross-validate` 跨组嵌合投票（2026-08-18 追加）

承接上节归因。深挖结论修正了初步判断：

* 归因到 master：逐链复现（`--kmer K..160` 子阶梯）显示连接在
  **k≥101 的全部链**存在（m101/m121/m128/m160），非 m160 独有；
  m128 链 quast 0 mis 只是 junction 恰在 contig 开头 146 bp，未达
  extensive misassembly 报告标准；
* junction 桥接 reads 计数 **9**（正常位置 16–21）→ 不是低支持噪声，
  而是**真实存在的低丰度菌株 A-B 相邻结构**（菌株结构变异）。
  单组内 count 阈值不可切、也不应切；唯一正确压制层是跨组一致性；
* 跨组 `olc --unitigs` 保留嵌合的机制：`filter_contained` 把
  "被嵌合 contig 前缀整条包含"的**其他组正确 contigs** 删掉，
  嵌合反而存活。

修复：`olc --unitigs --cross-validate`（默认关；模板
`7_merge_anchors.tera.sh` 已开启）。判定：一条 contig 的两端各自被
≥2 个其他文件的 contigs 覆盖（flank=min_overlap），且中部 junction
窗口（span=min_overlap/2）无任何其他文件 contig 横跨 → 删除该
contig（序列由其他组的分开 contigs 完整提供）。横跨判定先按来源
contig 合并对齐区间（精确重叠链在错配处断裂成多段，见 MG1655
236 kb contig 误删修复）。实现 `libs/olc/overlap.rs
drop_cross_chimeras`，4 个单测覆盖删/留/同文件/断裂链场景。

| Dataset | Assembly | # contigs | N50 | # mis | GF% | 对照 |
| ------- | -------- | --------: | ---: | ----: | --: | ---- |
| G37     | 7 组 + cross-validate | 14 | 121664 | **0** | 98.797 | 无 cv：1 mis/134724/98.774；GF 反升 |
| MG1655  | 5 组 + cross-validate | — | 118731 | **0** | 99.072 | 与基线逐指标相同（236 kb 正确 contig 保留） |

* 首版（区间不合并）曾误删 MG1655 anchor.2:anchor_1（236 kb）：
  组0 等价 contig 的精确重叠链在错配处断为 [0,48851)+[48757,236423)
  两段，单段横跨检查失败；按来源合并后恢复，N50 回到 118,731；
* cargo test 421 通过（含 4 个新增 cross 单测）；clippy/fmt 干净；
  模板 `7_merge_anchors.tera.sh` 接入 `--cross-validate`
  （单文件流程 covers 为空、自然无操作）。

#### 验证边界（重要，防过度信任）

* **仅在 mock 数据集（G37/MG1655 人工菌群）验证**，未在真实宏基因组
  验证。两数据集的其他组组装均足够完整、能横跨真连接，真实数据
  （覆盖不均、strain variation、碎片化）不保证此前提；
* 已知误删风险场景（均未实测）：
  1. 其他组只提供**碎片**：真连接两端各被 ≥2 组碎片覆盖、但无人
     横跨中部 → 真连接被误删（假阴性，断真 contig）；`min_groups≥2`
     缓解不消除；
  2. **长重复边界**：连接点处在比其他组 contig 更长的重复里，横跨
     必然失败 → 真连接被删（与基准中 82 kb 重复嵌合同构，机制上
     无法区分）；
  3. **3 组样本** 2:1 票太弱；基准均为 5–7 组；
* 规则设计依据是先验因果（低丰度菌株私有连接 ≠ 跨样本保守连接），
  参数锚定既有值（flank=min_overlap、span=min_overlap/2），无针对
  基准调优的自由量——但上述场景仍是开放风险；
* 后续验证要求：上真实宏基因组前，先做 cv 开/关 A/B（对比 N50/GF/
  contig 数；N50 显著下降即误删信号），确认无误删再固化模板默认。

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

anchr ena meta source.csv > ena_info.json
anchr ena manifest ena_info.json

tva to md ena_info.tsv --fmt

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

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 4583637 \
    --parallel 16 \
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
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

### dh5alpha: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

| Group             |  Mean | Median |  STDev | Pairs%/Orientation |
|-------------------|------:|-------:|-------:|-------------------:|
| R.genome.bbtools  | 470.7 |    346 | 2460.3 |             99.98% |
| R.tadpole.bbtools | 389.0 |    340 |  206.1 |             94.97% |
| R.genome.picard   | 394.8 |    346 |  208.3 |                 FR |
| R.tadpole.picard  | 389.0 |    340 |  205.8 |                 FR |

Table: statInsertSize

| K    | property              |          min |          max |
|------|-----------------------|-------------:|-------------:|
| R.21 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,422,294 bp |
|      | Genome Repeat Length  |   120,434 bp |   120,514 bp |
|      | Genome Unique Length  | 4,300,398 bp | 4,303,243 bp |
|      | Model Fit             |     97.6629% |     98.0363% |
|      | Read Error Rate       |              |    0.209622% |
|      | Kmer Cov              |              |        267.1 |
| R.51 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,441,231 bp |
|      | Genome Repeat Length  |    94,823 bp |    94,884 bp |
|      | Genome Unique Length  | 4,344,991 bp | 4,347,764 bp |
|      | Model Fit             |     97.9127% |      98.422% |
|      | Read Error Rate       |              |    0.131634% |
|      | Kmer Cov              |              |        185.6 |
| R.81 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,466,733 bp |
|      | Genome Repeat Length  |    87,251 bp |    87,306 bp |
|      | Genome Unique Length  | 4,378,088 bp | 4,380,823 bp |
|      | Model Fit             |     97.4061% |      98.803% |
|      | Read Error Rate       |              |    0.109848% |
|      | Kmer Cov              |              |        108.3 |

Table: statFastK

| chr       | chrLength |   size | coverage |
| --------- | --------: | -----: | -------: |
| NZ_CP017100 | 4583637 | 110586 |   0.0241 |
| all       |   4583637 | 110586 |   0.0241 |

Table: statRepetitive

| Name       |     N50 |     Sum |        # |
|------------|--------:|--------:|---------:|
| genome     | 4583637 |   4.58M |        1 |
| paralogs   |    1737 | 188.16K |      111 |
| repetitive |    1175 | 110.59K |      190 |
| Illumina.R |     125 |   1.47G | 11763308 |
| trim.R     |     125 |   1.37G | 10962178 |
| Q0L0       |     125 |   1.37G | 10962178 |
| Q25L60     |     125 |   1.25G | 10280852 |
| Q30L60     |     125 |   1.13G |  9405463 |

Table: statReads

| Name     | N50 |     Sum |        # |
|----------|----:|--------:|---------:|
| clumpify | 125 |   1.37G | 10970448 |
| highpass | 125 |   1.37G | 10966054 |
| trim     | 125 |   1.37G | 10962178 |
| filter   | 125 |   1.37G | 10962178 |
| R1       | 125 | 682.99M |  5481089 |
| R2       | 125 | 683.63M |  5481089 |
| Rs       |   0 |       0 |        0 |

Table: statTrimReads

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

| Name          | N50 |     Sum |        # |
|---------------|----:|--------:|---------:|
| clumped       | 125 |   1.37G | 10959360 |
| ecco          | 125 |   1.37G | 10959360 |
| eccc          | 125 |   1.37G | 10959360 |
| ecct          | 125 |   1.37G | 10952518 |
| extended      | 165 |    1.8G | 10952518 |
| merged.raw    | 343 |    1.1G |  3510900 |
| unmerged.raw  | 165 | 646.07M |  3930718 |
| unmerged.trim | 165 | 646.07M |  3930718 |
| M1            | 343 |   1.06G |  3403638 |
| U1            | 165 | 322.94M |  1965359 |
| U2            | 165 | 323.12M |  1965359 |
| Us            |   0 |       0 |        0 |
| M.cor         | 250 |   1.71G | 10737994 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
|--------------------|------:|-------:|------:|-------:|
| M.ihist.merge1.txt | 172.9 |    173 |  27.8 | 21.04% |
| M.ihist.merge.txt  | 312.5 |    310 |  87.7 | 64.11% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real | RunTime |
|----------|------:|-------:|---------:|-----:|------:|------:|---------:|--------:|
| Q0L0.R   | 298.2 |  262.0 |   12.13% | "87" | 4.58M |  4.6M |     1.00 | 0:02:07 |
| Q25L60.R | 273.7 |  253.4 |    7.43% | "87" | 4.58M | 4.53M |     0.99 | 0:01:57 |
| Q30L60.R | 246.5 |  232.2 |    5.80% | "87" | 4.58M | 4.52M |     0.99 | 0:01:44 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.983 |     57903 | 4.45M | 141 |    37.71K |     40 |   5 |   8.3 |   110 |
| Q0L0X40P001   |   40.0 |  0.982 |     54887 | 4.46M | 146 |    28.86K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P002   |   40.0 |  0.984 |     56551 | 4.46M | 140 |    38.42K |     40 |   5 |   8.3 |   110 |
| Q0L0X80P000   |   80.0 |   0.97 |     37978 | 4.45M | 197 |    21.17K |     80 |   9 |  17.7 |   214 |
| Q0L0X80P001   |   80.0 |  0.971 |     35821 | 4.45M | 211 |    26.49K |     80 |   9 |  17.7 |   214 |
| Q0L0X80P002   |   80.0 |  0.972 |     39975 | 4.45M | 200 |    21.03K |     80 |   9 |  17.7 |   214 |
| Q25L60X40P000 |   40.0 |  0.984 |     63101 | 4.46M | 132 |    28.01K |     40 |   6 |   7.3 |   116 |
| Q25L60X40P001 |   40.0 |  0.984 |     60855 | 4.46M | 128 |    25.15K |     40 |   6 |   7.3 |   116 |
| Q25L60X40P002 |   40.0 |  0.983 |     63658 | 4.45M | 129 |    37.22K |     40 |   5 |   8.3 |   110 |
| Q25L60X80P000 |   80.0 |  0.973 |     54848 | 4.46M | 157 |    16.31K |     80 |   9 |  17.7 |   214 |
| Q25L60X80P001 |   80.0 |  0.972 |     51779 | 4.46M | 157 |    19.03K |     80 |   9 |  17.7 |   214 |
| Q25L60X80P002 |   80.0 |  0.973 |     46203 | 4.44M | 158 |    18.13K |     80 |   9 |  17.7 |   214 |
| Q30L60X40P000 |   40.0 |  0.984 |     63579 | 4.46M | 131 |    26.39K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P001 |   40.0 |  0.985 |     63667 | 4.46M | 131 |    34.87K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P002 |   40.0 |  0.984 |     60836 | 4.46M | 132 |    26.21K |     40 |   6 |   7.3 |   116 |
| Q30L60X80P000 |   80.0 |  0.976 |     54900 | 4.46M | 154 |    22.07K |     80 |   9 |  17.7 |   214 |
| Q30L60X80P001 |   80.0 |  0.975 |     54879 | 4.45M | 156 |    17.76K |     80 |   9 |  17.7 |   214 |

Table: statUnitigsBcalm.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.972 |     73699 | 4.45M | 115 |    23.97K |     40 |   5 |   8.3 |   110 |
| MRX40P001 |   40.0 |  0.972 |     63662 | 4.45M | 119 |     17.4K |     40 |   5 |   8.3 |   110 |
| MRX40P002 |   40.0 |   0.97 |     73648 | 4.45M | 116 |    17.53K |     40 |   5 |   8.3 |   110 |
| MRX80P000 |   80.0 |  0.969 |     67221 | 4.45M | 118 |    16.65K |     80 |   8 |  18.7 |   208 |
| MRX80P001 |   80.0 |  0.968 |     73681 | 4.46M | 115 |    13.48K |     80 |   8 |  18.7 |   208 |
| MRX80P002 |   80.0 |  0.968 |     73658 | 4.45M | 115 |    18.68K |     80 |   8 |  18.7 |   208 |

Table: statMRUnitigsBcalm.md

| Name          | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.983 |     63651 | 4.45M | 121 |    32.58K |     40 |   5 |   8.3 |   110 |
| Q0L0X40P001   |   40.0 |  0.981 |     64143 | 4.46M | 118 |    18.77K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P002   |   40.0 |  0.981 |     64135 | 4.46M | 118 |    20.76K |     40 |   6 |   7.3 |   116 |
| Q0L0X80P000   |   80.0 |  0.973 |     73667 | 4.45M | 113 |    14.92K |     80 |   9 |  17.7 |   214 |
| Q0L0X80P001   |   80.0 |  0.973 |     67322 | 4.46M | 115 |    16.59K |     80 |   9 |  17.7 |   214 |
| Q0L0X80P002   |   80.0 |  0.974 |     67333 | 4.46M | 112 |    15.33K |     80 |   9 |  17.7 |   214 |
| Q25L60X40P000 |   40.0 |  0.982 |     63609 | 4.46M | 121 |    23.16K |     40 |   5 |   8.3 |   110 |
| Q25L60X40P001 |   40.0 |  0.983 |     64147 | 4.46M | 118 |    19.77K |     40 |   6 |   7.3 |   116 |
| Q25L60X40P002 |   40.0 |  0.982 |     64144 | 4.46M | 117 |    20.57K |     40 |   6 |   7.3 |   116 |
| Q25L60X80P000 |   80.0 |  0.974 |     67333 | 4.46M | 112 |    11.92K |     80 |   9 |  17.7 |   214 |
| Q25L60X80P001 |   80.0 |  0.973 |     63670 | 4.46M | 116 |    16.48K |     80 |   9 |  17.7 |   214 |
| Q25L60X80P002 |   80.0 |  0.974 |     67335 | 4.46M | 112 |    15.45K |     80 |   9 |  17.7 |   214 |
| Q30L60X40P000 |   40.0 |  0.985 |     63579 | 4.46M | 121 |    22.51K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P001 |   40.0 |  0.984 |     64139 | 4.46M | 117 |     25.3K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P002 |   40.0 |  0.983 |     67343 | 4.46M | 119 |    19.78K |     40 |   6 |   7.3 |   116 |
| Q30L60X80P000 |   80.0 |  0.977 |     64133 | 4.45M | 113 |    15.08K |     80 |   9 |  17.7 |   214 |
| Q30L60X80P001 |   80.0 |  0.974 |     67327 | 4.45M | 113 |    14.97K |     80 |   9 |  17.7 |   214 |

Table: statUnitigsBifrost.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |   0.97 |     73691 | 4.45M | 116 |    17.84K |     40 |   5 |   8.3 |   110 |
| MRX40P001 |   40.0 |  0.971 |     63653 | 4.46M | 120 |    17.86K |     40 |   5 |   8.3 |   110 |
| MRX40P002 |   40.0 |  0.969 |     73648 | 4.46M | 116 |    16.94K |     40 |   5 |   8.3 |   110 |
| MRX80P000 |   80.0 |  0.968 |     73694 | 4.45M | 117 |    16.15K |     80 |   8 |  18.7 |   208 |
| MRX80P001 |   80.0 |  0.968 |     73669 | 4.46M | 113 |    13.42K |     80 |   8 |  18.7 |   208 |
| MRX80P002 |   80.0 |  0.968 |     66255 | 4.45M | 117 |    17.03K |     80 |   8 |  18.7 |   208 |

Table: statMRUnitigsBifrost.md

| Name          | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|---------------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| Q0L0X40P000   |   40.0 |  0.972 |     28270 | 4.45M | 266 |    27.86K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P001   |   40.0 |  0.972 |     29053 | 4.46M | 283 |    23.31K |     40 |   6 |   7.3 |   116 |
| Q0L0X40P002   |   40.0 |  0.973 |     28264 | 4.46M | 267 |     26.8K |     40 |   6 |   7.3 |   116 |
| Q0L0X80P000   |   80.0 |  0.958 |     13252 | 4.41M | 533 |    46.59K |     79 |   9 |  17.3 |   212 |
| Q0L0X80P001   |   80.0 |   0.96 |     12411 | 4.42M | 552 |    44.08K |     79 |   9 |  17.3 |   212 |
| Q0L0X80P002   |   80.0 |  0.962 |     13206 | 4.43M | 520 |    41.43K |     80 |   9 |  17.7 |   214 |
| Q25L60X40P000 |   40.0 |  0.977 |     42440 | 4.46M | 186 |    19.29K |     40 |   6 |   7.3 |   116 |
| Q25L60X40P001 |   40.0 |  0.978 |     42093 | 4.46M | 187 |    21.82K |     40 |   6 |   7.3 |   116 |
| Q25L60X40P002 |   40.0 |  0.977 |     40745 | 4.46M | 198 |     26.4K |     40 |   6 |   7.3 |   116 |
| Q25L60X80P000 |   80.0 |  0.971 |     25042 | 4.45M | 305 |    23.03K |     80 |   9 |  17.7 |   214 |
| Q25L60X80P001 |   80.0 |   0.97 |     24364 | 4.45M | 317 |    32.85K |     79 |   9 |  17.3 |   212 |
| Q25L60X80P002 |   80.0 |   0.97 |     24816 | 4.45M | 313 |    26.94K |     80 |   9 |  17.7 |   214 |
| Q30L60X40P000 |   40.0 |  0.977 |     49612 | 4.46M | 176 |    19.24K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P001 |   40.0 |  0.979 |     46331 | 4.46M | 164 |    26.95K |     40 |   6 |   7.3 |   116 |
| Q30L60X40P002 |   40.0 |  0.979 |     44717 | 4.46M | 175 |    20.63K |     40 |   6 |   7.3 |   116 |
| Q30L60X80P000 |   80.0 |  0.973 |     30313 | 4.45M | 271 |    24.42K |     80 |   9 |  17.7 |   214 |
| Q30L60X80P001 |   80.0 |  0.973 |     28301 | 4.45M | 281 |    24.25K |     80 |   9 |  17.7 |   214 |

Table: statUnitigsSuperreads.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|-----------|-------:|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| MRX40P000 |   40.0 |  0.976 |     73680 | 4.45M | 116 |    25.76K |     40 |   5 |   8.3 |   110 |
| MRX40P001 |   40.0 |  0.975 |     63584 | 4.45M | 125 |    20.31K |     40 |   5 |   8.3 |   110 |
| MRX40P002 |   40.0 |  0.976 |     73659 | 4.46M | 119 |    20.13K |     40 |   5 |   8.3 |   110 |
| MRX80P000 |   80.0 |  0.975 |     67221 | 4.45M | 122 |    23.32K |     80 |   8 |  18.7 |   208 |
| MRX80P001 |   80.0 |  0.976 |     67271 | 4.46M | 116 |    18.38K |     80 |   8 |  18.7 |   208 |
| MRX80P002 |   80.0 |  0.975 |     64139 | 4.45M | 121 |    21.57K |     80 |   8 |  18.7 |   208 |

Table: statMRUnitigsSuperreads.md

| Name                          | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|-------------------------------|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| 7_merge_anchors               |  0.986 |     78598 | 4.46M | 110 |         0 |    263 |  22 |  65.7 |   658 |
| 7_merge_mr_unitigs_bcalm      |   0.99 |     73710 | 4.45M | 105 |         0 |    262 |  22 |  65.3 |   656 |
| 7_merge_mr_unitigs_bifrost    |  0.991 |     73706 | 4.45M | 104 |         0 |    263 |  22 |  65.7 |   658 |
| 7_merge_mr_unitigs_superreads |   0.99 |     73704 | 4.45M | 105 |         0 |    263 |  22 |  65.7 |   658 |
| 7_merge_unitigs_bcalm         |   0.99 |     73693 | 4.45M | 106 |         0 |    263 |  22 |  65.7 |   658 |
| 7_merge_unitigs_bifrost       |  0.989 |     78591 | 4.46M | 110 |         0 |    263 |  22 |  65.7 |   658 |
| 7_merge_unitigs_superreads    |   0.99 |     73706 | 4.45M | 103 |         0 |    263 |  22 |  65.7 |   658 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |   Sum |   # | SumOthers | median | MAD | lower | upper |
|--------------|-------:|----------:|------:|----:|----------:|-------:|----:|------:|------:|
| 8_spades     |  0.984 |    112448 | 4.47M |  76 |    18.41K |    263 |  22 |  65.7 |   658 |
| 8_mr_spades  |  0.989 |    132590 | 4.49M |  69 |    21.23K |    376 |  28 |  97.3 |   920 |
| 8_megahit    |  0.985 |     67322 | 4.46M | 117 |    25.03K |    263 |  22 |  65.7 |   658 |
| 8_mr_megahit |  0.992 |    132754 |  4.5M |  73 |    21.88K |    376 |  28 |  97.3 |   920 |

Table: statOtherAnchors.md

| Name                     |     N50 |     Sum |   # |
|--------------------------|--------:|--------:|----:|
| Genome                   | 4583637 |   4.58M |   1 |
| Paralogs                 |    1737 | 188.16K | 111 |
| repetitive               |    1175 | 110.59K | 190 |
| 7_merge_anchors.anchors  |   78598 |   4.46M | 110 |
| glue_anchors             |   85585 |   4.46M | 105 |
| fill_anchors             |  102321 |   4.46M |  92 |
| spades.contig            |  114710 |   4.52M | 171 |
| spades.scaffold          |  143522 |   4.52M | 163 |
| spades.non-contained     |  132337 |   4.49M |  78 |
| mr_spades.contig         |  178373 |   4.52M |  87 |
| mr_spades.scaffold       |  203812 |   4.52M |  84 |
| mr_spades.non-contained  |  178373 |   4.51M |  59 |
| megahit.contig           |   85613 |   4.51M | 175 |
| megahit.non-contained    |   85613 |   4.49M | 115 |
| mr_megahit.contig        |  133730 |   4.56M | 154 |
| mr_megahit.non-contained |  133730 |   4.52M |  70 |

Table: statFinal

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

Table: statBusco run_bacteria_odb10

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

Table: statBusco run_enterobacterales_odb10

### dh5alpha: multik 单次调用（--all-masters）全流程门禁（2026-08-18 追加）

承接 `g37/mg1655: olc --cross-validate 跨组嵌合投票`（见上）后的第三条基准
数据集验证。`asm multik --all-masters`（auto 阶梯按 reads N50 生成，
本数据 N50 250 → 31..160）单次调用，13 组 MR anchor 合并口径，
跨组 `olc --unitigs --cross-validate`：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 13 组 anchor + cross-validate + extend 跨 contig 护栏（现行） | 122 | 259019 | 4581642 |  82991 | 0 | 98.342 | 1.016 | 0.00 | 1.27 | 0.20 |
| 13 组 anchor + cross-validate（消除前，2 mis） | 105 | 259450 | 4619600 | 112966 | 2 | 98.848 | 1.019 | 0.00 | 2.71 | 0.24 |
| spades（旧运行） | 78 | 132337 | 4490000 | 112448 | — | — | — | — | — | — |
| mr_spades（旧运行） | 59 | 178373 | 4510000 | 132590 | — | — | — | — | — | — |
| mr_megahit（旧运行） | 70 | 133730 | 4520000 | 132754 | — | — | — | — | — | — |

要点（单组计时与资源）：
* 13 组 MR（MRX40P000-008 + MRX80P000-003）每组合计 multik 20–33 s
  （-p 8、2 组并发）、olc ~75–95 s、extend ~2 s、anchor ~2–3 s；
  峰值内存 8.4–14.9 GB/进程（40× 组 ~8.6–10.5 GB、80× 组
  ~10.8–14.9 GB），2 组并发峰值 ~30 GB（< 机器 88 GB 的 1/2）；
* 跨组 `olc --unitigs --cross-validate` 96.4 s / 820 MB，extend 1.4 s；
  QUAST 用 `quast.py -m 500 -r genome.fa --min-contig 200`。

#### 2 mis 归因（跨组保守重复区 relocation，已消除）

> 结论（2026-08-18）：**misassemblies 已归 0**。两条 contamination 根源于 `asm extend`
> 延伸 walk；由 extend 低覆盖缝检查（contig_4）与 extend 跨 contig 所有权护栏
> （contig_22）消除，详见 `notes/audit/audit-asm-relocation-chimeras.md` §3.4-3.5。
> G37 / MG1655 回归门禁 mis 均保持 0（GF -0.2~0.4 pp、N50 -1%~-11%）。

QUAST 报 2 条 relocation（contig_3 203,757 bp 与 contig_18 96,201 bp）：

| contig | length | junction (contig pos) | reloc 区间 (ref) | inconsistency | 组内存在 |
| ------ | -----: | --------------------: | ---------------- | ------------: | -------- |
| contig_3  | 203757 | 54973/54974 | 54973↔203757 | 1338 | 6/13 组 |
| contig_18 |  96201 |  9203/9204 | 1..9203 ↔ 9204..96201 | 48378 | 13/13 组 |

* **contig_18（13/13 组一致）**：junction 120 bp 序列在参考基因组
  NZ_CP017100 中有两处完全一致的同源拷贝（pos 1,206,312 与 1,254,690，
  相距 48,378 bp）——跨重复区 relocation。**所有 13 组 anchor 都保守地
  连出该 junction**，`--cross-validate` 的"两端各 ≥2 组覆盖且中部无人
  横跨"判定不成立（中部被其他组同样保守地横跨），故无法消除。这是
  基准自身的重复区结构决定，与 G37 MRX40P002 单组低丰度菌株连接
  （junction 仅在 1 组）不同；
* **contig_3（6/13 组存在）**：junction 120 bp 在参考中无同源拷贝
  （非重复区），但 6 组（≥2 组门槛）anchor 一致连出 → cross-validate
  判定为"跨组保守连接"而保留。疑为 reads 覆盖在该位点的一致性结构
  （真实菌株连接或参考 NZ_CP017100 与该菌株序列差异），机制与
  contig_18 不同、待专项归因；
* 单组 multik 输出与 08-17 版本字节级一致（`MRX40P000` 重跑 diff 为空），
  本记录可复现：`/tmp/dh5alpha_gate/`（run.sh 复刻 model_org.md 标准
 门禁链）。

## *Bacillus cereus* ATCC 10987

> GAGE-B MiSeq 数据集（100× 子集），reads 与参考**同株**（见 reference
> 节实测）。老流程基线（2025-06）来自 `results/gage_b.md`；现代流程门禁待跑。

### bcer: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/Bcer_100x/1_genome
cd ~/data/anchr/Bcer_100x/1_genome

cp ~/data/anchr/ref/Bcer/genome.fa .
cp ~/data/anchr/ref/Bcer/paralogs.fa .
cp ~/data/anchr/ref/Bcer/repetitive.fa .

```

* 参考构成：NC_003909（染色体 5,224,283 bp）+ NC_005707（质粒 pBc10987，
  208,369 bp）= 5,432,652 bp / 2 复制子；GC ~38%（低 GC）。
* **reads 与参考同株**：GAGE-B 数据集即 ATCC 10987。本地实测（50k 原始
  R1，2026-08-18）：reads 的 solid 31-mer 99.66% 出现在参考中、85.7%
  整条 perfect 回贴——菌株差异≈0，QUAST 判据可直接使用（对比 Mabs/Vcho
  的跨株情况，见下）。

### bcer: download

* Illumina（GAGE-B 100× 子集）

```shell
cd ~/data/anchr/Bcer_100x

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 4 -s 2 -c https://ccb.jhu.edu/gage_b/datasets/B_cereus_MiSeq.tar.gz

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

* 本地已就位：`~/data/anchr/Bcer_100x/2_illumina/` 有 R1.fq.gz + R2.fq.gz
  （2024-12 下载）；`1_genome/` 尚空，按 reference 节 cp 即可。

### bcer: template

* 现代流程模板（与 g37/dh5alpha 同口径；待跑）

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 5432652 \
    --parallel 24 \
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
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig bcalm" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

* 老流程模板（2025-06 GAGE-B 运行，基线数据见下节）：`--genome 5432652
  --fastqc --insertsize --kat --trim "--dedupe --tile --cutoff 5 --cutk 31"
  --qual "20 25 30" --len "60" --filter "adapter artifact" --quorum --merge
  --ecphase "1 2 3" --cov "40 50 60 all" --unitigger "superreads bcalm
  tadpole" --statp 2 --readl 250 --uscale 2 --lscale 3 --redo --extend`
  （完整命令见 `results/gage_b.md`）。

### bcer: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_100x

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

现代流程门禁**待跑**（本地数据已齐，可直接执行）；老流程 2025-06 运行见
`results/gage_b.md`。

### bcer: 老流程基线（2025-06 GAGE-B）

reads：2,080,000 对 × 250 bp（~100×，481.02 Mb）；接头污染极低
（#R.trim 匹配 0.06267%）。老流程 quorum（Q0L0）丢弃 12.98%。

Table: statReads（老流程）

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 5224283 | 5432652 |       2 |
| Paralogs    |    2295 |  220468 |     101 |
| repetitive |    2461 |  113050 |     173 |
| Illumina.R  |     251 | 481.02M | 2080000 |
| trim.R      |     250 | 404.36M | 1807384 |
| Q20L60      |     250 | 396.79M | 1758525 |
| Q25L60      |     250 | 379.57M | 1706208 |
| Q30L60      |     250 | 344.25M | 1611233 |

Table: statTrimReads（老流程）

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

Table: statMergeReads（老流程）

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

Table: statQuorum（老流程）

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   |  74.4 |   64.8 |   12.98% | "127" | 5.43M | 5.35M |     0.98 | 0:01'00'' |
| Q20L60.R |  73.0 |   64.7 |   11.48% | "127" | 5.43M | 5.35M |     0.98 | 0:00'52'' |
| Q25L60.R |  69.9 |   63.8 |    8.71% | "127" | 5.43M | 5.34M |     0.98 | 0:00'50'' |
| Q30L60.R |  63.4 |   59.7 |    5.75% | "127" | 5.43M | 5.34M |     0.98 | 0:00'49'' |

Table: statMRUnitigsSuperreads（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.62% |     37925 | 5.31M | 240 |        81 | 22.62K | 519 |   39.0 |  7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:34 |
| MRX40P001  |   40.0 |  97.54% |     40827 | 5.31M | 237 |        64 | 18.73K | 454 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:32 |
| MRX50P000  |   50.0 |  97.58% |     37545 | 5.31M | 249 |        90 | 22.64K | 547 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:36 |
| MRX60P000  |   60.0 |  97.53% |     37952 | 5.31M | 250 |        92 | 20.41K | 527 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:37 |
| MRXallP000 |   85.9 |  97.46% |     36340 | 5.33M | 253 |        90 |  18.9K | 537 |   84.0 | 15.0 |  13.0 | 258.0 | "31,41,51,61,71,81" |   0:02:00 |   0:00:35 |

Table: statMRUnitigsBcalm（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.81% |     39857 | 5.34M | 232 |        87 | 18.83K | 442 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:32 |
| MRX40P001  |   40.0 |  97.77% |     42781 | 5.31M | 229 |        86 | 18.97K | 412 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:31 |
| MRX50P000  |   50.0 |  97.44% |     39857 | 5.34M | 236 |        80 | 16.34K | 466 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:33 |
| MRX60P000  |   60.0 |  97.42% |     41648 | 5.38M | 238 |        98 | 17.86K | 463 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:35 |
| MRXallP000 |   85.9 |  97.38% |     39857 | 5.31M | 246 |        97 | 16.29K | 489 |   85.0 | 15.0 |  13.3 | 260.0 | "31,41,51,61,71,81" |   0:01:54 |   0:00:35 |

Table: statMRUnitigsTadpole（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.83% |     44440 | 5.32M | 217 |       107 |  18.4K | 404 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:32 |
| MRX40P001  |   40.0 |  97.77% |     44376 | 5.36M | 219 |        83 | 16.18K | 378 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:33 |
| MRX50P000  |   50.0 |  97.76% |     42799 | 5.34M | 224 |       131 | 19.98K | 444 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:32 |
| MRX60P000  |   60.0 |  97.77% |     42132 | 5.34M | 228 |       131 | 19.13K | 443 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:34 |
| MRXallP000 |   85.9 |  97.70% |     41648 | 5.32M | 239 |       120 | 18.37K | 484 |   85.0 | 15.0 |  13.3 | 260.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:35 |

Table: statMergeAnchors（老流程）

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.76% |     36737 |  5.3M | 254 |     41729 | 291.14K | 43 |   64.0 | 12.0 |   9.3 | 200.0 |   0:00:42 |
| 7_merge_mr_unitigs_bcalm      |  98.15% |     37087 | 5.21M | 252 |     41793 |  70.11K |  9 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:49 |
| 7_merge_mr_unitigs_superreads |  97.83% |     34002 | 5.15M | 259 |      1091 |    6.8K |  7 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:43 |
| 7_merge_mr_unitigs_tadpole    |  98.14% |     42769 | 5.21M | 232 |     41634 |  53.78K | 12 |   63.0 | 13.0 |   8.0 | 204.0 |   0:00:43 |
| 7_merge_unitigs_bcalm         |  98.04% |     32860 | 5.31M | 273 |     36517 |  63.09K | 19 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:49 |
| 7_merge_unitigs_superreads    |  98.10% |     32344 | 5.31M | 273 |     41729 |  64.41K | 24 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:46 |
| 7_merge_unitigs_tadpole       |  98.02% |     32710 | 5.31M | 277 |     62295 |  81.41K | 21 |   64.0 | 12.0 |   9.3 | 200.0 |   0:00:45 |

Table: statOtherAnchors（老流程）

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.96% |     49586 | 2.14M |  80 |       899 |  9.35K | 106 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:33 |
| 8_mr_spades  |  98.81% |     78041 | 5.35M | 130 |       217 | 10.11K | 212 |   85.0 | 15.0 |  13.3 | 260.0 |   0:00:36 |
| 8_megahit    |  98.69% |     37648 | 4.85M | 218 |       125 | 21.73K | 356 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:40 |
| 8_mr_megahit |  98.85% |     65403 | 5.36M | 153 |       451 | 13.53K | 265 |   85.0 | 15.0 |  13.3 | 260.0 |   0:00:36 |

Table: statFinal（老流程）

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 5224283 | 5432652 |   2 |
| Paralogs                 |    2295 |  220468 | 101 |
| repetitive              |    2461 |  113050 | 173 |
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

### bcer: 现代流程门禁（待跑）

* 数据已就位（本地 R1/R2 + ref），按 template/run 节命令执行即可；
* 关注点：低 GC（~38%）5.43 Mb 基因组（首个革兰氏阳性条目）下
  `multik --all-masters` + 跨组 `olc --unitigs --cross-validate` 的表现；
  老流程基线 fill N50 61,075 / 159 条（Sum 5,317,582 ≈ 参考 97.9%）；
* 复现记录建议放 `/tmp/bcer_gate/`（参照 dh5alpha_gate 的门禁链格式）。

## *Rhodobacter sphaeroides* 2.4.1

> GAGE-B MiSeq 数据集（100× 子集）。参考即 2.4.1（名义同株，无本地
> reads 实测）；结构复杂（2 染色体 + 5 质粒 = 7 复制子、重复 12.8%、
> GC ~68%）正是纳入理由——对"无 N 染色体"目标的多复制子压力测试。
> 老流程基线（2025-06）来自 `results/gage_b.md`；现代流程门禁待跑。

### rsph: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/Rsph_100x/1_genome
cd ~/data/anchr/Rsph_100x/1_genome

cp ~/data/anchr/ref/Rsph/genome.fa .
cp ~/data/anchr/ref/Rsph/paralogs.fa .
cp ~/data/anchr/ref/Rsph/repetitive.fa .

```

* 参考构成：NC_007493（染色体 1，3,188,524 bp）+ NC_007494（染色体 2，
  1,314,453 bp）+ 5 个质粒（NC_009007/NC_007488/NC_007489/NC_007490/
  NC_009008）= 4,602,977 bp / 7 复制子；GC ~68%（极高 GC），重复含量
  12.8%（GAGE-B 四菌株中最高）。

### rsph: download

* Illumina（GAGE-B 100× 子集）

```shell
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

* 本地暂无 reads（`~/data/anchr/Rsph_100x/` 为空），需先下载。

### rsph: template

* 现代流程模板（与 g37/dh5alpha 同口径；待跑）

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 4602977 \
    --parallel 24 \
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
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig bcalm" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

* 老流程模板（2025-06 GAGE-B 运行，基线数据见下节）：`--genome 4602977
  --fastqc --insertsize --kat --trim "--dedupe" --qual "20 25 30" --len
  "60" --filter "adapter artifact" --quorum --merge --ecphase "1 3" --cov
  "30 all" --unitigger "superreads bcalm tadpole" --statp 2 --readl 250
  --uscale 2 --lscale 3 --redo --extend`（完整命令见
  `results/gage_b.md`）。

### rsph: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_100x

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

现代流程门禁**待跑**（需先按 download 节取数）；老流程 2025-06 运行见
`results/gage_b.md`（QUAST 用 `--no-check`，多复制子参考）。

### rsph: 老流程基线（2025-06 GAGE-B）

reads：1,800,000 对 × 250 bp（~100×，451.8 Mb）；接头匹配 6.39%
（Reverse_adapter 4.58% 为主）。老流程 quorum（Q0L0）丢弃 11.10%。

Table: statReads（老流程）

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 3188524 | 4602977 |       7 |
| Paralogs    |    2337 |  146789 |      66 |
| repetitive |     572 |   57281 |     165 |
| Illumina.R  |     251 |  451.8M | 1800000 |
| trim.R      |     148 |  200.1M | 1452706 |
| Q20L60      |     148 | 193.66M | 1401466 |
| Q25L60      |     139 | 169.12M | 1304628 |
| Q30L60      |     119 | 125.02M | 1123194 |

Table: statTrimReads（老流程）

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

Table: statMergeReads（老流程）

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

Table: statQuorum（老流程）

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   |  43.5 |   38.7 |   11.10% | "39" |  4.6M | 4.55M |     0.99 | 0:00'37'' |
| Q20L60.R |  42.1 |   37.9 |    9.98% | "39" |  4.6M | 4.55M |     0.99 | 0:00'32'' |
| Q25L60.R |  36.8 |   34.9 |    5.03% | "35" |  4.6M | 4.54M |     0.99 | 0:00'29'' |
| Q30L60.R |  27.2 |   26.6 |    2.20% | "31" |  4.6M | 4.52M |     0.98 | 0:00'26'' |

Table: statMRUnitigsSuperreads（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.39% |     20508 | 4.35M | 364 |      4332 |  140.5K | 739 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRX30P001  |   30.0 |  97.40% |     21605 | 4.34M | 357 |      3817 | 154.53K | 743 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRXallP000 |   60.4 |  97.36% |     22279 | 4.34M | 345 |      5800 | 167.85K | 749 |   55.0 | 9.0 |   9.3 | 164.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:29 |

Table: statMRUnitigsBcalm（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.41% |     18774 | 4.35M | 404 |      5004 | 173.16K | 899 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:30 |
| MRX30P001  |   30.0 |  97.35% |     19071 | 4.34M | 389 |      6101 | 196.85K | 857 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:28 |
| MRXallP000 |   60.4 |  97.44% |     20747 | 4.33M | 353 |      6101 |  199.4K | 792 |   55.0 | 8.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:30 |

Table: statMRUnitigsTadpole（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.50% |     19948 | 4.34M | 379 |      5224 | 161.26K | 801 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:28 |
| MRX30P001  |   30.0 |  97.55% |     19862 | 4.34M | 372 |      5475 | 173.65K | 799 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:29 |
| MRXallP000 |   60.4 |  97.54% |     21368 | 4.33M | 344 |      6100 | 177.53K | 759 |   55.0 | 8.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:29 |

Table: statMergeAnchors（老流程）

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  91.17% |     36740 | 4.41M | 246 |      4689 | 360.25K | 133 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:33 |
| 7_merge_mr_unitigs_bcalm      |  92.62% |     21353 | 4.33M | 349 |      5224 | 203.81K |  65 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:35 |
| 7_merge_mr_unitigs_superreads |  92.48% |     22986 | 4.39M | 328 |      5487 | 188.38K |  62 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:34 |
| 7_merge_mr_unitigs_tadpole    |  92.68% |     22194 | 4.33M | 335 |      5819 | 199.44K |  57 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:35 |
| 7_merge_unitigs_bcalm         |  90.95% |     17117 | 4.33M | 424 |      5224 | 311.57K | 102 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:30 |
| 7_merge_unitigs_superreads    |  91.58% |     30174 | 4.35M | 271 |      4689 | 308.94K | 114 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:31 |
| 7_merge_unitigs_tadpole       |  91.26% |     22998 | 4.34M | 319 |      5018 | 335.15K | 115 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:30 |

Table: statOtherAnchors（老流程）

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  99.20% |     45176 | 1.94M |  83 |      7794 |  86.43K | 125 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:29 |
| 8_mr_spades  |  99.13% |     30528 | 3.93M | 216 |     16123 |  151.9K | 304 |   55.0 |  8.0 |  10.3 | 158.0 |   0:00:31 |
| 8_megahit    |  98.54% |     37964 |  4.1M | 230 |      5525 | 153.99K | 387 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:28 |
| 8_mr_megahit |  99.34% |     26589 | 4.38M | 287 |     16123 | 164.89K | 523 |   55.0 |  8.0 |  10.3 | 158.0 |   0:00:31 |

Table: statFinal（老流程）

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 3188524 | 4602977 |    7 |
| Paralogs                 |    2337 |  146789 |   66 |
| repetitive              |     572 |   57281 |  165 |
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

### rsph: 现代流程门禁（待跑）

* 需先下载 reads（ccb.jhu.edu），再按 template/run 节命令执行；
* 关注点：7 复制子（2 染色体 + 5 质粒）下 GF/dup 统计会被质粒放大，
  QUAST 用 `--no-check`（老流程即如此），mis 判定看染色体级骨架；
  `--cross-validate` 跨组投票对高重复（12.8%）菌株是核心验证；
* 老流程基线 fill N50 48,535 / 183 条（Sum 4,406,444 ≈ 参考 95.7%），
  anchors Mapped% 仅 91.17%（GAGE-B 四菌株最低，重复/质粒干扰）；
* 复现记录建议放 `/tmp/rsph_gate/`（参照 dh5alpha_gate 的门禁链格式）。
