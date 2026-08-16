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

模板 multik 分支自 2026-08-16 起自动包含：主 K `31..121 128`（6_ MR 链；
4_ 链 `31..91`，bcalm 分支 `31..121`）、multik 气泡合并（默认
`--merge-similar 0.95 --merge-len 20`）、olc `--min-contig-len 200`（仅
multik 分支，bcalm/unitig 分支保持 1000）、`asm extend --min-len 1000`
（短于 1000 bp 的 contig 不延伸，避免重复区嵌合）。详细机制与中间实验见
`notes/benchmarks/g37-megahit-spades.md` §7/§8/§10。

7 组 MR（MRX40P000-004 + MRX80P000-001）QUAST：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_mr_multik（旧，K31..81） | 17 | 114442 | 565518 |  55049 | 0 | 96.963 | 1.003 | 0.00 | 77.46 | 19.14 |
| merge_mr_multik（新链） | 18 | 236494 | 582604 | 83756 | 0 | 98.983 | 1.002 | 0.00 | 226.03 | 66.37 |
| mr_spades（旧运行） |  2 | 580506 | 580632 | 580506 | 0 | 100.000 | 1.001 | 0.00 | 300.77 | 91.64 |
| mr_megahit（旧运行） | 44 | 319186 | 599025 | 319186 | 3 | 99.893 | 1.025 | 0.00 | 311.14 | 88.39 |

要点：
* 新链 N50 55.0K→**83.8K**（+52%）、GF 96.963→**98.983%**、**0 mis 保持**，
  并跨复制起点接出 236.5K 最大 contig；
* mm/100k 上升（77→226）是 QUAST 对参考的口径问题：新覆盖集中在低复杂
  区，reads 与 NC_000908 在那里差异 ~5%，consensus 忠实反映 reads
  （reads-vs-contigs 一致率 99.96%）；bwa 口径下真缺口仅 ~1.6K bp；
* QUAST minimap 会漏对齐低复杂度区短 contig，报告的未覆盖 bp 偏大。

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
`MRX80P000/P001`），当前 multik 链（K31..121+128 + 气泡合并 + extend
`--min-len 1000` + `--min-contig-len 200`）QUAST：

| Assembly | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| 新 multik 链（5 组） | 88 | 268332 | 4573304 | 123988 | 0 | 98.364 | 1.002 | 0.00 | 0.22 | 0.09 |
| 旧 multik 链（同 5 组，K31..81） | 107 | — | — | 95478 | 0 | 97.61 | ~1.00 | — | — | — |
| merge_mr_multik（旧运行，全组） | 317 | 73037 | 4527602 | 23594 | 5 | 96.49 | 1.011 | 0.00 | 0.24 | 0.02 |
| spades（旧运行） | 125 | 224028 | 4572988 | 125607 | 0 | 98.47 | 1.000 | 0.00 | 0.74 | 0.17 |
| mr_spades（旧运行） | 152 | 284843 | 4588125 | 148607 | 0 | 98.64 | 1.002 | 0.00 | 0.92 | 0.17 |
| megahit（旧运行） | 273 | 175838 | 4588105 | 43891 | 92 | 98.26 | 1.004 | 0.00 | 3.73 | 0.50 |
| mr_megahit（旧运行） | 138 | 311797 | 4611104 | 126312 | 2 | 98.89 | 1.004 | 0.00 | 2.13 | 0.26 |

要点（详细分析见 `notes/benchmarks/mg1655-process-compare.md`）：
* 新链 N50 95.5K→**124.0K**（+30%）、GF 97.61→**98.364%**、**0 mis**；
  同 5 组输入口径下已超过 megahit（82.8K / 1 mis）并与 spades 持平，
  仅 mr_spades（148.6K，全量输入）更高；
* **extend 必须 `--min-len 1000`**：extend 短碎片会把 1.2 Mb 处重复元件
  拷贝接成嵌合体（238 bp 碎片被长成 1,238 bp relocation，3 mis）；加门槛
  后 0 mis（G37 同步验证）；
* 处理过程中顺带修复 `asm olc` consensus 的 O(n²×L) 性能瓶颈（种子索引
  预筛，输出逐位一致），MG1655 单组 9 主池从数小时降到 3 分钟。
